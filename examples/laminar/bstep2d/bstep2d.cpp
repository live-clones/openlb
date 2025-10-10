/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006, 2007, 2012, 2025 Jonas Latt, Mathias J. Krause,
 *  Louis Kronberg, Christian Vorwerk, Bastian Schäffauer, Yuji (Sam) Shimojima
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

/* bstep2d.cpp:
 * The implementation of a backward facing step. It is furthermore
 * shown how to use checkpointing to save the state of the
 * simulation regularly.
 * The geometry of the step is based on the experiment described in
 * [Armaly, B.F., Durst, F., Pereira, J. C. F. and Schönung, B. Experimental
 * and theoretical investigation of backward-facing step flow. 1983.
 * J. Fluid Mech., vol. 127, pp. 473-496, DOI: 10.1017/S0022112083002839]
 */

#include <olb.h>
using namespace olb;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<NavierStokes, Lattice<double, descriptors::D2Q9<>>>;

namespace olb::parameters {
  struct PHYS_LENGTH_OF_STEP : public descriptors::FIELD_BASE<1> {};
  struct PHYS_HEIGHT_OF_STEP : public descriptors::FIELD_BASE<1> {};
  struct PHYS_ORIGIN_OF_STEP : public descriptors::FIELD_BASE<1> {};
} // namespace olb::parameters

/// @brief Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters)
{
  using T = MyCase::value_t;
  // setup channel
  const Vector extendChannel = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector originChannel = parameters.get<parameters::ORIGIN>();
  std::shared_ptr<IndicatorF2D<T>> channel = std::make_shared<IndicatorCuboid2D<T>>(extendChannel, originChannel);
  // setup step
  const T lengthStep      = parameters.get<parameters::PHYS_LENGTH_OF_STEP>(); // length of step in meter
  const T heightStep      = parameters.get<parameters::PHYS_HEIGHT_OF_STEP>(); // height of step in meter
  const Vector originStep = parameters.get<parameters::PHYS_ORIGIN_OF_STEP>(); // origin of step in meter

  const Vector                     extendStep(lengthStep, heightStep);
  std::shared_ptr<IndicatorF2D<T>> step = std::make_shared<IndicatorCuboid2D<T>>(extendStep, originStep);

  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();

  Mesh<T, MyCase::d> mesh(*(channel - step), physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

/// @brief Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
/// @note The material numbers will be used to assign physics to lattice nodes
void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;
  using T          = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry  = myCase.getGeometry();
  // Parameters for the simulation setup

  // setup channel
  const Vector extendChannel = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector originChannel = parameters.get<parameters::ORIGIN>();
  std::shared_ptr<IndicatorF2D<T>> channel = std::make_shared<IndicatorCuboid2D<T>>(extendChannel, originChannel);

  // setup step
  const T lengthStep = parameters.get<parameters::PHYS_LENGTH_OF_STEP>(); // length of step in meter
  const T heightStep = parameters.get<parameters::PHYS_HEIGHT_OF_STEP>(); // height of step in meter

  const Vector                     extendStep(lengthStep, heightStep);
  const Vector                     originStep(0, 0);
  std::shared_ptr<IndicatorF2D<T>> step = std::make_shared<IndicatorCuboid2D<T>>(extendStep, originStep);

  const T physDeltaX = parameters.get<parameters::PHYS_CHAR_LENGTH>() / parameters.get<parameters::RESOLUTION>();
  // material numbers from zero to 2 inside geometry defined by indicator
  geometry.rename(0, 2, channel - step);
  geometry.rename(2, 1, {1, 1});
  const T      lengthChannel = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T      heightChannel = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const T      heightInlet   = heightChannel - heightStep;
  Vector<T, 2> extendBC_out((T)0.0 + (T)1.0 * physDeltaX, heightChannel);
  Vector<T, 2> extendBC_in((T)0.0, heightInlet);
  Vector<T, 2> originBC_out(lengthChannel - (T)1.0 * physDeltaX, 0);
  Vector<T, 2> originBC_in((T)0.0, heightStep);

  IndicatorCuboid2D<T> inflow(extendBC_in, originBC_in);
  // Set material number for inflow
  geometry.rename(2, 3, 1, inflow);

  IndicatorCuboid2D<T> outflow(extendBC_out, originBC_out);
  // Set material number for outflow
  geometry.rename(2, 4, 1, outflow);

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  // Removes all not needed boundary voxels inside the surface
  geometry.innerClean();
  geometry.checkForErrors();
  geometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

/// @brief Set lattice dynamics
/// @param myCase The Case instance which keeps the simulation data
void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T = MyCase::value_t;

  auto& geometry   = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes{});
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  {
    using namespace olb::parameters;
    lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR>>(
        parameters.get<RESOLUTION>(),       // resolution
        parameters.get<LATTICE_RELAXATION_TIME>(),  // relaxation time
        parameters.get<PHYS_CHAR_LENGTH>(), // charPhysLength: reference length of simulation geometry
        parameters.get<PHYS_CHAR_VELOCITY>(), // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
        parameters.get<PHYS_CHAR_VISCOSITY>(), // physViscosity: physical kinematic viscosity in __m^2 / s__
        parameters.get<PHYS_CHAR_DENSITY>()    // physDensity: physical density in __kg / m^3__
    );
  }
  const auto& converter = lattice.getUnitConverter();

  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("bstep2d");

  auto bulkIndicator = geometry.getMaterialIndicator({1, 3, 4});

  // Material=1 -->bulk dynamics
  // Material=3 -->bulk dynamics (inflow)
  // Material=4 -->bulk dynamics (outflow)
  lattice.defineDynamics<BGKdynamics<T, DESCRIPTOR>>(bulkIndicator);
  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(lattice, geometry, 2);

  //if boundary conditions are chosen to be local
  boundary::set<boundary::LocalVelocity>(lattice, geometry, 3);
  boundary::set<boundary::LocalPressure>(lattice, geometry, 4);

  //if boundary conditions are chosen to be interpolated
  // boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 3);
  // boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 4);

  clout << "Prepare Lattice ... OK" << std::endl;
}

/// Set initial condition for primal variables (velocity and density)
/// @param myCase The Case instance which keeps the simulation data
/// @note Be careful: initial values have to be set using lattice units
void setInitialValues(MyCase& myCase)
{
  OstreamManager clout(std::cout, "Initialization");
  clout << "Initialize lattice ..." << std::endl;
  using T         = MyCase::value_t;
  auto& lattice  = myCase.getLattice(NavierStokes {});
  auto& geometry = myCase.getGeometry();

  // Initial conditions
  AnalyticalConst2D<T, T>    ux(0.);
  AnalyticalConst2D<T, T>    uy(0.);
  AnalyticalConst2D<T, T>    rho(1.);
  AnalyticalComposed2D<T, T> u(ux, uy);

  //Initialize all values of distribution functions to their local equilibrium
  auto bulkIndicator = geometry.getMaterialIndicator({1, 3, 4});
  lattice.defineRhoU(bulkIndicator, rho, u);
  lattice.iniEquilibrium(bulkIndicator, rho, u);

  const T omega = lattice.getUnitConverter().getLatticeRelaxationFrequency();
  lattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  lattice.initialize();
  clout << "Initialize lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setTemporalValues(MyCase& myCase, std::size_t iT)
{
  OstreamManager clout(std::cout, "setTemporalValues");
  using T = MyCase::value_t;

  auto&   lattice    = myCase.getLattice(NavierStokes {});
  auto&   converter  = lattice.getUnitConverter();
  auto&   geometry   = myCase.getGeometry();
  auto&   parameters = myCase.getParameters();
  const T maxPhysT   = parameters.get<parameters::MAX_PHYS_T>();
  const T physStartT = parameters.get<parameters::PHYS_START_T>();
  const T iTUpdate   = parameters.get<parameters::PHYS_BOUNDARY_VALUE_UPDATE_T>();

  // time for smooth start-up
  std::size_t iTmaxStart = converter.getLatticeTime(physStartT);

  if (iT % iTupdate == 0 && iT <= iTmaxStart) {
    // Smooth start curve, sinus
    // SinusStartScale<T,std::size_t> StartScale(iTmaxStart, (T)1);
    // Smooth start curve, polynomial
    PolynomialStartScale<T, std::size_t> StartScale(iTmaxStart, T(1));
    // Creates and sets the Poiseuille inflow profile using functors
    std::size_t iTvec[1] = {iT};
    T           frac[1]  = {};
    StartScale(frac, iTvec);
    T               maxVelocity   = converter.getCharLatticeVelocity() * (T)3.0 / (T)2.0 * frac[0];
    T               distance2Wall = converter.getPhysDeltaX() / (T)2.0;
    Poiseuille2D<T> poiseuilleU(geometry, 3, maxVelocity, distance2Wall);
    // define lattice speed on inflow
    lattice.defineU(geometry, 3, poiseuilleU);

    lattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
        ProcessingContext::Simulation
    );
  }
}

void getResults(MyCase& myCase, std::size_t iT, util::Timer<MyCase::value_t> &timer)
{
  OstreamManager clout(std::cout, "getResults");
  using T = MyCase::value_t;
  auto&          parameters    = myCase.getParameters();
  const T        heightChannel = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const T        heightStep    = parameters.get<parameters::PHYS_HEIGHT_OF_STEP>(); // height of step in meter
  const T        lengthStep    = parameters.get<parameters::PHYS_LENGTH_OF_STEP>(); // length of step in meter
  const T        heightInlet   = heightChannel - heightStep;
  auto&          lattice       = myCase.getLattice(NavierStokes{});
  auto&          converter     = lattice.getUnitConverter();
  auto&          geometry      = myCase.getGeometry();

  // instantiate reusable functors
  SuperPlaneIntegralFluxVelocity2D<T> velocityFlux(lattice, converter, geometry,
                                                   {lengthStep / (T)2.0, heightInlet / (T)2.0}, {(T)0.0, (T)1.0});

  SuperPlaneIntegralFluxPressure2D<T> pressureFlux(lattice, converter, geometry,
                                                   {lengthStep / (T)2.0, heightInlet / (T)2.0}, {(T)0.0, (T)1.0});

  SuperVTMwriter2D<T>                 vtmWriter("bstep2d");

  const int vtkIter  = converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T>());
  const int statIter = converter.getLatticeTime(parameters.get<parameters::PHYS_STAT_ITER_T>());

  if (iT == 0) {
    // Writes geometry, cuboid no. and rank no. to file system
    SuperLatticeCuboid2D cuboid(lattice);
    SuperLatticeRank2D   rank(lattice);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();
  }

  // Writes every 0.1 simulated
  if (iT % statIter == 0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    velocityFlux.print();
    pressureFlux.print();

    // write to terminal
    timer.update(iT);
    timer.printStep();
    // Lattice statistics console output
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
  }

  if (iT % vtkIter == 0) {
    SuperLatticePhysVelocity2D velocity(lattice, converter);
    SuperLatticePhysPressure2D pressure(lattice, converter);
    vtmWriter.addFunctor(velocity);
    vtmWriter.addFunctor(pressure);
    // write vtk to file system
    vtmWriter.write(iT);
    using T = MyCase::value_t_of<NavierStokes>;
    SuperEuklidNorm2D     normVel(velocity);
    BlockReduction2D2D<T> planeReduction(normVel, 1200, BlockDataSyncMode::ReduceOnly);
    // write output as JPEG
    heatmap::plotParam<T> jpeg_Param;
    jpeg_Param.maxValue       = converter.getCharPhysVelocity() * 3. / 2.;
    jpeg_Param.minValue       = 0.0;
    jpeg_Param.fullScreenPlot = true;
    heatmap::write(planeReduction, iT, jpeg_Param);
  }
}

void simulate(MyCase& myCase)
{
  OstreamManager clout(std::cout, "simulate");

  using T             = MyCase::value_t;
  auto&   parameters  = myCase.getParameters();
  auto&   lattice     = myCase.getLattice(NavierStokes {});
  const T maxPhysT    = parameters.get<parameters::MAX_PHYS_T>();
  util::Timer<T> timer(lattice.getUnitConverter().getLatticeTime(maxPhysT),
                       myCase.getGeometry().getStatistics().getNvoxel());

  clout << "starting simulation..." << std::endl;
  timer.start();

  for (std::size_t iT = 0; iT < lattice.getUnitConverter().getLatticeTime(maxPhysT); ++iT) {

    setTemporalValues(myCase, iT);

    lattice.collideAndStream();

    getResults(myCase, iT, timer);
  }

  lattice.setProcessingContext(ProcessingContext::Evaluation);
  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);
  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION>(60);
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.518);
    myCaseParameters.set<PHYS_LENGTH_OF_STEP>(0.2);
    myCaseParameters.set<PHYS_HEIGHT_OF_STEP>(4.9e-3);
    myCaseParameters.set<PHYS_CHAR_VELOCITY>(1.0);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(5.2e-6);
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1.0);
    myCaseParameters.set<DOMAIN_EXTENT>({0.7, 0.0101});
    myCaseParameters.set<ORIGIN>({0., 0.});
    myCaseParameters.set<PHYS_ORIGIN_OF_STEP>({0., 0.});
    myCaseParameters.set<MAX_PHYS_T>(2.0);
    myCaseParameters.set<PHYS_START_T>([&]{return
      0.2 * myCaseParameters.get<MAX_PHYS_T>();
    });
    myCaseParameters.set<PHYS_BOUNDARY_VALUE_UPDATE_T>(1);
    myCaseParameters.set<PHYS_CHAR_LENGTH>([&]{return 2.0 *
      (myCaseParameters.get<DOMAIN_EXTENT>()[1] - myCaseParameters.get<PHYS_HEIGHT_OF_STEP>());
    });
    myCaseParameters.set<PHYS_CHAR_LENGTH>([&]{return
      parameters.get<PHYS_CHAR_LENGTH>() / parameters.get<RESOLUTION>()};
    );
  }
  myCaseParameters.fromCLI(argc, argv);

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);
  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase);

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);

  return 0;
}
