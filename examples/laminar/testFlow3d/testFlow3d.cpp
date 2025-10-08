/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006-2025 Mathias J. Krause, Fabian Klemens,
 *  Julius Jessberger, Shota Ito
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


/** @file testFlow3d.cpp
 * @brief In this example, a 3d flow in a cube is simulated.
 *
 * It shows the basic structure of an OpenLB simulation and should be suitable for beginners.
 *
 * As every OpenLB simulation, it consists of eight steps:
 * @li Step 1: Declarations
 * @li Step 2: Initialization
 * @li Step 3: Create mesh
 * @li Step 4: Create case
 * @li Step 5: Prepare geometry
 * @li Step 6: Prepare lattice
 * @li Step 7: Set initial values
 * @li Step 8: Simulate
 *
 * For more information, we refer to the OpenLB user guide.
 */

/// Include OpenLB library and load analytical solution
#include <olb.h>
#include "analyticalSolutionTestFlow3D.h"

using namespace olb;
using namespace olb::names;

/// @brief Step 1: Declare simulation structure.
/// Model name and lattice type are collected in a Case class
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q19<descriptors::FORCE>>
>;

/// @brief Step 3: Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector origin{0, 0, 0};
  IndicatorCuboid3D<T> cuboid(extent, origin);

  const T physDeltaX = extent[0] / parameters.get<parameters::RESOLUTION>();
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

/// @brief Step 5: Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
/// @note The material numbers are used to assign physics to lattice nodes
void prepareGeometry(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;
  auto& geometry = myCase.getGeometry();

  /// @li Set material numbers
  geometry.rename(0, 2);
  geometry.rename(2, 1, {1,1,1});
  geometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

/// @brief Step 6: Set lattice dynamics and boundary conditions
/// @param myCase The Case instance which keeps the simulation data
void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& parameters = myCase.getParameters();

  /// @li Set up a unit converter with the characteristic physical units
  lattice.setUnitConverter<UnitConverterFromResolutionAndLatticeVelocity<MyCase::value_t,
                                                                         MyCase::descriptor_t_of<NavierStokes>>>(
    parameters.get<parameters::RESOLUTION>(),               // resolution
    parameters.get<parameters::LATTICE_CHAR_VELOCITY>(),    // charLatticeVelocity
    parameters.get<parameters::DOMAIN_EXTENT>()[0],         // charPhysLength: reference length of simulation geometry in [m]
    parameters.get<parameters::PHYS_CHAR_VELOCITY>(),       // charPhysVelocity: highest expected velocity during simulation in [m/s]
    parameters.get<parameters::PHYS_CHAR_VISCOSITY>(),      // physViscosity: physical kinematic viscosity in [m^2/s]
    parameters.get<parameters::PHYS_CHAR_DENSITY>()         // physDensity: physical density [kg/m^3]
  );
  auto& converter = lattice.getUnitConverter();
  converter.print();

  /// @li Material=1 --> bulk dynamics
  lattice.template defineDynamics<ForcedBGKdynamics>(geometry.getMaterialIndicator({1}));
  /// @li Material=2,3 --> velocity boundary
  boundary::set<boundary::LocalVelocity>(lattice, geometry.getMaterialIndicator({2}));
  /// @li Set lattice relaxation frequency
  lattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
}

// Step 7: Set initial values for primal variables (e.g. velocity, density) and additional fields
/// @param myCase The Case instance which keeps the simulation data
/// @note Initial values have to be set using lattice units
void setInitialValues(MyCase& myCase)
{
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();

  /// @li Initialize density to one everywhere
  /// @li Initialize velocity to be tangential at the lid and zero in the bulk and at the walls
  AnalyticalConst3D<T,T> rhoF(1);
  const Vector<T,3> u{0,0,0};
  AnalyticalConst3D<T,T> uF(u);

  /// @li Initialize populations to equilibrium state
  lattice.defineRhoU(geometry.getMaterialIndicator({1,2}), rhoF, uF);
  lattice.iniEquilibrium(geometry.getMaterialIndicator({1,2}), rhoF, uF);

  /// @li Set force field
  ForceTestFlow3D<T,T,MyCase::descriptor_t> forceF(converter);
  const T latticeScaling(converter.getConversionFactorMass() / converter.getConversionFactorForce());
  AnalyticalScaled3D<T,T> scaledForceF(forceF, latticeScaling);  // conversion to lattice units
  lattice.template defineField<descriptors::FORCE>(geometry.getMaterialIndicator({1}), scaledForceF);

  /// @li Initialize lattice
  lattice.initialize();
}

/// Step 8.1: Update boundary values at times (and additional fields, if needed)
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
/// In order to start the simulation slowly, we start with zero boundary velocity and gradually increase it
/// @note Boundary values have to be set using lattice units
void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& converter = lattice.getUnitConverter();
  const T physStartT = myCase.getParameters().get<parameters::MAX_PHYS_T>() * (2.0 / 3.0);

  const std::size_t itStart = converter.getLatticeTime(physStartT);
  if (iT <= itStart) {
    /// @li Compute scaling factor
    PolynomialStartScale<T,T> startScaleF(itStart, T(1));
    const T iTvec[1] = {T(iT)};
    T frac[1] = {};
    startScaleF(frac, iTvec);

    /// @li Take analytical velocity solution, scale it to lattice units, set the boundary data
    VelocityTestFlow3D<T,T,MyCase::descriptor_t> velocityF(converter);
    AnalyticalScaled3D<T,T> uBoundaryStartF(velocityF, frac[0] / converter.getConversionFactorVelocity());
    lattice.defineU(geometry, 2, uBoundaryStartF);

    /// @li Communicate the new boundary velocity to GPU (if needed)
    lattice.template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }
}

/// Step 8.3: Compute simulation results at times
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();
  const T physMaxT = myCase.getParameters().get<parameters::MAX_PHYS_T>();
  const std::size_t iTvtk = converter.getLatticeTime(physMaxT/5);
  const std::size_t iTlog = converter.getLatticeTime(physMaxT/5);

  SuperVTMwriter3D<T> vtmWriter("testFlow3d");
  SuperLatticePhysVelocity3D velocityF(lattice, converter);
  SuperLatticePhysPressure3D pressureF(lattice, converter);
  vtmWriter.addFunctor(velocityF);
  vtmWriter.addFunctor(pressureF);

  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  // Writes the VTK files
  if (iT%iTvtk == 0 && iT > 0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }

  /// @li Print some (numerical and computational) statistics
  if (iT%iTlog == 0) {
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
    timer.print(iT);
  }
}

/// Step 8: Execute simulation
/// @param myCase The Case instance which keeps the simulation data
/// Run time loop
void simulate(MyCase& myCase)
{
  using T = MyCase::value_t;
  const T physMaxT = myCase.getParameters().get<parameters::MAX_PHYS_T>();
  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(physMaxT);
  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT)
  {
    /// @li Step 8.1: Update the Boundary Values and Fields at Times
    setTemporalValues(myCase, iT);

    /// @li Step 8.2: Collide and Stream Execution
    myCase.getLattice(NavierStokes{}).collideAndStream();

    /// @li Stripe off density offset due to Dirichlet boundary conditions
    myCase.getLattice(NavierStokes{}).stripeOffDensityOffset(
      myCase.getLattice(NavierStokes{}).getStatistics().getAverageRho() - T{1});

    /// @li Step 8.3: Computation and Output of the Results
    getResults(myCase, timer, iT);
  }

  /// @li Evaluate timer
  timer.stop();
  timer.printSummary();
}


/// Steps 2-8: Set up and run a simulation
int main(int argc, char* argv[])
{
  /// @par Contents
  /// @li Step 2: Initialization
  initialize(&argc, &argv);
  // Step 2.1: Set parameters
  MyCase::ParametersD myCaseParametersD;
  {
    using namespace parameters;
    myCaseParametersD.set<DOMAIN_EXTENT        >({1.0, 1.0, 1.0});
    myCaseParametersD.set<PHYS_CHAR_VELOCITY   >(            1.0);
    myCaseParametersD.set<PHYS_CHAR_VISCOSITY  >(            0.1);
    myCaseParametersD.set<PHYS_CHAR_DENSITY    >(            1.0);
    myCaseParametersD.set<MAX_PHYS_T           >(            6.0);
    myCaseParametersD.set<RESOLUTION           >(             11);
    myCaseParametersD.set<LATTICE_CHAR_VELOCITY>(           0.07);
  }

  /// @li Step 3: Create Mesh
  Mesh mesh = createMesh(myCaseParametersD);

  /// @li Step 4: Create Case
  MyCase myCase(myCaseParametersD, mesh);

  /// @li Step 5: Prepare Geometry
  prepareGeometry(myCase);

  /// @li Step 6: Prepare Lattice
  prepareLattice(myCase);

  /// @li Step 7: Definition of Initial and Boundary Values and Fields
  setInitialValues(myCase);

  /// @li Step 8: Simulate
  simulate(myCase);
}
