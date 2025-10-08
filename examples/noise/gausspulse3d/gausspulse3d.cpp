/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Philipp Spelten
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

/* gausspulse3d.cpp:
 * Benchmark from C. Bogey and C. Bailly, “Three-dimensional non-reﬂective boundary conditions for acoustic simulations: far ﬁeld formulation and validation test cases,” vol. 88, 2002.
 * Method from H. Xu and P. Sagaut, “Analysis of the absorbing layers for the weakly-compressible lattice Boltzmann methods,” Journal of Computational Physics, vol. 245, pp. 14–42, Jul. 2013, doi: 10.1016/j.jcp.2013.02.051.
 * Use the option --BOUNDARY_CONDITION to vary how the far field is modeled:
 *  0=eternal (3x domain size to capture the vanishing pulse as ideal reference)
 *  1=periodic,
 *  2=local equilibrium BC,
 *  3=damping (with periodic boundaries)
 */

#include "olb.h"
#include "../noiseauxiliary.h"

using namespace olb;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q19<>>
>;

enum class BoundaryCondition: int {
  eternal   = 0,
  periodic  = 1,
  local     = 2,
  damping   = 3
};

namespace olb::parameters {
  // DOMAIN
  struct LATT_CHAR_VELOCITY : public descriptors::FIELD_BASE<1> { };    // characteristic lattice velocity ("Mach number")
  struct CORE_EXTENT        : public descriptors::FIELD_BASE<0,1> { };  // fluid domain size [m]
  // PULSE
  struct AMPLITUDE          : public descriptors::FIELD_BASE<1> { };    // pressure amplitude of Gaussian pulse
  struct ALPHA              : public descriptors::FIELD_BASE<1> { };    // distribution factor of Gaussian pulse
  // BOUNDARY CONDITIONS
  struct BOUNDARY_CONDITION : public descriptors::TYPED_FIELD_BASE<BoundaryCondition,1> { };  // counter, see beginning of file
  struct ETERNALSCALE       : public descriptors::FIELD_BASE<1> { };               // amout to extend domain by for 'eternal' case
  struct DAMPING_DEPTH_LU   : public descriptors::TYPED_FIELD_BASE<size_t,1> { };  // number of points for sponge layer
  struct DAMPING_STRENGTH   : public descriptors::FIELD_BASE<1> { };               // maximum damping strength of sponge layer
  // TIMING AND OUTPUTS
  struct T_TABLE            : public descriptors::FIELD_BASE<1> { };  // tabular outputs every T_TABLE seconds (for line plot)
  struct T_GRAPHICAL_OUTPUT : public descriptors::FIELD_BASE<1> { };  // graphical outputs every T_GRAPHICAL seconds
  struct T_LOG              : public descriptors::FIELD_BASE<1> { };  // timer outputs every T_LOG seconds
  struct DO_GRAPHICAL_OUTPUT: public descriptors::TYPED_FIELD_BASE<bool,1> { };    // set to do images
}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  OstreamManager clout(std::cout, "createMesh");
  clout << std::endl << "createMesh ..." << std::endl;
  using T = MyCase::value_t;

  Vector            extentFluid     = parameters.get<parameters::CORE_EXTENT>();
  BoundaryCondition boundarytype    = parameters.get<parameters::BOUNDARY_CONDITION>();
  const T           physDeltaX      = parameters.get<parameters::PHYS_DELTA_X>();
  T                 dampingDepthPU  = physDeltaX * parameters.get<parameters::DAMPING_DEPTH_LU>();
  T                 eternalscale    = parameters.get<parameters::ETERNALSCALE>();

  Vector extentDomain = extentFluid;
  switch ( boundarytype ) {
    case BoundaryCondition::eternal:  extentDomain *= eternalscale;       break; // extend the domain
    case BoundaryCondition::periodic:                                     break; // reference domain size
    case BoundaryCondition::local:    extentDomain += 2 * physDeltaX;     break; // add one layer in each direction
    case BoundaryCondition::damping:  extentDomain += 2 * dampingDepthPU; break; // add boundary layer
  }
  clout << "Fluid Domain = " << extentFluid[0] << "^3; Simulation Domain = " << extentDomain[0] << "^3" << std::endl;
  parameters.set<parameters::DOMAIN_EXTENT>( {extentDomain[0],extentDomain[1],extentDomain[2]} );

  const Vector originDomain = -0.5*extentDomain;
  IndicatorCuboid3D<T> cuboid( extentDomain, originDomain );

  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  if ( boundarytype == BoundaryCondition::local ) {
    mesh.getCuboidDecomposition().setPeriodicity( {false, false, false} );
  } else {
    mesh.getCuboidDecomposition().setPeriodicity( {true, true, true} );
  }
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry( MyCase& myCase ) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << std::endl << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto&             geometry          = myCase.getGeometry();
  auto&             parameters        = myCase.getParameters();
  BoundaryCondition boundaryCondition = parameters.get<parameters::BOUNDARY_CONDITION>();
  const T           physDeltaX2       = parameters.get<parameters::PHYS_DELTA_X>() / 2.;
  const Vector      extentFluid       = parameters.get<parameters::CORE_EXTENT>();
  const Vector      originFluid       = -0.5*extentFluid;
  IndicatorCuboid3D<T> domainFluid( extentFluid, originFluid );

  // all nodes to temporary type
  geometry.rename(0, 2);
  Vector origin = geometry.getStatistics().getMinPhysR(2);
  Vector extent = geometry.getStatistics().getMaxPhysR(2) - geometry.getStatistics().getMinPhysR(2);

  if ( boundaryCondition == BoundaryCondition::local ) {
    // set fluid material, keeping the outside layer on 2
    geometry.rename(2, 1, {1, 1, 1});

    // Set material number for inflow
    origin        -= physDeltaX2;
    extent        += 2 * physDeltaX2;
    extent[0]     = 2 * physDeltaX2;
    IndicatorCuboid3D<T> inflow(extent, origin);
    geometry.rename(2, 4, 1, inflow);

    // Set material number for outflow
    origin[0] = geometry.getStatistics().getMaxPhysR(2)[0] - physDeltaX2;
    IndicatorCuboid3D<T> outflow(extent, origin);
    geometry.rename(2, 5, 1, outflow);

    origin    = geometry.getStatistics().getMinPhysR(2) - physDeltaX2;
    extent    = geometry.getStatistics().getMaxPhysR(2) - geometry.getStatistics().getMinPhysR(2) + 2 * physDeltaX2;
    extent[1] = 2 * physDeltaX2;
    IndicatorCuboid3D<T> bottom(extent, origin);
    geometry.rename(2, 6, 1, bottom);
    origin[1] = geometry.getStatistics().getMaxPhysR(2)[1] - physDeltaX2;
    IndicatorCuboid3D<T> top(extent, origin);
    geometry.rename(2, 6, 1, top);

    origin    = geometry.getStatistics().getMinPhysR(2) - physDeltaX2;
    extent    = geometry.getStatistics().getMaxPhysR(2) - geometry.getStatistics().getMinPhysR(2) + 2 * physDeltaX2;
    extent[2] = 2 * physDeltaX2;
    IndicatorCuboid3D<T> front(extent, origin);
    geometry.rename(2, 6, 1, front);
    origin[2] = geometry.getStatistics().getMaxPhysR(2)[2] - physDeltaX2;
    IndicatorCuboid3D<T> back(extent, origin);
    geometry.rename(2, 6, 1, back);

  } else if ( boundaryCondition == BoundaryCondition::damping ) {

    extent -= 6 * physDeltaX2;
    origin += 3 * physDeltaX2;
    IndicatorCuboid3D<T> spongeOutside(extent, origin);
    geometry.rename(2, 3, spongeOutside);
    geometry.rename(3, 1, domainFluid);
    geometry.rename(2, 1);

  } else {

    geometry.rename(2, 1);

  }

  geometry.checkForErrors();
  geometry.print();
  clout
      << "Materials:" << std::endl
      << "0 = no Material (default, should be immediately renamed to Check)" << std::endl
      << "1 = fluid (except damping and eternal, then it's the far field)" << std::endl
      << "2 = check (should be empty, but corners may have issues with the normal; check is therefore treated as inflow"
      << std::endl
      << "3 = sponge" << std::endl
      << "4 = inflow (local)" << std::endl
      << "5 = outflow (local)" << std::endl
      << "6 = around domain (top/bottom/front/back)" << std::endl;

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( MyCase& myCase ) {
  OstreamManager clout(std::cout, "prepareLattice");
  clout << std::endl << "prepareLattice ..." << std::endl;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});

  const T physExtentX       = parameters.get<parameters::CORE_EXTENT>()[0];
  const T physCharLength    = physExtentX;
  const T physDeltaX        = parameters.get<parameters::PHYS_DELTA_X>();
  const T latticeCharVelocity = parameters.get<parameters::LATT_CHAR_VELOCITY>();
  const T physCharVelocity  = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T physCharViscosity = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T physCharDensity   = parameters.get<parameters::PHYS_CHAR_DENSITY>();
  const T physDeltaT        = latticeCharVelocity / physCharVelocity * physDeltaX;

  // Set up a unit converter with the characteristic physical units
  lattice.setUnitConverter(
    physDeltaX,         // physDeltaX: spacing between two lattice cells in [m]
    physDeltaT,         // physDeltaT: time step in [s]
    physCharLength,     // charPhysLength: reference length of simulation geometry in [m]
    physCharVelocity,   // physCharVelocity: highest expected velocity during simulation in [m/s]
    physCharViscosity,  // physCharViscosity: physical kinematic viscosity in [m^2/s]
    physCharDensity     // physCharDensity: physical density [kg/m^3]
  );
  auto& converter = lattice.getUnitConverter();
  converter.print();

  auto bulkIndicator = geometry.getMaterialIndicator({1,3,4,5,6});
  BoundaryCondition boundaryCondition = parameters.get<parameters::BOUNDARY_CONDITION>();
  switch ( boundaryCondition ) {
    case BoundaryCondition::eternal:
    case BoundaryCondition::periodic:
      lattice.defineDynamics<BGKdynamics>( geometry, 1 );
      break;
    case BoundaryCondition::local:
      lattice.defineDynamics<BGKdynamics>( geometry, 1 );
      boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 4);
      boundary::set<boundary::LocalPressure>(lattice, geometry, 5);
      boundary::set<boundary::LocalVelocity>(lattice, geometry, 6);
      break;
    case BoundaryCondition::damping:
      lattice.defineDynamics<
        SpongeLayerDynamics<
          T,DESCRIPTOR,momenta::BulkTuple,equilibria::SecondOrder>
        >( geometry, 3 );
      lattice.defineDynamics<BGKdynamics>( geometry, 1 );
      break;
  }
  lattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  // Define far field observable fields
  AnalyticalConst3D<T,T>  uxFar(converter.getCharLatticeVelocity());
  AnalyticalConst3D<T,T>  uyFar(0.);
  AnalyticalConst3D<T,T>  uzFar(0.);
  AnalyticalConst3D<T,T>  rhoFar(1.);
  lattice.defineField<descriptors::UX>( bulkIndicator, uxFar );
  lattice.defineField<descriptors::UY>( bulkIndicator, uyFar );
  lattice.defineField<descriptors::UZ>( bulkIndicator, uzFar );
  lattice.defineField<descriptors::DENSITY>( bulkIndicator, rhoFar );

  // define and output damping layer parameter
  const T       dampingDepthPU  = parameters.get<parameters::DAMPING_DEPTH_LU>() * converter.getPhysDeltaX();
  const Vector  domainLengths   = parameters.get<parameters::DOMAIN_EXTENT>();
  const T       dampingStrength = parameters.get<parameters::DAMPING_STRENGTH>();
  DampingTerm<3,T,DESCRIPTOR> sigma(dampingDepthPU, domainLengths, dampingStrength);
  lattice.defineField<descriptors::DAMPING>(bulkIndicator, sigma);
  // === alternatively, define sigma as a global constant
  // AnalyticalConst<3,T,T> sigma( dampingStrength );

  clout << "prepareLattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  OstreamManager clout(std::cout, "setInitialValues");
  clout << std::endl << "setInitialValues ..." << std::endl;
  using T = MyCase::value_t;
  auto& lattice     = myCase.getLattice(NavierStokes{});
  auto& geometry    = myCase.getGeometry();
  auto& parameters  = myCase.getParameters();
  auto& converter   = lattice.getUnitConverter();

  T amplitude = parameters.get<parameters::AMPLITUDE>();
  T alpha     = parameters.get<parameters::ALPHA>();

  // Initial velocity and density
  AnalyticalConst3D<T,T>  u(converter.getCharLatticeVelocity(), 0., 0.);
  AcousticPulse<3,T>      densityProfile(1., amplitude, alpha);
  size_t res  = converter.getResolution();
  T physDx    = converter.getPhysDeltaX();
  linePlot<3,T>( densityProfile, res, physDx, "pulse_diag", "density [LU]", diagonal2d);
  linePlot<3,T>( densityProfile, res, physDx, "pulse_hline", "density [LU]", horizontal);

  /// Initialize populations to equilibrium state
  auto domain = geometry.getMaterialIndicator({1,3,4,5,6});
  lattice.iniEquilibrium( domain, densityProfile, u);
  lattice.defineRhoU(     domain, densityProfile, u);
  lattice.initialize();

  clout << "setInitialValues ... OK" << std::endl;
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{ }

// write data to file system
void setPlotData( MyCase& myCase,
                  util::Timer<MyCase::value_t>& timer,
                  std::size_t iT,
                  Gnuplot<MyCase::value_t>& gplot_l2_abs,
                  MyCase::value_t Lp0 ) {
  using T = MyCase::value_t_of<NavierStokes>;
  auto& lattice     = myCase.getLattice(NavierStokes{});
  auto& geometry    = myCase.getGeometry();
  auto& parameters  = myCase.getParameters();
  auto& converter   = lattice.getUnitConverter();

  const size_t iTout = converter.getLatticeTime(parameters.get<parameters::T_TABLE>());
  const Vector extent = parameters.get<parameters::CORE_EXTENT>();
  const Vector origin = -0.5*extent;  // === calculate output intervals
  SuperIndicatorFfromIndicatorF3D<T> domainFluid( new IndicatorCuboid3D<T>( extent, origin ), geometry );

  if (iT % iTout == 0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    T Lpi = L2Norm<3,T, MyCase::descriptor_t_of<NavierStokes>>(lattice, converter, domainFluid);
    gplot_l2_abs.setData(T(iT), Lpi / Lp0);
  }
}

void getGraphicalResults( MyCase& myCase,
                          util::Timer<MyCase::value_t>& timer,
                          std::size_t iT ) {

  /// Write vtk plots every 0.3 seconds (of phys. simulation time)
  using T = MyCase::value_t;
  auto& parameters      = myCase.getParameters();
  auto& lattice         = myCase.getLattice(NavierStokes{});
  auto& converter       = lattice.getUnitConverter();
  const T amplitude     = parameters.get<parameters::AMPLITUDE>();
  const std::string name= "gausspulse3d";
  const size_t iTgraph  = converter.getLatticeTime(parameters.get<parameters::T_GRAPHICAL_OUTPUT>());

  SuperVTMwriter3D<T> vtmWriter(name);
  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  if ( iT % iTgraph == 0 ) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    lattice.scheduleBackgroundOutputVTK([&, name, iT](auto task) {
      SuperVTMwriter3D<T>        vtmWriter(name);
      SuperLatticePhysVelocity3D velocityF(lattice, converter);
      SuperLatticePhysPressure3D pressureF(lattice, converter);
      vtmWriter.addFunctor(velocityF);
      vtmWriter.addFunctor(pressureF);
      task(vtmWriter, iT);
    });

    // output pressure image
    SuperLatticePhysPressure3D pressure(lattice, converter);
    BlockReduction3D2D<T>      pressureReduction(pressure, Vector<T,3>({0, 0, 1}));
    heatmap::plotParam<T>      jpeg_ParamP;
    jpeg_ParamP.maxValue       = converter.getPhysPressure(+amplitude / 200);
    jpeg_ParamP.minValue       = converter.getPhysPressure(-amplitude / 200);
    jpeg_ParamP.colour         = "rainbow";
    jpeg_ParamP.fullScreenPlot = true;
    heatmap::write(pressureReduction, iT, jpeg_ParamP);

    std::stringstream ss;
    ss << std::setw(4) << std::setfill('0') << iT;
    T                          dist        = converter.getPhysDeltaX();
    T                          ndatapoints = converter.getResolution(); // number of data points on line
    AnalyticalFfromSuperF3D<T> pressure_interpolation(pressure, true, true);
    T                          pmin(converter.getPhysPressure(-amplitude / 50));
    T                          pmax(converter.getPhysPressure(+amplitude / 50));
    linePlot<3,T>(pressure_interpolation, ndatapoints, dist,
                  "pressure_hline_" + ss.str(), "pressure [PU]",
                  horizontal, false, true, pmin, pmax);
    linePlot<3,T>(pressure_interpolation, ndatapoints, dist,
                  "pressure_vline_" + ss.str(), "pressure [PU]",
                  vertical, false, true, pmin, pmax);
    linePlot<3,T>(pressure_interpolation, ndatapoints, dist,
                  "pressure_diagonal_" + ss.str(), "pressure [PU]",
                  diagonal2d, false, true, pmin, pmax);
  }
}

void simulate(MyCase& myCase) {
  using T = MyCase::value_t;
  auto&         parameters  = myCase.getParameters();
  auto&         geometry    = myCase.getGeometry();
  auto&         lattice     = myCase.getLattice(NavierStokes{});
  auto&         converter   = lattice.getUnitConverter();
  const size_t  iTmax       = converter.getLatticeTime(parameters.get<parameters::MAX_PHYS_T>());
  const size_t  iTlog       = converter.getLatticeTime(parameters.get<parameters::T_LOG>());

  // === Initialize pressure L2 norm plot
  const Vector extent = parameters.get<parameters::CORE_EXTENT>();
  const Vector origin = -0.5*extent;  // === calculate output intervals
  SuperIndicatorFfromIndicatorF3D<T> domainFluid( new IndicatorCuboid3D<T>( extent, origin ), geometry );
  Gnuplot<T> gplot_l2_abs("l2_absolute");
  gplot_l2_abs.setLabel("time []", "absolute L2 norm []");
  T Lp0 = L2Norm<3,T>( lattice, converter, domainFluid );
  const bool doGraphicalOutput = parameters.get<parameters::DO_GRAPHICAL_OUTPUT>();

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    myCase.getLattice(NavierStokes{}).collideAndStream();

    /// === Step 8.3: Computation and Output of the Results ===
    setPlotData(myCase, timer, iT, gplot_l2_abs, Lp0);
    if ( doGraphicalOutput ) getGraphicalResults(myCase, timer, iT);

    /// === Step 8.4: Print some (numerical and computational) statistics ===
    if ( iT%iTlog == 0 ) {
      lattice.getStatistics().print(iT, converter.getPhysTime(iT));
      timer.print(iT);
    }
  }

  gplot_l2_abs.setYrange(1e-5, 1);
  gplot_l2_abs.setLogScale(2);
  gplot_l2_abs.writePNG(-1, -1, "gplot_l2_abs");

  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");
  clout << "Starting gausspulse3d ..." << std::endl;

  /// === Step 2a: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    // DOMAIN
    myCaseParameters.set<LATT_CHAR_VELOCITY >(        0.1 );
    myCaseParameters.set<CORE_EXTENT        >( {1.,1.,1.} );
    // PULSE
    myCaseParameters.set<AMPLITUDE          >(       1e-3 );
    myCaseParameters.set<ALPHA              >( log(2.) / (1. / 20. * 1. / 20.) );
    // BOUNDARY CONDITIONS
    myCaseParameters.set<BOUNDARY_CONDITION >( BoundaryCondition::damping );
    myCaseParameters.set<DAMPING_STRENGTH   >(         1. );
    myCaseParameters.set<DAMPING_DEPTH_LU   >(         20 );
    myCaseParameters.set<ETERNALSCALE       >(         3. );
    // TIMING AND OUTPUTS
    myCaseParameters.set<DO_GRAPHICAL_OUTPUT>(      false );
    myCaseParameters.set<T_GRAPHICAL_OUTPUT >(        .05 );
    myCaseParameters.set<T_LOG              >(        .05 );
    myCaseParameters.set<T_TABLE            >(       .025 );
    // STANDARD PARAMETERS
    myCaseParameters.set<MAX_PHYS_T         >(        .25 );
    myCaseParameters.set<PHYS_DELTA_T       >(       1e-3 );
    myCaseParameters.set<PHYS_CHAR_VELOCITY >(        1.0 );
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(      0.001 );
    myCaseParameters.set<PHYS_CHAR_DENSITY  >(        1.0 );
    myCaseParameters.set<RESOLUTION         >(        101 );
    myCaseParameters.set<PHYS_DELTA_X       >(
      myCaseParameters.get<CORE_EXTENT>()[0]
      / myCaseParameters.get<RESOLUTION>());
  }
  myCaseParameters.fromCLI(argc, argv);
  // === Step 2b: set output directory depending on input values ===
  // === Step 2c: set boundarytype depending on input values ===
  using T = MyCase::value_t;
  BoundaryCondition boundaryCondition = myCaseParameters.get<parameters::BOUNDARY_CONDITION>();
  T                 eternalscale      = myCaseParameters.get<parameters::ETERNALSCALE>();
  T                 amplitude         = myCaseParameters.get<parameters::AMPLITUDE>();
  Vector            coreExtent        = myCaseParameters.get<parameters::CORE_EXTENT>();
  size_t            res               = myCaseParameters.get<parameters::RESOLUTION>();
  size_t            dampingDepthLU    = myCaseParameters.get<parameters::DAMPING_DEPTH_LU>();
  T                 dampingStrength   = myCaseParameters.get<parameters::DAMPING_STRENGTH>();

  std::stringstream outdir_mod;
  outdir_mod << "./tmp";
  switch ( boundaryCondition ) {
    case BoundaryCondition::eternal:
      std::cout << "Boundary condition is solved by just extending the domain to " << eternalscale << " times" << std::endl;
      outdir_mod << "_eternal"; break;
    case BoundaryCondition::periodic:
      clout << "Boundary condition type specified to periodic." << std::endl;
      outdir_mod << "_periodic"; break;
    case BoundaryCondition::local:
      clout << "Boundary condition type specified to local." << std::endl;
      outdir_mod << "_local"; break;
    case BoundaryCondition::damping:
      clout << "Boundary condition type specified to damping." << std::endl;
      outdir_mod << "_damping"; break;
  }
  if (amplitude != 1e-3)    outdir_mod << "_a" << amplitude;
  if (coreExtent[0] != 1.)  outdir_mod << "_l" << coreExtent[0];
  if (res != 101)           outdir_mod << "_res" << res;
  if (boundaryCondition == BoundaryCondition::eternal)
                            outdir_mod << "_scale" << eternalscale;
  if (boundaryCondition == BoundaryCondition::damping)
                            outdir_mod << "_bd" << dampingDepthLU << "x" << dampingStrength;

  singleton::directories().setOutputDir(outdir_mod.str() + "/");
  clout << "Output directory set to " << outdir_mod.str() << std::endl;

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

  clout << "gausspulse3d ... done" << std::endl;
}