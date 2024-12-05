/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2013 Mathias J. Krause, Thomas Henn, Tim Dornieden
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

/* du93_3d.cpp:
 * TODO: change description
 */


#include "olb3D.h"
#include "olb3D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = D3Q19<>;

// T* uAverage = NULL;

#define BOUZIDI

// Parameters for the simulation setup
// const int N = 10;        // resolution of the model
// const T Re = 20.;       // Reynolds number
// const T maxPhysT = 50.; // max. simulation time in s, SI unit

// Stores data from stl file in geometry in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter, //std::shared_ptr<IndicatorF3D<T>> domain,
                      STLreader<T>& stlReader, SuperGeometry<T,3>& superGeometry, bool withTrailingEdge )
{
  OstreamManager clout( std::cout,"prepareGeometry" );

  superGeometry.rename( 0,2 );

  Vector<T,3> origin = superGeometry.getStatistics().getMinPhysR( 2 );
  Vector<T,3> extend = superGeometry.getStatistics().getMaxPhysR( 2 ) - superGeometry.getStatistics().getMinPhysR( 2 );
  origin[0] += converter.getConversionFactorLength();
  extend[0] -= 2*converter.getConversionFactorLength();
  IndicatorCuboid3D<T> domain( extend, origin );

  // all nodes except x-boundaries to fluid
  superGeometry.rename( 2, 1, domain );
  // du93 profile
  superGeometry.rename( 1,5,stlReader );
  if ( withTrailingEdge ) {
    // trailing edge area
    origin[0] = 0.7;
    extend[0] = 0.3;
    IndicatorCuboid3D<T> trailingEdge( extend, origin );
    superGeometry.rename( 5, 6, trailingEdge );
  }

  // Set material number for inflow
  origin[0] = superGeometry.getStatistics().getMinPhysR( 2 )[0]-converter.getConversionFactorLength();
  extend[0] = 2*converter.getConversionFactorLength();
  IndicatorCuboid3D<T> inflow( extend,origin );
  superGeometry.rename( 2,3,1,inflow );

  // Set material number for outflow
  origin[0] = superGeometry.getStatistics().getMaxPhysR( 2 )[0]-converter.getConversionFactorLength();
  extend[0] = 2*converter.getConversionFactorLength();
  IndicatorCuboid3D<T> outflow( extend,origin );
  superGeometry.rename( 2,4,1,outflow );

  // Removes all not needed boundary voxels outside the surface
  // superGeometry.clean(5);
  superGeometry.checkForErrors();
  clout << "Prepare Geometry ..." << std::endl
        << "Materials:" << std::endl
        << "1 = Fluid" << std::endl
        << "2 = Check (should be empty!)" << std::endl
        << "3 = Inflow" << std::endl
        << "4 = Outflow" << std::endl
        << "5 = Airfoil Bounce Back or Bouzidi" << std::endl
        << "6 = Airfoil porous trailing edge" << std::endl;
  superGeometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( SuperLattice<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     STLreader<T>& stlReader,
                     SuperGeometry<T,3>& superGeometry,
                     bool withTrailingEdge,
                     T Kin )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  auto bulkIndicator = superGeometry.getMaterialIndicator({1});
  sLattice.defineDynamics<BGKdynamics>(bulkIndicator);

  // Material=2 -->bounce back
  setBounceBackBoundary(sLattice, superGeometry, 2);

  // Material=3 -->inlet
  setLocalVelocityBoundary(sLattice, omega, superGeometry, 3);
  // Material=4 -->outlet
  setInterpolatedConvectionBoundary(sLattice, omega, superGeometry, 4);

  // Material=5 -->bouzidi / bounce back
  #ifdef BOUZIDI
  setBouzidiBoundary<T,DESCRIPTOR>(sLattice, superGeometry, 5, stlReader);
  #else
  setBounceBackBoundary(sLattice, superGeometry, 5);
  #endif

  if ( withTrailingEdge ) {
    // Material=6 --> porous media
    sLattice.defineDynamics<PorousBGKdynamics>(superGeometry, 6);  // velocity will always be multiplied by d
    T tau = converter.getLatticeRelaxationTime();
    T nu = (tau-0.5)/3.;
    T h = converter.getPhysDeltaX();
    T d = 1. - (h*h*nu*tau/Kin);
    clout << "Lattice Porosity: " << d << "(1 = permeable , 0 = not permeable)" << std::endl;
    clout << "Kmin: " << h*h*nu*tau << std::endl;
    if (Kin < h*h*nu*tau) {
      clout << "WARNING: Chosen K is too small!" << std::endl;
      exit(1);
    }
    AnalyticalConst3D<T,T> porosity(d);
    sLattice.defineField<POROSITY>(superGeometry.getMaterialIndicator({0,1,2,3,4,5,6}), porosity);
  }

  // Initial conditions
  AnalyticalConst3D<T,T> rhoF( 1 );
  Vector<T,3> velocityV;
  AnalyticalConst3D<T,T> uF(velocityV);

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, rhoF, uF );
  sLattice.iniEquilibrium( bulkIndicator, rhoF, uF );

  sLattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( SuperLattice<T, DESCRIPTOR>& sLattice,
                        UnitConverter<T,DESCRIPTOR> const& converter,
                        size_t iT,
                        SuperGeometry<T,3>& superGeometry,
                        size_t iTmaxStart,
                        T ux_max_LU,
                        T maxPhysU )
{
  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up#
  int iTupdate = int(iTmaxStart/50/10)*10;  // --> about 50 updates, rounded to 10 iterations

  if ( iT%iTupdate == 0 && iT <= iTmaxStart ) {
    // Smooth start curve, sinus
    // SinusStartScale<T,int> StartScale(iTmaxStart, T( ux_max_LU ));

    // Smooth start curve, polynomial
    PolynomialStartScale<T, size_t> StartScale( iTmaxStart, T( ux_max_LU ) );

    // Creates and sets the Poiseuille inflow profile using functors
    size_t iTvec[1] = { iT };
    T ux_i[1] = {};
    StartScale( ux_i,iTvec );
    AnalyticalConst3D<T,T> ux( ux_i[0] );
    AnalyticalConst3D<T,T> uy( 0. );
    AnalyticalConst3D<T,T> uz( 0. );
    AnalyticalComposed3D<T,T> u( ux, uy, uz );
    sLattice.defineU( superGeometry, 3, u );

    clout << "startup step=" << iT << "/" << iTmaxStart << "; t=" << converter.getPhysTime(iT)
          << "; ux_LU=" << ux_i[0] << "; ux_max_LU=" << ux_max_LU
          << "; ux_PU=" << converter.getPhysVelocity( ux_i[0] ) << "; ux_max_LU=" << maxPhysU << std::endl;

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }
}

void vtkOutput( SuperVTMwriter3D<T> vtmWriter, int iT,
                SuperLatticePhysVelocity3D<T, DESCRIPTOR>& velocity,
                SuperLatticePhysPressure3D<T, DESCRIPTOR>& pressure,
                SuperDiscretizationF3D<T>& discretization ) {
  vtmWriter.write( iT );

  {
    SuperEuklidNorm3D<T> normVel( velocity );
    BlockReduction3D2D<T> planeReduction( normVel, Vector<T,3>({0, 0, 1}) );
    // write output as JPEG
    heatmap::write(planeReduction, iT);
  }

  {
    BlockReduction3D2D<T> planeReduction( discretization, Vector<T,3>({0, 0, 1}) );
    heatmap::plotParam<T> jpeg_scale;
    jpeg_scale.colour = "blackbody";
    jpeg_scale.name = "quality";
    heatmap::write( planeReduction, iT, jpeg_scale );
  }
}

// Computes the pressure drop between the voxels before and after the airfoil
void getResults( SuperLattice<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer,
                 STLreader<T>& stlReader,
                 size_t iTmax, size_t iTmaxStart, size_t iTvtk )
{

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter3D<T> vtmWriter( "du93_3d" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  SuperLatticeYplus3D<T, DESCRIPTOR> yPlus( sLattice, converter, superGeometry, stlReader, 5 );
  SuperLatticeRefinementMetricKnudsen3D<T, DESCRIPTOR> quality( sLattice, converter );
  SuperRoundingF3D<T, T> roundedQuality ( quality, RoundingMode::NearestInteger );
  SuperDiscretizationF3D<T> discretization ( roundedQuality, 0., 2. );

  vtmWriter.addFunctor( quality );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );
  vtmWriter.addFunctor( yPlus );

  const size_t statIter = iTmaxStart;
  const size_t iTcheck  = 500;
  bool lastIteration = false;

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }

  if ( iT%iTcheck == 0 ) {
    if ( sLattice.getStatistics().getAverageRho() > 2. ) {
      lastIteration = true;
    }
  }

  // Writes output on the console
  if ( iT%statIter == 0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Drag, lift, pressure drop
    AnalyticalFfromSuperF3D<T> intpolatePressure( pressure, true );
    SuperLatticePhysDrag3D<T,DESCRIPTOR> drag( sLattice, superGeometry.getMaterialIndicator({5,6}), converter );

    olb::Vector<T, 3> point1V = superGeometry.getStatistics().getCenterPhysR( 5 );
    olb::Vector<T, 3> point2V = superGeometry.getStatistics().getCenterPhysR( 5 );
    T point1[3] = {};
    T point2[3] = {};
    for ( int i = 0; i<3; i++ ) {
      point1[i] = point1V[i];
      point2[i] = point2V[i];
    }
    point1[0] = std::min(superGeometry.getStatistics().getMinPhysR( 5 )[0], superGeometry.getStatistics().getMinPhysR( 6 )[0]) - converter.getConversionFactorLength();
    point2[0] = std::max(superGeometry.getStatistics().getMaxPhysR( 5 )[0], superGeometry.getStatistics().getMaxPhysR( 6 )[0]) + converter.getConversionFactorLength();

    T p1, p2;
    intpolatePressure( &p1,point1 );
    intpolatePressure( &p2,point2 );

    clout << "pressure1=" << p1;
    clout << "; pressure2=" << p2;

    T pressureDrop = p1-p2;
    clout << "; pressureDrop=" << pressureDrop;

    T dragA[3];
    int input1[0];
    drag( dragA, input1 );
    clout << "; drag=" << dragA[0] << "; lift=" << dragA[1] << std::endl;

    int input[4] = {};
    SuperMax3D<T> yPlusMaxF( yPlus, superGeometry, 1 );
    T yPlusMax[1];
    yPlusMaxF( yPlusMax,input );
    clout << "yPlusMax=" << yPlusMax[0] << std::endl;

    if ( p1 != p1 ) lastIteration = true;

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }

  // Writes the vtk files
  if ( iT%iTvtk == 0 || lastIteration ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    vtkOutput( vtmWriter, iT, velocity, pressure, discretization );
    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }

  if ( lastIteration ) {
      timer.stop();
      timer.printSummary();
      exit(1);
  }
}

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  //clout.setMultiOutput(true);
  CLIreader args(argc, argv);
  std::string outdir        = args.getValueOrFallback<std::string>( "--outdir", "./tmp" );
  const T lengthDomain      = args.getValueOrFallback( "--lx",            6);   // in m
  const T heightDomain      = args.getValueOrFallback( "--ly",            3);   // in m
  const T depthDomain       = args.getValueOrFallback( "--lz",            1.2); // in m
  const size_t res          = args.getValueOrFallback( "--res",           50);  // voxel/m (dx_LU/m)
  T maxPhysT                = args.getValueOrFallback( "--tmax",          10);  // in s
  size_t iTmax              = args.getValueOrFallback( "--imax",          0);
  size_t nout               = args.getValueOrFallback( "--nout",          5);   // minimum number of vtk outputs
  size_t iout               = args.getValueOrFallback( "--iout",          0);   // iterations for vtk outputs
  T tout                    = args.getValueOrFallback( "--tout",          0);   // timestep for vtk outputs
  T maxPhysU                = args.getValueOrFallback( "--umax",          1);   // in m/s
  T Re                      = args.getValueOrFallback( "--Re",            0);
  T tau                     = args.getValueOrFallback( "--tau",           0);   // previously tau=0.53 fixed
  const T Kin               = args.getValueOrFallback( "--permeability",  1e-8);
  const T tMaxInit          = args.getValueOrFallback( "--tmaxinit",      2);
  T charL                   = args.getValueOrFallback( "--charL",         1);
  const bool withTrailingEdge = !args.contains("--no-porous");

  const T cs_LU             = 1 / std::sqrt(3.0);
  T charV                   = 2 * maxPhysU;
  T viscosity;
  if ( tau == 0 ) {
    if ( Re == 0 ) {
      viscosity             = 1.48e-5; // 1e-2;
      Re                    = charV * charL / viscosity;
    } else viscosity        = charV * charL / Re;
    tau                     = viscosity / (cs_LU*cs_LU) + 0.5;
  }
  else {
    if ( Re == 0 ) Re       = 200;
    viscosity               = charV * charL / Re;
  }

  std::stringstream outdir_mod;
  outdir_mod << outdir;
  if ( !withTrailingEdge ) outdir_mod << "_noPorous";
  else outdir_mod << "_" << Kin << "porous";
  outdir_mod << "_u" << charV << "_Re" << Re << "_" << lengthDomain << "x" << heightDomain << "x" << depthDomain << "_res" << res;

  singleton::directories().setOutputDir( outdir_mod.str()+"/" );  // set output directory
  if ( withTrailingEdge ) clout << "Calculating with trailing edge" << std::endl;
  else clout << "Calculating without trailing edge" << std::endl;

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    (size_t)  res,            // resolution: number of voxels per charPhysL
    (T)       tau,            // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)       charL,          // charPhysLength: reference length of simulation geometry
    (T)       charV,          // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)       viscosity,      // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)       1.204           // physDensity: physical density in __kg / m^3__; air at 101.325 kPa and 20 °C
  );
  T umax_LU                 = converter.getLatticeVelocity( maxPhysU );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("du93_3d");

  // === 2nd Step: Prepare Geometry ===

  // Instantiation of the STLreader class
  // file name, voxel size in meter, stl unit in meter, outer voxel no., inner voxel no.
  STLreader<T> stlReader( "extruded_shape.stl", converter.getConversionFactorLength() );
  IndicatorLayer3D<T> extendedDomain( stlReader, converter.getConversionFactorLength() );

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 7;
#endif
  // setup domain
  Vector<T,3> originDomain( -lengthDomain/3, -heightDomain/2, 0.1 );  //
  Vector<T,3> extendDomain( lengthDomain, heightDomain, depthDomain );  // size of the domain
  std::shared_ptr<IndicatorF3D<T>> domain = std::make_shared<IndicatorCuboid3D<T>>( extendDomain, originDomain );

  CuboidGeometry3D<T> cuboidGeometry(*(domain),
                                     converter.getConversionFactorLength(),
                                     noOfCuboids
                                     );
  cuboidGeometry.setPeriodicity(false, true, true);

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,3> superGeometry( cuboidGeometry, loadBalancer );

  prepareGeometry( converter, stlReader, superGeometry, withTrailingEdge );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  //prepareLattice and set boundaryCondition
  prepareLattice( sLattice, converter, stlReader, superGeometry, withTrailingEdge, Kin );

  // === 3a-rd Step: calculate iterations from input ===
  // iTmax depends on maximum physical time. If iTmax is provided in command line, it is an upper bound
  if ( iTmax == 0 ) iTmax = converter.getLatticeTime( maxPhysT );
  else iTmax = std::min( iTmax, converter.getLatticeTime( maxPhysT ) );
  size_t iTmaxStart = converter.getLatticeTime( tMaxInit );
  // nout is the minimum number of vtk outputs --> take max between nout and nout derived from iout or tout
  size_t nout_from_iout = 0, nout_from_tout = 0;
  if ( iout != 0 ) { nout_from_iout = size_t( iTmax / iout ); nout = std::max( nout, nout_from_iout ); }
  if ( tout != 0 ) { nout_from_tout = size_t( iTmax / tout ); nout = std::max( nout, nout_from_tout ); }
  size_t iTvtk = size_t( iTmax / nout );

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( iTmax, superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for (size_t iT = 0; iT < iTmax; ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( sLattice, converter, iT, superGeometry, iTmaxStart, umax_LU, maxPhysU );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer, stlReader, iTmax, iTmaxStart, iTvtk );
  }

  timer.stop();
  timer.printSummary();
}
