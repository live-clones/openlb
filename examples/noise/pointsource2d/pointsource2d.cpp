/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006, 2007, 2012 Jonas Latt, Mathias J. Krause,
 *  Louis Kronberg, Christian Vorwerk, Bastian Schäffauer
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

/* pointsource2d.cpp:
 * The implementation of a backward facing step. It is furthermore
 * shown how to use checkpointing to save the state of the
 * simulation regularly.
 * The geometry of the step is based on the experiment described in
 * [Armaly, B.F., Durst, F., Pereira, J. C. F. and Schönung, B. Experimental
 * and theoretical investigation of backward-facing step flow. 1983.
 * J. Fluid Mech., vol. 127, pp. 473-496, DOI: 10.1017/S0022112083002839]
 */

#include "olb.h"
//#include "olb2D.hh"
#include <iostream>
#include <sstream>
#include <iomanip>


using namespace olb;
using namespace olb::descriptors;
using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = D2Q9<>;

using BulkDynamics = BGKdynamics<T,DESCRIPTOR>;

// Parameters for the simulation setup
const T lengthDomain      = 1;                           // length of domain in meter
const T heightDomain      = lengthDomain; // 0.0101;                        // height of domain in meter
const T charL             = lengthDomain; //2 * heightInlet;               // characteristic length
const size_t N            = 960 + 1;                            // resolution of the model
const size_t maxN         = 1e10; // maximum iteration number
const T maxPhysT          = lengthDomain;  // / 2 / 343.;  // time for sound wave to hit boundary
const T frequency         = 10; // Hz
const T relaxationTime    = 0.518;
const T radiusPointSource = lengthDomain / N;
const T amplitude         = 2e-3; // like Bocanegra // 20 * 1e-6; // hearing threshold

// Stores geometry information in form of material numbers
void prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperGeometry<T,2>& superGeometry,
                     std::shared_ptr<IndicatorF2D<T>> domain,
                     std::shared_ptr<IndicatorF2D<T>> pointsource
                     )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  // defining pressure point source
  superGeometry.rename( 0, 2 ) ;//, domain-pointsource );
  superGeometry.rename( 2, 3, pointsource );
  superGeometry.rename( 2, 1 );
  superGeometry.communicate();

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperLattice<T,DESCRIPTOR>& sLattice,
                     SuperGeometry<T,2>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  auto bulkIndicator = superGeometry.getMaterialIndicator({1,3});
  sLattice.defineDynamics<BulkDynamics>(bulkIndicator);
  // Material=2 -->bulk dynamics (pressure source)
//  setLocalPressureBoundary(sLattice, omega, superGeometry, 3);
  //if boundary conditions are chosen to be interpolated
//  setInterpolatedVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 3);

  // Initial conditions
  AnalyticalConst2D<T,T> ux( 0. );
  AnalyticalConst2D<T,T> uy( 0. );
  AnalyticalConst2D<T,T> rho( 1. );
  AnalyticalComposed2D<T,T> u( ux,uy );

  //Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, rho, u );
  AnalyticalConst2D<T,T> rho2( 2. );
  auto bulkIndicator2 = superGeometry.getMaterialIndicator({3});
  sLattice.defineRhoU( bulkIndicator2, rho2, u );
//  sLattice.iniEquilibrium( bulkIndicator, rho, u );

//  AnalyticalConst2D<T,T> rho2( 2. );
//  sLattice.defineRhoU( superGeometry, 1, rho2, u );
//  sLattice.defineRhoU( superGeometry, 2, rho2, u );
  sLattice.iniEquilibrium( bulkIndicator, rho, u );
  sLattice.iniEquilibrium( bulkIndicator2, rho, u );

  sLattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a sinusoidal pressure
void setBoundaryValues( UnitConverter<T,DESCRIPTOR> const& converter,
                        SuperLattice<T,DESCRIPTOR>& sLattice, int iT,
                        SuperGeometry<T,2>& superGeometry )
{
  OstreamManager clout( std::cout,"setBoundaryValues" );
//  int iTmaxStart = converter.getLatticeTime( maxPhysT );

  T rho = converter.getLatticeDensityFromPhysPressure(
          amplitude * std::sin(converter.getPhysTime(iT) * frequency * 2 * 3.14)
      );
  AnalyticalConst2D<T,T> rho2d(rho);
  AnalyticalConst2D<T,T> ux( 0. );
  AnalyticalConst2D<T,T> uy( 0. );
  AnalyticalComposed2D<T,T> u2d( ux,uy );

  auto bulkIndicator = superGeometry.getMaterialIndicator({3});
  sLattice.defineRhoU(bulkIndicator, rho2d, u2d);
//  sLattice.iniEquilibrium( bulkIndicator, rho2d, u2d );

  sLattice.setProcessingContext<
    Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
    ProcessingContext::Simulation);

//  sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
}

// write data to termimal and file system
void getResults( SuperLattice<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer,
                 SuperPlaneIntegralFluxVelocity2D<T>& velocityFlux,
                 SuperPlaneIntegralFluxPressure2D<T>& pressureFlux )
{
  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter2D<T> vtmWriter( "pointsource2d" );

  if ( iT==0 ) {
    // Writes geometry, cuboid no. and rank no. to file system
//    SuperLatticeGeometry2D<T,DESCRIPTOR> geometry( sLattice, superGeometry );
//    SuperLatticeCuboid2D<T,DESCRIPTOR> cuboid( sLattice );
//    SuperLatticeRank2D<T,DESCRIPTOR> rank( sLattice );
//    vtmWriter.write( geometry );
//    vtmWriter.write( cuboid );
//    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }


  // Writes to console 10 times
  if ( iT%converter.getLatticeTime( maxPhysT / 10 )==0 ) {
//    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    velocityFlux.print();
    pressureFlux.print();
    // write to terminal
    timer.update( iT );
    timer.printStep();
    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }

  // often get LU density plot along x
  if ( iT%100 == 0 ) {
    SuperLatticeDensity2D<T, DESCRIPTOR> density(sLattice);

    std::stringstream ss;
    ss << std::setw(4) << std::setfill('0') << iT;
    Gnuplot<T> gplot_centerline( "pressure_centerline_" + ss.str());
    Gnuplot<T> gplot_diagonal( "pressure_diagonal_" + ss.str());
    gplot_centerline.setLabel("distance [LU]", "density [LU]");
    gplot_diagonal.setLabel("distance [LU]", "density [LU]");
    // Vectors for simulated solution
    T densities_centerline[1] = {T()};
    T densities_diagonal[1] = {T()};
    T dist = converter.getPhysDeltaX();
    T ndatapoints = converter.getResolution(); // number of data points on line
    CSV<T> csvWriterConcentration;
    //save concentration along the middle of the PRF
    for (int n = int(ndatapoints/2); n <= ndatapoints; n++) {
      T input_centerline[2] =  {n*dist, lengthDomain/2};
      AnalyticalFfromSuperF2D<T> interpolation_centerline( density, true, true );
      interpolation_centerline(densities_centerline, input_centerline);
      csvWriterConcentration.writeDataFile(input_centerline[0], densities_centerline[0], "simulation" , 16);
      gplot_centerline.setData(input_centerline[0], densities_centerline[0]);

      T input_diagonal[2] =  {n*dist, n*dist};
      AnalyticalFfromSuperF2D<T> interpolation_diagonal( density, true, true );
      interpolation_diagonal(densities_diagonal, input_diagonal);
      csvWriterConcentration.writeDataFile(input_diagonal[0], densities_diagonal[0], "simulation" , 16);
      gplot_diagonal.setData(input_diagonal[0], densities_diagonal[0]);
    }
    T ymin(converter.getLatticeDensityFromPhysPressure(-amplitude));
    T ymax(converter.getLatticeDensityFromPhysPressure(+amplitude));
    gplot_centerline.setYrange(ymin, ymax);
    gplot_diagonal.setYrange(ymin, ymax);
    // plot is generated
    gplot_centerline.writePNG(-1, -1, "densities_centerline");
//    gplot_centerline.writePDF("densities_centerline");
    gplot_diagonal.writePNG(-1, -1, "densities_diagonal");
//    gplot_diagonal.writePDF("densities_diagonal");
  }


  // get VTK and images
  if ( iT%converter.getLatticeTime( maxPhysT / 5 )==0 ) {
    SuperLatticePhysVelocity2D<T,DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure( sLattice, converter );
    // VTK
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    vtmWriter.write( iT );

    // pressure image
    BlockReduction2D2D<T> planeReductionP( pressure, 1200, BlockDataSyncMode::ReduceOnly );
    heatmap::plotParam<T> jpeg_ParamP;
    jpeg_ParamP.maxValue = amplitude/10;
    jpeg_ParamP.minValue = 0.;
    jpeg_ParamP.colour = "rainbow";
    jpeg_ParamP.fullScreenPlot = true;
    heatmap::write(planeReductionP, iT, jpeg_ParamP);
  }

  // Saves lattice data
  if ( iT%converter.getLatticeTime( maxPhysT )==0 && iT>0 ) {
    //clout << "Checkpointing the system at t=" << iT << std::endl;
    //sLattice.save( "pointsource2d.checkpoint" );
    // The data can be reloaded using
    //     sLattice.load("pointsource2d.checkpoint");
  }
}


int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );  // set output directory
  OstreamManager clout( std::cout, "main" );

//  CommandLineParser parser(argc, argv);
//  parser.setArgument<size_t>("N", "number of grid points", N);
//  parser.setArgument<T>("tau", "relaxation time", relaxationTime);
//  parser.setArgument<T>("f", "frequency", frequency);
//  parser.setArgument<T>("A", "amplitude", amplitude);
//  parser.setArgument<T>("tmax", "physical simulation time", maxPhysT);
//  parser.setArgument<T>("nmax", "maximum iterations", maxN);

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> converter(
    (size_t) N, //parser.getArgument<size_t>("N"),      // resolution
    (T)      relaxationTime, //parser.getArgument<T>("tau"),         // relaxation time
    (T)      charL,             // charPhysLength: reference length of simulation geometry
    (T)      1,                 // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)      1.48e-5,           // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)      1.293              // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("pointsource2d");

  // === 2nd Step: Prepare Geometry ===

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 6;
#endif

  // setup domain
  Vector<T,2> originDomain( 0, 0 );  //
  Vector<T,2> extendDomain( lengthDomain, heightDomain );  // size of the domain
  std::shared_ptr<IndicatorF2D<T>> domain = std::make_shared<IndicatorCuboid2D<T>>( extendDomain, originDomain );

  Vector<T,2> originSource( lengthDomain/2., heightDomain/2.);
  std::shared_ptr<IndicatorF2D<T>> pointsource = std::make_shared<IndicatorCircle2D<T>>( originSource, radiusPointSource );

  CuboidGeometry2D<T> cuboidGeometry( *(domain-pointsource), converter.getConversionFactorLength(), noOfCuboids );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,2> superGeometry( cuboidGeometry, loadBalancer );

  prepareGeometry( converter, superGeometry, domain, pointsource );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T,DESCRIPTOR> sLattice( superGeometry );

  //prepare Lattice and set boundaryConditions
  prepareLattice( converter, sLattice, superGeometry );

  // instantiate reusable functors
  SuperPlaneIntegralFluxVelocity2D<T> velocityFlux( sLattice,
      converter,
      superGeometry,
  {lengthDomain/T(2),  heightDomain / T(2)},
  {0.,  1.} );

  SuperPlaneIntegralFluxPressure2D<T> pressureFlux( sLattice,
      converter,
      superGeometry,
  {lengthDomain/T(2),  heightDomain / T(2) },
  {0.,  1.} );

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( std::size_t iT = 0; iT < std::min( converter.getLatticeTime( maxPhysT ), maxN ); ++iT ) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( converter, sLattice, iT, superGeometry );
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer, velocityFlux, pressureFlux );
  }

  timer.stop();
  timer.printSummary();
}
