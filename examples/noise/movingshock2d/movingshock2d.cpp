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

/* movingshock2d.cpp:
 * The implementation of a backward facing step. It is furthermore
 * shown how to use checkpointing to save the state of the
 * simulation regularly.
 * The geometry of the step is based on the experiment described in
 * [Armaly, B.F., Durst, F., Pereira, J. C. F. and Schönung, B. Experimental
 * and theoretical investigation of backward-facing step flow. 1983.
 * J. Fluid Mech., vol. 127, pp. 473-496, DOI: 10.1017/S0022112083002839]
 */

#include "olb2D.h"
#include "olb2D.hh"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <getopt.h>

using namespace olb;
using namespace olb::descriptors;
using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = D2Q9<>;

using BulkDynamics = BGKdynamics<T,DESCRIPTOR>;

// Parameters for the simulation setup
//const T lengthDomain      = 1;                           // length of domain in meter
//const T heightDomain      = lengthDomain; // 0.0101;                        // height of domain in meter
//const T charL             = 1; //2 * heightInlet;               // characteristic length
//const size_t N            = 100 + 1;                            // resolution of the model
//const size_t maxN         = 1000; // maximum iteration number
//const size_t nout         = 10;  // frequency of images output
//const T maxPhysT          = 10; // / 2 / 343.;  // time for sound wave to hit boundary
//const T physTout          = maxPhysT/10; // time between output
//const T relaxationTime    = 0.55;
//const T amplitude         = 1e-3; // like Bocanegra // 20 * 1e-6; // hearing threshold
//const T Ma                = 0.5; // maximum stable Ma: .3
//const T targetMa          = 0.5;
//const T charV             = targetMa/Ma; // characteristic velocity
const T b                 = 3;
const T alpha             = log(2)/(b*b);
//const T rho_old    = 1.;//.293;
typedef enum {periodic, bounceBack, local, interpolated} BoundaryType;
BoundaryType boundarytype = interpolated;

// Stores geometry information in form of material numbers
void prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperGeometry<T,2>& superGeometry,
                     std::shared_ptr<IndicatorF2D<T>> domain
                     )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  // defining pressure point source
  if (boundarytype == periodic) {  // periodic
    superGeometry.rename(0, 1);
  } else {
    superGeometry.rename( 0, 2 ) ;
    superGeometry.rename( 2, 1, {1,1} );
    superGeometry.rename( 2, 3 );

  }
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
void prepareLattice(UnitConverter<T,DESCRIPTOR> const& converter,
                    SuperLattice<T,DESCRIPTOR>& sLattice,
                    SuperGeometry<T,2>& superGeometry,
                    T rho0,
                    T Ma,
                    T amplitude
                    )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  auto bulkIndicator = superGeometry.getMaterialIndicator({1});
  if ( boundarytype != periodic ) {
    bulkIndicator = superGeometry.getMaterialIndicator({1, 3});
  }
  sLattice.defineDynamics<BulkDynamics>(bulkIndicator);

  if (boundarytype == bounceBack) {
    setBounceBackBoundary(sLattice, superGeometry, 3);
  }
  else {
    if (boundarytype == local) {
      setLocalPressureBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 3);
      setLocalVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 3);
    }
    else {
      setInterpolatedPressureBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 3);
      setInterpolatedVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 3);
    }
  }

  // Initial conditions
//  bulkIndicator = superGeometry.getMaterialIndicator({0, 1});
  AnalyticalConst2D<T,T> *ux;
  if ( boundarytype != bounceBack ) {
    ux = new AnalyticalConst2D<T,T>( Ma );
  } else {
    ux = new AnalyticalConst2D<T,T>( 0. );
  }
  AnalyticalConst2D<T,T> uy( 0. );
  AnalyticalComposed2D<T,T> u( *ux, uy );
  AcousticPulse<2,T> pressureProfile( rho0, amplitude, alpha );
  //Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, pressureProfile, u );
  sLattice.iniEquilibrium( bulkIndicator, pressureProfile, u );

  sLattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a sinusoidal pressure
void setBoundaryValues( UnitConverter<T,DESCRIPTOR> const& converter,
                        SuperLattice<T,DESCRIPTOR>& sLattice, int iT,
                        SuperGeometry<T,2>& superGeometry, T rho0, T Ma )
{

  if ( boundarytype != periodic and boundarytype != bounceBack) {
    OstreamManager clout( std::cout,"setBoundaryValues" );

    AnalyticalConst2D<T,T> ux( Ma );
    AnalyticalConst2D<T,T> uy( 0. );
    AnalyticalComposed2D<T,T> u( ux,uy );
    AnalyticalConst2D<T,T> rho( rho0 );
    sLattice.defineRhoU( superGeometry, 3, rho, u );

    sLattice.iniEquilibrium( superGeometry.getMaterialIndicator({3}), rho, u );

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }

}

T L2Norm(SuperLattice<T,DESCRIPTOR>& sLattice,
        SuperGeometry<T,2>& superGeometry,
        int iT,
        UnitConverter<T,DESCRIPTOR> const& converter) {
  OstreamManager clout(std::cout, "error");

  T result[2] = {T(), T()};
  int tmp[] = {int()};

  SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure2D(sLattice, converter);
  AnalyticalConst<2,T,T> rho02D( 0. );

  auto indicatorF = superGeometry.getMaterialIndicator({1});
  SuperAbsoluteErrorL2Norm2D<T> absPressureErrorL2Norm(sLattice,
                                                       pressure2D,
                                                       rho02D,
                                                       indicatorF);

  absPressureErrorL2Norm(result, tmp);
  T l2_abs = result[0];

  return l2_abs;
}

// write data to termimal and file system
void getResults(SuperLattice<T,DESCRIPTOR>& sLattice,
                UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer,
                SuperPlaneIntegralFluxVelocity2D<T>& velocityFlux,
                SuperPlaneIntegralFluxPressure2D<T>& pressureFlux,
                T rho0,
                T amplitude,
                Gnuplot<T>& gplot_l2_abs, T Lp0,
                size_t iout,
                T tmax,
                size_t imax
                )
{
#ifdef PLATFORM_GPU_CUDA
  gpu::cuda::device::synchronize();
#endif
  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter2D<T> vtmWriter( "movingshock2d" );

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
  if ( iT%iout==0 ) {
//    sLattice.setProcessingContext(ProcessingContext::Evaluation);
//    velocityFlux.print();
//    pressureFlux.print();
    // write to terminal
    timer.update( iT );
    timer.printStep();
    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }

  if ( iT%int(iout/10)==0 ) {
    gplot_l2_abs.setData(T(iT), L2Norm(sLattice, superGeometry, iT,
                                       converter) / Lp0 );
  }

  // often get LU density plot along x
  if ( iT%iout == 0 ) {
    SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure( sLattice, converter );

    std::stringstream ss;
    ss << std::setw(4) << std::setfill('0') << iT;
    Gnuplot<T> gplot_hline( "pressure_hline_" + ss.str());
    Gnuplot<T> gplot_vline( "pressure_vline_" + ss.str());
    Gnuplot<T> gplot_diagonal( "pressure_diagonal_" + ss.str());
    gplot_hline.setLabel("distance [m]", "density [LU]");
    gplot_vline.setLabel("distance [m]", "density [LU]");
    gplot_diagonal.setLabel("distance [m]", "density [LU]");
    // Vectors for simulated solution
    T densities_hline[1] = {T()};
    T densities_vline[1] = {T()};
    T densities_diagonal[1] = {T()};
    T dist = converter.getPhysDeltaX();
    T ndatapoints = converter.getResolution(); // number of data points on line
//    CSV<T> csvWriterConcentration;
    //save concentration along the middle of the PRF
    for (int n = 0; n <= int(ndatapoints/2); n++) {
      T input_hline[2] =  {n*dist, 0};
      AnalyticalFfromSuperF2D<T> interpolation_hline( pressure, true, true );
      interpolation_hline(densities_hline, input_hline);
//      csvWriterConcentration.writeDataFile(input_hline[0], densities_hline[0], "simulation" , 16);
      gplot_hline.setData(input_hline[0], densities_hline[0]);

      T input_vline[2] =  {0, n*dist};
      AnalyticalFfromSuperF2D<T> interpolation_vline( pressure, true, true );
      interpolation_vline(densities_vline, input_vline);
      gplot_vline.setData(input_vline[1], densities_vline[0]);

      T input_diagonal[2] =  {n*dist, n*dist};
      AnalyticalFfromSuperF2D<T> interpolation_diagonal( pressure, true, true );
      interpolation_diagonal(densities_diagonal, input_diagonal);
//      csvWriterConcentration.writeDataFile(input_diagonal[0], densities_diagonal[0], "simulation" , 16);
      gplot_diagonal.setData(input_diagonal[0], densities_diagonal[0]);
    }
    T ymin(-amplitude/100);
    T ymax(+amplitude/100);
    gplot_hline.setYrange(ymin, ymax);
    gplot_vline.setYrange(ymin, ymax);
    gplot_diagonal.setYrange(ymin, ymax);
    // plot is generated
    gplot_hline.writePNG(-1, -1, "pressure_hline");
    gplot_vline.writePNG(-1, -1, "pressure_vline");
//    gplot_hline.writePDF("densities_hline");
    gplot_diagonal.writePNG(-1, -1, "pressure_diagonal");
//    gplot_diagonal.writePDF("densities_diagonal");
  }

  // get VTK and images
  if ( iT%(iout*2)==0) {
    SuperLatticePhysVelocity2D<T,DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure( sLattice, converter );
    // VTK
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    vtmWriter.write( iT );

    // pressure image
    BlockReduction2D2D<T> planeReductionP( pressure, 600, BlockDataSyncMode::ReduceOnly );
    heatmap::plotParam<T> jpeg_ParamP;
    jpeg_ParamP.maxValue = converter.getPhysPressure(+amplitude/50);
    jpeg_ParamP.minValue = converter.getPhysPressure(-amplitude/50);
    jpeg_ParamP.colour = "rainbow";
    jpeg_ParamP.fullScreenPlot = true;
    heatmap::write(planeReductionP, iT, jpeg_ParamP);
  }

  // Saves lattice data
  if ( iT%imax==0 && iT>0 ) {
    //clout << "Checkpointing the system at t=" << iT << std::endl;
    //sLattice.save( "movingshock2d.checkpoint" );
    // The data can be reloaded using
    //     sLattice.load("movingshock2d.checkpoint");
  }
}


int main( int argc, char* argv[], char *envp[] )
{
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  if ( singleton::mpi().isMainProcessor() ) {
    std::copy(argv, argv + argc, std::ostream_iterator<char *>(std::cout, " "));
    printf("\n");
  }
  std::copy(argv, argv + argc, std::ostream_iterator<char *>(std::cout, " "));
  printf("\n");
  singleton::directories().setOutputDir( "./tmp/" );  // set output directory
  OstreamManager clout( std::cout, "main" );

  // === get Command Line Arguments
  static const struct option long_options[] =
      {
          { "rho0",         required_argument, 0,  1  },
          { "res",          required_argument, 0,  2  },
          { "Ma",           required_argument, 0,  3  },
          { "targetMa",     required_argument, 0,  4  },
          { "amplitude",    required_argument, 0,  5  },
          { "length",       required_argument, 0, 'l' },
          { "tau",          required_argument, 0, 't' },
          { "imax",         required_argument, 0,  6  },
          { "tmax",         required_argument, 0,  7  },
          { "iout",         required_argument, 0,  8  },
          { "nout",         required_argument, 0,  9  },
          0
  };
  // initialize arguments
  int res       = 0;
  size_t imax   = 0;
  size_t iout   = 0;
  size_t nout   = 0;
  T tmax        = 0.;
  T rho0        = 0.;
  T Ma          = 0.;
  T targetMa       = 0.;
  T amplitude   = 0.;
  T tau         = 0.;
  T lengthDomain = 0.;

  while (1)
  {
    int index = -1;
    struct option * opt = 0;
    int result = getopt_long(argc, argv, "abc:d", long_options, &index);
    if (result == -1) break; /* end of list */
    switch (result)
    {
    case 1: /* same as index==0 */
      rho0 = T(atof(optarg));
      clout << "'rho0' was specified to " << optarg << std::endl;
      break;
    case 2: /* same as index==1 */
      res = atoi(optarg);
      clout << "'res' was specified to " << optarg << std::endl;
      break;
    case 3: /* same as index==1 */
      Ma = T(atof(optarg));
      clout << "'Ma' was specified to " << optarg << std::endl;
      break;
    case 4: /* same as index==1 */
      targetMa = T(atof(optarg));
      clout << "'targetMa' was specified to " << optarg << std::endl;
      break;
    case 5: /* same as index==1 */
      amplitude = T(atof(optarg));
      clout << "'amplitude' was specified to " << optarg << std::endl;
      break;
    case 6: /* same as index==1 */
      imax = atoi(optarg);
      clout << "'imax' was specified to " << optarg << std::endl;
      break;
    case 7: /* same as index==1 */
      tmax = T(atof(optarg));
      clout << "'tmax' was specified to " << optarg << std::endl;
      break;
    case 8: /* same as index==1 */
      iout = atoi(optarg);
      clout << "'iout' was specified to " << optarg << std::endl;
      break;
    case 9: /* same as index==1 */
      nout = atoi(optarg);
      clout << "'nout' was specified to " << optarg << std::endl;
      break;
    case 'l': /* same as index==2 */
      lengthDomain = T(atof(optarg));
      clout << "'l'/'length' was specified to " << optarg << std::endl;
      break;
    case 't': /* same as index==3 */
      tau = T(atof(optarg));
      clout << "'t'/'tau' was specified to " << optarg << std::endl;
      break;
    case 0: /* all parameter that do not appear in the optstring */
      opt = (struct option *)&(long_options[index]);
      clout << "'" << opt->name << "' was specified to ";
      if (opt->has_arg == required_argument) {
        clout << " to " << optarg;
      }
      clout << ".\n";
      break;
    default: /* unknown */
      break;
    }
  }
  /* print all other parameters */
  while (optind < argc)
  {
    clout << "other parameter: <" << argv[optind++] << ">\n";
  }
  /* default values */
  if ( res == 0 )           res             = 401;
  if ( rho0 == 0 )          rho0            = 1.;
  if ( imax == 0 )          imax            = 10000;
  if ( nout == 0 )          nout            = 20;
  if ( tmax == 0 )          tmax            = 2.;  // should be t=100 [LU]
  if ( Ma == 0 )            Ma              = 0.3;
  if ( targetMa == 0 )      targetMa        = 0.3;
  T charV                 = Ma/targetMa;
  if ( amplitude == 0 )     amplitude       = 0.001;
  if ( tau == 0 )           tau             = 0.518;
  if ( lengthDomain == 0 )  lengthDomain    = 1;
  const T heightDomain    = lengthDomain;

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> converter(
    (size_t) res,       // resolution
    (T)      tau,       // relaxation time
    (T)      lengthDomain,         // charPhysLength: reference length of simulation geometry
    (T)      charV,     // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)      1.48e-5,   // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)      rho0       // physDensity: physical density in __kg / m^3__
  );

  imax = std::min( converter.getLatticeTime( tmax ), imax );
  if ( iout == 0 )          iout            = imax / nout;

  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("movingshock2d");

  // === 2nd Step: Prepare Geometry ===

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 6;
#endif

  // setup domain
  T lengthDomainLU = 1.;
  T heightDomainLU = 1.;
  Vector<T,2> originDomain( - lengthDomain / 2, - heightDomain / 2 );  //
  Vector<T,2> extendDomain( lengthDomain, heightDomain );  // size of the domain
  std::shared_ptr<IndicatorF2D<T>> domain = std::make_shared<IndicatorCuboid2D<T>>( extendDomain, originDomain );

  CuboidGeometry2D<T> cuboidGeometry( *(domain), converter.getConversionFactorLength(), noOfCuboids );
  if ( boundarytype == periodic ) cuboidGeometry.setPeriodicity(true, true);

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,2> superGeometry( cuboidGeometry, loadBalancer );

  prepareGeometry( converter, superGeometry, domain );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T,DESCRIPTOR> sLattice( superGeometry );

  //prepare Lattice and set boundaryConditions
  prepareLattice( converter, sLattice, superGeometry, rho0, Ma, amplitude );

  // instantiate reusable functors
  SuperPlaneIntegralFluxVelocity2D<T> velocityFlux(
      sLattice,
      converter,
      superGeometry,
      {lengthDomain/T(2),  heightDomain/T(2)},
      {0.,  1.}
 );

  SuperPlaneIntegralFluxPressure2D<T> pressureFlux(
      sLattice,
      converter,
      superGeometry,
      {lengthDomain/T(2),  heightDomain/T(2) },
      {0.,  1.}
 );

  Gnuplot<T> gplot_l2_abs("l2_absolute");//, Gnuplot<T>::LOGLOG, Gnuplot<T>::OFF);
  gplot_l2_abs.setLabel("time []", "absolute L2 norm []");
  T Lp0 = L2Norm( sLattice, superGeometry, 0, converter );

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer(converter.getLatticeTime( tmax ),
                       superGeometry.getStatistics().getNvoxel()
                       );
  timer.start();

  for ( size_t iT = 0; iT < imax; ++iT ) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( converter, sLattice, iT, superGeometry, rho0, Ma );
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults(sLattice, converter, iT, superGeometry, timer,
               velocityFlux, pressureFlux,
               rho0, amplitude,
               gplot_l2_abs, Lp0,
               iout, tmax, imax);
    // test density failure
    if (sLattice.getStatistics().getAverageRho() > 2.) {
      clout << "breaking simulation loop because density is too high: "
            << sLattice.getStatistics().getAverageRho();
      break;
    }
    if (sLattice.getStatistics().getAverageEnergy() > 2.) {
      clout << "breaking simulation loop because energy is too high: "
            << sLattice.getStatistics().getAverageEnergy();
      break;
    }
  }

  gplot_l2_abs.setYrange(1e-3, 1);
  gplot_l2_abs.setLogScale(2);
  gplot_l2_abs.writePNG(-1, -1, "gplot_l2_abs");
  timer.stop();
  timer.printSummary();
}
