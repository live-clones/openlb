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

/* movingshock3d.cpp:
 * The implementation of a backward facing step. It is furthermore
 * shown how to use checkpointing to save the state of the
 * simulation regularly.
 * The geometry of the step is based on the experiment described in
 * [Armaly, B.F., Durst, F., Pereira, J. C. F. and Schönung, B. Experimental
 * and theoretical investigation of backward-facing step flow. 1983.
 * J. Fluid Mech., vol. 127, pp. 473-496, DOI: 10.1017/S0022112083002839]
 */

#include "olb3D.h"
#include "olb3D.hh"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <getopt.h>

using namespace olb;
using namespace olb::descriptors;
using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = D3Q27<>;

using BulkDynamics = BGKdynamics<T,DESCRIPTOR>;
typedef enum {periodic, local, damping} BoundaryType;

// Stores geometry information in form of material numbers
void prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperGeometry<T,3>& superGeometry,
                     std::shared_ptr<IndicatorF3D<T>> domain,
                     BoundaryType boundarytype,
                     int res,
                     bool debug,
                     int boundary_depth,
                     Vector<T,3> domain_lengths
                     )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl
        << "Materials:" << std::endl
        << "1 = Fluid" << std::endl
        << "2 = Check (should be empty)" << std::endl
        << "3 = Outer Boundaries (far field)" << std::endl
        << "4 = Cylinder" << std::endl;

  superGeometry.rename( 0, 2 );   // all nodes to temporary type

  // cylinder
  const T radiusCylinder = std::max(domain_lengths[1]/20., converter.getPhysDeltaX()*4);
  Vector<T,3> center1( -domain_lengths[0]/4.+converter.getPhysDeltaX()*3./2., converter.getPhysDeltaX()*3./2., -domain_lengths[2]/2.);
  Vector<T,3> center2( -domain_lengths[0]/4.+converter.getPhysDeltaX()*3./2., converter.getPhysDeltaX()*3./2., +domain_lengths[2]/2.);
  std::shared_ptr<IndicatorF3D<T>> cylinder = std::make_shared<IndicatorCylinder3D<T>>( center1, center2, radiusCylinder );

  // fluid domain
  if ( boundarytype == periodic ) {
    superGeometry.rename( 2, 4, cylinder );
    superGeometry.rename( 2, 1 );   // all remaining are fluid
  } else if ( boundarytype == damping ) {
    superGeometry.rename( 2, 4, cylinder );
    T bd_pu = converter.getPhysLength(boundary_depth);
    bd_pu = converter.getPhysLength(boundary_depth);
    Vector<T,3> extend( domain_lengths[0]-2*bd_pu, domain_lengths[1]-2*bd_pu, domain_lengths[2]/2);//-2*bd_pu );
    Vector<T,3> origin( -domain_lengths[0]/2+bd_pu, -domain_lengths[1]/2+bd_pu, -domain_lengths[2]/2);//+bd_pu );
    std::shared_ptr<IndicatorF3D<T>> fluid_domain = std::make_shared<IndicatorCuboid3D<T>>( extend, origin );
    superGeometry.rename( 2, 1, fluid_domain );
    superGeometry.rename( 2, 3 );                 // remaining to far-field boundaries
  } else if ( boundarytype == local) {
    superGeometry.rename( 2, 1, {1,1,0} );        // all fluid except outer layer
    superGeometry.rename( 1, 4, cylinder );
    superGeometry.rename( 2, 4, cylinder );
    superGeometry.rename( 2, 3 );                 // remaining to outer boundary
  }

  superGeometry.communicate();

  // Removes all not needed boundary voxels outside the surface
  // superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  // superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(UnitConverter<T,DESCRIPTOR> const& converter,
                    SuperLattice<T,DESCRIPTOR>& sLattice,
                    SuperGeometry<T,3>& superGeometry,
                    T rho0,
                    T Ma,
                    BoundaryType boundarytype,
                    int boundary_depth,
                    T damping_strength,
                    Vector<T,3> domain_lengths
                    )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 --> bulk dynamics
  auto bulkIndicator = superGeometry.getMaterialIndicator({1});
  sLattice.defineDynamics<BulkDynamics>(bulkIndicator);

  // Material=3 --> dynamics depend..
  clout << "boundarytype=" << boundarytype << std::endl;
  if ( boundarytype == damping ) {
    sLattice.defineDynamics<BulkDynamics>(superGeometry.getMaterialIndicator({3}));
  }
  if ( boundarytype == local ) {
    setLocalPressureBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 3);
    setLocalVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 3);
  }
  else if ( boundarytype == damping ) {
    clout << "Setting damping boundary" << std::endl;
    bulkIndicator = superGeometry.getMaterialIndicator({1, 3});
    setDampingBoundary<T,DESCRIPTOR>(sLattice, superGeometry, 3);
  }

  // Material=4 --> cylinder simple bounce back
  setBounceBackBoundary<T,DESCRIPTOR>(sLattice, superGeometry, 4);

  // Initial conditions
  AnalyticalConst3D<T,T> *ux, *uy, *uz;
  if ( boundarytype == local ) Ma = 0;
  ux = new AnalyticalConst3D<T,T>( Ma );
  uy = new AnalyticalConst3D<T,T>( 0. );
  uz = new AnalyticalConst3D<T,T>( 0. );
  AnalyticalComposed3D<T,T> u( *ux, *uy, *uz );

  AnalyticalConst3D<T,T> rho(rho0);
  //Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( superGeometry.getMaterialIndicator({1, 3}), rho, u );
  sLattice.iniEquilibrium( superGeometry.getMaterialIndicator({1, 3}), rho, u );

  // Define fields
  sLattice.setParameter<descriptors::OMEGA>( omega );
  AnalyticalConst3D<T,T> rhoField( rho0 );
  sLattice.defineField<descriptors::DENSITY>( superGeometry, 3, rhoField );
  AnalyticalConst3D<T,T> uxField( Ma );
  sLattice.defineField<descriptors::UX>( superGeometry, 3, uxField );
  AnalyticalConst3D<T,T> uyField( 0. );
  sLattice.defineField<descriptors::UY>( superGeometry, 3, uyField );
  AnalyticalConst3D<T,T> uzField( 0. );
  sLattice.defineField<descriptors::UZ>( superGeometry, 3, uzField );  

  DampingTerm<3,T,DESCRIPTOR> sigma( converter, boundary_depth, domain_lengths, damping_strength );
  sLattice.defineField<descriptors::DAMPING>( superGeometry, 3, sigma );
  
  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a sinusoidal pressure
void setBoundaryValues( UnitConverter<T,DESCRIPTOR> const& converter,
                        SuperLattice<T,DESCRIPTOR>& sLattice, int iT,
                        SuperGeometry<T,3>& superGeometry,
                        T rho0, T Ma,
                        BoundaryType boundarytype
                        )
{
  if ( boundarytype == local ) {
    OstreamManager clout( std::cout,"setBoundaryValues" );
    int iTmaxStart(1000);
    SinusStartScale<T,int> startScale( iTmaxStart, Ma );
    int iTvec[1]= {iT};
    T velocity[1]= {};
    startScale( velocity,iTvec );
    AnalyticalConst3D<T,T> ux( velocity[0] );
    AnalyticalConst3D<T,T> uy( 0. );
    AnalyticalConst3D<T,T> uz( 0. );
    AnalyticalComposed3D<T,T> u( ux, uy, uz );
    AnalyticalConst3D<T,T> rho( rho0 );

    sLattice.defineRhoU( superGeometry, 3, rho, u );
    sLattice.iniEquilibrium( superGeometry.getMaterialIndicator({3}), rho, u );

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  } else if ( boundarytype == damping ) {
    // OstreamManager clout( std::cout,"setBoundaryValues" );
    // int iTmaxStart(1000);
    // SinusStartScale<T,int> startScale( iTmaxStart, Ma );
    // int iTvec[1]= {iT};
    // T velocity[1]= {};
    // startScale( velocity,iTvec );
    // AnalyticalConst3D<T,T> ux( velocity[0] );
    // AnalyticalConst3D<T,T> uy( 0. );
    // AnalyticalConst3D<T,T> uz( 0. );
    // AnalyticalComposed3D<T,T> u( ux, uy, uz );
    // AnalyticalConst3D<T,T> rho( rho0 );

    // sLattice.defineRhoU( superGeometry, 3, rho, u );

    // sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }
}

T L2Norm(SuperLattice<T,DESCRIPTOR>& sLattice,
        SuperGeometry<T,3>& superGeometry,
        int iT,
        UnitConverter<T,DESCRIPTOR> const& converter) {
  OstreamManager clout(std::cout, "L2Norm");

  T result[3] = {T(), T(), T()};
  int tmp[] = {int()};

  SuperLatticePhysPressure3D<T,DESCRIPTOR> pressure3D(sLattice, converter);
  AnalyticalConst<3,T,T> rho03D( 0. );

  auto indicatorF = superGeometry.getMaterialIndicator({1});
  SuperAbsoluteErrorL2Norm3D <T> absPressureErrorL2Norm(pressure3D, rho03D, indicatorF);

  absPressureErrorL2Norm(result, tmp);
  T l2_abs = result[0];

  return l2_abs;
}

// write data to termimal and file system
void getResults(SuperLattice<T,DESCRIPTOR>& sLattice,
                UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer,
                SuperPlaneIntegralFluxVelocity3D<T>& velocityFlux,
                SuperPlaneIntegralFluxPressure3D<T>& pressureFlux,
                T rho0,
                Gnuplot<T>& gplot_l2_abs, T Lp0,
                size_t iout, T tmax, size_t imax
                )
{
#ifdef PLATFORM_GPU_CUDA
  gpu::cuda::device::synchronize();
#endif
  OstreamManager clout( std::cout,"getResults" );
  // clout << "getResults startet at iT=" << iT << ", iout=" << iout << std::endl;
  SuperVTMwriter3D<T> vtmWriter( "cylinder3dnoise" );

  if ( iT==0 ) {
    // Writes geometry, cuboid no. and rank no. to file system
//    SuperLatticeGeometry3D<T,DESCRIPTOR> geometry( sLattice, superGeometry );
//    SuperLatticeCuboid3D<T,DESCRIPTOR> cuboid( sLattice );
//    SuperLatticeRank3D<T,DESCRIPTOR> rank( sLattice );
//    vtmWriter.write( geometry );
//    vtmWriter.write( cuboid );
//    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Writes to console 10 times
  if ( iT%iout==0 ) {
    // sLattice.setProcessingContext(ProcessingContext::Evaluation);
    // velocityFlux.print();
    // pressureFlux.print();
    // write to terminal
    timer.update( iT );
    timer.printStep();
    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }

  if ( iT%std::max(int(iout/10), 1)==0 ) {
    gplot_l2_abs.setData(T(iT), L2Norm(sLattice, superGeometry, iT,
                                       converter) / Lp0 );
  }

  // often get LU density plot along x
  if ( iT%iout == 0 ) {
    SuperLatticePhysPressure3D<T,DESCRIPTOR> pressure( sLattice, converter );
    SuperLatticePhysVelocity3D<T,DESCRIPTOR> velocity( sLattice, converter );

    std::stringstream ss;
    ss << std::setw(4) << std::setfill('0') << iT;
    Gnuplot<T> gplot_hline( "pressure_hline_" + ss.str());
    Gnuplot<T> gplot_vline( "pressure_vline_" + ss.str());
    Gnuplot<T> gplot_diagonal( "pressure_diagonal_" + ss.str());
    Gnuplot<T> gplot_udiagonal( "velocity_diagonal_" + ss.str());
    gplot_hline.setLabel("distance [m]", "density [LU]");
    gplot_vline.setLabel("distance [m]", "density [LU]");
    gplot_diagonal.setLabel("distance [m]", "density [LU]");
    gplot_udiagonal.setLabel("distance [m]", "velocity [LU]");
    // Vectors for simulated solution
    T densities_hline[1] = {T()};
    T densities_vline[1] = {T()};
    T densities_diagonal[1] = {T()};
    T velocities_diagonal[1] = {T()};
    T dist = converter.getPhysDeltaX();
    T ndatapoints = converter.getResolution(); // number of data points on line
    // CSV<T> csvWriterConcentration;
    // save concentration along the middle of the PRF
    for (int n = -int(ndatapoints/2); n <= int(ndatapoints/2); n++) {
      T input_hline[3] =  {n*dist, 0, 0};
      AnalyticalFfromSuperF3D<T> interpolation_hline( pressure, true, true );
      interpolation_hline(densities_hline, input_hline);
    // csvWriterConcentration.writeDataFile(input_hline[0], densities_hline[0], "simulation" , 16);
      gplot_hline.setData(input_hline[0], densities_hline[0]);

      T input_vline[3] =  {0, n*dist, 0};
      AnalyticalFfromSuperF3D<T> interpolation_vline( pressure, true, true );
      interpolation_vline(densities_vline, input_vline);
      gplot_vline.setData(input_vline[1], densities_vline[0]);

      T input_diagonal[3] =  {n*dist, n*dist, n*dist};
      AnalyticalFfromSuperF3D<T> interpolation_diagonal( pressure, true, true );
      interpolation_diagonal(densities_diagonal, input_diagonal);
      // csvWriterConcentration.writeDataFile(input_diagonal[0], densities_diagonal[0], "simulation" , 16);
      gplot_diagonal.setData(input_diagonal[0], densities_diagonal[0]);
      AnalyticalFfromSuperF3D<T> interpolation_udiagonal( velocity, true, true );
      interpolation_udiagonal(velocities_diagonal, input_diagonal);
      gplot_udiagonal.setData(input_diagonal[0], velocities_diagonal[0]);
    }

    // T ymin(converter.getPhysPressure(-amplitude/200));
    // T ymax(converter.getPhysPressure(+amplitude/200));
    // gplot_hline.setYrange(ymin, ymax);
    // gplot_vline.setYrange(ymin, ymax);
    // gplot_diagonal.setYrange(ymin, ymax);

    // plot is generated
    gplot_hline.writePNG(-1, -1, "pressure_hline");
    gplot_vline.writePNG(-1, -1, "pressure_vline");
    // gplot_hline.writePDF("densities_hline");
    gplot_diagonal.writePNG(-1, -1, "pressure_diagonal");
    gplot_udiagonal.writePNG(-1, -1, "velocity_diagonal");
    // gplot_diagonal.writePDF("densities_diagonal");

    // get VTK and images
    if ( iT%(iout*2)==0) {
      clout << "vtmWriter startet" << std::endl;
      SuperLatticePhysVelocity3D<T,DESCRIPTOR> velocity( sLattice, converter );
      SuperLatticePhysPressure3D<T,DESCRIPTOR> pressure( sLattice, converter );
      // VTK
      vtmWriter.addFunctor( velocity );
      vtmWriter.addFunctor( pressure );
      vtmWriter.write( iT );

      // pressure image
      BlockReduction3D2D<T> pressureReduction( pressure, Vector<T,3>({0, 0, 1}) );
      heatmap::plotParam<T> jpeg_ParamP;
      // jpeg_ParamP.maxValue = converter.getPhysPressure(+amplitude/200);
      // jpeg_ParamP.minValue = converter.getPhysPressure(-amplitude/200);
      jpeg_ParamP.colour = "rainbow";
      jpeg_ParamP.fullScreenPlot = true;
      heatmap::write(pressureReduction, iT, jpeg_ParamP);

      // velocity image
      SuperEuklidNorm3D<T> normVel( velocity );
      // BlockReduction3D2D<T> planeReduction( normVel, origin, u, v, 600, BlockDataSyncMode::ReduceOnly );
      BlockReduction3D2D<T> planeReduction( normVel, Vector<T,3>({0, 0, 1}) );
      heatmap::plotParam<T> plotParam;
      jpeg_ParamP.maxValue = converter.getCharPhysVelocity()*1.1;
      jpeg_ParamP.minValue = converter.getCharPhysVelocity()*0.9;
      jpeg_ParamP.colour = "rainbow";
      jpeg_ParamP.fullScreenPlot = true;
      heatmap::write(planeReduction, iT, plotParam);
    }
  }

  // Saves lattice data
  if ( iT%imax==0 && iT>0 ) {
    //clout << "Checkpointing the system at t=" << iT << std::endl;
    //sLattice.save( "movingshock3d.checkpoint" );
    // The data can be reloaded using
    //     sLattice.load("movingshock3d.checkpoint");
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
  CLIreader args(argc, argv);
  const std::string outdir = args.getValueOrFallback<std::string>( "--outdir", "./tmp/" );
  singleton::directories().setOutputDir( outdir );  // set output directory
  OstreamManager clout( std::cout, "main" );

  // === get Command Line Arguments
  const int res = args.getValueOrFallback( "--res", 81);
  const T rho0 = args.getValueOrFallback( "--rho0", 1.);
  const T Ma = args.getValueOrFallback( "--Ma", 0.1);
  T Re = args.getValueOrFallback( "--Re", 0.);
  const T targetMa = args.getValueOrFallback( "--targetMa", 0.1);
  const T lengthDomain = args.getValueOrFallback( "--lengthDomain", 2.);
  T tau = args.getValueOrFallback( "--tau", 0.);
  size_t imax = args.getValueOrFallback( "--imax", 10000);
  T tmax = args.getValueOrFallback( "--tmax", 1.);
  int iout = args.getValueOrFallback( "--iout", 0);
  int nout = args.getValueOrFallback( "--nout", 20);
  const int boundary_condition = args.getValueOrFallback( "--boundary_condition", 3);
  const bool debug = args.contains("--debug");
  const int boundary_depth = args.getValueOrFallback( "--boundary_depth", 10);
  const T damping_strength = args.getValueOrFallback( "--damping_strength", 1.);

  /* default values */
  BoundaryType boundarytype;
  switch ( boundary_condition ) {
    case 1: boundarytype = periodic; clout << "Boundary condition type specified to periodic." << std::endl; break;
    case 2: boundarytype = local; clout << "Boundary condition type specified to local." << std::endl; break;
    case 3: boundarytype = damping; clout << "Boundary condition type specified to damping." << std::endl; break;
    default: boundarytype = damping; clout << "Boundary condition type not specified. Default to damping." << std::endl; break;
  }
  T heightDomain            = 0.5*lengthDomain;
  T depthDomain             = 0.25*lengthDomain;

  const T charL             = lengthDomain;
  T charV                   = std::max(Ma*1.1, 0.1);  // /targetMa;
  const T cs_LU             = 1 / std::sqrt(3.0);
  T viscosity;
  if ( tau == 0 ) {
    if ( Re == 0 ) {
      viscosity             = 1.48e-5; // 1e-2;
      Re                    = charV * charL / viscosity;
    } else viscosity        = charV * charL / Re;
    tau                     = viscosity / (cs_LU*cs_LU) + 0.5;
  }
  else {
    if ( Re == 0 ) Re       = 20000;
    viscosity               = charV * charL / Re;
  }
  if ( debug ) {
    Re                      = 200;
    charV                   = 0.25;
    heightDomain            = lengthDomain;
    depthDomain             = lengthDomain;
    viscosity               = charV * charL / Re;
    tau                     = viscosity / (cs_LU*cs_LU) + 0.5;
  }

  clout << "Input to unitconverter: res=" << res << ", tau=" << tau << ", lengthDomain=" << lengthDomain << ", charV=" << charV << ", viscosity=" << viscosity << ", rho0=" << rho0 << std::endl;
  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> converter(
    (size_t) res,               // resolution
    (T)      tau,               // relaxation time
    (T)      lengthDomain,      // charPhysLength: reference length of simulation geometry
    (T)      charV,             // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)      viscosity,         // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)      rho0               // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("movingshock3d");

  // === 2nd Step: Prepare Geometry ===

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 6;
#endif

  // setup domain
  Vector<T,3> domain_lengths = {lengthDomain, heightDomain, depthDomain};
  Vector<T,3> originDomain( - lengthDomain/2, - heightDomain/2, - depthDomain/2 );  //
  Vector<T,3> extendDomain( lengthDomain, heightDomain, depthDomain );  // size of the domain
  std::shared_ptr<IndicatorF3D<T>> domain =
      std::make_shared<IndicatorCuboid3D<T>>( extendDomain, originDomain );

  CuboidGeometry3D<T> cuboidGeometry(*(domain),
                                     converter.getConversionFactorLength(),
                                     noOfCuboids
                                     );
  if ( boundarytype == periodic || boundarytype == damping ) {
    cuboidGeometry.setPeriodicity(true, true, true);
  } else {
    cuboidGeometry.setPeriodicity(false, false, false);
  }

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,3> superGeometry( cuboidGeometry, loadBalancer );

  T boundary_depth_pu = converter.getPhysLength(boundary_depth);
  clout << "Setup: debug=" << debug << "; boundary_depth=" << boundary_depth << "; bd_depth_pu=" << boundary_depth_pu << "; damping_strength=" << damping_strength << "; overlap=" << superGeometry.getOverlap() << std::endl;
  prepareGeometry( converter, superGeometry, domain, boundarytype, res, debug, boundary_depth, domain_lengths );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T,DESCRIPTOR> sLattice( superGeometry, debug );

  //prepare Lattice and set boundaryConditions
  prepareLattice( converter, sLattice, superGeometry, rho0, Ma, boundarytype, boundary_depth, damping_strength, domain_lengths );

  // instantiate reusable functors
  SuperPlaneIntegralFluxVelocity3D<T> velocityFlux(
      sLattice,
      converter,
      superGeometry,
      { lengthDomain/T(2), heightDomain/T(2), depthDomain/T(2) },
      Vector<T, 3>( 0., 1., 1. )
  );

  SuperPlaneIntegralFluxPressure3D<T> pressureFlux(
      sLattice,
      converter,
      superGeometry,
      { lengthDomain/T(2), heightDomain/T(2), depthDomain/T(2) },
      Vector<T, 3>( 0., 1., 1. )
  );

  Gnuplot<T> gplot_l2_abs("l2_absolute");//, Gnuplot<T>::LOGLOG, Gnuplot<T>::OFF);
  gplot_l2_abs.setLabel("time []", "absolute L2 norm []");
  T Lp0 = L2Norm( sLattice, superGeometry, 0, converter );

  clout << "tmax=" << tmax << "[PU]=" << converter.getLatticeTime( tmax ) << "[LU], while imax=" << imax << " as input." << std::endl;
  imax = std::min( converter.getLatticeTime( tmax ), imax );
  tmax = converter.getPhysTime( imax );
  if ( debug ) {            imax=100;         iout=10;  }
  if ( iout == 0 )          iout            = imax / nout;
  clout << "tmax=" << tmax << ", while imax=" << imax << " recalculated. iout=" << iout << std::endl;

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer(converter.getLatticeTime( tmax ),
                       superGeometry.getStatistics().getNvoxel()
                       );
  timer.start();

  size_t iT = 0;
  while ( iT < imax ) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( converter, sLattice, iT, superGeometry, rho0, Ma, boundarytype );
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults(sLattice, converter, iT, superGeometry, timer,
               velocityFlux, pressureFlux,
               rho0,
               gplot_l2_abs, Lp0,
               iout, tmax, imax
               );
    // test density failure
    if (sLattice.getStatistics().getAverageRho() > 2.) {
      clout << "breaking simulation loop because density is too high: "
            << sLattice.getStatistics().getAverageRho() << std::endl;
      break;
    }
    // if (sLattice.getStatistics().getAverageEnergy() > 2.) {
    //   clout << "breaking simulation loop because energy is too high: "
    //         << sLattice.getStatistics().getAverageEnergy() << std::endl;
    //   break;
    // }
    iT++;
  }
  clout << "Simulation stopped after " << iT << " iterations." << std::endl;

  gplot_l2_abs.setYrange(1e-3, 1);
  gplot_l2_abs.setLogScale(2);
  gplot_l2_abs.writePNG(-1, -1, "gplot_l2_abs");
  timer.stop();
  timer.printSummary();
}
