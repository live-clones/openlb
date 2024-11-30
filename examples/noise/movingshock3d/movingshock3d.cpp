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
 * Method from H. Xu and P. Sagaut, “Analysis of the absorbing layers for the weakly-compressible lattice Boltzmann methods,” Journal of Computational Physics, vol. 245, pp. 14–42, Jul. 2013, doi: 10.1016/j.jcp.2013.02.051.
 * Benchmark from C. Bogey and C. Bailly, “Three-dimensional non-reﬂective boundary conditions for acoustic simulations: far ﬁeld formulation and validation test cases,” vol. 88, 2002.
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
using DESCRIPTOR = D3Q19<>;
using BulkDynamics = KBCdynamics<T,DESCRIPTOR>;
// using DESCRIPTOR = D3Q27<>;
// using BulkDynamics = BGKdynamics<T,DESCRIPTOR>;

const int noMat   = 0;
const int checMat = 2;
const int dampMat = 1;
const int fluiMat = 3;
const int inflMat = 4;
const int outfMat = 5;
const unsigned int ndim = 3;

typedef enum {periodic, bounceBack, local, damping, dampingAndLocal} BoundaryType;
typedef enum {shock, pointsource} SourceType;

// Stores geometry information in form of material numbers
void prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperGeometry<T,ndim>& superGeometry,
                     std::shared_ptr<IndicatorF<T,ndim>> domain,
                     BoundaryType boundarytype,
                     SourceType source,
                     int res,
                     bool debug,
                     int boundary_depth
                     )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl
        << "Materials:" << std::endl
        << noMat   << " = no Material (default, should be immediately renamed to Check)" << std::endl
        << fluiMat << " = Fluid" << std::endl
        << checMat << " = Check (should be empty)" << std::endl
        << dampMat << " = far field" << std::endl
        << inflMat << " = inflow" << std::endl
        << outfMat << " = outflow" << std::endl;

  superGeometry.rename( noMat, checMat );   // all nodes to temporary type

  // if ( boundarytype == dampingAndLocal ) {

  //   superGeometry.rename( checMat, noMat, {1,0,0} );  // temporarily set inside except outlet/inlet to 0
  //   superGeometry.rename( noMat, checMat );  // set inside except outlet/inlet back to noMat
  // }

  if ( boundarytype == dampingAndLocal ) {
    boundary_depth += 1;
  }
  T bd_pu = converter.getPhysLength(boundary_depth);
  Vector<T,ndim> originFluid = superGeometry.getStatistics().getMinPhysR( checMat ) + bd_pu;
  Vector<T,ndim> extendFluid = superGeometry.getStatistics().getMaxPhysR( checMat ) - superGeometry.getStatistics().getMinPhysR( checMat ) - 2*bd_pu;
  std::shared_ptr<IndicatorF<T,ndim>> fluid_domain = std::make_shared<IndicatorCuboid<T,ndim>>( extendFluid, originFluid );
  if ( boundarytype != dampingAndLocal ) superGeometry.rename( noMat, checMat, fluid_domain );
  superGeometry.getStatistics().print();
  switch ( boundarytype ) {
    case periodic: superGeometry.rename( noMat, fluiMat ); superGeometry.clean(); break; // all nodes are fluidbreak;
    case bounceBack: superGeometry.rename( noMat, dampMat, fluiMat, domain ); superGeometry.clean(); break;
    case local: superGeometry.rename( noMat, dampMat, fluiMat, domain ); superGeometry.clean(); break;
    case damping: superGeometry.rename( noMat, dampMat ); break;  // all remaining nodes are far-field
    case dampingAndLocal:
      // fluid domain including damping layer to damping layer
      Vector<T,ndim> origin = superGeometry.getStatistics().getMinPhysR( checMat );
      Vector<T,ndim> extend = superGeometry.getStatistics().getMaxPhysR( checMat ) - superGeometry.getStatistics().getMinPhysR( checMat );
      origin[0] += converter.getConversionFactorLength();
      extend[0] -= 2*converter.getConversionFactorLength();
      std::shared_ptr<IndicatorF<T,ndim>> dampingField = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      // clout << "checMat=" << checMat << " to dampMat=" << dampMat << " from " << origin << " to " << extend  << std::endl;
      superGeometry.rename( checMat, dampMat, dampingField );  // all remaining nodes are far-field, except...
      // superGeometry.getStatistics().print();

      // fluid domain part to fluid
      // clout << "dampMat=" << dampMat << " to fluidMat=" << fluiMat << " from " << originFluid << " to " << extendFluid  << std::endl;
      superGeometry.rename( dampMat, fluiMat, fluid_domain );  // all nodes except boundary length
      // superGeometry.getStatistics().print();
  
      // Set material number for inflow
      origin[0] -= 2*converter.getConversionFactorLength();
      extend[0] = 2*converter.getConversionFactorLength();
      std::shared_ptr<IndicatorF<T,ndim>> inflow = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      // clout << "checMat=" << checMat << " to inflMat=" << inflMat << " pointing to dampMat=" << dampMat << " from " << origin << " to " << extend  << std::endl;
      superGeometry.rename( checMat, inflMat, dampMat, inflow );  // inflow to bc, pointing to porous
      // superGeometry.getStatistics().print();
      #
      // Set material number for outflow
      origin[0] = superGeometry.getStatistics().getMaxPhysR( checMat )[0]-converter.getConversionFactorLength();
      std::shared_ptr<IndicatorF<T,ndim>> outflow = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      // superGeometry.rename( noMat, noMat, outflow );  // outflow back to temporary
      // clout << "checMat=" << checMat << " to outfMat=" << outfMat << " pointing to dampMat=" << dampMat << " from " << origin << " to " << extend  << std::endl;
      superGeometry.rename( checMat, outfMat, dampMat, outflow );  // outflow to bc, pointing to porous
      // superGeometry.getStatistics().print();
      break;
  }

  superGeometry.communicate();
  // Removes all not needed boundary voxels inside the surface
  // superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.updateStatistics();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(UnitConverter<T,DESCRIPTOR> const& converter,
                    SuperLattice<T,DESCRIPTOR>& sLattice,
                    SuperGeometry<T,ndim>& superGeometry,
                    T rho0,
                    T Ma,
                    T amplitude,
                    T alpha,
                    BoundaryType boundarytype,
                    SourceType source,
                    int boundary_depth,
                    Vector<T,ndim> domain_lengths,
                    T damping_strength
                    )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=fluiMat --> bulk dynamics
  auto bulkIndicator = superGeometry.getMaterialIndicator({fluiMat});
  sLattice.defineDynamics<BulkDynamics>(bulkIndicator);

  // Material=dampMat --> dynamics depend..
  switch ( boundarytype ) {
    case periodic: sLattice.defineDynamics<BulkDynamics>(superGeometry.getMaterialIndicator({dampMat})); break;
    case local:
      setLocalPressureBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, dampMat);
      setLocalVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, dampMat);
      break;
    case bounceBack: setBounceBackBoundary<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator({dampMat})); break;
    case damping:
      bulkIndicator = superGeometry.getMaterialIndicator({fluiMat, dampMat});
      setDampingBoundary<T,DESCRIPTOR>(sLattice, superGeometry, dampMat);
      break;
    case dampingAndLocal:
      domain_lengths -= 2*converter.getConversionFactorLength();
      bulkIndicator = superGeometry.getMaterialIndicator({fluiMat, dampMat});
      setDampingBoundary<T,DESCRIPTOR>(sLattice, superGeometry, dampMat);
      // Material=inflMat -->inlet
      setLocalVelocityBoundary(sLattice, omega, superGeometry, inflMat);
      // Material=outfMat -->outlet
      setInterpolatedConvectionBoundary(sLattice, omega, superGeometry, outfMat);
  }

  // Initial conditions
  AnalyticalConst<ndim,T,T> ux = AnalyticalConst<ndim,T,T>( Ma );
  AnalyticalConst<ndim,T,T> uy = AnalyticalConst<ndim,T,T>( 0. );
  AnalyticalConst<ndim,T,T> uz = AnalyticalConst<ndim,T,T>( 0. );
  AnalyticalComposed<ndim,T,T> u( ux, uy, uz );

  T ndatapoints = converter.getResolution();
  T dist = converter.getPhysDeltaX();
  if ( source == shock ) {
    AcousticPulse<3,T> pressureProfile( rho0, amplitude, alpha );
    
    Gnuplot<T> gplot_hline_p( "shock_hline");
    Gnuplot<T> gplot_diag_p( "shock_diag");
    gplot_hline_p.setLabel("distance [m]", "density [LU]");
    gplot_diag_p.setLabel("distance [m]", "density [LU]");
    for (int n = 0; n <= int(ndatapoints/2); n++) {
      T input_hline[3] =  {n*dist, 0, 0};
      T output_hline[3];
      pressureProfile(output_hline, input_hline);
      gplot_hline_p.setData(input_hline[0], output_hline[0]);
      T input_diag[3] =  {n*dist, n*dist, 0};
      T output_diag[3];
      pressureProfile(output_diag, input_diag);
      gplot_diag_p.setData(input_diag[0], output_diag[0]);
    }
    gplot_hline_p.writePNG(-1, -1, "shock_hline");
    gplot_diag_p.writePNG(-1, -1, "shock_diag");

    //Initialize all values of distribution functions to their local equilibrium
    sLattice.defineRhoU( bulkIndicator, pressureProfile, u );
    sLattice.iniEquilibrium( bulkIndicator, pressureProfile, u );
  } else if ( source == pointsource ) {
    AnalyticalConst<ndim,T,T> rho(rho0);
    //Initialize all values of distribution functions to their local equilibrium
    sLattice.defineRhoU( superGeometry.getMaterialIndicator({3}), rho, u );
    sLattice.iniEquilibrium( superGeometry.getMaterialIndicator({3}), rho, u );
  }

  sLattice.setParameter<descriptors::OMEGA>( omega );
  // Define fields; could also use setParameter as for omega, but then I could not use a flow profile
  sLattice.defineField<descriptors::UX>( bulkIndicator, ux );
  sLattice.defineField<descriptors::UY>( bulkIndicator, uy );
  sLattice.defineField<descriptors::UZ>( bulkIndicator, uz );
  AnalyticalConst<ndim,T,T> rhoField( rho0 );
  sLattice.defineField<descriptors::DENSITY>( bulkIndicator, rhoField );

  // output damping layer parameter
  DampingTerm<3,T,DESCRIPTOR> sigma_plot( converter, boundary_depth, domain_lengths );
  Gnuplot<T> gplot_hline( "sigma_hline");
  Gnuplot<T> gplot_diag( "sigma_diag");
  gplot_hline.setLabel("distance [m]", "density [LU]");
  gplot_diag.setLabel("distance [m]", "density [LU]");
  for (int n = 0; n <= int(ndatapoints/2); n++) {
    T input_hline[3] =  {n*dist, 0, 0};
    T output_hline[3];
    sigma_plot(output_hline, input_hline);
    gplot_hline.setData(input_hline[0], output_hline[0]);
    T input_diag[3] =  {n*dist, n*dist, 0};
    T output_diag[3];
    sigma_plot(output_diag, input_diag);
    gplot_diag.setData(input_diag[0], output_diag[0]);
  }
  gplot_hline.writePNG(-1, -1, "sigma_hline");
  gplot_diag.writePNG(-1, -1, "sigma_diag");

  DampingTerm<3,T,DESCRIPTOR> sigma( converter, boundary_depth, domain_lengths, damping_strength );
  sLattice.defineField<descriptors::DAMPING>( superGeometry.getMaterialIndicator({ dampMat }), sigma );
  
  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a sinusoidal pressure
void setBoundaryValues( UnitConverter<T,DESCRIPTOR> const& converter,
                        SuperLattice<T,DESCRIPTOR>& sLattice, int iT,
                        SuperGeometry<T,3>& superGeometry,
                        T rho0, T Ma,
                        BoundaryType boundarytype,
                        SourceType source,
                        T amplitude )
{
  OstreamManager clout( std::cout,"setBoundaryValues" );
  AnalyticalConst<ndim,T,T> ux( Ma );
  AnalyticalConst<ndim,T,T> uy( 0. );
  AnalyticalConst<ndim,T,T> uz( 0. );
  AnalyticalComposed<ndim,T,T> u( ux, uy, uz );
  AnalyticalConst<ndim,T,T> rho( rho0 );

  auto farFieldIndicator = superGeometry.getMaterialIndicator({dampMat});
  if ( boundarytype == dampingAndLocal ) farFieldIndicator = superGeometry.getMaterialIndicator({inflMat,outfMat});
  sLattice.defineRhoU( farFieldIndicator, rho, u );
  sLattice.iniEquilibrium( farFieldIndicator, rho, u );

  sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
}

T L2Norm(SuperLattice<T,DESCRIPTOR>& sLattice,
         SuperGeometry<T,ndim>& superGeometry,
         int iT,
         UnitConverter<T,DESCRIPTOR> const& converter) {
  OstreamManager clout(std::cout, "L2Norm");

  T result[ndim];
  int tmp[] = {int()};

  SuperLatticePhysPressure3D<T,DESCRIPTOR> pressurenD(sLattice, converter);
  AnalyticalConst<ndim,T,T> rho0nD( 0. );

  auto indicatorF = superGeometry.getMaterialIndicator({fluiMat});
  SuperAbsoluteErrorL2Norm3D<T> absPressureErrorL2Norm(pressurenD, rho0nD, indicatorF);

  absPressureErrorL2Norm(result, tmp);
  T l2_abs = result[0];

  return l2_abs;
}

// write data to termimal and file system
void getResults(SuperLattice<T,DESCRIPTOR>& sLattice,
                UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                SuperGeometry<T,ndim>& superGeometry, util::Timer<T>& timer,
                SuperPlaneIntegralFluxVelocity3D<T>& velocityFlux,
                SuperPlaneIntegralFluxPressure3D<T>& pressureFlux,
                T rho0, T amplitude,
                Gnuplot<T>& gplot_l2_abs, T Lp0,
                size_t iout, T tmax, size_t imax,
                SourceType sourcetype
                )
{
#ifdef PLATFORM_GPU_CUDA
  gpu::cuda::device::synchronize();
#endif
  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter3D<T> vtmWriter( "movingshock3d" );

  if ( iT==0 ) {
    // Writes geometry, cuboid no. and rank no. to file system
    SuperLatticeGeometry3D<T,DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid3D<T,DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T,DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Writes to console 10 times
  if ( iT%iout==0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);  // important for synchronization (?) on GPU
    // velocityFlux.print();
    // pressureFlux.print();
    // write to terminal
    timer.update( iT );
    timer.printStep();
    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    if ( iT%std::max(int(iout/10), 1)==0 ) {
      gplot_l2_abs.setData(T(iT), L2Norm(sLattice, superGeometry, iT,
                                        converter) / Lp0 );
    }

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
    T ymin(converter.getPhysPressure(-amplitude/200));
    T ymax(converter.getPhysPressure(+amplitude/200));
    gplot_hline.setYrange(ymin, ymax);
    gplot_vline.setYrange(ymin, ymax);
    gplot_diagonal.setYrange(ymin, ymax);
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
      jpeg_ParamP.maxValue = converter.getPhysPressure(+amplitude/200);
      jpeg_ParamP.minValue = converter.getPhysPressure(-amplitude/200);
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
  // if ( singleton::mpi().isMainProcessor() ) {
  //   std::copy(argv, argv + argc, std::ostream_iterator<char *>(std::cout, " "));
  //   printf("\n");
  // }
  // std::copy(argv, argv + argc, std::ostream_iterator<char *>(std::cout, " "));
  // printf("\n");
  CLIreader args(argc, argv);
  std::string outdir = args.getValueOrFallback<std::string>( "--outdir", "./tmp" );

  // === get Command Line Arguments
  const int res = args.getValueOrFallback( "--res", 81);
  const T rho0 = args.getValueOrFallback( "--rho0", 1.);
  const T Ma = args.getValueOrFallback( "--Ma", 0.1);
  T Re = args.getValueOrFallback( "--Re", 0.);
  const T lengthDomain = args.getValueOrFallback( "--lengthDomain", 2.);
  const T heightDomain = args.getValueOrFallback( "--heightDomain", 4.);
  const T depthDomain = args.getValueOrFallback( "--depthDomain", 4.);
  T tau = args.getValueOrFallback( "--tau", 0.);
  size_t imax = args.getValueOrFallback( "--imax", 300);
  T tmax = args.getValueOrFallback( "--tmax", 1.);
  int iout = args.getValueOrFallback( "--iout", 0);
  int nout = args.getValueOrFallback( "--nout", 20);
  const int boundary_condition = args.getValueOrFallback( "--boundary_condition", 3);
  const int source_type = args.getValueOrFallback( "--source_type", 1);
  const bool debug = args.contains("--debug");
  const size_t boundary_depth = args.getValueOrFallback( "--boundary_depth", 10);
  const T amplitude = args.getValueOrFallback( "--amplitude", 0.1);
  const T damping_strength = args.getValueOrFallback( "--damping_strength", 1.);
  size_t overlap = args.getValueOrFallback( "--overlap", 3);

  // initialize arguments
  const T b_shock           = 3;
  const T alpha_shock       = log(2)/(b_shock*b_shock);

  // managing outputs
  std::stringstream outdir_mod;
  outdir_mod << outdir;
  switch ( source_type ) {
    case 1: outdir_mod << "_shock"; break;
    case 2: outdir_mod << "_pointsource"; break;
    default: outdir_mod << "_shock"; break;
  }
  switch ( boundary_condition ) {
    case 1: outdir_mod << "_periodic"; break;
    case 2: outdir_mod << "_local"; break;
    case 3: outdir_mod << "_damping"; break;
    case 4: outdir_mod << "_dampingAndLocal"; break;
    default: outdir_mod << "_damping"; break;
  }
  outdir_mod << "_" << Ma << "_Re" << Re << "_a" << amplitude << "_" << lengthDomain << "x" << heightDomain << "x" << depthDomain << "_res" << res << "_overlap" << overlap << "_bd" << boundary_depth << "x" << damping_strength;

  singleton::directories().setOutputDir( outdir_mod.str()+"/" );  // set output directory
  OstreamManager clout( std::cout, "main" );
  clout << "outdir set to " << outdir_mod.str() << std::endl;

  /* default values */
  BoundaryType boundarytype;
  switch ( boundary_condition ) {
    case 1: boundarytype = periodic; clout << "Boundary condition type specified to periodic." << std::endl; break;
    case 2: boundarytype = local; clout << "Boundary condition type specified to local." << std::endl; break;
    case 3: boundarytype = damping; clout << "Boundary condition type specified to damping." << std::endl; break;
    case 4: boundarytype = dampingAndLocal; clout << "Boundary condition type specified to damping with local outside bc." << std::endl; break;
    default: boundarytype = damping; clout << "Boundary condition type not specified. Default to damping." << std::endl; break;
  }

  SourceType source;
  switch ( source_type ) {
    case 1: source = shock; clout << "Source type specified to shock." << std::endl; break;
    case 2: source = pointsource; clout << "Source type specified to point source." << std::endl; break;
    default: source = shock; clout << "Source type not specified. Default to shock." << std::endl; break;
  }

  // determining Reynolds regime (incl. viscosity and relaxation time)
  const T charL             = lengthDomain;
  T charV                   = std::max(Ma*1.1, 0.1);  // /targetMa;
  const T cs_LU             = 1 / std::sqrt(3.0);
  T viscosity;
  if ( tau == 0 ) {
    if ( Re == 0 ) {
      viscosity             = 1.48e-5; // 1e-2;
      Re                    = charV * charL / viscosity;
    } else {
      viscosity        = charV * charL / Re;
    }
    tau                     = viscosity / (cs_LU*cs_LU) + 0.5;
  }
  else {
    if ( Re == 0 ) Re       = 20000;
    viscosity               = charV * charL / Re;
  }
  if ( debug ) {
    Re                      = 200;
    charV                   = 0.25;
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

  if ( noOfCuboids > 1 && ( boundarytype == damping || boundarytype == dampingAndLocal) ) {
    clout << "overlap=" << overlap << "; boundary_depth_LU=" << boundary_depth << "; setting overlap to >=boundary_depth." << std::endl;
    overlap = std::max(overlap, boundary_depth);
  }

  // setup domain
  Vector<T,ndim> domain_lengths = {lengthDomain, heightDomain, depthDomain};
  Vector<T,ndim> originDomain( - lengthDomain/2, - heightDomain/2, - depthDomain/2 );  //
  Vector<T,ndim> extendDomain( lengthDomain, heightDomain, depthDomain );  // size of the domain
  std::shared_ptr<IndicatorF<T,ndim>> domain = std::make_shared<IndicatorCuboid<T,ndim>>( extendDomain, originDomain );

  CuboidGeometry<T,ndim> cuboidGeometry(*(domain),
                                        converter.getConversionFactorLength(),
                                        noOfCuboids
                                        );
  if ( boundarytype == dampingAndLocal ) cuboidGeometry.setPeriodicity(false,true,true);
  else cuboidGeometry.setPeriodicity(true, true, true);

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,ndim> superGeometry( cuboidGeometry, loadBalancer, overlap );

  T boundary_depth_pu = converter.getPhysLength(boundary_depth);
  clout << "Setup: debug=" << debug << "; boundary_depth=" << boundary_depth << "; bd_depth_pu=" << boundary_depth_pu << "; overlap=" << superGeometry.getOverlap() << std::endl;
  prepareGeometry( converter, superGeometry, domain, boundarytype, source, res, debug, boundary_depth );
  
  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T,DESCRIPTOR> sLattice( superGeometry, debug );

  //prepare Lattice and set boundaryConditions
  prepareLattice( converter, sLattice, superGeometry, rho0, Ma, amplitude, alpha_shock, boundarytype, source, boundary_depth, domain_lengths, damping_strength );

  // instantiate reusable functors
  SuperPlaneIntegralFluxVelocity3D<T> velocityFlux(
      sLattice,
      converter,
      superGeometry,
      { lengthDomain/T(2), heightDomain/T(2), depthDomain/T(2) },
      Vector<T, ndim>( 0., 1., 1. )
  );

  SuperPlaneIntegralFluxPressure3D<T> pressureFlux(
      sLattice,
      converter,
      superGeometry,
      { lengthDomain/T(2), heightDomain/T(2), depthDomain/T(2) },
      Vector<T, ndim>( 0., 1., 1. )
  );

  Gnuplot<T> gplot_l2_abs("l2_absolute");//, Gnuplot<T>::LOGLOG, Gnuplot<T>::OFF);
  gplot_l2_abs.setLabel("time []", "absolute L2 norm []");
  T Lp0 = L2Norm( sLattice, superGeometry, 0, converter );

  SuperVTMwriter<T,ndim> vtmWriter_init( "movingshock3d_init" );
  vtmWriter_init.createMasterFile();
  SuperLatticePhysField3D<T,DESCRIPTOR,DAMPING> damping( sLattice, 1. );
  damping.getName() = "dampingField";
  SuperLatticePhysField3D<T,DESCRIPTOR,UX> ux( sLattice, 1. );
  ux.getName() = "uxField";
  SuperLatticePhysField3D<T,DESCRIPTOR,DENSITY> density( sLattice, 1. );
  density.getName() = "densityField";
  vtmWriter_init.addFunctor( damping );
  vtmWriter_init.addFunctor( ux );
  vtmWriter_init.addFunctor( density );
  vtmWriter_init.write( 0 );

  clout << "tmax=" << tmax << "[PU]=" << converter.getLatticeTime( tmax ) << "[PU], while imax=" << imax << " as input." << std::endl;
  imax = std::min( converter.getLatticeTime( tmax ), imax );
  tmax = converter.getPhysTime( imax );
  clout << "tmax=" << tmax << ", while imax=" << imax << " recalculated." << std::endl;
  clout << "nout=" << nout << "; iout=" << iout << " as input." << std::endl;
  nout = std::max(1, nout);
  if ( iout == 0 ) iout = int(imax / nout);
  iout = std::min(iout, int(imax / nout));
  clout << "nout=" << nout << "; iout=" << iout << " recalculated." << std::endl;

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer(converter.getLatticeTime( tmax ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  size_t iT = 0;
  while ( iT < imax ) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    if ( boundarytype == local || boundarytype == dampingAndLocal ) {
      setBoundaryValues( converter, sLattice, iT, superGeometry, rho0, Ma, boundarytype, source, amplitude );
    }
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults(sLattice, converter, iT, superGeometry, timer,
               velocityFlux, pressureFlux,
               rho0, amplitude,
               gplot_l2_abs, Lp0,
               iout, tmax, imax,
               source );

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
