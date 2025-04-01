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
 * Use the option --boundaryCondition to vary how the far field is modeled:
 *  0=eternal (2x domain size to capture the vanishing shock as ideal reference)
 *  1=periodic,
 *  2=local,
 *  3=damping (with periodic boundaries),
 *  4=dampingAndLocal,
 *  5=dampingAndZeroGradient, (zero gradient boundary around flow)
 *  6=localAndZeroGrad (zero gradient boundary around flow)
 *  7=localAndConvection (local velocity around, interpolated convection outlet
 *  8=interpolated (interpolated velocity and pressure)
 *  9=interpolatedAndZeroGrad (same but zero gradient around flow)
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
// using DESCRIPTOR = D3Q27<>;
using BulkDynamics = BGKdynamics<T,DESCRIPTOR>;
// using BulkDynamics = KBCdynamics<T,DESCRIPTOR>;
// using BulkDynamics = olb::dynamics::Tuple<T,DESCRIPTOR,momenta::BulkTuple,equilibria::Incompressible,collision::TRT>;

const int noMat   = 0;
const int dampMat = 1;
const int checMat = 2;
const int fluiMat = 3;
const int inflMat = 4;
const int outfMat = 5;
const int rounMat = 6;
const int ndim = 3;  // a few things (e.g. SuperSum3D cannot be adapted to 2D, but this should help speed it up)

typedef enum {eternal, periodic, local, damping, dampingAndLocal, dampingAndZeroGrad, localAndZeroGrad, localAndConvection, interpolated, interpolatedAndZeroGrad} BoundaryType;
typedef enum {horizontal, vertical, diagonal2d, diagonal3d} SamplingDirection;

T* uAverage = NULL;

// Stores geometry information in form of material numbers
void prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperGeometry<T,ndim>& superGeometry,
                     IndicatorF3D<T>& domain,
                     IndicatorF3D<T>& fluidDomain,
                     BoundaryType boundarytype,
                     int res,
                     int boundaryDepth
                     )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << std::endl << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( noMat, checMat );   // all nodes to temporary type

  T dx = converter.getConversionFactorLength();
  Vector<T,ndim> origin, extend, origin_d, extend_d;
  
  switch ( boundarytype ) {
    // eternal and damping: fluiMat is the actual fluid; periodic: dampMat is the fluid
    case eternal:
    case damping:
      superGeometry.rename( checMat, fluiMat, fluidDomain );  // just to be able to get uAverage for inside
    case periodic:
      superGeometry.rename( checMat, dampMat );
      break;
    
    case local:
    case localAndZeroGrad:
    case localAndConvection:
    case interpolated:
    case interpolatedAndZeroGrad:
      // rename fluid area (must be 1 for normals)
      superGeometry.rename( checMat, dampMat, fluidDomain );
    case dampingAndLocal:
    case dampingAndZeroGrad: {
      if ( boundarytype == dampingAndLocal || boundarytype == dampingAndZeroGrad ) {
        // rename damping areas
        origin = superGeometry.getStatistics().getMinPhysR( checMat ) + dx;
        extend = superGeometry.getStatistics().getMaxPhysR( checMat ) - superGeometry.getStatistics().getMinPhysR( checMat ) - 2*dx;
        IndicatorCuboid3D<T> dampingField( extend, origin );
        superGeometry.rename( checMat, dampMat, dampingField );
        // rename fluid area within damping area
        superGeometry.rename( dampMat, fluiMat, fluidDomain );
      }
  
      // Set material number for inflow
      origin = superGeometry.getStatistics().getMinPhysR( checMat ) - dx;
      extend = superGeometry.getStatistics().getMaxPhysR( checMat ) - superGeometry.getStatistics().getMinPhysR( checMat ) + 2*dx;
      extend[0] = 2*dx;
      IndicatorCuboid3D<T> inflow( extend, origin );
      superGeometry.rename( checMat, inflMat, dampMat, inflow );
      
      // Set material number for outflow
      origin[0] = superGeometry.getStatistics().getMaxPhysR( checMat )[0] - dx;
      IndicatorCuboid3D<T> outflow( extend, origin );
      superGeometry.rename( checMat, outfMat, dampMat, outflow );

      // Set remaining to outside
      origin = superGeometry.getStatistics().getMinPhysR( checMat ) - dx;
      extend = superGeometry.getStatistics().getMaxPhysR( checMat ) - superGeometry.getStatistics().getMinPhysR( checMat ) + 2*dx;
      extend[1] = 2*dx;
      IndicatorCuboid3D<T> bottom( extend, origin );
      superGeometry.rename( checMat, rounMat, dampMat, bottom );

      origin[1] = superGeometry.getStatistics().getMaxPhysR( checMat )[1] - dx;
      IndicatorCuboid3D<T> top( extend, origin );
      superGeometry.rename( checMat, rounMat, dampMat, top );

      origin = superGeometry.getStatistics().getMinPhysR( checMat ) - dx;
      extend = superGeometry.getStatistics().getMaxPhysR( checMat ) - superGeometry.getStatistics().getMinPhysR( checMat ) + 2*dx;
      extend[2] = 2*dx;
      IndicatorCuboid3D<T> front( extend, origin );
      superGeometry.rename( checMat, rounMat, dampMat, front );

      origin[2] = superGeometry.getStatistics().getMaxPhysR( checMat )[2] - dx;
      IndicatorCuboid3D<T> back( extend, origin );
      superGeometry.rename( checMat, rounMat, dampMat, back );
      break;
    }
  }

  superGeometry.communicate();
  // superGeometry.clean() deletes some nodes although it shouldn't - this is probably to do with the assumption that 1 is fluid and only needs one layer next to it (rather than 20 for pml)
  superGeometry.checkForErrors();
  superGeometry.updateStatistics();
  clout << "Materials:" << std::endl
        << noMat   << " = no Material (default, should be immediately renamed to Check)" << std::endl
        << dampMat << " = far field (must be 1 so outside boundaries can calculate their normal to the damping area)" << std::endl
        << checMat << " = check (should be empty, but corners may have issues with the normal; check is therefore treated as inflow" << std::endl
        << fluiMat << " = fluid (sometimes must be set to 1 to work with normals of boundary conditions)" << std::endl
        << inflMat << " = inflow" << std::endl
        << outfMat << " = outflow" << std::endl
        << rounMat << " = around domain" << std::endl;
  superGeometry.getStatistics().print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void setUAverage( SuperLattice<T,DESCRIPTOR>& sLattice,
                  SuperGeometry<T,ndim>& superGeometry,
                  int fluidMaterial
                  )
{
  sLattice.setProcessingContext(ProcessingContext::Evaluation);  // important for synchronization (?) on GPU

  // update outflow boundary value (adaptive convection boundary for smaller domains)
  SuperLatticeVelocity3D<T, DESCRIPTOR> velocity( sLattice );
  SuperSum3D<T> sum( velocity, superGeometry, fluidMaterial );
  int input[1];
  T output[3];
  sum(output, input);
  *uAverage = output[0] / superGeometry.getStatistics().getNvoxel( fluidMaterial );

  sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
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
                    int boundaryDepth,
                    Vector<T,ndim> domainLengths,
                    T dampingStrength
                    )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << std::endl << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();
  T tmp = T();
  uAverage = &tmp;  //*uAverage = 0.0;
  if ( ( boundarytype == local ) || ( boundarytype == dampingAndLocal ) ) {
    setUAverage(sLattice, superGeometry, boundarytype);    
  }

  // Material=fluiMat --> bulk dynamics
  std::unique_ptr<SuperIndicatorF<T,ndim>> bulkIndicator = superGeometry.getMaterialIndicator( {dampMat, fluiMat, checMat} );
  sLattice.defineDynamics<BulkDynamics>( bulkIndicator );

  std::unique_ptr<SuperIndicatorF<T,ndim>> inletIndicator = superGeometry.getMaterialIndicator( {inflMat} );
  std::unique_ptr<SuperIndicatorF<T,ndim>> aroundIndicator = superGeometry.getMaterialIndicator( {rounMat} );
  switch ( boundarytype ) {
    case eternal:
    case periodic:
      break;

    case local:
    case localAndZeroGrad:
      boundary::set<boundary::LocalVelocity>( sLattice, inletIndicator );
      if ( boundarytype == local ) boundary::set<boundary::LocalVelocity>( sLattice, aroundIndicator );
      if ( boundarytype == localAndZeroGrad ) setZeroGradientBoundary<T,DESCRIPTOR>( sLattice, aroundIndicator );      
      boundary::set<boundary::LocalPressure>( sLattice, superGeometry, outfMat );
      break;
    
    case localAndConvection:
      boundary::set<boundary::LocalVelocity>( sLattice, inletIndicator );
      boundary::set<boundary::LocalVelocity>( sLattice, aroundIndicator );
      boundary::set<boundary::InterpolatedConvection>( sLattice, superGeometry, outfMat );
      break;

    case interpolated:
    case interpolatedAndZeroGrad:
      boundary::set<boundary::InterpolatedVelocity>( sLattice, inletIndicator );
      if ( boundarytype == interpolated ) boundary::set<boundary::LocalVelocity>( sLattice, aroundIndicator );
      if ( boundarytype == interpolatedAndZeroGrad ) setZeroGradientBoundary<T,DESCRIPTOR>( sLattice, aroundIndicator );
      boundary::set<boundary::InterpolatedPressure>( sLattice, superGeometry, outfMat );
      break;
      
    case dampingAndLocal:
    case dampingAndZeroGrad:
      domainLengths -= 2*converter.getConversionFactorLength();
      boundary::set<boundary::LocalVelocity>( sLattice, inletIndicator );
      boundary::set<boundary::LocalPressure>( sLattice, superGeometry, outfMat );
      if ( boundarytype == dampingAndLocal ) boundary::set<boundary::LocalVelocity>( sLattice, aroundIndicator );
      if ( boundarytype == dampingAndZeroGrad ) setZeroGradientBoundary<T,DESCRIPTOR>( sLattice, aroundIndicator );
    case damping:
      boundary::set<boundary::PerfectlyMatchedLayer>( sLattice, superGeometry, dampMat );
      break;
  }

  // Initial conditions
  AnalyticalConst<ndim,T,T> ux = AnalyticalConst<ndim,T,T>( Ma );
  AnalyticalConst<ndim,T,T> uy = AnalyticalConst<ndim,T,T>( 0. );
  AnalyticalConst<ndim,T,T> uz = AnalyticalConst<ndim,T,T>( 0. );
  AnalyticalComposed<ndim,T,T> u( ux, uy, uz );

  T ndatapoints = converter.getResolution();
  T dist = converter.getPhysDeltaX();
  bulkIndicator = superGeometry.getMaterialIndicator({1,2,3,4,5,6,7,8,9,10});  // got all indicators to define RhoU and fields globally
  AcousticPulse<3,T> densityProfile( rho0, amplitude, alpha );
  linePlot<ndim,T>( densityProfile, ndatapoints, dist, "shock_diag", "density [LU]", diagonal2d );
  linePlot<ndim,T>( densityProfile, ndatapoints, dist, "shock_hline", "density [LU]", horizontal );

  //Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, densityProfile, u );
  sLattice.iniEquilibrium( bulkIndicator, densityProfile, u );

  sLattice.setParameter<descriptors::OMEGA>( omega );

  if ( boundarytype == damping || boundarytype == dampingAndLocal ) {
    // Define fields; could also use setParameter as for omega, but then I could not use a flow profile
    sLattice.defineField<descriptors::UX>( bulkIndicator, ux );
    sLattice.defineField<descriptors::UY>( bulkIndicator, uy );
    sLattice.defineField<descriptors::UZ>( bulkIndicator, uz );
    AnalyticalConst<ndim,T,T> rhoField( rho0 );
    sLattice.defineField<descriptors::DENSITY>( bulkIndicator, rhoField );
    // output damping layer parameter
    DampingTerm<3,T,DESCRIPTOR> sigma_plot( converter, boundaryDepth, domainLengths );
    linePlot<ndim,T>( sigma_plot, ndatapoints, dist, "sigma_hline", "sigma", horizontal );
    linePlot<ndim,T>( sigma_plot, ndatapoints, dist, "sigma_diag", "sigma", diagonal2d );
    DampingTerm<3,T,DESCRIPTOR> sigma( converter, boundaryDepth, domainLengths, dampingStrength );
    sLattice.defineField<descriptors::DAMPING>( bulkIndicator, sigma );
  }
  
  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Sets fixed far field average values
void setFarFieldValues( SuperLattice<T,DESCRIPTOR>& sLattice,
                        SuperGeometry<T,3>& superGeometry,
                        T rho0,
                        T farFieldVelocity,
                        BoundaryType boundarytype )
{
  OstreamManager clout( std::cout,"setFarFieldValues" );
  AnalyticalConst<ndim,T,T> ux( farFieldVelocity );
  AnalyticalConst<ndim,T,T> uy( 0. );
  AnalyticalConst<ndim,T,T> uz( 0. );
  AnalyticalComposed<ndim,T,T> u( ux, uy, uz );
  AnalyticalConst<ndim,T,T> rho( rho0 );

  auto farFieldIndicator = superGeometry.getMaterialIndicator( {inflMat, outfMat, rounMat,checMat} );
  sLattice.defineRhoU( farFieldIndicator, rho, u );
  sLattice.iniEquilibrium( farFieldIndicator, rho, u );
}

T L2Norm(SuperLattice<T,DESCRIPTOR>& sLattice,
         SuperGeometry<T,ndim>& superGeometry,
         UnitConverter<T,DESCRIPTOR> const& converter,
         int fluidMaterial
         ) {
  OstreamManager clout(std::cout, "L2Norm");
  SuperLatticePhysPressure3D<T,DESCRIPTOR> pressurenD( sLattice, converter );
  AnalyticalConst<ndim,T,T> rho0nD( 0. );
  SuperAbsoluteErrorL2Norm3D<T> absPressureErrorL2Norm( pressurenD, rho0nD, superGeometry.getMaterialIndicator( fluidMaterial ) );
  T result[ndim];
  int tmp[] = {int()};
  absPressureErrorL2Norm( result, tmp );
  return result[0];
}

// write data to termimal and file system
void getResults(SuperLattice<T,DESCRIPTOR>& sLattice,
                UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                SuperGeometry<T,ndim>& superGeometry, util::Timer<T>& timer,
                T rho0, T amplitude,
                Gnuplot<T>& gplot_l2_abs, T Lp0,
                size_t iTvtk, size_t imax,
                size_t iTplot,
                BoundaryType boundarytype,
                int fluidMaterial
                )
{
#ifdef PLATFORM_GPU_CUDA
  gpu::cuda::device::synchronize();
#endif
  OstreamManager clout( std::cout,"getResults" );

  std::string name = "movingshock3d";

  if ( iT==0 ) {
    SuperVTMwriter3D<T> vtmWriter( name );
    // Writes geometry, cuboid no. and rank no. to file system
    SuperLatticeCuboid3D<T,DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T,DESCRIPTOR> rank( sLattice );
    SuperGeometryF<T,DESCRIPTOR::d> geometryF(superGeometry);
    vtmWriter.write(geometryF);
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );

  // Writes to plot 100 times
  if ( iT%iTplot==0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);  // important for synchronization (?) on GPU

    if ( iT%( iTplot*5 ) == 0 ) {      
      // write to terminal
      timer.update( iT );
      timer.printStep();

      // Lattice statistics console output
      sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
    }

    gplot_l2_abs.setData(T(iT), L2Norm( sLattice, superGeometry, converter, fluidMaterial ) / Lp0 );

    if ( iT%( iTplot*10 ) == 0 ) {
      std::stringstream ss;
      ss << std::setw(4) << std::setfill('0') << iT;
      T dist = converter.getPhysDeltaX();
      T ndatapoints = converter.getResolution(); // number of data points on line
      AnalyticalFfromSuperF3D<T> pressure_interpolation( pressure, true, true );
      T pmin(converter.getPhysPressure(-amplitude/50));
      T pmax(converter.getPhysPressure(+amplitude/50));
      linePlot<ndim,T>( pressure_interpolation, ndatapoints, dist, "pressure_hline_" + ss.str(), "pressure [PU]", horizontal, false, true, pmin, pmax );
      linePlot<ndim,T>( pressure_interpolation, ndatapoints, dist, "pressure_vline_" + ss.str(), "pressure [PU]", vertical, false, true, pmin, pmax );
      linePlot<ndim,T>( pressure_interpolation, ndatapoints, dist, "pressure_diagonal_" + ss.str(), "pressure [PU]", diagonal2d, false, true, pmin, pmax );
    }

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }

  if ( ( boundarytype == dampingAndLocal || boundarytype == local ) && ( iT%50 == 0 ) ) {
    setUAverage( sLattice, superGeometry, fluidMaterial );
    clout << "uAverage(" << iT << ")= " << *uAverage << std::endl;
  }

  // get VTK and images
  if ( iT%iTvtk==0 || sLattice.getStatistics().getAverageRho() > 2. ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);  // important for synchronization on GPU (or is sync above enough?)
    sLattice.scheduleBackgroundOutputVTK([&,name,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriter(name);
      SuperLatticePhysVelocity3D velocityF(sLattice, converter);
      SuperLatticePhysPressure3D pressureF(sLattice, converter);
      vtmWriter.addFunctor(velocityF);
      vtmWriter.addFunctor(pressureF);
      task(vtmWriter, iT);
    });

    // pressure image
    BlockReduction3D2D<T> pressureReduction( pressure, Vector<T,3>({0, 0, 1}) );
    heatmap::plotParam<T> jpeg_ParamP;
    jpeg_ParamP.maxValue = converter.getPhysPressure(+amplitude/200);
    jpeg_ParamP.minValue = converter.getPhysPressure(-amplitude/200);
    jpeg_ParamP.colour = "rainbow";
    jpeg_ParamP.fullScreenPlot = true;
    heatmap::write(pressureReduction, iT, jpeg_ParamP);

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }
}


int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  CLIreader args(argc, argv);

  // === get Command Line Arguments
  std::string outdir            = args.getValueOrFallback<std::string>( "--outdir",            "" );
  const int res                 = args.getValueOrFallback( "--res",               101 );
  const T rho0                  = args.getValueOrFallback( "--rho0",              1.  );
  const T Ma                    = args.getValueOrFallback( "--Ma",                0.1 );  // can be quite high
  const T charV                 = args.getValueOrFallback( "--charV",             1.  );
  T Re                          = args.getValueOrFallback( "--Re",                0.  );
  const T lx                    = args.getValueOrFallback( "--lx",                1.  );
  const T ly                    = args.getValueOrFallback( "--ly",                1.  );
  const T lz                    = args.getValueOrFallback( "--lz",                1.  );
  T tau                         = args.getValueOrFallback( "--tau",               0.  );
  T maxPhysT                    = args.getValueOrFallback( "--tmax",              0.  );
  size_t maxLatticeT            = args.getValueOrFallback( "--iTmax",             0   );
  size_t nout                   = args.getValueOrFallback( "--nout",              5   );  // minimum number of vtk outputs
  size_t iTout                  = args.getValueOrFallback( "--iTout",             0   );  // iterations for vtk outputs
  T tout                        = args.getValueOrFallback( "--tout",              0.  );  // timestep for vtk outputs
  size_t nplot                  = args.getValueOrFallback( "--nplot",             100 );  // minimum number of plot points // number of plots will be 1/10th
  const int boundaryCondition   = args.getValueOrFallback( "--boundaryCondition", 3   );
  const size_t boundaryDepth    = args.getValueOrFallback( "--boundaryDepth",     20  );
  const T amplitude             = args.getValueOrFallback( "--amplitude",         1e-3 );  // physical pressure amplitude
  const T dampingStrength       = args.getValueOrFallback( "--dampingStrength",   1.  );
  size_t overlap                = args.getValueOrFallback( "--overlap",           3   );
  const T bShock                = args.getValueOrFallback( "--bShock",            1./20.);
  const T aShock                = args.getValueOrFallback( "--aShock",            log(2.)/(bShock*bShock) );
  const int eternalscale        = args.getValueOrFallback( "--eternalscale",      2   );

  // managing outputs
  std::stringstream outdir_mod;
  if ( outdir == "" ) {
    outdir_mod << "./tmp";
    switch ( boundaryCondition ) {
      case 0: outdir_mod << "_eternal"; break;
      case 1: outdir_mod << "_periodic"; break;
      case 2: outdir_mod << "_local"; break;
      case 3: outdir_mod << "_damping"; break;
      case 4: outdir_mod << "_dampingAndLocal"; break;
      case 5: outdir_mod << "_dampingAndZeroGrad"; break;
      case 6: outdir_mod << "_localAndZeroGrad"; break;
      case 7: outdir_mod << "_localAndConvection"; break;
      case 8: outdir_mod << "_interpolated"; break;
      case 9: outdir_mod << "_interpolatedAndZeroGrad"; break;
    }
    if ( Ma != .1 ) outdir_mod << "_Ma" << Ma;
    if ( Re != 0 ) outdir_mod << "_Re" << Re;
    if ( amplitude != 1e-3 ) outdir_mod << "_a" << amplitude;
    outdir_mod  << "_" << lx << "x" << ly << "x" << lz << "_res" << res;
    if ( boundaryCondition == 0 ) outdir_mod << "x" << eternalscale;
    if ( boundaryCondition == 3 || boundaryCondition == 4 ) outdir_mod << "_bd" << boundaryDepth << "_bs" << dampingStrength;
  } else {
    outdir_mod << outdir;
  }
  singleton::directories().setOutputDir( outdir_mod.str()+"/" );  // set output directory
  // Create a TeeBuffer that writes to both the console and the file
  std::ofstream fileStream( outdir_mod.str() + "/output.txt" );
  DoubleBuffer doubleBuffer( std::cout.rdbuf(), fileStream.rdbuf() );
  std::streambuf* originalCoutBuffer = std::cout.rdbuf( &doubleBuffer );
  OstreamManager clout( std::cout, "main" );
  clout << std::endl << "outdir set to " << outdir_mod.str() << std::endl;
  clout << "b=" << bShock << "; alpha=" << aShock << std::endl;

  /* default values */
  BoundaryType boundarytype;
  switch ( boundaryCondition ) {
    case 0:   boundarytype = eternal;         clout << "Boundary condition is solved by just extending the domain to " << eternalscale << " times" << std::endl; break;
    case 1:   boundarytype = periodic;        clout << "Boundary condition type specified to periodic."                       << std::endl; break;
    case 2:   boundarytype = local;           clout << "Boundary condition type specified to local."                          << std::endl; break;
    case 3:   boundarytype = damping;         clout << "Boundary condition type specified to damping."                        << std::endl; break;
    case 4:   boundarytype = dampingAndLocal; clout << "Boundary condition type specified to damping with local outside bc."  << std::endl; break;
    case 5:   boundarytype = dampingAndZeroGrad; clout << "Boundary condition type specified to damping with zero gradient outside bc." << std::endl; break;
    case 6:   boundarytype = localAndZeroGrad; clout << "Boundary condition type specified to local with zero gradient outside bc." << std::endl; break;
    case 7:   boundarytype = localAndConvection; clout << "Boundary condition type specified to local with interpolated convection outflow." << std::endl; break;
    case 8:   boundarytype = interpolated;    clout << "Boundary condition type specified to interpolated." << std::endl; break;
    case 9:   boundarytype = interpolatedAndZeroGrad; clout << "Boundary condition type specified to interpolated with zero gradient outside bc." << std::endl; break;
    default:  boundarytype = periodic;        clout << "Boundary condition type specified to periodic."                       << std::endl; break;
  }

  // determining Reynolds regime (incl. viscosity and relaxation time)
  T farFieldVelocity = Ma;
  const T charL             = lx;
  const T cs_LU             = 1 / std::sqrt( 3.0 );
  T viscosity;
  if ( tau == 0 ) {
    if ( Re == 0 ) {
      viscosity             = 1.48e-5;
      Re                    = charV * charL / viscosity;
    } else {
      viscosity             = charV * charL / Re;
    }
    tau                     = viscosity / (cs_LU*cs_LU) + 0.5;
  }
  else {
    if ( Re == 0 ) Re       = 20000;
    viscosity               = charV * charL / Re;
  }

  clout << std::endl << "Input to unitconverter: res=" << res << ", tau=" << tau
  << ", charL=" << lx << ", charV=" << charV << ", farFieldVelocity=" << farFieldVelocity << ", viscosity=" << viscosity
        << ", physDensity=1.204" << std::endl;
  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> converter(
    (size_t) res,         // resolution
    (T)      tau,         // relaxation time
    (T)      charL,       // charPhysLength: reference length of simulation geometry
    (T)      charV,       // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)      viscosity,   // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)      rho0         // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write( "movingshock3d" );
  
  // === setup domain
  // changing domain size depending on how much area is lost to the damping boundary layer
  T boundaryDepth_pu = converter.getPhysLength( boundaryDepth );
  T lengthDomain, heightDomain, depthDomain;
  switch ( boundarytype ) {
    case eternal:         // extend the domain all around
      lengthDomain  = eternalscale*lx;
      heightDomain  = eternalscale*ly;
      depthDomain   = eternalscale*lz;
      break;
    case periodic:        // reference domain size
      lengthDomain  = lx;
      heightDomain  = ly;
      depthDomain   = lz;
      break;
    case local:           // add one layer in each direction
    case localAndZeroGrad:
    case localAndConvection:
    case interpolated:
    case interpolatedAndZeroGrad:
      lengthDomain  = lx + 2*converter.getConversionFactorLength();
      heightDomain  = ly + 2*converter.getConversionFactorLength();
      depthDomain   = lz + 2*converter.getConversionFactorLength();
      break;
    case damping:         // add boundary layer
      lengthDomain  = lx + 2*boundaryDepth_pu;
      heightDomain  = ly + 2*boundaryDepth_pu;
      depthDomain   = lz + 2*boundaryDepth_pu;
      break;
    case dampingAndLocal: // add one layer in flow direction and boundary layer all round
    case dampingAndZeroGrad:
      lengthDomain  = lx + 2*boundaryDepth_pu + 2*converter.getConversionFactorLength();
      heightDomain  = ly + 2*boundaryDepth_pu + 2*converter.getConversionFactorLength();
      depthDomain   = lz + 2*boundaryDepth_pu + 2*converter.getConversionFactorLength();
      break;
  }
  clout << "Fluid Domain = " << lx << "x" << ly << "x" << lz << std::endl;
  clout << "Simulation Domain = " << lengthDomain << "x" << heightDomain << "x" << depthDomain << std::endl;
  Vector<T,ndim> domainLengths = { lengthDomain, heightDomain, depthDomain };
  Vector<T,ndim> originDomain( -lengthDomain/2, -heightDomain/2, -depthDomain/2 );  //
  Vector<T,ndim> extendDomain( lengthDomain, heightDomain, depthDomain );  // size of the domain
  IndicatorCuboid3D<T> domain( extendDomain, originDomain );
  Vector<T,3> fluidExtend = { lx, ly, lz };
  Vector<T,3> fluidOrigin = { -lx/2, -ly/2, -lz/2 };
  IndicatorCuboid3D<T> fluidDomain( fluidExtend, fluidOrigin );
  // Instantiation of a cuboidGeometry with weights
  #ifdef PARALLEL_MODE_MPI
    const int noOfCuboids = singleton::mpi().getSize();
  #else
    const int noOfCuboids = 6;
  #endif
  CuboidGeometry<T,ndim> cuboidGeometry( domain, converter.getConversionFactorLength(), noOfCuboids );
  switch ( boundarytype ) {
    case eternal:
    case periodic:
    case damping:
      cuboidGeometry.setPeriodicity({true, true, true});    break;
    case local:
    case localAndZeroGrad:
    case localAndConvection:
    case interpolated:
    case interpolatedAndZeroGrad:
    case dampingAndLocal:
    case dampingAndZeroGrad:
      cuboidGeometry.setPeriodicity({false, false, false}); break;
  }
  
  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  if ( noOfCuboids > 1 && ( boundarytype == damping || boundarytype == dampingAndLocal || boundarytype == dampingAndZeroGrad ) ) {
    clout << "overlap=" << overlap << "; boundaryDepth_LU=" << boundaryDepth << "; setting overlap to >=boundaryDepth." << std::endl;
    overlap = std::max( overlap, boundaryDepth );
  }
  SuperGeometry<T,ndim> superGeometry( cuboidGeometry, loadBalancer, overlap );

  clout << std:: endl << "Domain setup: boundaryDepthLU=" << boundaryDepth << "; boundaryDepthPU=" << boundaryDepth_pu << "; overlapLU=" << superGeometry.getOverlap() << std::endl;
  prepareGeometry( converter, superGeometry, domain, fluidDomain, boundarytype, res, boundaryDepth );
  int fluidMaterial = fluiMat;
  switch ( boundarytype ) {  // for correct normals, fluid material number must be 1
    case periodic:
    case local:
    case localAndConvection:
    case localAndZeroGrad:
    case interpolated:
    case interpolatedAndZeroGrad:
      fluidMaterial = dampMat; break;
    default: break;
  }
  clout << "Actual fluid material number is " << fluidMaterial << "; number of fluid voxels: " << superGeometry.getStatistics().getNvoxel( fluidMaterial ) << std::endl;
  
  SuperLattice<T,DESCRIPTOR> sLattice( superGeometry );

  //prepare Lattice and set boundaryConditions
  prepareLattice( converter, sLattice, superGeometry, rho0, Ma, amplitude, aShock, boundarytype, boundaryDepth, domainLengths, dampingStrength );

  // Initialize pressure L2 norm plot
  Gnuplot<T> gplot_l2_abs("l2_absolute");
  gplot_l2_abs.setLabel("time []", "absolute L2 norm []");
  T Lp0 = L2Norm( sLattice, superGeometry, converter, fluidMaterial );

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

  // === 3a-rd Step: calculate iterations from input ===
  // maxLatticeT depends on maximum physical time. If maxLatticeT is provided in command line, it is an upper bound
  if ( maxPhysT == 0. ) maxPhysT = .015/Ma;  // expecting cs = 343 m/s, so this should be enough for 5 revolutions through the domain
  if ( maxLatticeT == 0 ) maxLatticeT = converter.getLatticeTime( maxPhysT );
  else maxLatticeT = std::min( maxLatticeT, converter.getLatticeTime( maxPhysT ) );
  // nout is the minimum number of vtk outputs --> take max between nout and nout derived from iTout or tout
  size_t nout_from_iTout = 0, nout_from_tout = 0;
  if ( iTout != 0 ) { nout_from_iTout = size_t( maxLatticeT / iTout ); nout = std::max( nout, nout_from_iTout ); }
  if ( tout != 0 ) { nout_from_tout = size_t( maxLatticeT / tout ); nout = std::max( nout, nout_from_tout ); }
  size_t iTvtk    = std::max(int( maxLatticeT / nout ), 1);
  size_t iTplot   = std::max(int( maxLatticeT / nplot ), 1);

  clout << std::endl << "Timing setup:" << std::endl
        << "maxLatticeT=" << maxLatticeT << std::endl
        << "maxPhysT=" << maxPhysT << std::endl
        << "dt=" << converter.getPhysDeltaT() << std::endl
        << "iTout=" << iTout << std::endl
        << "tout=" << tout << std::endl
        << "iTvtk=" << iTvtk << std::endl
        << "iTplot=" << iTplot << std::endl;

  // === 4th Step: Main Loop with Timer ===
  clout << std::endl << "starting simulation..." << std::endl;
  util::Timer<T> timer( maxLatticeT, superGeometry.getStatistics().getNvoxel() );
  timer.start();

  size_t iT = 0;
  while ( iT < maxLatticeT ) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    if ( boundarytype == local || boundarytype == dampingAndLocal || boundarytype == dampingAndZeroGrad ) {
      setFarFieldValues( sLattice, superGeometry, rho0, farFieldVelocity, boundarytype );
    }
    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer,
                rho0, amplitude,
                gplot_l2_abs, Lp0,
                iTvtk, maxLatticeT, iTplot,
                boundarytype, fluidMaterial );
    iT++;
  }
  clout << "Simulation stopped after " << iT << "/" << maxLatticeT << " iterations." << std::endl;

  // output the development of pressure L2 norm over time
  gplot_l2_abs.setYrange(1e-3, 1);
  gplot_l2_abs.setLogScale(2);
  gplot_l2_abs.writePNG(-1, -1, "gplot_l2_abs");

  timer.stop();
  timer.printSummary();
  std::cout.rdbuf(originalCoutBuffer);
}
