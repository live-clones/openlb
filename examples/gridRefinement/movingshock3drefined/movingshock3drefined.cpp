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
 * Method from H. Xu and P. Sagaut, “Analxsis of the absorbing layers for the weaklx-compressible lattice Boltzmann methods,” Journal of Computational Physics, vol. 245, pp. 14–42, Jul. 2013, doi: 10.1016/j.jcp.2013.02.051.
 * Benchmark from C. Bogey and C. Baillx, “Three-dimensional non-reﬂective boundary conditions for acoustic simulations: far ﬁeld formulation and validation test cases,” vol. 88, 2002.
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
// using DESCRIPTOR = D3Q19<>;
using DESCRIPTOR = D3Q27<>;
using BulkDynamics = BGKdynamics<T,DESCRIPTOR>;
// using BulkDynamics = KBCdynamics<T,DESCRIPTOR>;
// using BulkDynamics = olb::dynamics::Tuple<T,DESCRIPTOR,momenta::BulkTuple,equilibria::Incompressible,collision::TRT>;

const int noMat   = 0;
const int fluiMat = 1;
const int checMat = 2;
const int farfMat = 3;
const int ndim = 3;  // a few things (e.g. SuperSum3D cannot be adapted to 2D, but this should help speed it up)

typedef enum {eternal, periodic, stretchedLayer, dampingAndStretched} BoundaryType;

T* uAverage = NULL;

// Stores geometry information in form of material numbers
void prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperGeometry<T,ndim>& sGeometry,
                     IndicatorF3D<T>& domain,
                     IndicatorF3D<T>& fluidDomain,
                     BoundaryType boundarytype,
                     int res,
                     int boundaryDepth
                     )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << std::endl << "Prepare Geometry ..." << std::endl;

  sGeometry.rename( noMat, checMat );   // all nodes to temporary type  
  sGeometry.rename( checMat, fluiMat, fluidDomain );
  sGeometry.rename( checMat, farfMat );

  sGeometry.communicate();
  // sGeometry.clean() deletes some nodes although it shouldn't - this is probablx to do with the assumption that 1 is fluid and onlx needs one layer next to it (rather than 20 for pml)
  sGeometry.checkForErrors();
  sGeometry.updateStatistics();
  clout << "Materials:" << std::endl
        << noMat   << " = no Material (default, should be immediately renamed to Check)" << std::endl
        << checMat << " = check (should be empty, but corners may have issues with the normal; check is therefore treated as inflow" << std::endl
        << fluiMat << " = fluid (sometimes must be set to 1 to work with normals of boundary conditions)" << std::endl
        << farfMat << " = far field (must be 1 so outside boundaries can calculate their normal to the damping area)" << std::endl;
  sGeometry.getStatistics().print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// write data to termimal and file system
void getResults(SuperLattice<T,DESCRIPTOR>& sLattice,
                UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                SuperGeometry<T,ndim>& sGeometry, util::Timer<T>& timer,
                T amplitude,
                Gnuplot<T>& gplot_l2_abs, T Lp0,
                size_t iTvtk,
                size_t iTplot,
                std::string layername="L0",
                bool finestLayer=false
                )
{
#ifdef PLATFORM_GPU_CUDA
  gpu::cuda::device::synchronize();
#endif
  OstreamManager clout( std::cout,"getResults" );

  std::string name = "movingshock3d" + layername;

  if ( iT==0 ) {
    SuperVTMwriter3D<T> vtmWriter( name, 0 );
    // Writes geometry, cuboid no. and rank no. to file system
    SuperLatticeCuboid3D<T,DESCRIPTOR> cuboid( sLattice );
    cuboid.getName() += layername;
    vtmWriter.write( cuboid );
    SuperLatticeRank3D<T,DESCRIPTOR> rank( sLattice );
    rank.getName() += layername;
    vtmWriter.write( rank );
    SuperGeometryF<T,DESCRIPTOR::d> geometryF(sGeometry);
    geometryF.getName() += layername;
    vtmWriter.write(geometryF);
    SuperLatticePhysField3D<T,DESCRIPTOR,DAMPING> damping( sLattice, 1. );
    damping.getName() = "dampingField" + layername;
    vtmWriter.write( damping );
    vtmWriter.createMasterFile();
  }
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );

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
    jpeg_ParamP.name += layername;
    heatmap::write(pressureReduction, iT, jpeg_ParamP);

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }

  // Writes to plot 100 times
  if ( finestLayer && iT%iTplot==0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);  // important for synchronization (?) on GPU

    if ( iT%( iTplot*5 ) == 0 ) {      
      // write to terminal
      timer.update( iT );
      timer.printStep();

      // Lattice statistics console output
      sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
    }

    gplot_l2_abs.setData(T(iT), L2Norm<ndim,T,DESCRIPTOR>( sLattice, sGeometry, converter, fluiMat ) / Lp0 );

    if ( iT%( iTplot*10 ) == 0 ) {
      std::stringstream ss;
      ss << std::setw(4) << std::setfill('0') << iT;
      T dist = converter.getPhysDeltaX();
      T ndatapoints = converter.getResolution(); // number of data points on line
      AnalyticalFfromSuperF3D<T> pressure_interpolation( pressure, true, true );
      T pmin(converter.getPhysPressure(-amplitude/50));
      T pmax(converter.getPhysPressure(+amplitude/50));
      linePlot<ndim,T>( pressure_interpolation, ndatapoints, dist, "pressure_hline_" + ss.str(), "pressure [PU]" + layername, horizontal, false, true, pmin, pmax );
      linePlot<ndim,T>( pressure_interpolation, ndatapoints, dist, "pressure_vline_" + ss.str(), "pressure [PU]" + layername, vertical, false, true, pmin, pmax );
      linePlot<ndim,T>( pressure_interpolation, ndatapoints, dist, "pressure_diagonal_" + ss.str(), "pressure [PU]" + layername, diagonal2d, false, true, pmin, pmax );
    }

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
  size_t res                    = args.getValueOrFallback( "--res",               101 );
  const T rho0                  = args.getValueOrFallback( "--rho0",              1.  );
  const T Ma                    = args.getValueOrFallback( "--Ma",                0.1 );  // can be quite high
  const T charV                 = args.getValueOrFallback( "--charV",             1.  );
  T Re                          = args.getValueOrFallback( "--Re",                0.  );
  const T lx                    = args.getValueOrFallback( "--lx",                1.  );
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
      case 2: outdir_mod << "_stretchedLayer"; break;
      case 3: outdir_mod << "_dampingAndStretched"; break;
      default: throw std::range_error("Invalid boundary type input!"); break;
    }
    if ( Ma != .1 ) outdir_mod << "_Ma" << Ma;
    if ( Re != 0 ) outdir_mod << "_Re" << Re;
    if ( amplitude != 1e-3 ) outdir_mod << "_a" << amplitude;
    outdir_mod  << "_" << lx << "^3" << "_res" << res;
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
    case 0:   boundarytype = eternal;         clout << "Boundary condition is solved by just extending the domain to "   << eternalscale << " times" << std::endl; break;
    case 2:   boundarytype = stretchedLayer;  clout << "Boundary condition speficied to periodic with 2 stretched layers outside." << std::endl; break;
    case 3:   boundarytype = dampingAndStretched; clout << "Boundary condition specified to damping with 1 additional stretched layer." << std::endl; break;
    case 1:   
    default:  boundarytype = periodic;        clout << "Boundary condition specified to periodic."                       << std::endl; break;
  }

  // determining Reynolds regime (incl. viscosity and relaxation time)
  const T charL             = lx;
  const T cs_LU             = 1 / std::sqrt( 3.0 );
  T refDeg = 1.;  // refinement degree
  if ( boundarytype == dampingAndStretched ) refDeg = 2.;
  else if ( boundarytype == stretchedLayer ) refDeg = 4.;
  T viscosity;
  viscosity   = 1.48e-5;
  viscosity  *= refDeg;
  res        /= refDeg;
  Re          = charV * charL / viscosity;
  T tau       = viscosity / (cs_LU*cs_LU) + 0.5;
  clout << std::endl << "Input to unitconverter: res=" << res << ", tau=" << tau
        << ", charL=" << charL << ", charV=" << charV
        << ", viscosity=" << viscosity
        << ", physDensity=1.0" << std::endl;
  UnitConverter<T, DESCRIPTOR> converterL0( (charL/res),  // dx
                                            (tau - 0.5) / descriptors::invCs2<T,DESCRIPTOR>() * util::pow(charL/res,2) / viscosity,  // dt
                                            charL, charV, viscosity, rho0, 0.);
  // Prints the converterL0 log as console output
  converterL0.print();
  // Writes the converterL0 log in a file
  converterL0.write( "movingshock3d" );
  
  // === setup domain
  // changing domain size depending on how much area is lost to the damping boundary layer
  T bdPU = converterL0.getPhysLength( boundaryDepth );
  T lengthDomain = lx;
  T scaling = 1.6;
  switch ( boundarytype ) {
    case eternal:             lengthDomain = eternalscale*lx;     break;  // extend the domain all around
    case stretchedLayer:      lengthDomain = lx*scaling*scaling;  break;  // twice add 12.5 % of domain extend in each direction
    case dampingAndStretched: lengthDomain = (lx+2*bdPU)*scaling; break;  // add boundary layer, then add 12.5 % of domain extend
    case periodic: default: break;  // reference domain size
  }
  clout << "Fluid Domain = " << lx << "^" << ndim << std::endl;
  clout << "Simulation Domain = " << lengthDomain << "^" << ndim << std::endl;
  Vector<T,ndim> domainLengths = { lengthDomain };
  Vector<T,ndim> originDomain( 0. );  // -lengthDomain/2.
  Vector<T,ndim> extendDomain( lengthDomain );  // size of the domain
  IndicatorCuboid3D<T> domain( extendDomain, originDomain );
  Vector<T,3> fluidExtend( lx );
  Vector<T,3> fluidOrigin( (lengthDomain-lx)*.5 ); // -lx/2.
  IndicatorCuboid3D<T> fluidDomain( fluidExtend, fluidOrigin );
  // Instantiation of a cGeometryL0 with weights
  CuboidGeometry<T,ndim> cGeometryL0( domain, converterL0.getConversionFactorLength() );
  cGeometryL0.setPeriodicity({true, true, true});
  auto cGeometryL1 = cGeometryL0;
  auto cGeometryL2 = cGeometryL0;
  Vector<T,3> extendRefinement, originRefinement;
  if ( boundarytype == stretchedLayer) {
    // cGeometryL0.splitFractional(0, 0, {0.1,0.8,0.1});  // 0,1,2
    // cGeometryL0.splitFractional(1, 1, {0.1,0.8,0.1});  // 2 back to 1; split 1 into 2,3,4
    // cGeometryL0.splitFractional(3, 2, {0.1,0.8,0.1});  // 4 back to 3; split 3 into 3,4,5
    // cGeometryL0.splitFractional(5, 0, {0.1,0.8,0.1});
    // cGeometryL0.splitFractional(7, 1, {0.1,0.8,0.1});
    // cGeometryL0.splitFractional(9, 2, {0.1,0.8,0.1});
    clout << "fluidExtend=" << fluidExtend << ",fluidOrigin" << fluidOrigin << ",lengthDomain" << lengthDomain << ",lx" << lx << std::endl;
    extendRefinement = {lx*scaling, lx*scaling, lx*scaling};
    originRefinement = {(lengthDomain-lx*scaling)/2.,(lengthDomain-lx*scaling)/2.,(lengthDomain-lx*scaling)/2.};  // {-scaling/2., -scaling/2., -scaling/2.};
    IndicatorCuboid3D<T> toBeRefinedI(extendRefinement, originRefinement);
    cGeometryL1 = cGeometryL0;
    cGeometryL1.remove(toBeRefinedI);
    cGeometryL1.shrink(toBeRefinedI);
    cGeometryL1.refine(2);
    cGeometryL1.print();
    IndicatorCuboid3D<T> toBeRefinedII(fluidExtend, fluidOrigin);
    cGeometryL2 = cGeometryL1;
    cGeometryL2.remove(toBeRefinedII);
    cGeometryL2.shrink(toBeRefinedII);
    cGeometryL2.refine(2);
    cGeometryL2.print();
    // Adjust weights for balancing
    for (int iC=0; iC < cGeometryL0.size(); ++iC) {
      auto& cuboid = cGeometryL0.get(iC);
      cuboid.setWeight(cuboid.getLatticeVolume());
      auto origin = cGeometryL0.get(iC).getOrigin();
      for (int jC=0; jC < cGeometryL1.size(); ++jC) {
        if (cGeometryL1.get(jC).getOrigin() == origin) {
          cuboid.setWeight(cuboid.getWeight() + 2*cGeometryL1.get(jC).getLatticeVolume());
        }
      }
      for (int jC=0; jC < cGeometryL2.size(); ++jC) {
        if (cGeometryL2.get(jC).getOrigin() == origin) {
          cuboid.setWeight(cuboid.getWeight() + 4*cGeometryL2.get(jC).getLatticeVolume());
        }
      }
    }
  } else if ( boundarytype == dampingAndStretched ) {
    // cGeometryL0.splitFractional(0, 0, {0.1,0.8,0.1});  // 0,1,2
    // cGeometryL0.splitFractional(1, 1, {0.1,0.8,0.1});  // 2 back to 1; split 1 into 2,3,4
    // cGeometryL0.splitFractional(3, 2, {0.1,0.8,0.1});  // 4 back to 3; split 3 into 3,4,5
    extendRefinement = {lengthDomain/scaling, lengthDomain/scaling, lengthDomain/scaling};
    originRefinement = {(1.-scaling)*lx,(1.-scaling)*lx,(1.-scaling)*lx};  // {lengthDomain*((scaling-1.)/2.), lengthDomain*((scaling-1.)/2.), lengthDomain*((scaling-1.)/2.)};
    IndicatorCuboid3D<T> toBeRefinedI(extendRefinement, originRefinement);
    cGeometryL1 = cGeometryL0;
    cGeometryL1.remove(toBeRefinedI);
    cGeometryL1.shrink(toBeRefinedI);
    cGeometryL1.refine(2);
    cGeometryL1.print();
    // Adjust weights for balancing
    for (int iC=0; iC < cGeometryL0.size(); ++iC) {
      auto& cuboid = cGeometryL0.get(iC);
      cuboid.setWeight(cuboid.getLatticeVolume());
      auto origin = cGeometryL0.get(iC).getOrigin();
      for (int jC=0; jC < cGeometryL1.size(); ++jC) {
        if (cGeometryL1.get(jC).getOrigin() == origin) {
          cuboid.setWeight(cuboid.getWeight() + 2*cGeometryL1.get(jC).getLatticeVolume());
        }
      }
    }
  }
  
  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cGeometryL0 );

  // Instantiation of a sGeometry
  if ( boundarytype == dampingAndStretched ) {
    clout << "overlap=" << overlap << "; boundaryDepth_LU=" << boundaryDepth << "; setting overlap to >=boundaryDepth." << std::endl;
    overlap = std::max( overlap, boundaryDepth );
  }
  SuperGeometry<T,ndim> sGeometryL0( cGeometryL0, loadBalancer, overlap );

  clout << std:: endl << "Domain setup: boundaryDepthLU=" << boundaryDepth << "; boundaryDepthPU=" << bdPU << "; overlapLU=" << sGeometryL0.getOverlap() << std::endl;
  prepareGeometry( converterL0, sGeometryL0, domain, fluidDomain, boundarytype, res, boundaryDepth );
  
  RefinedLoadBalancer<T,ndim> loadBalancerL1(cGeometryL0, loadBalancer, cGeometryL1);

  refDeg = 1.;  // refinement degree
  if ( boundarytype == stretchedLayer ) refDeg = 2.;
  viscosity /= 2;
  res *= 2;
  // UnitConverter<T, DESCRIPTOR> converterL1 ( (charL/res),  // dx
  //                                            (tau - 0.5) / descriptors::invCs2<T,DESCRIPTOR>() * util::pow(charL/res,2) / viscosity,  // dt
  //                                            charL, charV, viscosity, rho0, 0.);
  UnitConverter<T, DESCRIPTOR> converterL1 = convectivelyRefineUnitConverter(converterL0, 2);
  converterL1.print();
  SuperGeometry<T,ndim> sGeometryL1(cGeometryL1, loadBalancerL1);
  prepareGeometry(converterL1, sGeometryL1, domain, fluidDomain, boundarytype, res, boundaryDepth);

  RefinedLoadBalancer<T,ndim> loadBalancerL2(cGeometryL1, loadBalancerL1, cGeometryL2);
  refDeg = 1.;  // refinement degree
  viscosity /= 2;
  res *= 2;
  // UnitConverter<T, DESCRIPTOR> converterL2( (charL/res),  // dx
  //                                           (tau - 0.5) / descriptors::invCs2<T,DESCRIPTOR>() * util::pow(charL/res,2) / viscosity,  // dt
  //                                           charL, charV, viscosity, rho0, 0.);
  UnitConverter<T, DESCRIPTOR> converterL2 = convectivelyRefineUnitConverter(converterL1, 2);
  converterL2.print();
  SuperGeometry<T,ndim> sGeometryL2(cGeometryL2, loadBalancerL2);
  prepareGeometry(converterL2, sGeometryL2, domain, fluidDomain, boundarytype, res, boundaryDepth);

  SuperVTMwriter<T,ndim> vtmWriterGeometry( "geometry_init", 0 );
  SuperGeometryF<T,DESCRIPTOR::d> geometryFL0(sGeometryL0);
  geometryFL0.getName() += "L0";
  vtmWriterGeometry.write( geometryFL0 );
  SuperGeometryF<T,DESCRIPTOR::d> geometryFL1(sGeometryL1);
  geometryFL1.getName() += "L1";
  vtmWriterGeometry.write( geometryFL1 );
  SuperGeometryF<T,DESCRIPTOR::d> geometryFL2(sGeometryL2);
  geometryFL2.getName() += "L2";
  vtmWriterGeometry.write( geometryFL2 );

  // Initial conditions
  Vector<T,ndim> x0 = sGeometryL0.getStatistics().getCenterPhysR(fluiMat);
  AcousticPulse<ndim,T> densityProfile( rho0, amplitude, aShock, x0 );
  AnalyticalConst<ndim,T,T> ux = AnalyticalConst<ndim,T,T>( Ma );
  AnalyticalConst<ndim,T,T> uy = AnalyticalConst<ndim,T,T>( 0. );
  AnalyticalConst<ndim,T,T> uz = AnalyticalConst<ndim,T,T>( 0. );
  AnalyticalComposed<ndim,T,T> u( ux, uy, uz );

  // === prepare Lattice and set boundaryConditions
  SuperLattice<T,DESCRIPTOR> sLatticeL0( cGeometryL0, loadBalancer, 3, converterL0 );
  std::unique_ptr<SuperIndicatorF<T,ndim>> bulkIndicator = sGeometryL0.getMaterialIndicator( {farfMat, fluiMat, checMat} );
  sLatticeL0.defineDynamics<BulkDynamics>( bulkIndicator );
  sLatticeL0.defineRhoU( bulkIndicator, densityProfile, u );
  sLatticeL0.iniEquilibrium( bulkIndicator, densityProfile, u );
  sLatticeL0.setParameter<descriptors::OMEGA>( converterL0.getLatticeRelaxationFrequency() );
  sLatticeL0.initialize();

  SuperLattice<T,DESCRIPTOR> sLatticeL1(cGeometryL1, loadBalancerL1, 3, converterL1);
  std::unique_ptr<SuperIndicatorF<T,ndim>> bulkIndicatorL1 = sGeometryL1.getMaterialIndicator( {farfMat, fluiMat, checMat} );
  sLatticeL1.defineDynamics<BulkDynamics>( bulkIndicatorL1 );
  if ( boundarytype == dampingAndStretched ) boundary::set<boundary::SpongeLayer>( sLatticeL1, sGeometryL2, farfMat );
  sLatticeL1.defineRhoU( bulkIndicatorL1, densityProfile, u );
  sLatticeL1.iniEquilibrium( bulkIndicatorL1, densityProfile, u );
  sLatticeL1.setParameter<descriptors::OMEGA>( converterL1.getLatticeRelaxationFrequency() );
  sLatticeL1.initialize();

  SuperLattice<T,DESCRIPTOR> sLatticeL2(cGeometryL2, loadBalancerL2, 3, converterL2);
  std::unique_ptr<SuperIndicatorF<T,ndim>> bulkIndicatorL2 = sGeometryL2.getMaterialIndicator( {farfMat, fluiMat, checMat} );
  sLatticeL2.defineDynamics<BulkDynamics>( bulkIndicatorL2 );
  sLatticeL2.defineRhoU( bulkIndicatorL2, densityProfile, u );
  sLatticeL2.iniEquilibrium( bulkIndicatorL2, densityProfile, u );
  sLatticeL2.setParameter<descriptors::OMEGA>( converterL2.getLatticeRelaxationFrequency() );
  sLatticeL2.initialize();

  if ( boundarytype == dampingAndStretched ) {
    // Define fields; could also use setParameter as for omega, but then I could not use a flow profile
    sLatticeL1.defineField<descriptors::UX>( bulkIndicatorL1, ux );
    sLatticeL1.defineField<descriptors::UY>( bulkIndicatorL1, uy );
    sLatticeL1.defineField<descriptors::UZ>( bulkIndicatorL1, uz );
    AnalyticalConst<ndim,T,T> rhoField( rho0 );
    sLatticeL1.defineField<descriptors::DENSITY>( bulkIndicatorL1, rhoField );
    // output damping layer parameter
    T ndatapoints = converterL0.getResolution();
    T dist = converterL0.getPhysDeltaX();
    DampingTerm<ndim,T,DESCRIPTOR> sigma_plot( converterL0, boundaryDepth, domainLengths );
    linePlot<ndim,T>( sigma_plot, ndatapoints, dist, "sigma_hline", "sigma", horizontal );
    linePlot<ndim,T>( sigma_plot, ndatapoints, dist, "sigma_diag", "sigma", diagonal2d );
    DampingTerm<ndim,T,DESCRIPTOR> sigma( converterL0, boundaryDepth, domainLengths, dampingStrength );
    sLatticeL1.defineField<descriptors::DAMPING>( bulkIndicatorL1, sigma );
  }

  auto coarseToFineL1 = refinement::lagrava::makeCoarseToFineCoupler( sLatticeL0, sGeometryL0, sLatticeL1, sGeometryL1 );
  auto fineToCoarseL1 = refinement::lagrava::makeFineToCoarseCoupler( sLatticeL0, sGeometryL0, sLatticeL1, sGeometryL1 );
  auto coarseToFineL2 = refinement::lagrava::makeCoarseToFineCoupler( sLatticeL1, sGeometryL1, sLatticeL2, sGeometryL2 );
  auto fineToCoarseL2 = refinement::lagrava::makeFineToCoarseCoupler( sLatticeL1, sGeometryL1, sLatticeL2, sGeometryL2 );
  if ( boundarytype == dampingAndStretched ) {
    coarseToFineL1->initialize_prev(meta::id<refinement::lagrava::FullTimeCoarseToFineO>{});
  } else if ( boundarytype == stretchedLayer ) {
    coarseToFineL1->initialize_prev(meta::id<refinement::lagrava::FullTimeCoarseToFineO>{});
    coarseToFineL2->initialize_prev(meta::id<refinement::lagrava::FullTimeCoarseToFineO>{});
  }

  clout << "Number of fluid voxels: " << sGeometryL0.getStatistics().getNvoxel( fluiMat ) << std::endl;
  // Initialize pressure L2 norm plot
  Gnuplot<T> gplot_l2_abs("l2_absolute");
  gplot_l2_abs.setLabel("time []", "absolute L2 norm []");
  T Lp0;
  if ( boundarytype == stretchedLayer ) Lp0 = L2Norm<ndim,T,DESCRIPTOR>( sLatticeL2, sGeometryL2, converterL2, fluiMat );
  else if ( boundarytype == dampingAndStretched ) Lp0 = L2Norm<ndim,T,DESCRIPTOR>( sLatticeL1, sGeometryL1, converterL1, fluiMat );
  else Lp0 = L2Norm<ndim,T,DESCRIPTOR>( sLatticeL0, sGeometryL0, converterL0, fluiMat );

  // === 3a-rd Step: calculate iterations from input ===
  // maxLatticeT depends on maximum physical time. If maxLatticeT is provided in command line, it is an upper bound
  if ( maxPhysT == 0. ) maxPhysT = .015/Ma;  // expecting cs = 343 m/s, so this should be enough for 5 revolutions through the domain
  if ( maxLatticeT == 0. ) maxLatticeT = converterL0.getLatticeTime( maxPhysT );
  else maxLatticeT = std::min( maxLatticeT, converterL0.getLatticeTime( maxPhysT ) );
  // nout is the minimum number of vtk outputs --> take max between nout and nout derived from iTout or tout
  size_t nout_from_iTout = 0, nout_from_tout = 0;
  if ( iTout != 0 ) { nout_from_iTout = size_t( maxLatticeT / iTout ); nout = std::max( nout, nout_from_iTout ); }
  if ( tout != 0 ) { nout_from_tout = size_t( maxLatticeT / tout ); nout = std::max( nout, nout_from_tout ); }
  size_t iTvtk    = std::max(int( maxLatticeT / nout ), 1);
  size_t iTplot   = std::max(int( maxLatticeT / nplot ), 1);

  clout << std::endl << "Timing setup:" << std::endl
        << "maxLatticeT=" << maxLatticeT << std::endl
        << "maxPhysT=" << maxPhysT << std::endl
        << "dt=" << converterL0.getPhysDeltaT() << std::endl
        << "iTout=" << iTout << std::endl
        << "tout=" << tout << std::endl
        << "iTvtk=" << iTvtk << std::endl
        << "iTplot=" << iTplot << std::endl;

  // === 4th Step: Main Loop with Timer ===
  clout << std::endl << "starting simulation..." << std::endl;
  size_t Nvoxel = sGeometryL0.getStatistics().getNvoxel();
  switch ( boundarytype ) {
    case dampingAndStretched: Nvoxel = sGeometryL0.getStatistics().getNvoxel() + 2*sGeometryL1.getStatistics().getNvoxel();
    case stretchedLayer:      Nvoxel = sGeometryL0.getStatistics().getNvoxel() + 2*sGeometryL1.getStatistics().getNvoxel() + 4*sGeometryL2.getStatistics().getNvoxel();
    default: Nvoxel = sGeometryL0.getStatistics().getNvoxel();
  }
  util::Timer<T> timer( maxLatticeT, Nvoxel );
  timer.start();

  size_t iT = 0;
  while ( iT < maxLatticeT ) {
    sLatticeL0.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
    // === 6th Step: Collide and Stream Execution ===
    sLatticeL0.collideAndStream();
    if ( boundarytype == stretchedLayer || boundarytype == dampingAndStretched ) {
      sLatticeL1.collideAndStream();
      coarseToFineL1->apply(meta::id<refinement::lagrava::HalfTimeCoarseToFineO>{});
      if ( boundarytype == stretchedLayer ) {
        sLatticeL2.collideAndStream();
        coarseToFineL2->apply(meta::id<refinement::lagrava::HalfTimeCoarseToFineO>{});
        sLatticeL2.collideAndStream();
        coarseToFineL2->apply(meta::id<refinement::lagrava::FullTimeCoarseToFineO>{});
        fineToCoarseL2->apply(meta::id<refinement::lagrava::FineToCoarseO>{});
      }
      sLatticeL1.collideAndStream();
      coarseToFineL1->apply(meta::id<refinement::lagrava::FullTimeCoarseToFineO>{});
      if ( boundarytype == stretchedLayer ) {
        sLatticeL2.collideAndStream();
        coarseToFineL2->apply(meta::id<refinement::lagrava::HalfTimeCoarseToFineO>{});
        sLatticeL2.collideAndStream();
        coarseToFineL2->apply(meta::id<refinement::lagrava::FullTimeCoarseToFineO>{});
        fineToCoarseL2->apply(meta::id<refinement::lagrava::FineToCoarseO>{});
      }
      fineToCoarseL1->apply(meta::id<refinement::lagrava::FineToCoarseO>{});
    }
    // === 7th Step: Computation and Output of the Results ===
    if ( boundarytype == stretchedLayer ) {
      getResults( sLatticeL2, converterL2, iT, sGeometryL2, timer, amplitude, gplot_l2_abs, Lp0, iTvtk, iTplot, "L2", true );
      getResults( sLatticeL1, converterL1, iT, sGeometryL1, timer, amplitude, gplot_l2_abs, Lp0, iTvtk, iTplot, "L1" );
      getResults( sLatticeL0, converterL0, iT, sGeometryL0, timer, amplitude, gplot_l2_abs, Lp0, iTvtk, iTplot, "L0" );
    } else if ( boundarytype == dampingAndStretched ) {
      getResults( sLatticeL1, converterL1, iT, sGeometryL1, timer, amplitude, gplot_l2_abs, Lp0, iTvtk, iTplot, "L1", true );
      getResults( sLatticeL0, converterL0, iT, sGeometryL0, timer, amplitude, gplot_l2_abs, Lp0, iTvtk, iTplot, "L0" );
    } else {
      getResults( sLatticeL0, converterL0, iT, sGeometryL0, timer, amplitude, gplot_l2_abs, Lp0, iTvtk, iTplot, "L0", true );
    }    
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
