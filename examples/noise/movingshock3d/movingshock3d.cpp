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
// using DESCRIPTOR = D3Q19<>;
// using BulkDynamics = KBCdynamics<T,DESCRIPTOR>;
using DESCRIPTOR = D3Q27<>;
// using BulkDynamics = BGKdynamics<T,DESCRIPTOR>;
using BulkDynamics = olb::dynamics::Tuple<T,
                                          DESCRIPTOR,
                                          momenta::BulkTuple,
                                          equilibria::Incompressible,
                                          collision::TRT>;

const int noMat   = 0;
const int dampMat = 1;
const int checMat = 2;
const int fluiMat = 3;
const int inflMat = 4;
const int outfMat = 5;
const int bottMat = 6;
const int topfMat = 7;
const int fronMat = 8;
const int backMat = 9;
const unsigned int ndim = 3;

typedef enum {periodic, local, damping, dampingAndLocal} BoundaryType;
typedef enum {shock, pointsource} SourceType;
typedef enum {horizontal, vertical, diagonal2d, diagonal3d} SamplingDirection;

T* uAverage = NULL;

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
        << dampMat << " = far field" << std::endl
        << checMat << " = check (should be empty)" << std::endl
        << fluiMat << " = fluid" << std::endl
        << inflMat << " = inflow" << std::endl
        << outfMat << " = outflow" << std::endl;

  superGeometry.rename( noMat, checMat );   // all nodes to temporary type
  
  clout << "Geometry before different BC types:" << std::endl;
  clout << "minPhysR=" << superGeometry.getStatistics().getMinPhysR( checMat ) << "; maxPhysR=" << superGeometry.getStatistics().getMaxPhysR( checMat ) << std::endl;
  superGeometry.getStatistics().print();

  T bd_pu, dx, dx2;
  dx = converter.getConversionFactorLength();
  dx2 = converter.getConversionFactorLength() * T(0.5);
  clout << "dx=" << dx << "; dx/2=" << dx2 << "; 1 m_LU = " << converter.getPhysLength(1) << " m_PU" << std::endl;
  std::shared_ptr<IndicatorF<T,ndim>> fluid_domain, dampingField, inflow, outflow, bottomflow, topflow, backflow, frontflow;
  Vector<T,ndim> origin, extend, origin_f, extend_f, origin_d, extend_d;
  SuperVTMwriter3D<T> vtmWriter( "movingshock3d_local" );
  SuperLatticeGeometry3D<T,DESCRIPTOR> geometry( superGeometry );
  geometry.getName() = "geometry0";
  vtmWriter.write( geometry );
  
  switch ( boundarytype ) {
    case periodic:
      superGeometry.rename( checMat, dampMat );
      break; // all nodes are fluid

    case local:
      // fluid domain including damping layer to damping layer
      origin_d = superGeometry.getStatistics().getMinPhysR( checMat ) - dx2;
      extend_d = superGeometry.getStatistics().getMaxPhysR( checMat ) - superGeometry.getStatistics().getMinPhysR( checMat ) + dx;
      origin_f = origin_d;
      extend_f = extend_d;
      origin_f[0] += dx;
      origin_f[1] += dx;
      extend_f[0] -= T(2)*dx;
      extend_f[1] -= T(2)*dx;
      clout << "fluid_domain: origin=" << origin_f << "; origin+extend=" << origin_f+extend_f << std::endl;
      fluid_domain = std::make_shared<IndicatorCuboid<T,ndim>>( extend_f, origin_f );
      superGeometry.rename( checMat, dampMat, fluid_domain );
      geometry.getName() = "geometry1";
      vtmWriter.write( geometry );
      
      // Set material number for inflow
      origin = origin_d;
      extend = extend_d;
      extend[0] = dx;
      clout << "inflow: origin=" << origin << "; origin+extend=" << origin+extend << std::endl;
      inflow = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      superGeometry.rename( checMat, inflMat, dampMat, inflow );  // inflow to bc, pointing to porous
      geometry.getName() = "geometry2";
      vtmWriter.write( geometry );
      
      // // Set material number for outflow
      origin[0] = origin_d[0] + extend_d[0] - dx;
      clout << "outflow: origin=" << origin << "; origin+extend=" << origin+extend << std::endl;
      outflow = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      superGeometry.rename( checMat, outfMat, dampMat, outflow );  // outflow to bc, pointing to porous
      geometry.getName() = "geometry3";
      vtmWriter.write( geometry );

      // superGeometry.rename( checMat, inflMat, dampMat, domain );  // inflow to bc, pointing to porous
      // geometry.getName() = "geometry4";
      // vtmWriter.write( geometry );
      // superGeometry.rename( checMat, inflMat, dampMat, domain );  // inflow to bc, pointing to porous
      // geometry.getName() = "geometry5";
      // vtmWriter.write( geometry );
      // superGeometry.rename( checMat, inflMat, dampMat, domain );  // inflow to bc, pointing to porous
      // geometry.getName() = "geometry6";
      // vtmWriter.write( geometry );

      // Set material number for bottomflow
      origin = origin_d;
      extend = extend_d;
      extend[1] = dx;
      clout << "bottomflow: origin=" << origin << "; origin+extend=" << origin+extend << std::endl;
      bottomflow = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      superGeometry.rename( checMat, bottMat, dampMat, bottomflow );
      superGeometry.getStatistics().print();
      geometry.getName() = "geometry4";
      vtmWriter.write( geometry );

      // Set material number for topflow
      origin[1] = origin_d[1] + extend_d[1] - dx;
      clout << "topflow: origin=" << origin << "; origin+extend=" << origin+extend << std::endl;
      topflow = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      superGeometry.rename( checMat, topfMat, dampMat, topflow );
      superGeometry.getStatistics().print();
      geometry.getName() = "geometry5";
      vtmWriter.write( geometry );

      // Set material number for frontflow
      origin = origin_d;
      extend = extend_d;
      extend[2] = dx;
      clout << "frontflow: origin=" << origin << "; origin+extend=" << origin+extend << std::endl;
      frontflow = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      superGeometry.rename( checMat, fronMat, dampMat, frontflow );
      superGeometry.getStatistics().print();
      geometry.getName() = "geometry6";
      vtmWriter.write( geometry );

      // Set material number for backflow
      origin[2] = origin_d[2] + extend_d[2] - dx;
      clout << "backflow: origin=" << origin << "; origin+extend=" << origin+extend << std::endl;
      backflow = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      superGeometry.rename( checMat, backMat, dampMat, backflow );
      superGeometry.getStatistics().print();
      geometry.getName() = "geometry7";
      vtmWriter.write( geometry );

      // all remaining to outsMat
      // superGeometry.rename( checMat, inflMat );
      break;

    case damping:
      bd_pu = converter.getPhysLength(boundary_depth);
      origin = superGeometry.getStatistics().getMinPhysR( checMat ) + bd_pu;
      extend = superGeometry.getStatistics().getMaxPhysR( checMat ) - superGeometry.getStatistics().getMinPhysR( checMat ) - 2*bd_pu;
      fluid_domain = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );

      superGeometry.rename( checMat, fluiMat, fluid_domain );
      superGeometry.rename( checMat, dampMat );
      break;  // all remaining nodes are far-field

    case dampingAndLocal:
      // fluid domain including damping layer to damping layer
      origin = superGeometry.getStatistics().getMinPhysR( checMat );
      extend = superGeometry.getStatistics().getMaxPhysR( checMat ) - superGeometry.getStatistics().getMinPhysR( checMat );
      origin[0] += converter.getConversionFactorLength();
      extend[0] -= 2*converter.getConversionFactorLength();
      dampingField = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      superGeometry.rename( checMat, dampMat, dampingField );  // all remaining nodes are far-field, except...

      // fluid domain part to fluid
      bd_pu = converter.getPhysLength(boundary_depth);
      origin = superGeometry.getStatistics().getMinPhysR( checMat ) + bd_pu;
      extend = superGeometry.getStatistics().getMaxPhysR( checMat ) - superGeometry.getStatistics().getMinPhysR( checMat ) - 2*bd_pu;
      fluid_domain = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      superGeometry.rename( dampMat, fluiMat, fluid_domain );  // all nodes except boundary length
  
      // Set material number for inflow
      origin = superGeometry.getStatistics().getMinPhysR( checMat );
      extend = superGeometry.getStatistics().getMaxPhysR( checMat ) - superGeometry.getStatistics().getMinPhysR( checMat );
      origin[0] -= converter.getConversionFactorLength();
      extend[0] = 2*converter.getConversionFactorLength();
      inflow = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      superGeometry.rename( checMat, inflMat, dampMat, inflow );  // inflow to bc, pointing to porous
      
      // Set material number for outflow
      origin[0] = superGeometry.getStatistics().getMaxPhysR( checMat )[0]-converter.getConversionFactorLength();
      outflow = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      superGeometry.rename( checMat, outfMat, dampMat, outflow );  // outflow to bc, pointing to porous
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

void plotSamplings( AnalyticalF3D<T,T>& data, size_t ndatapoints, T dist,
                    std::string title, std::string ylabel,
                    SamplingDirection direction, bool halfDomain = true,
                    bool setRange = false, T ymin=0, T ymax=0 )
{
  Gnuplot<T> gplot( title );
  gplot.setLabel( "distance [m]", ylabel );
  int nmin = 0;
  if ( !halfDomain ) nmin = -int( ndatapoints/2 );
  for ( int n = nmin; n <= int(ndatapoints/2); n++ ) {
    T input[ndim] = {0,0,0};
    T distance = 0;
    switch ( direction ) {
      case horizontal:  input[0] = n*dist;                                        distance = n*dist;              break;
      case vertical:    input[1] = n*dist;                                        distance = n*dist;              break;
      case diagonal2d:  input[0] = n*dist; input[1] = n*dist;                     distance = n*dist*std::sqrt(2); break;
      case diagonal3d:  input[0] = n*dist; input[1] = n*dist; input[2] = n*dist;  distance = n*dist*std::sqrt(3); break;
    }
    T output[1];
    data( output, input );
    gplot.setData( distance, output[0] );
  }
  if ( setRange ) gplot.setYrange( ymin, ymax );
  gplot.writePNG( -1, -1, title );
}

void setUAverage( SuperLattice<T,DESCRIPTOR>& sLattice,
                  SuperGeometry<T,ndim>& superGeometry,
                  BoundaryType boundarytype
                  )
{
  sLattice.setProcessingContext(ProcessingContext::Evaluation);  // important for synchronization (?) on GPU

  // update outflow boundary value (adaptive convection boundary for smaller domains)
  SuperLatticeVelocity3D<T, DESCRIPTOR> velocity( sLattice );
  int fluidMaterial = dampMat;
  if ( boundarytype == dampingAndLocal ) fluidMaterial = fluiMat;
  std::unique_ptr<SuperIndicatorF<T,ndim>> fluidIndicator = superGeometry.getMaterialIndicator({fluidMaterial});
  SuperSum3D<T> sum( velocity, fluidIndicator );
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
                    SourceType source,
                    int boundary_depth,
                    Vector<T,ndim> domain_lengths,
                    T damping_strength
                    )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();
  T tmp = T();
  uAverage = &tmp;  //*uAverage = 0.0;
  if ( ( boundarytype == local ) || ( boundarytype == dampingAndLocal ) ) {
    setUAverage(sLattice, superGeometry, boundarytype);    
  }

  // Material=fluiMat --> bulk dynamics
  std::unique_ptr<SuperIndicatorF<T,ndim>> bulkIndicator;

  // Material=dampMat --> dynamics depend..
  switch ( boundarytype ) {
    case periodic:
      bulkIndicator = superGeometry.getMaterialIndicator({dampMat});
      sLattice.defineDynamics<BulkDynamics>(bulkIndicator);
      break;
    case local:
      sLattice.defineDynamics<BulkDynamics>(superGeometry, dampMat);
      // domain_lengths -= 2*converter.getConversionFactorLength();
      bulkIndicator = superGeometry.getMaterialIndicator({dampMat, inflMat, outfMat});
      // Material=inflMat -->inlet
      setLocalVelocityBoundary(sLattice, omega, superGeometry, inflMat);
      // Material=outfMat -->outlet
      // setInterpolatedConvectionBoundary(sLattice, omega, superGeometry, outfMat, uAverage);
      setLocalPressureBoundary(sLattice, omega, superGeometry, outfMat);
      // Material=outsMat -->zero gradient
      setLocalVelocityBoundary(sLattice, omega, superGeometry, bottMat);
      setLocalVelocityBoundary(sLattice, omega, superGeometry, topfMat);
      // setZeroGradientBoundary(sLattice, superGeometry, bottMat);
      // setZeroGradientBoundary(sLattice, superGeometry, topfMat);
      setZeroGradientBoundary(sLattice, superGeometry, backMat);
      setZeroGradientBoundary(sLattice, superGeometry, fronMat);
      break;
    case damping:
      sLattice.defineDynamics<BulkDynamics>(superGeometry, fluiMat);
      bulkIndicator = superGeometry.getMaterialIndicator({fluiMat, dampMat});
      setDampingBoundary<T,DESCRIPTOR>(sLattice, superGeometry, dampMat);
      break;
    case dampingAndLocal:
      sLattice.defineDynamics<BulkDynamics>(superGeometry, fluiMat);
      domain_lengths -= 2*converter.getConversionFactorLength();
      bulkIndicator = superGeometry.getMaterialIndicator({fluiMat, dampMat, inflMat, outfMat});
      setDampingBoundary<T,DESCRIPTOR>(sLattice, superGeometry, dampMat);
      // Material=inflMat -->inlet
      setLocalVelocityBoundary(sLattice, omega, superGeometry, inflMat);
      // Material=outfMat -->outlet
      setInterpolatedConvectionBoundary(sLattice, omega, superGeometry, outfMat, uAverage);
      break;
  }

  // Initial conditions
  AnalyticalConst<ndim,T,T> ux = AnalyticalConst<ndim,T,T>( Ma );
  AnalyticalConst<ndim,T,T> uy = AnalyticalConst<ndim,T,T>( 0. );
  AnalyticalConst<ndim,T,T> uz = AnalyticalConst<ndim,T,T>( 0. );
  AnalyticalComposed<ndim,T,T> u( ux, uy, uz );

  T ndatapoints = converter.getResolution();
  T dist = converter.getPhysDeltaX();
  bulkIndicator = superGeometry.getMaterialIndicator({1,2,3,4,5});  // got all indicators to define RhoU and fields globally
  if ( source == shock ) {
    AcousticPulse<3,T> pressureProfile( rho0, amplitude, alpha );
    plotSamplings( pressureProfile, ndatapoints, dist, "shock_diag", "density [LU]", diagonal2d );
    plotSamplings( pressureProfile, ndatapoints, dist, "shock_hline", "density [LU]", horizontal );

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
  plotSamplings( sigma_plot, ndatapoints, dist, "sigma_hline", "sigma", horizontal );
  plotSamplings( sigma_plot, ndatapoints, dist, "sigma_diag", "sigma", diagonal2d );

  DampingTerm<3,T,DESCRIPTOR> sigma( converter, boundary_depth, domain_lengths, damping_strength );
  sLattice.defineField<descriptors::DAMPING>( bulkIndicator, sigma );
  
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

  auto farFieldIndicator = superGeometry.getMaterialIndicator( {inflMat, outfMat} );
  sLattice.defineRhoU( farFieldIndicator, rho, u );
  sLattice.iniEquilibrium( farFieldIndicator, rho, u );

  sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
}

T L2Norm(SuperLattice<T,DESCRIPTOR>& sLattice,
         SuperGeometry<T,ndim>& superGeometry,
         int iT,
         UnitConverter<T,DESCRIPTOR> const& converter,
         int fluidMaterial ) {
  OstreamManager clout(std::cout, "L2Norm");

  T result[ndim];
  int tmp[] = {int()};

  SuperLatticePhysPressure3D<T,DESCRIPTOR> pressurenD(sLattice, converter);
  AnalyticalConst<ndim,T,T> rho0nD( 0. );

  std::unique_ptr<SuperIndicatorF<T,ndim>> indicatorF = superGeometry.getMaterialIndicator({fluidMaterial});
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
  SuperVTMwriter3D<T> vtmWriter( "movingshock3d" );

  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );

  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

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

  // Writes to plot 100 times
  if ( iT%iTplot==0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);  // important for synchronization (?) on GPU

    if ( iT%( iTplot*5 ) == 0 ) {
      if ( iT%( iTplot*10 ) == 0 ) {
        velocityFlux.print();
        pressureFlux.print();
      }
      
      // write to terminal
      timer.update( iT );
      timer.printStep();

      // Lattice statistics console output
      sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
    }

    gplot_l2_abs.setData(T(iT), L2Norm( sLattice, superGeometry, iT,
                                        converter, fluidMaterial ) / Lp0 );

    if ( iT%( iTplot*20 ) == 0 ) {
      std::stringstream ss;
      ss << std::setw(4) << std::setfill('0') << iT;
      T dist = converter.getPhysDeltaX();
      T ndatapoints = converter.getResolution(); // number of data points on line
      AnalyticalFfromSuperF3D<T> pressure_interpolation( pressure, true, true );
      T pmin(converter.getPhysPressure(-amplitude/200));
      T pmax(converter.getPhysPressure(+amplitude/200));
      plotSamplings( pressure_interpolation, ndatapoints, dist, "pressure_hline_" + ss.str(), "pressure [PU]", horizontal, false, true, pmin, pmax );
      plotSamplings( pressure_interpolation, ndatapoints, dist, "pressure_vline_" + ss.str(), "pressure [PU]", vertical, false, true, pmin, pmax );
      plotSamplings( pressure_interpolation, ndatapoints, dist, "pressure_diagonal_" + ss.str(), "pressure [PU]", diagonal2d, false, true, pmin, pmax );
      // AnalyticalFfromSuperF3D<T> velocity_interpolation( velocity, true, true );
      // plotSamplings( velocity_interpolation, ndatapoints, dist, "velocity_diagonal_" + ss.str(), "velocity [LU]", diagonal2d, false );
    }

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }

  if ( ( boundarytype == dampingAndLocal || boundarytype == local ) && ( iT%50 == 0 ) ) {
    setUAverage(sLattice, superGeometry, boundarytype);
  }

  // get VTK and images
  if ( iT%iTvtk==0 || sLattice.getStatistics().getAverageRho() > 2. ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);  // important for synchronization (?) on GPU

    clout << "vtmWriter startet" << std::endl;
    // VTK
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
    // jpeg_ParamP.maxValue = converter.getCharPhysVelocity()*1.1;
    // jpeg_ParamP.minValue = 0;
    jpeg_ParamP.colour = "rainbow";
    jpeg_ParamP.fullScreenPlot = true;
    heatmap::write(planeReduction, iT, plotParam);

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }

  // Saves lattice data
  if ( iT%imax==0 && iT>0 ) {
    //clout << "Checkpointing the system at t=" << iT << std::endl;
    //sLattice.save( "movingshock3d.checkpoint" );
    // The data can be reloaded using
    //     sLattice.load("movingshock3d.checkpoint");
  }
}


int main( int argc, char* argv[] )
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
  const int res                 = args.getValueOrFallback( "--res",           141  );
  const T rho0                  = args.getValueOrFallback( "--rho0",          1   );
  const T Ma                    = args.getValueOrFallback( "--Ma",            0.25 );  // can be quite high
  const T charV                 = args.getValueOrFallback( "--charV",         1 );
  T Re                          = args.getValueOrFallback( "--Re",            0   );
  T lengthDomain                = args.getValueOrFallback( "--lengthDomain",  2   );
  T heightDomain                = args.getValueOrFallback( "--heightDomain",  4   );
  T depthDomain                 = args.getValueOrFallback( "--depthDomain",   4   );
  T tau                         = args.getValueOrFallback( "--tau",           0   );
  T maxPhysT                    = args.getValueOrFallback( "--tmax",          1   );
  size_t iTmax                  = args.getValueOrFallback( "--imax",          0   );
  size_t nout                   = args.getValueOrFallback( "--nout",          5   );  // minimum number of vtk outputs
  size_t iout                   = args.getValueOrFallback( "--iout",          0   );  // iterations for vtk outputs
  T tout                        = args.getValueOrFallback( "--tout",          0   );  // timestep for vtk outputs
  size_t nplot                  = args.getValueOrFallback( "--nplot",         100 );  // minimum number of plot points // number of plots will be 1/10th
  const int boundary_condition  = args.getValueOrFallback( "--boundary_condition", 3 );
  const int source_type         = args.getValueOrFallback( "--source_type",   1   );
  const bool debug              = args.contains("--debug");
  const size_t boundary_depth   = args.getValueOrFallback( "--boundary_depth", 20 );
  const T amplitude             = args.getValueOrFallback( "--amplitude",     1e-6 );  // physical pressure amplitude
  const T damping_strength      = args.getValueOrFallback( "--damping_strength", 1.);
  size_t overlap                = args.getValueOrFallback( "--overlap",       3   );

  // initialize arguments
  const T b_shock           = 3;
  const T alpha_shock       = log(2)/(b_shock*b_shock);

  // managing outputs
  std::stringstream outdir_mod;
  outdir_mod << outdir;
  switch ( source_type ) {
    case 1: outdir_mod << "_shock"; break;
    case 2: outdir_mod << "_pointsource"; break;
  }
  switch ( boundary_condition ) {
    case 1: outdir_mod << "_periodic"; break;
    case 2: outdir_mod << "_local"; break;
    case 3: outdir_mod << "_damping"; break;
    case 4: outdir_mod << "_dampingAndLocal"; break;
  }
  outdir_mod << "_Ma" << Ma << "_Re" << Re << "_a" << amplitude << "_" << lengthDomain << "x" << heightDomain << "x" << depthDomain << "_res" << res << "_overlap" << overlap << "_bd" << boundary_depth << "x" << damping_strength;

  singleton::directories().setOutputDir( outdir_mod.str()+"/" );  // set output directory
  std::ofstream fileStream(outdir_mod.str()+"/output.txt");
  // Create a TeeBuffer that writes to both the console and the file
  DoubleBuffer doubleBuffer(std::cout.rdbuf(), fileStream.rdbuf());
  std::streambuf* originalCoutBuffer = std::cout.rdbuf(&doubleBuffer);
  OstreamManager clout( std::cout, "main" );
  clout << "outdir set to " << outdir_mod.str() << std::endl;

  /* default values */
  BoundaryType boundarytype;
  switch ( boundary_condition ) {
    case 1:   boundarytype = periodic;        clout << "Boundary condition type specified to periodic."                       << std::endl; break;
    case 2:   boundarytype = local;           clout << "Boundary condition type specified to local."                          << std::endl; break;
    case 3:   boundarytype = damping;         clout << "Boundary condition type specified to damping."                        << std::endl; break;
    case 4:   boundarytype = dampingAndLocal; clout << "Boundary condition type specified to damping with local outside bc."  << std::endl; break;
    default:  boundarytype = damping;         clout << "Boundary condition type not specified. Default to damping."           << std::endl; break;
  }

  SourceType source = shock;
  switch ( source_type ) {
    case 1:   source = shock;       clout << "Source type specified to shock."              << std::endl; break;
    case 2:   source = pointsource; clout << "Source type specified to point source."       << std::endl; break;
    default:  source = shock;       clout << "Source type not specified. Default to Shock." << std::endl; break;
  }

  // determining Reynolds regime (incl. viscosity and relaxation time)
  const T charL             = lengthDomain;
  const T cs_LU             = 1 / std::sqrt(3.0);
  T viscosity;
  if ( tau == 0 ) {
    if ( Re == 0 ) {
      viscosity             = 1.48e-5; // 1e-2;
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

  clout << "Input to unitconverter: res=" << res << ", tau=" << tau << ", lengthDomain=" << lengthDomain << ", charV=" << charV << ", viscosity=" << viscosity << ", physDensity=1.204" << std::endl;
  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> converter(
    (size_t) res,               // resolution
    (T)      tau,               // relaxation time
    (T)      charL,             // charPhysLength: reference length of simulation geometry
    (T)      charV,             // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)      viscosity,         // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)      1.204              // physDensity: physical density in __kg / m^3__
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

  // changing domain size depending on how much area is lost to the damping boundary layer
  T boundary_depth_pu = converter.getPhysLength(boundary_depth);
  switch ( boundarytype )
  {
    case periodic:
      lengthDomain  -= 2*boundary_depth_pu;
      heightDomain  -= 2*boundary_depth_pu;
      depthDomain   -= 2*boundary_depth_pu;
      break;
    case local:           // add one layer in each direction but remove damping layer depth
      lengthDomain  -= 2*boundary_depth_pu;
      heightDomain  -= 2*boundary_depth_pu;
      depthDomain   -= 2*boundary_depth_pu;
      lengthDomain  += 2*converter.getConversionFactorLength();
      heightDomain  += 2*converter.getConversionFactorLength();
      depthDomain   += 2*converter.getConversionFactorLength();
      break;
    case damping:         // reference domain size
      break;
    case dampingAndLocal: // add one layer in flow direction
      lengthDomain  += 2*converter.getConversionFactorLength();
      break;
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

  switch ( boundarytype ) {
    case periodic:        cuboidGeometry.setPeriodicity(true, true, true);    break;
    case local:           cuboidGeometry.setPeriodicity(false, false, false); break;
    case damping:         cuboidGeometry.setPeriodicity(true, true, true);    break;
    case dampingAndLocal: cuboidGeometry.setPeriodicity(false, true, true);   break;
  }
  
  int fluidMaterial = fluiMat;
  if ( ( boundarytype == periodic ) || ( boundarytype == local ) ) fluidMaterial = 1;

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,ndim> superGeometry( cuboidGeometry, loadBalancer, overlap );

  clout << "Setup: debug=" << debug << "; boundary_depth=" << boundary_depth << "; bd_depth_pu=" << boundary_depth_pu << "; overlap=" << superGeometry.getOverlap() << std::endl;
  prepareGeometry( converter, superGeometry, domain, boundarytype, source, res, debug, boundary_depth );
  clout << "Number of fluid voxels: " << superGeometry.getStatistics().getNvoxel( fluidMaterial ) << std::endl;
  
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
      Vector<T, ndim>( 0., 0., 1. )
  );

  SuperPlaneIntegralFluxPressure3D<T> pressureFlux(
      sLattice,
      converter,
      superGeometry,
      { lengthDomain/T(2), heightDomain/T(2), depthDomain/T(2) },
      Vector<T, ndim>( 0., 0., 1. )
  );

  Gnuplot<T> gplot_l2_abs("l2_absolute");//, Gnuplot<T>::LOGLOG, Gnuplot<T>::OFF);
  gplot_l2_abs.setLabel("time []", "absolute L2 norm []");
  T Lp0 = L2Norm( sLattice, superGeometry, 0, converter, fluidMaterial );

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
  // iTmax depends on maximum physical time. If iTmax is provided in command line, it is an upper bound
  if ( iTmax == 0 ) iTmax = converter.getLatticeTime( maxPhysT );
  else iTmax = std::min( iTmax, converter.getLatticeTime( maxPhysT ) );
  // nout is the minimum number of vtk outputs --> take max between nout and nout derived from iout or tout
  size_t nout_from_iout = 0, nout_from_tout = 0;
  if ( iout != 0 ) { nout_from_iout = size_t( iTmax / iout ); nout = std::max( nout, nout_from_iout ); }
  if ( tout != 0 ) { nout_from_tout = size_t( iTmax / tout ); nout = std::max( nout, nout_from_tout ); }
  size_t iTvtk    = std::max(int( iTmax / nout ), 1);
  size_t iTplot   = std::max(int( iTmax / nplot ), 1);

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( iTmax, superGeometry.getStatistics().getNvoxel() );
  timer.start();

  size_t iT = 0;
  while ( iT < iTmax ) {
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
               iTvtk, iTmax, iTplot,
               boundarytype, fluidMaterial );

    // test density failure
    if ( sLattice.getStatistics().getAverageRho() > 2. ) {
      clout << "breaking simulation loop because density is too high: "
            << sLattice.getStatistics().getAverageRho() << std::endl;
      break;
    }
    iT++;
  }
  
  clout << "Simulation stopped after " << iT << " iterations." << std::endl;

  gplot_l2_abs.setYrange(1e-3, 1);
  gplot_l2_abs.setLogScale(2);
  gplot_l2_abs.writePNG(-1, -1, "gplot_l2_abs");
  timer.stop();
  timer.printSummary();
  std::cout.rdbuf(originalCoutBuffer);
}
