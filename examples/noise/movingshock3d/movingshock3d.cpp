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
const int sourMat = 10;
const unsigned int ndim = 3;

const int eternalscale = 4;

typedef enum {periodic, local, damping, dampingAndLocal, eternal} BoundaryType;
typedef enum {shock, point, movingpoint} SourceType;
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
                     int boundaryDepth,
                     Vector<T,3> fluidExtend,
                     Vector<T,3> fluidOrigin
                     )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << std::endl
        << "Prepare Geometry ... Materials:" << std::endl
        << noMat   << " = no Material (default, should be immediately renamed to Check)" << std::endl
        << dampMat << " = far field" << std::endl
        << checMat << " = check (should be empty)" << std::endl
        << fluiMat << " = fluid" << std::endl
        << inflMat << " = inflow" << std::endl
        << outfMat << " = outflow" << std::endl
        << sourMat << " = point source" << std::endl;

  superGeometry.rename( noMat, checMat );   // all nodes to temporary type

  T dx = converter.getConversionFactorLength();
  T dx2 = dx * T(0.5);
  std::shared_ptr<IndicatorF<T,ndim>> dampingField, inflow, outflow, bottomflow, topflow, backflow, frontflow;
  std::shared_ptr<IndicatorF<T,ndim>> fluidDomain = std::make_shared<IndicatorCuboid<T,ndim>>( fluidExtend, fluidOrigin );
  Vector<T,ndim> origin, extend, origin_d, extend_d;

  if ( source == point || source == movingpoint ) {
    Vector<T,3> originSource( 0., 0., 0. );
    std::shared_ptr<IndicatorF3D<T>> pointSourceIndicator = std::make_shared<IndicatorSphere3D<T>>( originSource, dx2 );
    superGeometry.rename( checMat, sourMat, pointSourceIndicator );
  }
  
  switch ( boundarytype ) {
    case eternal:

      clout << "Renaming " << checMat << " to " << dampMat << " in fluidDomain" << std::endl;
      superGeometry.rename( checMat, dampMat, fluidDomain );
      superGeometry.getStatistics().print();

      clout << "Renaming " << checMat << " to " << fluiMat << " everywhere else" << std::endl;
      superGeometry.rename( checMat, fluiMat);
      superGeometry.getStatistics().print();
      break;

    case periodic:

      clout << "Renaming " << checMat << " to " << fluiMat << " everywhere" << std::endl;
      superGeometry.rename( checMat, dampMat );
      superGeometry.getStatistics().print();
      break; // all nodes are fluid

    case local:

      // fluid domain
      clout << "Renaming " << checMat << " to " << dampMat << " in fluidDomain" << std::endl;
      superGeometry.rename( checMat, dampMat, fluidDomain );
      superGeometry.getStatistics().print();
      
      // Set material number for inflow
      origin = superGeometry.getStatistics().getMinPhysR( checMat ) - converter.getConversionFactorLength();
      extend = superGeometry.getStatistics().getMaxPhysR( checMat ) - superGeometry.getStatistics().getMinPhysR( checMat ) + 2*converter.getConversionFactorLength();
      extend[0] = 2*converter.getConversionFactorLength();
      inflow = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      clout << "inflow: origin=" << origin << "; origin+extend=" << origin+extend << std::endl;
      clout << "Renaming " << checMat << " to " << inflMat << " pointing to " << dampMat << " in inflow" << std::endl;
      superGeometry.rename( checMat, inflMat, dampMat, inflow );
      superGeometry.getStatistics().print();
      
      // // Set material number for outflow
      origin[0] = superGeometry.getStatistics().getMaxPhysR( checMat )[0] - converter.getConversionFactorLength();
      outflow = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      clout << "outflow: origin=" << origin << "; origin+extend=" << origin+extend << std::endl;
      clout << "Renaming " << checMat << " to " << outfMat << " pointing to " << dampMat << " in outflow" << std::endl;
      superGeometry.rename( checMat, outfMat, dampMat, outflow );  // outflow to bc, pointing to porous
      superGeometry.getStatistics().print();
      
      // Set material number for bottomflow
      origin = superGeometry.getStatistics().getMinPhysR( checMat ) - converter.getConversionFactorLength();
      extend = superGeometry.getStatistics().getMaxPhysR( checMat ) - superGeometry.getStatistics().getMinPhysR( checMat ) + 2*converter.getConversionFactorLength();
      extend[1] = 2*converter.getConversionFactorLength();
      bottomflow = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      clout << "bottomflow: origin=" << origin << "; origin+extend=" << origin+extend << std::endl;
      clout << "Renaming " << checMat << " to " << bottMat << " pointing to " << dampMat << " in bottomflow" << std::endl;
      superGeometry.rename( checMat, bottMat, dampMat, bottomflow );
      superGeometry.getStatistics().print();

      // Set material number for topflow
      origin[1] = superGeometry.getStatistics().getMaxPhysR( checMat )[1] - converter.getConversionFactorLength();
      topflow = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      clout << "topflow: origin=" << origin << "; origin+extend=" << origin+extend << std::endl;
      clout << "Renaming " << checMat << " to " << topfMat << " pointing to " << dampMat << " in topflow" << std::endl;
      superGeometry.rename( checMat, topfMat, dampMat, topflow );
      superGeometry.getStatistics().print();

      // Set material number for frontflow
      origin = superGeometry.getStatistics().getMinPhysR( checMat ) - converter.getConversionFactorLength();
      extend = superGeometry.getStatistics().getMaxPhysR( checMat ) - superGeometry.getStatistics().getMinPhysR( checMat ) + 2*converter.getConversionFactorLength();
      extend[2] = 2*converter.getConversionFactorLength();
      frontflow = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      clout << "frontflow: origin=" << origin << "; origin+extend=" << origin+extend << std::endl;
      clout << "Renaming " << checMat << " to " << fronMat << " pointing to " << dampMat << " in frontflow" << std::endl;
      superGeometry.rename( checMat, fronMat, dampMat, frontflow );
      superGeometry.getStatistics().print();

      // Set material number for backflow
      origin[2] = superGeometry.getStatistics().getMaxPhysR( checMat )[2] - converter.getConversionFactorLength();
      backflow = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      clout << "Renaming " << checMat << " to " << backMat << " pointing to " << dampMat << " in backflow" << std::endl;
      superGeometry.rename( checMat, backMat, dampMat, backflow );
      superGeometry.getStatistics().print();
      break;

    case damping:

      clout << "Renaming " << checMat << " to " << fluiMat << " in fluidDomain" << std::endl;
      superGeometry.rename( checMat, fluiMat, fluidDomain );
      superGeometry.getStatistics().print();
      
      clout << "Renaming " << checMat << " to " << dampMat << " everywhere else" << std::endl;
      superGeometry.rename( checMat, dampMat );
      superGeometry.getStatistics().print();
      break;  // all remaining nodes are far-field

    case dampingAndLocal:

      // fluid domain including damping layer to damping layer
      origin = superGeometry.getStatistics().getMinPhysR( checMat );
      extend = superGeometry.getStatistics().getMaxPhysR( checMat ) - superGeometry.getStatistics().getMinPhysR( checMat );
      origin[0] += converter.getConversionFactorLength();
      extend[0] -= 2*converter.getConversionFactorLength();
      dampingField = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      clout << "Renaming " << checMat << " to " << dampMat << " in dampingField" << std::endl;
      superGeometry.rename( checMat, dampMat, dampingField );  // all remaining nodes are far-field, except...
      superGeometry.getStatistics().print();
      
      clout << "Renaming " << dampMat << " to " << fluiMat << " in fluidDomain" << std::endl;
      superGeometry.rename( dampMat, fluiMat, fluidDomain );  // all nodes except boundary length
      superGeometry.getStatistics().print();
  
      // Set material number for inflow
      origin = superGeometry.getStatistics().getMinPhysR( checMat ) - converter.getConversionFactorLength();
      extend = superGeometry.getStatistics().getMaxPhysR( checMat ) - superGeometry.getStatistics().getMinPhysR( checMat ) + 2*converter.getConversionFactorLength();
      extend[0] = 2*converter.getConversionFactorLength();
      inflow = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      clout << "Renaming " << checMat << " to " << inflMat << " pointing to " << dampMat << " in inflow" << std::endl;
      superGeometry.rename( checMat, inflMat, dampMat, inflow );  // inflow to bc, pointing to porous
      superGeometry.getStatistics().print();
      
      // Set material number for outflow
      origin[0] = superGeometry.getStatistics().getMaxPhysR( checMat )[0] - converter.getConversionFactorLength();
      outflow = std::make_shared<IndicatorCuboid<T,ndim>>( extend, origin );
      clout << "Renaming " << checMat << " to " << outfMat << " pointing to " << dampMat << " in outflow" << std::endl;
      superGeometry.rename( checMat, outfMat, dampMat, outflow );  // outflow to bc, pointing to porous
      superGeometry.getStatistics().print();
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
                    int boundaryDepth,
                    Vector<T,ndim> domain_lengths,
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
  std::unique_ptr<SuperIndicatorF<T,ndim>> bulkIndicator;

  // Material=dampMat --> dynamics depend..
  switch ( boundarytype ) {
    case eternal:
      bulkIndicator = superGeometry.getMaterialIndicator({dampMat,fluiMat});
      sLattice.defineDynamics<BulkDynamics>(bulkIndicator);
      break;

    case periodic:
      bulkIndicator = superGeometry.getMaterialIndicator({dampMat});
      sLattice.defineDynamics<BulkDynamics>(bulkIndicator);
      break;

    case local:
      sLattice.defineDynamics<BulkDynamics>(superGeometry, dampMat);
      bulkIndicator = superGeometry.getMaterialIndicator({dampMat, inflMat, outfMat, backMat, fronMat, topfMat, bottMat, checMat});
      // boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, inflMat);
      boundary::set<boundary::LocalVelocity>(sLattice, superGeometry.getMaterialIndicator({inflMat, outfMat, backMat, fronMat, topfMat, bottMat, checMat}));
      // setZeroGradientBoundary( sLattice, superGeometry, fronMat );

      // setInterpolatedConvectionBoundary(sLattice, omega, superGeometry, outfMat, uAverage);
      boundary::set<boundary::LocalPressure>(sLattice, superGeometry, outfMat);
      // setLocalVelocityBoundary(sLattice, omega, superGeometry, outfMat);
      // setZeroGradientBoundary(sLattice, superGeometry, outfMat);

      break;

    case damping:
      sLattice.defineDynamics<BulkDynamics>( superGeometry, fluiMat );
      bulkIndicator = superGeometry.getMaterialIndicator( {fluiMat, dampMat} );
      boundary::set<boundary::PerfectlyMatchedLayer>( sLattice, superGeometry, dampMat );
      break;
      
    case dampingAndLocal:
      sLattice.defineDynamics<BulkDynamics>(superGeometry, fluiMat);
      domain_lengths -= 2*converter.getConversionFactorLength();
      bulkIndicator = superGeometry.getMaterialIndicator( {fluiMat, dampMat, inflMat, outfMat} );
      boundary::set<boundary::PerfectlyMatchedLayer>( sLattice, superGeometry, dampMat );
      boundary::set<boundary::LocalVelocity>( sLattice, superGeometry, inflMat );
      boundary::set<boundary::LocalPressure>( sLattice, superGeometry, outfMat );
      // boundary::set<boundary::InterpolatedConvection>(sLattice, superGeometry, outfMat, uAverage);
      break;
  }

  if ( source == point || source == movingpoint ) {
    sLattice.defineDynamics<BulkDynamics>( superGeometry, sourMat );
  }

  // Initial conditions
  AnalyticalConst<ndim,T,T> *ux;
  if ( source == point )    ux = new AnalyticalConst<ndim,T,T>( 0. );
  else                      ux = new AnalyticalConst<ndim,T,T>( Ma );
  AnalyticalConst<ndim,T,T> uy = AnalyticalConst<ndim,T,T>( 0. );
  AnalyticalConst<ndim,T,T> uz = AnalyticalConst<ndim,T,T>( 0. );
  AnalyticalComposed<ndim,T,T> u( *ux, uy, uz );

  T ndatapoints = converter.getResolution();
  T dist = converter.getPhysDeltaX();
  bulkIndicator = superGeometry.getMaterialIndicator({1,2,3,4,5,6,7,8,9,10});  // got all indicators to define RhoU and fields globally
  if ( source == shock ) {
    AcousticPulse<3,T> densityProfile( rho0, amplitude, alpha );
    plotSamplings( densityProfile, ndatapoints, dist, "shock_diag", "density [LU]", diagonal2d );
    plotSamplings( densityProfile, ndatapoints, dist, "shock_hline", "density [LU]", horizontal );
    //Initialize all values of distribution functions to their local equilibrium
    sLattice.defineRhoU( bulkIndicator, densityProfile, u );
    sLattice.iniEquilibrium( bulkIndicator, densityProfile, u );
  } else if ( source == point || source == movingpoint ) {
    AnalyticalConst<ndim,T,T> rho( rho0 );
    //Initialize all values of distribution functions to their local equilibrium
    sLattice.defineRhoU( bulkIndicator, rho, u );
    sLattice.iniEquilibrium( bulkIndicator, rho, u );
  }

  sLattice.setParameter<descriptors::OMEGA>( omega );

  if ( boundarytype == damping || boundarytype == dampingAndLocal ) {
    // Define fields; could also use setParameter as for omega, but then I could not use a flow profile
    sLattice.defineField<descriptors::UX>( bulkIndicator, *ux );
    sLattice.defineField<descriptors::UY>( bulkIndicator, uy );
    sLattice.defineField<descriptors::UZ>( bulkIndicator, uz );
    AnalyticalConst<ndim,T,T> rhoField( rho0 );
    sLattice.defineField<descriptors::DENSITY>( bulkIndicator, rhoField );
    // output damping layer parameter
    DampingTerm<3,T,DESCRIPTOR> sigma_plot( converter, boundaryDepth, domain_lengths );
    plotSamplings( sigma_plot, ndatapoints, dist, "sigma_hline", "sigma", horizontal );
    plotSamplings( sigma_plot, ndatapoints, dist, "sigma_diag", "sigma", diagonal2d );
    DampingTerm<3,T,DESCRIPTOR> sigma( converter, boundaryDepth, domain_lengths, dampingStrength );
    sLattice.defineField<descriptors::DAMPING>( bulkIndicator, sigma );
  }
  
  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Sets fixed far field average values
void setFarFieldValues( UnitConverter<T,DESCRIPTOR> const& converter,
                        SuperLattice<T,DESCRIPTOR>& sLattice, int iT,
                        SuperGeometry<T,3>& superGeometry,
                        T rho0, T farFieldVelocity,
                        BoundaryType boundarytype )
{
  OstreamManager clout( std::cout,"setFarFieldValues" );
  AnalyticalConst<ndim,T,T> ux( farFieldVelocity );
  AnalyticalConst<ndim,T,T> uy( 0. );
  AnalyticalConst<ndim,T,T> uz( 0. );
  AnalyticalComposed<ndim,T,T> u( ux, uy, uz );
  AnalyticalConst<ndim,T,T> rho( rho0 );

  auto farFieldIndicator = superGeometry.getMaterialIndicator( {inflMat, outfMat, topfMat, bottMat, fronMat, backMat} );
  sLattice.defineRhoU( farFieldIndicator, rho, u );
  sLattice.iniEquilibrium( farFieldIndicator, rho, u );
}

// Generates a sinusoidal pressure
void setSourceValues( UnitConverter<T,DESCRIPTOR> const& converter,
                        SuperLattice<T,DESCRIPTOR>& sLattice, int iT,
                        SuperGeometry<T,3>& superGeometry,
                        T rho0, T farFieldVelocity,
                        SourceType source,
                        T amplitude, T frequency )
{
  OstreamManager clout( std::cout,"setSourceValues" );
  AnalyticalConst<ndim,T,T> ux( farFieldVelocity );
  AnalyticalConst<ndim,T,T> uy( 0. );
  AnalyticalConst<ndim,T,T> uz( 0. );
  AnalyticalComposed<ndim,T,T> u( ux, uy, uz );
  T rho = converter.getLatticeDensityFromPhysPressure(
    amplitude * std::sin(converter.getPhysTime(iT) * frequency * 2 * 3.14)
  );
  clout << "rho=" << rho << std::endl;
  AnalyticalConst<ndim,T,T> rho3d( rho );

  auto pointSourceIndicator = superGeometry.getMaterialIndicator( {sourMat} );
  sLattice.defineRhoU( pointSourceIndicator, rho3d, u );
  if ( source == movingpoint ) sLattice.iniEquilibrium( pointSourceIndicator, rho3d, u );
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
    SuperLatticeCuboid3D<T,DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T,DESCRIPTOR> rank( sLattice );
    SuperGeometryF<T,DESCRIPTOR::d> geometryF(superGeometry);
    vtmWriter.write(geometryF);
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

    if ( iT%( iTplot*10 ) == 0 ) {
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
    clout << "uAverage(" << iT << ")= " << *uAverage << std::endl;
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
    clout << "vtmWriter finished" << std::endl;

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

  // === get Command Line Arguments
  std::string outdir            = args.getValueOrFallback<std::string>( "--outdir",            "" );
  const int res                 = args.getValueOrFallback( "--res",               200 );
  const T rho0                  = args.getValueOrFallback( "--rho0",              1.  );
  const T Ma                    = args.getValueOrFallback( "--Ma",                0.1 );  // can be quite high
  const T charV                 = args.getValueOrFallback( "--charV",             1.  );
  T Re                          = args.getValueOrFallback( "--Re",                0.  );
  const T lx                    = args.getValueOrFallback( "--lx",                1.  );
  const T ly                    = args.getValueOrFallback( "--ly",                1.  );
  const T lz                    = args.getValueOrFallback( "--lz",                1.  );
  T tau                         = args.getValueOrFallback( "--tau",               0.  );
  T maxPhysT                    = args.getValueOrFallback( "--tmax",              0.3 );
  size_t maxLatticeT            = args.getValueOrFallback( "--iTmax",             0   );
  size_t nout                   = args.getValueOrFallback( "--nout",              5   );  // minimum number of vtk outputs
  size_t iTout                  = args.getValueOrFallback( "--iTout",             0   );  // iterations for vtk outputs
  T tout                        = args.getValueOrFallback( "--tout",              0.  );  // timestep for vtk outputs
  size_t nplot                  = args.getValueOrFallback( "--nplot",             100 );  // minimum number of plot points // number of plots will be 1/10th
  const int boundaryCondition   = args.getValueOrFallback( "--boundaryCondition", 3   );
  const int sourceType          = args.getValueOrFallback( "--sourceType",        1   );
  const size_t boundaryDepth    = args.getValueOrFallback( "--boundaryDepth",     20  );
  const T amplitude             = args.getValueOrFallback( "--amplitude",         1e-3 );  // physical pressure amplitude
  const T frequency             = args.getValueOrFallback( "--frequency",         100.);  // 1/s
  const T dampingStrength       = args.getValueOrFallback( "--dampingStrength",   1.  );
  size_t overlap                = args.getValueOrFallback( "--overlap",           3   );
  const bool debug              = args.contains("--debug");

  // initialize arguments
  const T b_shock           = 1./20.;
  const T alpha_shock       = log(2.)/(b_shock*b_shock);

  // managing outputs
  std::stringstream outdir_mod;
  if ( outdir == "" ) {
    outdir_mod << "./tmp_";
    switch ( sourceType ) {
      case 1: outdir_mod << "shock"; break;
      case 2: outdir_mod << "point"; break;
      case 3: outdir_mod << "movingpoint"; break;
    }
    switch ( boundaryCondition ) {
      case 1: outdir_mod << "_periodic"; break;
      case 2: outdir_mod << "_local"; break;
      case 3: outdir_mod << "_damping"; break;
      case 4: outdir_mod << "_dampingAndLocal"; break;
      case 5: outdir_mod << "_eternal"; break;
    }
    outdir_mod  << "_Ma" << Ma << "_Re" << Re << "_a" << amplitude
                << "_" << lx << "x" << ly << "x" << lz
                << "_res" << res;
    if ( boundaryCondition == 3 || boundaryCondition == 4 ) outdir_mod << "_bd" << boundaryDepth << "x" << dampingStrength;
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
  clout << "b=" << b_shock << "; alpha=" << alpha_shock << std::endl;

  /* default values */
  BoundaryType boundarytype;
  switch ( boundaryCondition ) {
    case 1:   boundarytype = periodic;        clout << "Boundary condition type specified to periodic."                       << std::endl; break;
    case 2:   boundarytype = local;           clout << "Boundary condition type specified to local."                          << std::endl; break;
    case 3:   boundarytype = damping;         clout << "Boundary condition type specified to damping."                        << std::endl; break;
    case 4:   boundarytype = dampingAndLocal; clout << "Boundary condition type specified to damping with local outside bc."  << std::endl; break;
    case 5:   boundarytype = eternal;         clout << "Boundary condition is solved by just extending the domain to " << eternalscale << " times" << std::endl; break;
  }

  SourceType source = shock;
  T farFieldVelocity = Ma;
  switch ( sourceType ) {
    case 1:   source = shock;       clout << "Source type specified to shock."              << std::endl; break;
    case 2:   source = point;       clout << "Source type specified to point source."       << std::endl;
              farFieldVelocity = 0.; break;
    case 3:   source = movingpoint; clout << "Source type specified to moving point source."<< std::endl; break;
  }

  // determining Reynolds regime (incl. viscosity and relaxation time)
  const T charL             = lx;
  const T cs_LU             = 1 / std::sqrt( 3.0 );
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

  clout << std::endl << "Input to unitconverter: res=" << res << ", tau=" << tau
        << ", charL=" << lx << ", charV=" << charV << ", viscosity=" << viscosity
        << ", physDensity=1.204" << std::endl;
  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> converter(
    (size_t) res,               // resolution
    (T)      tau,               // relaxation time
    (T)      charL,             // charPhysLength: reference length of simulation geometry
    (T)      charV,             // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)      viscosity,         // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)      1.//.204           // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write( "movingshock3d" );

  // === 2nd Step: Prepare Geometry ===

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 6;
#endif

  if ( noOfCuboids > 1 && ( boundarytype == damping || boundarytype == dampingAndLocal) ) {
    clout << "overlap=" << overlap << "; boundaryDepth_LU=" << boundaryDepth << "; setting overlap to >=boundaryDepth." << std::endl;
    overlap = std::max( overlap, boundaryDepth );
  }

  // changing domain size depending on how much area is lost to the damping boundary layer
  T boundaryDepth_pu = converter.getPhysLength( boundaryDepth );
  T lengthDomain, heightDomain, depthDomain;
  switch ( boundarytype )
  {
    case periodic:        // reference domain size
      lengthDomain  = lx;
      heightDomain  = ly;
      depthDomain   = lz;
      break;
    case local:           // add one layer in each direction
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
      lengthDomain  = lx + 2*boundaryDepth_pu + 2*converter.getConversionFactorLength();
      heightDomain  = ly + 2*boundaryDepth_pu;
      depthDomain   = lz + 2*boundaryDepth_pu;
      break;
    case eternal:         // extend the domain all around
      lengthDomain  = eternalscale*lx;
      heightDomain  = eternalscale*ly;
      depthDomain   = eternalscale*lz;
      break;
  }
  clout << "Fluid Domain = " << lx << "x" << ly << "x" << lz << std::endl;
  clout << "Simul Domain = " << lengthDomain << "x" << heightDomain << "x" << depthDomain << std::endl;
  
  // setup domain
  Vector<T,ndim> domain_lengths = { lengthDomain, heightDomain, depthDomain };
  Vector<T,ndim> originDomain( -lengthDomain/2, -heightDomain/2, -depthDomain/2 );  //
  Vector<T,ndim> extendDomain( lengthDomain, heightDomain, depthDomain );  // size of the domain
  std::shared_ptr<IndicatorF<T,ndim>> domain = std::make_shared<IndicatorCuboid<T,ndim>>( extendDomain, originDomain );

  CuboidGeometry<T,ndim> cuboidGeometry(  *(domain),
                                          converter.getConversionFactorLength(),
                                          noOfCuboids
                                          );

  switch ( boundarytype ) {
    case eternal:         cuboidGeometry.setPeriodicity({true, true, true});    break;
    case periodic:        cuboidGeometry.setPeriodicity({true, true, true});    break;
    case local:           cuboidGeometry.setPeriodicity({false, false, false}); break;
    case damping:         cuboidGeometry.setPeriodicity({true, true, true});    break;
    case dampingAndLocal: cuboidGeometry.setPeriodicity({false, true, true});   break;
  }
  
  int fluidMaterial = fluiMat;
  if ( ( boundarytype == periodic ) || ( boundarytype == local || ( boundarytype == eternal ) ) ) fluidMaterial = 1;
  Vector<T,3> fluidExtend = { lx, ly, lz };
  Vector<T,3> fluidOrigin = { -lx/2, -ly/2, -lz/2 };

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,ndim> superGeometry( cuboidGeometry, loadBalancer, overlap );

  clout << std:: endl << "Setup: debug=" << debug << "; boundaryDepthLU=" << boundaryDepth
        << "; boundaryDepthPU=" << boundaryDepth_pu << "; overlapLU=" << superGeometry.getOverlap() << std::endl;
  prepareGeometry( converter, superGeometry, domain, boundarytype, source, res, debug, boundaryDepth, fluidExtend, fluidOrigin );
  clout << "Number of fluid voxels: " << superGeometry.getStatistics().getNvoxel( fluidMaterial ) << std::endl;
  
  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T,DESCRIPTOR> sLattice( superGeometry );

  //prepare Lattice and set boundaryConditions
  prepareLattice( converter, sLattice, superGeometry, rho0, Ma, amplitude, alpha_shock, boundarytype, source, boundaryDepth, domain_lengths, dampingStrength );

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
  // maxLatticeT depends on maximum physical time. If maxLatticeT is provided in command line, it is an upper bound
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
    if ( boundarytype == local || boundarytype == dampingAndLocal ) {
      setFarFieldValues( converter, sLattice, iT, superGeometry, rho0, farFieldVelocity, boundarytype );
    }
    if ( source == point || source == movingpoint ) {
      setSourceValues( converter, sLattice, iT, superGeometry, rho0, farFieldVelocity, source, amplitude, frequency );
    }
    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults(sLattice, converter, iT, superGeometry, timer,
               velocityFlux, pressureFlux,
               rho0, amplitude,
               gplot_l2_abs, Lp0,
               iTvtk, maxLatticeT, iTplot,
               boundarytype, fluidMaterial );

    // test density failure
    if ( sLattice.getStatistics().getAverageRho() > 2. ) {
      clout << "breaking simulation loop because density is too high: "
            << sLattice.getStatistics().getAverageRho() << std::endl;
      break;
    }
    iT++;
  }
  
  clout << "Simulation stopped after " << iT << "/" << maxLatticeT << " iterations." << std::endl;

  gplot_l2_abs.setYrange(1e-3, 1);
  gplot_l2_abs.setLogScale(2);
  gplot_l2_abs.writePNG(-1, -1, "gplot_l2_abs");

  timer.stop();
  timer.printSummary();
  std::cout.rdbuf(originalCoutBuffer);
}
