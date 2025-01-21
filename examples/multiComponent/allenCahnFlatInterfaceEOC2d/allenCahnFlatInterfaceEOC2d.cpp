/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Tim Bingert
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

/** allenCahnFlatInterfaceEOC2d.cpp
 * In this example a Young-Laplace test or a flat interface test are
 * performed. A rectangular domain of fluid 2 is immersed in fluid 1.
 * A diffusive interface forms with the profile of a hyperbolic tangent
 * whose accuracy is measured for multiple resolutions in order to test
 * the models experimental order of convergence. The equilibrium pressure
 * is also investigated in a similar way. The grid refinement can be
 * performed with either a constant or a decreasing Cahn number.
 *
 * This example shows the simplest application of the hybrid phase field
 * Allen-Cahn model with periodic boundaries, based on:
 *
 * Liu, Xi, Zhenhua Chai, and Baochang Shi. "Improved hybrid Allen-Cahn
 * phase-field-based lattice Boltzmann method for incompressible two-phase
 * flows." Physical Review E 107.3 (2023): 035308.
 */
#include <fstream>
#include <iostream>
#include "olb2D.h"
#include "olb2D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;
using NSDESCRIPTOR = D2Q9<RHO,NABLARHO,FORCE,EXTERNAL_FORCE,TAU_EFF,STATISTIC>;
using ACDESCRIPTOR = D2Q9<FORCE,SOURCE,SOURCE_OLD,VELOCITY,OLD_PHIU,STATISTIC,PSI,NORMGRADPSI,SCALAR,PSI0,TOP,BOTTOM>;
using NSBulkDynamics = MPIncBGKdynamics<T,NSDESCRIPTOR>;
using ACBulkDynamics = AllenCahnBGKdynamics<T,ACDESCRIPTOR>;
using Coupling = AllenCahnPostProcessor;
using Helper = AllenCahnNonLocalHelper;

// Parameters for the simulation setup
const int Nx = 2;                                        // domain resolution x
int Ny = 100;                                            // domain resolution y
int phaseLength = Ny/2;                                  // [lattice units]
const T charPhysLength = 100e-6;                         // charPhysLength [physical units]
const T Re = 0.;                                         // definition: Reynolds number of continuous phase
T sigma = 0.64;                                          // surface tension [lattice units]
const T surfaceTension = 0.072;                          // surface tension [physical units]
const T viscosityH2O = 9e-7;                             // physViscosity H2O liquid [physical units]
const T DeltaRho = 40.;
const T tau_l = 0.55;                                    // relaxation time Water lattice [lattice units]
const T tau_v = 17.*(tau_l-0.5)+0.5;                     // relaxation time Air lattice [lattice units]
const T tau_mobil = 0.8;                                 // relaxation time for interface mobility [lattice units]
T w = 5.;
std::vector<T> rhos = {0.03, 25.};
const int maxIter  = 10000000;
const int vtkIter  = 20000;
const int statIter = 20000;
const bool Cahn_const = true;

void prepareGeometry( SuperGeometry<T,2>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;
  superGeometry.rename( 0,1 );
  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

template <typename STAGE>
void signedDistanceFunction( SuperLattice<T,ACDESCRIPTOR>& sLatticeAC,
                           SuperGeometry<T,2>& superGeometry, T w )
{
  int i = 0;
  T max = 1.;
  while ( max <= 3*w ) {
    int in[2];
    T out = 0;
    sLatticeAC.getCommunicator(STAGE{}).communicate();
    sLatticeAC.executePostProcessors(STAGE{});
    SuperLatticeExternalScalarField2D<T, ACDESCRIPTOR, PSI> psi( sLatticeAC );
    SuperEuklidNorm2D<T, ACDESCRIPTOR> normPsi(psi);
    SuperSum2D<T> Psi_total_(normPsi, superGeometry, 1);
    SuperMax2D<T,T> Max_psi_(psi, superGeometry, 1);
    Psi_total_(&out, in);
    Max_psi_(&max, in);
    i++;
  }
}

template <typename SuperLatticeCoupling>
void prepareLattice( SuperLattice<T,NSDESCRIPTOR>& sLatticeNS,
                     SuperLattice<T,ACDESCRIPTOR>& sLatticeAC,
                     SuperLatticeCoupling& coupling,
                     UnitConverter<T,NSDESCRIPTOR> const& converter,
                     SuperGeometry<T,2>& superGeometry, int Ny, int phaseLength, T w, T sigma )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  // define lattice Dynamics
  sLatticeNS.defineDynamics<NoDynamics>(superGeometry, 0);
  sLatticeNS.defineDynamics<NSBulkDynamics>(superGeometry, 1);
  sLatticeAC.defineDynamics<NoDynamics>(superGeometry, 0);
  sLatticeAC.defineDynamics<ACBulkDynamics>(superGeometry, 1);

  // bulk initial conditions
  Vector<T,2> u(0., 0.);
  AnalyticalConst2D<T,T> zeroVelocity( u );
  Vector<T,2> u0(0., 0.);
  AnalyticalConst2D<T,T> zeroVelocity0( u0 );

  AnalyticalConst2D<T,T> one ( 1. );
  AnalyticalConst2D<T,T> zero ( 0. );
  AnalyticalConst2D<T,T> point ( 0. );
  AnalyticalConst2D<T,T> rhov ( rhos[0] );
  AnalyticalConst2D<T,T> rhol ( rhos[1] );
  AnalyticalConst2D<T,T> tauv ( tau_v );
  AnalyticalConst2D<T,T> taul ( tau_l );

  SmoothIndicatorFactoredCuboid2D<T,T> interfaceAC( {Nx/2., Ny/2.}, 0, phaseLength+1, w/2, 0, {0,0}, 0, -1. );
  AnalyticalIdentity2D<T,T> phi( one + interfaceAC );
  AnalyticalIdentity2D<T,T> phiU( zeroVelocity0 );
  AnalyticalIdentity2D<T,T> rho( rhov + (rhol-rhov)*phi );
  AnalyticalIdentity2D<T,T> tau( tauv + (taul-tauv)*phi );
  AnalyticalIdentity2D<T,T> pressure( one );

  sLatticeNS.defineField<descriptors::RHO>( superGeometry, 1, rho );
  sLatticeNS.defineField<descriptors::TAU_EFF>( superGeometry, 1, tau );
  sLatticeAC.defineField<descriptors::OLD_PHIU>( superGeometry, 1, phiU );

  sLatticeNS.defineRhoU( superGeometry, 1, pressure, zeroVelocity );
  sLatticeNS.iniEquilibrium( superGeometry, 1, pressure, zeroVelocity );
  sLatticeAC.defineRhoU( superGeometry, 1, phi, zeroVelocity );
  sLatticeAC.iniEquilibrium( superGeometry, 1, phi, zeroVelocity );

  sLatticeAC.addPostProcessor<stage::PostPostProcess>(meta::id<Helper>{});
  coupling.template setParameter<Coupling::SIGMA>(sigma);
  coupling.template setParameter<Coupling::W>(w);
  coupling.template setParameter<Coupling::TAUS>({tau_v,tau_l,tau_mobil});
  coupling.template setParameter<Coupling::RHOS>(rhos);

  {
    auto& communicator = sLatticeNS.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(1);
    communicator.exchangeRequests();
  }

  sLatticeAC.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());
  sLatticeAC.addPostProcessor<stage::PreCoupling>(meta::id<initialPsi>{});
  sLatticeAC.addPostProcessor<stage::IterativePostProcess>(meta::id<normGradPsi>{});
  sLatticeAC.addPostProcessor<stage::IterativePostProcess>(meta::id<psiEvolve>{});
  sLatticeAC.setParameter<psiEvolve::DELTAT>(0.5);
  sLatticeNS.setParameter<descriptors::OMEGA>( 1./tau_l );
  sLatticeAC.setParameter<descriptors::OMEGA>( 1./tau_mobil );
  sLatticeAC.setParameter<descriptors::INTERFACE_WIDTH>( w );

  {
    auto& communicator = sLatticeAC.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(1);
    communicator.requestField<STATISTIC>();
    communicator.requestField<PSI>();
    communicator.exchangeRequests();
  }
  {
    auto& communicator = sLatticeAC.getCommunicator(stage::IterativePostProcess());
    communicator.requestOverlap(1);
    communicator.requestField<PSI>();
    communicator.exchangeRequests();
  }

  sLatticeAC.executePostProcessors(stage::PreCoupling());
  signedDistanceFunction<stage::IterativePostProcess>(sLatticeAC,superGeometry,w);
  sLatticeNS.initialize();
  sLatticeAC.initialize();
  sLatticeNS.getCommunicator(stage::PreCoupling()).communicate();
  sLatticeAC.getCommunicator(stage::PreCoupling()).communicate();

  clout << "Prepare Lattice ... OK" << std::endl;
}

std::vector<T> error( SuperGeometry<T,2>& superGeometry,
         SuperLattice<T, ACDESCRIPTOR>& sLatticeAC,
         SuperLattice<T, NSDESCRIPTOR>& sLatticeNS,
         int Nx, int Ny, int phaseLength, T w )
{
  OstreamManager clout( std::cout,"error" );
  //T phiL1AbsError, phiL1RelError, phiL2AbsError, phiL2RelError, phiLinfAbsError, phiLinfRelError;
  T phiL2RelError, pL2RelError, uL2AbsError, phiL2AbsError;
  int tmp[]= { };
  T result[2]= { };

  AnalyticalConst2D<T,T> one ( 1. );
  AnalyticalConst2D<T,T> zero ( 0. );
  Vector<T,2> u(0., 0.);
  AnalyticalConst2D<T,T> zeroVelocity( u );
  SmoothIndicatorFactoredCuboid2D<T,T> interfaceAC_diffuse( {Nx/2., Ny/2.}, 0, phaseLength+1, w/2, 0, {0,0}, 0, -1. );
  AnalyticalIdentity2D<T,T> phi0_diffuse( one + interfaceAC_diffuse );
  Vector<T,2> extend( 2., Ny/4.-1 );
  Vector<T,2> origin1( 0., 0.);
  Vector<T,2> origin2( 0., Ny*3./4.+1);
  IndicatorCuboid2D<T> ind1( extend, origin1 );
  IndicatorCuboid2D<T> ind2( extend, origin2 );
  SmoothIndicatorCuboid2D<T,T> interfaceAC_sharp1( ind1, T(0) );
  SmoothIndicatorCuboid2D<T,T> interfaceAC_sharp2( ind2, T(0) );
  AnalyticalIdentity2D<T,T> phi0_sharp( zero + interfaceAC_sharp1 + interfaceAC_sharp2 );

  SuperLatticeDensity2D<T, ACDESCRIPTOR> phi( sLatticeAC );
  SuperLatticeDensity2D<T, NSDESCRIPTOR> p_hydro( sLatticeNS );
  SuperLatticeVelocity2D<T, NSDESCRIPTOR> velocity( sLatticeNS );
  auto indicatorF = superGeometry.getMaterialIndicator(1);

  //for order parameter
  if (Cahn_const) {
    SuperAbsoluteErrorL1Norm2D<T> absPhiErrorNormL1(phi, phi0_diffuse, indicatorF);
    absPhiErrorNormL1(result, tmp);
    //clout << "phi-L1-error(abs)=" << result[0];
    //phiL1AbsError = result[0];

    SuperRelativeErrorL1Norm2D<T> relPhiErrorNormL1(phi, phi0_diffuse, indicatorF);
    relPhiErrorNormL1(result, tmp);
    //clout << "; phi-L1-error(rel)=" << result[0] << std::endl;
    //phiL1RelError = result[0];

    SuperAbsoluteErrorL2Norm2D<T> absPhiErrorNormL2(phi, phi0_diffuse, indicatorF);
    absPhiErrorNormL2(result, tmp);
    //clout << "phi-L2-error(abs)=" << result[0];
    //phiL2AbsError = result[0];

    SuperRelativeErrorL2Norm2D<T> relPhiErrorNormL2(phi, phi0_diffuse, indicatorF);
    relPhiErrorNormL2(result, tmp);
    //clout << "; phi-L2-error(rel)=" << result[0] << std::endl;
    phiL2RelError = result[0];

    SuperAbsoluteErrorLinfNorm2D<T> absPhiErrorNormLinf(phi, phi0_diffuse, indicatorF);
    absPhiErrorNormLinf(result, tmp);
    //clout << "phi-Linf-error(abs)=" << result[0];
    //phiLinfAbsError = result[0];

    SuperRelativeErrorLinfNorm2D<T> relPhiErrorNormLinf(phi, phi0_diffuse, indicatorF);
    relPhiErrorNormLinf(result, tmp);
    //clout << "; phi-Linf-error(rel)=" << result[0] << std::endl;
    //phiLinfRelError = result[0];
  }
  else {
    SuperAbsoluteErrorL1Norm2D<T> absPhiErrorNormL1(phi, phi0_sharp, indicatorF);
    absPhiErrorNormL1(result, tmp);
    //clout << "phi-L1-error(abs)=" << result[0];
    //phiL1AbsError = result[0];

    SuperRelativeErrorL1Norm2D<T> relPhiErrorNormL1(phi, phi0_sharp, indicatorF);
    relPhiErrorNormL1(result, tmp);
    //clout << "; phi-L1-error(rel)=" << result[0] << std::endl;
    //phiL1RelError = result[0];

    SuperAbsoluteErrorL2Norm2D<T> absPhiErrorNormL2(phi, phi0_sharp, indicatorF);
    absPhiErrorNormL2(result, tmp);
    //clout << "phi-L2-error(abs)=" << result[0];
    phiL2AbsError = result[0];

    SuperRelativeErrorL2Norm2D<T> relPhiErrorNormL2(phi, phi0_sharp, indicatorF);
    relPhiErrorNormL2(result, tmp);
    //clout << "; phi-L2-error(rel)=" << result[0] << std::endl;
    //phiL2RelError = result[0];

    SuperAbsoluteErrorLinfNorm2D<T> absPhiErrorNormLinf(phi, phi0_sharp, indicatorF);
    absPhiErrorNormLinf(result, tmp);
    //clout << "phi-Linf-error(abs)=" << result[0];
    //phiLinfAbsError = result[0];

    SuperRelativeErrorLinfNorm2D<T> relPhiErrorNormLinf(phi, phi0_sharp, indicatorF);
    relPhiErrorNormLinf(result, tmp);
    //clout << "; phi-Linf-error(rel)=" << result[0] << std::endl;
    //phiLinfRelError = result[0];
  }
  //for pressure
  SuperAbsoluteErrorL1Norm2D<T> absPErrorNormL1(p_hydro, one, indicatorF);
  absPErrorNormL1(result, tmp);
  //clout << "p-L1-error(abs)=" << result[0];
  //pL1AbsError = result[0];

  SuperRelativeErrorL1Norm2D<T> relPErrorNormL1(p_hydro, one, indicatorF);
  relPErrorNormL1(result, tmp);
  //clout << "; p-L1-error(rel)=" << result[0] << std::endl;
  //pL1RelError = result[0];

  SuperAbsoluteErrorL2Norm2D<T> absPErrorNormL2(p_hydro, one, indicatorF);
  absPErrorNormL2(result, tmp);
  //clout << "p-L2-error(abs)=" << result[0];
  //pL2AbsError = result[0];

  SuperRelativeErrorL2Norm2D<T> relPErrorNormL2(p_hydro, one, indicatorF);
  relPErrorNormL2(result, tmp);
  //clout << "; p-L2-error(rel)=" << result[0] << std::endl;
  pL2RelError = result[0];

  SuperAbsoluteErrorLinfNorm2D<T> absPErrorNormLinf(p_hydro, one, indicatorF);
  absPErrorNormLinf(result, tmp);
  //clout << "p-Linf-error(abs)=" << result[0];
  //pLinfAbsError = result[0];

  SuperRelativeErrorLinfNorm2D<T> relPErrorNormLinf(p_hydro, one, indicatorF);
  relPErrorNormLinf(result, tmp);
  //clout << "; p-Linf-error(rel)=" << result[0] << std::endl;
  //pLinfRelError = result[0];

  //for velocity
  SuperAbsoluteErrorL1Norm2D<T> absUErrorNormL1(velocity, zeroVelocity, indicatorF);
  absUErrorNormL1(result, tmp);
  //clout << "u-L1-error(abs)=" << result[0];
  //uL1AbsError = result[0];

  SuperRelativeErrorL1Norm2D<T> relUErrorNormL1(velocity, zeroVelocity, indicatorF);
  relUErrorNormL1(result, tmp);
  //clout << "; u-L1-error(rel)=" << result[0] << std::endl;
  //uL1RelError = result[0];

  SuperAbsoluteErrorL2Norm2D<T> absUErrorNormL2(velocity, zeroVelocity, indicatorF);
  absUErrorNormL2(result, tmp);
  //clout << "u-L2-error(abs)=" << result[0];
  uL2AbsError = result[0];

  SuperRelativeErrorL2Norm2D<T> relUErrorNormL2(velocity, zeroVelocity, indicatorF);
  relUErrorNormL2(result, tmp);
  //clout << "; u-L2-error(rel)=" << result[0] << std::endl;
  //uL2RelError = result[0];

  SuperAbsoluteErrorLinfNorm2D<T> absUErrorNormLinf(velocity, zeroVelocity, indicatorF);
  absUErrorNormLinf(result, tmp);
  //clout << "u-Linf-error(abs)=" << result[0];
  //uLinfAbsError = result[0];

  SuperRelativeErrorLinfNorm2D<T> relUErrorNormLinf(velocity, zeroVelocity, indicatorF);
  relUErrorNormLinf(result, tmp);
  //clout << "; u-Linf-error(rel)=" << result[0] << std::endl;
  //uLinfRelError = result[0];

  if (Cahn_const) return {phiL2RelError,pL2RelError,uL2AbsError};
  else return {phiL2AbsError,pL2RelError,uL2AbsError};
}

T getResults( SuperLattice<T,NSDESCRIPTOR>& sLatticeNS,
              SuperLattice<T,ACDESCRIPTOR>& sLatticeAC,
              int iT, SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer,
              UnitConverter<T,NSDESCRIPTOR> converter, int Ny )
{
  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter2D<T> vtmWriter( "allenCahnFlatInterfaceEOC2d" );
  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, NSDESCRIPTOR> cuboid( sLatticeNS );
    SuperLatticeRank2D<T, NSDESCRIPTOR> rank( sLatticeNS );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Get statistics
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();
    sLatticeNS.getStatistics().print( iT, converter.getPhysTime(iT) );
    sLatticeAC.getStatistics().print( iT, converter.getPhysTime(iT) );
  }

  T phiInterpol_total = 0.;
  // Writes the VTK files
  if ( iT%vtkIter==0 ) {
    SuperLatticeDensity2D<T, NSDESCRIPTOR> p_hydro_L( sLatticeNS );
    AnalyticalConst2D<T,T> ConversionPressure_( converter.getConversionFactorPressure() );
    SuperLatticeFfromAnalyticalF2D<T, NSDESCRIPTOR> ConversionPressure(ConversionPressure_, sLatticeNS);
    SuperIdentity2D<T,T> p_hydro( p_hydro_L * ConversionPressure );
    p_hydro.getName() = "p_hydro";
    SuperLatticeDensity2D<T, ACDESCRIPTOR> phi( sLatticeAC );
    phi.getName() = "phi";
    SuperLatticeExternalScalarField2D<T, ACDESCRIPTOR, PSI> psi( sLatticeAC );
    psi.getName() = "psi";
    SuperLatticeExternalScalarField2D<T, NSDESCRIPTOR, RHO> rho_L( sLatticeNS );
    AnalyticalConst2D<T,T> ConversionDensity_( converter.getConversionFactorDensity() );
    SuperLatticeFfromAnalyticalF2D<T, NSDESCRIPTOR> ConversionDensity(ConversionDensity_, sLatticeNS);
    SuperIdentity2D<T,T> rho( rho_L * ConversionDensity );
    rho.getName() = "rho";
    SuperLatticeVelocity2D<T, NSDESCRIPTOR> velocity_L( sLatticeNS );
    AnalyticalConst2D<T,T> ConversionVelocity_( converter.getConversionFactorVelocity() );
    SuperLatticeFfromAnalyticalF2D<T, NSDESCRIPTOR> ConversionVelocity(ConversionVelocity_, sLatticeNS);
    SuperIdentity2D<T,T> velocity( velocity_L * ConversionVelocity );
    velocity.getName() = "u";

    vtmWriter.addFunctor( p_hydro );
    vtmWriter.addFunctor( phi );
    vtmWriter.addFunctor( psi );
    vtmWriter.addFunctor( rho );
    vtmWriter.addFunctor( velocity );
    vtmWriter.write( iT );

    AnalyticalFfromSuperF2D<T,T> interpolPhi( phi, true, 1);
    int count_phi =0;
    for ( int i=0; i<Ny/2; ++i ) {
      T position[2] = { 1, static_cast<T>(i) };
      T phiInterpol = 0.;
      interpolPhi(&phiInterpol, position);
      //clout << "phi_int:" << phiInterpol << std::endl;
      if(phiInterpol <= 0.75 && phiInterpol > 0.25){
        count_phi++;
        phiInterpol_total += phiInterpol;
      }
    }
    //clout << "Count the phi between 0.25 and 0.75: " << count_phi << std::endl;
    phiInterpol_total /= count_phi;
  }
  return phiInterpol_total;
}

T helper( SuperLattice<T,ACDESCRIPTOR>& sLatticeAC,
          SuperGeometry<T,2>& superGeometry )
{
  sLatticeAC.executePostProcessors(stage::PostPostProcess());
  int in[2];
  T out[1];
  SuperLatticeExternalScalarField2D<T, ACDESCRIPTOR, TOP> top( sLatticeAC );
  SuperLatticeExternalScalarField2D<T, ACDESCRIPTOR, BOTTOM> bottom( sLatticeAC );
  SuperAverage2D<T> top_total_(top, superGeometry, 1);
  SuperAverage2D<T> bottom_total_(bottom, superGeometry, 1);
  top_total_(out, in);
  T top_total = out[0];
  bottom_total_(out, in);
  T bottom_total = out[0];
  return top_total/bottom_total;
}

std::vector<T> simulate( int Ny, int phaseLength, T w, T sigma )
{
  OstreamManager clout( std::cout,"main" );
  // === 1st Step: Initialization ===
  UnitConverterFromResolutionAndRelaxationTime<T,NSDESCRIPTOR> converter(
    int   {Ny},                     // resolution
    (T)   tau_l,                    // lattice relaxation time
    (T)   charPhysLength,           // charPhysLength: reference length of simulation geometry
    (T)   Re/charPhysLength*viscosityH2O,   // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   viscosityH2O,             // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   DeltaRho                  // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  converter.print();
  // === 2nd Step: Prepare Geometry ===
  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidGeometry2D<T> cGeometry({0, 0}, 1, {Nx, Ny}, noOfCuboids );
  // set periodic boundaries to the domain
  cGeometry.setPeriodicity({ true, true });
  // Instantiation of loadbalancer
  HeuristicLoadBalancer<T> loadBalancer( cGeometry );
  loadBalancer.print();
  // Instantiation of superGeometry
  SuperGeometry<T,2> superGeometry( cGeometry,loadBalancer,2 );
  prepareGeometry( superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T,NSDESCRIPTOR> sLatticeNS( superGeometry );
  SuperLattice<T,ACDESCRIPTOR> sLatticeAC( superGeometry );
  SuperLatticeCoupling coupling(
      AllenCahnPostProcessor{},
      names::NavierStokes{}, sLatticeNS,
      names::Component1{}, sLatticeAC);
  prepareLattice( sLatticeNS, sLatticeAC, coupling, converter, superGeometry, Ny, phaseLength, w, sigma );

  // === 4th Step: Main Loop with Timer ===
  int iT = 0;
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( maxIter, superGeometry.getStatistics().getNvoxel() );
  timer.start();
  std::vector<T> output = {static_cast<T>(Ny),0,0,0};
  T phi_total_old = 1.;
  T min_u_err = 1.;
  int count = 0;

  T tobo = helper(sLatticeAC,superGeometry);
  coupling.template setParameter<Coupling::NONLOCALITY>(tobo);
  coupling.execute();

  for ( iT=0; iT<=maxIter; ++iT ) {
    // Collide and stream (and coupling) execution
    sLatticeNS.collideAndStream();
    sLatticeAC.collideAndStream();

    sLatticeAC.executePostProcessors(stage::PreCoupling());
    signedDistanceFunction<stage::IterativePostProcess>(sLatticeAC,superGeometry,w);
    sLatticeNS.getCommunicator(stage::PreCoupling()).communicate();
    sLatticeAC.getCommunicator(stage::PreCoupling()).communicate();

    tobo = helper(sLatticeAC,superGeometry);
    coupling.template setParameter<Coupling::NONLOCALITY>(tobo);
    coupling.execute();

    // Computation and output of the results
    T results = getResults( sLatticeNS, sLatticeAC, iT, superGeometry, timer, converter, Ny );
    if ( std::isnan( sLatticeNS.getStatistics().getAverageEnergy() ) ) {
      break;
    }
    T phi_total = results;
    std::vector<T> err = error(superGeometry, sLatticeAC, sLatticeNS, Nx, Ny, phaseLength, w);
    T u_err = err[2];
    if(iT % (vtkIter) == 0) {
      clout << "current error phi: " << err[0] << std::endl;
      clout << "current error p: " << err[1] << std::endl;
      clout << "current error u: " << err[2] << std::endl;
      if(min_u_err > u_err){
        min_u_err = u_err;
        count = 0;
      }
      else count++;
      clout << "counter: " << count << std::endl;
      if(abs(phi_total - phi_total_old)/abs(phi_total) < 1e-9/* && count >= 5*/){
        clout << "Converged at: " << abs(phi_total - phi_total_old)/abs(phi_total) << std::endl;
        if (Cahn_const) output[1] = err[0];
        else output[1] = err[0]/(Nx*Ny);
        output[2] = err[1];
        output[3] = 1e-14;
        break;
      }
      phi_total_old = phi_total;
    }
    if(iT*converter.getPhysDeltaT() >= 0.1) break;
  }
  timer.stop();
  timer.printSummary();
  return output;
}

int main( int argc, char *argv[] )
{
  olbInit( &argc, &argv );
  if (argc > 1) {
    Ny = atof(argv[1]);
  }
  std::ofstream outfile;
  outfile.open ("EOC.dat");

  //Update the geometry, lattice, and other dependent variables for each Ny
  for (size_t i=0; i<4; ++i) {
    singleton::directories().setOutputDir("./tmp/sub" + std::to_string(i) + "/");
    std::vector<T> out = simulate(Ny*util::pow(2,i),
                                  phaseLength*util::pow(2,i),
                                  w*util::pow(2,i*(Cahn_const)),
                                  sigma/util::pow(2,i));
    outfile << out[0] << "," << out[1] << "," << out[2] << "," << out[3] << "\n";
  }
  outfile.close();
}
