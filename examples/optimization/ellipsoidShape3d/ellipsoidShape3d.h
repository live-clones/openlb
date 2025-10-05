/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Felix Schuhmann, Shota Ito
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

// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include "olb.h"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::opti;

using S = FLOATING_POINT_TYPE;
using U = util::ADf<S,2>;
using DESCRIPTOR = D3Q19<>;

const int N = 20;                    // resolution of the model
const S Re = 50;                     // Reynolds number
const S heightChannel = 0.4;         // height of the outer channel
const S lengthChannel = 1.2;         // length of the outer channel
const S epsilon = 5.;                // in cell layers
const S maxPhysT = 50.;              // Max. simulation time in s, SI unit
const S rampStartT = 30.;
const S rampUpdateT = 0.01;
const S latticeRelaxationTime = 0.55;
const S charPhysViscosity = 0.001;
const S physDensity = 1.0;

const S centerEllipsoidX = 0.5;
const S centerChannelYZ = heightChannel / 2.0;

template<typename T>
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter, IndicatorF3D<T>& channel,
                      SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0, 1, channel );

  T delta = converter.getPhysDeltaX()/2.;
  auto min = channel.getMin();
  auto max = channel.getMax();
  Vector<T,3> origin { min[0]-delta, min[1]-delta, min[2]-delta};
  Vector<T,3> extend { delta * 2.0, max[1] + delta*2.0, max[2] + delta*2. };
  IndicatorCuboid3D<T> inflow(extend, origin);
  superGeometry.rename(1,3,inflow);

  Vector<T,3> origin2 { max[0]-delta, min[1]-delta, min[2]-delta};
  Vector<T,3> extend2 { delta * 2.0, max[1] + delta*2.0, max[2] + delta*2. };
  IndicatorCuboid3D<T> outflow(extend2, origin2);
  superGeometry.rename(1,4,outflow);

  superGeometry.clean();
  superGeometry.checkForErrors();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

template<typename T>
void prepareLattice( SuperLattice<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperGeometry<T,3>& superGeometry,
                     std::shared_ptr<AnalyticalF3D<T,T>> smoothEllipsoid)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Bulk dynamics for all materials
  auto bulkIndicator = superGeometry.getMaterialIndicator( { 1,3,4,5 } );
  sLattice.template defineDynamics<PorousBGKdynamics<T,DESCRIPTOR>>(bulkIndicator);

  // boundary conditions
  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);

  // Set required porosities
  std::shared_ptr<AnalyticalF3D<T,T>> one = std::make_shared<AnalyticalConst3D<T,T>>(1.);
  std::shared_ptr<AnalyticalF3D<T,T>> ellipsoidField = one - smoothEllipsoid;

  // Initialize porosity
  sLattice.template defineField<POROSITY>( bulkIndicator, *ellipsoidField );

  // Initialize rho and u
  AnalyticalConst3D<T,T> rhoF(1);
  Vector<T,3> u_0 {0.0, 0.0, 0.0};
  AnalyticalConst3D<T,T> u(u_0);

  sLattice.defineRhoU(bulkIndicator, rhoF, u);
  sLattice.iniEquilibrium(bulkIndicator, rhoF, u);

  sLattice.template setParameter<descriptors::OMEGA>(omega);

  sLattice.initialize();
  clout << "Prepare Lattice ... OK" << std::endl;
}

template<typename T>
void setBoundaryValues( SuperLattice<T,DESCRIPTOR>& sLattice,
                        UnitConverter<T,DESCRIPTOR> const& converter,
                        SuperGeometry<T,3>& superGeometry, int iT )
{
  int iTmaxStart = converter.getLatticeTime( rampStartT );
  int iTupdate = converter.getLatticeTime( rampUpdateT );

  if (iT%iTupdate == 0 && iT <= iTmaxStart) {
    PolynomialStartScale<T,int> startScale(iTmaxStart, T(1));
    int iTvec[1] = {iT};
    T frac[1] = {};
    startScale(frac, iTvec);
    Vector<T,3> u_free {frac[0]*converter.getCharLatticeVelocity(), 0.0, 0.0};
    AnalyticalConst3D<T,T> u(u_free);
    sLattice.defineU(superGeometry, 3, u);

    sLattice.template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

template<typename T>
T objective_details( SuperLattice<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperGeometry<T,3>& superGeometry )
{
  OstreamManager clout( std::cout,"objective_details" );
  T result = 0.0;
  Vector<T,3> center( centerEllipsoidX, centerChannelYZ, centerChannelYZ );
  IndicatorCuboid3D<T> integrationDomain(1.0, 0.4, 0.4, center);
  SuperIndicatorFfromIndicatorF3D<T> iDomain(integrationDomain, superGeometry);

  // Viscous Dissipation
  SuperLatticePhysDissipation3D<T,DESCRIPTOR> viscous_dissipation( sLattice, converter );
  SuperIntegral3D<T,T> dissipationIntegral1(viscous_dissipation, iDomain);
  T vDissipation[1] = {0.};
  int input1[1];
  dissipationIntegral1( vDissipation, input1 );
  clout << "Viscous dissipation = " << vDissipation[0] << std::endl;
  result = vDissipation[0];

  // Porous Dissipation
  SuperLatticeFfromCallableF<T,DESCRIPTOR> porous_Dissipation(sLattice, [&](T* output, auto cell){
    T uTemp[DESCRIPTOR::d];
    cell.computeU(uTemp);
    const T porosity = cell.template getField<descriptors::POROSITY>();

    const T invPermeability = projection::porosityToInvPermeability(porosity, converter);
    const T uNormSq = util::euklidN2(uTemp, DESCRIPTOR::d);
    output[0] = converter.getPhysViscosity() * invPermeability * uNormSq;
  });
  SuperIntegral3D<T,T> dissipationIntegral2(porous_Dissipation, iDomain);
  T pDissipation[1] = {0.};
  dissipationIntegral2( pDissipation, input1 );
  clout << "Porous dissipation = " << pDissipation[0] << std::endl;
  result += pDissipation[0];
  return result;
}

// Computes and visualizes Lattice fields, computes dissipation
template<typename T>
void getResults( SuperLattice<T, DESCRIPTOR>& sLattice,
                 const UnitConverter<T,DESCRIPTOR>& converter, int iT,
                 SuperGeometry<T,3>& superGeometry,
                 util::Timer<T>& timer)
{
  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter3D<T> vtmWriter( "ellipsoid3dOpti" );
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T,DESCRIPTOR> pressure( sLattice, converter );
  SuperLatticePorosity3D<T,DESCRIPTOR> porosity( sLattice );
  SuperLatticePhysDissipation3D<T,DESCRIPTOR> viscous_dissipation( sLattice, converter );
  SuperLatticeFfromCallableF<T,DESCRIPTOR> porous_dissipation(sLattice, [&](T* output, auto cell){
    T uTemp[DESCRIPTOR::d];
    cell.computeU(uTemp);
    const T porosity = cell.template getField<descriptors::POROSITY>();

    const T invPermeability = projection::porosityToInvPermeability(porosity, converter);
    const T uNormSq = util::euklidN2(uTemp, DESCRIPTOR::d);
    output[0] = converter.getPhysViscosity() * invPermeability * uNormSq;
  });

  viscous_dissipation.getName() = "viscous dissipation";
  porous_dissipation.getName() = "porous dissipation";

  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );
  vtmWriter.addFunctor( porosity );
  vtmWriter.addFunctor( viscous_dissipation );
  vtmWriter.addFunctor( porous_dissipation );

  const int vtkIter  = converter.getLatticeTime( maxPhysT / 20 );
  const int statIter = converter.getLatticeTime( maxPhysT / 10 );

  if ( iT==0 ) {
    vtmWriter.createMasterFile();
  }

  // Writes output on the console
  if ( iT%statIter == 0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
    objective_details(sLattice, converter, superGeometry);
  }

  // Writes the vtk files, currently only for last simulation
  if ( iT%vtkIter == 0 ) {
    vtmWriter.write( iT );
  }
}

// Main simulation
template<typename T> // If getAllResults = false, only the objective is evaluated, no lattice fields will be computed and visualized
T simulateEllipsoid3D( std::vector<T> controls )
{
  OstreamManager clout( std::cout,"simulateEllipsoid" );

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> converter(
    int {N},                                      // resolution: number of voxels per charPhysL
    (T)   latticeRelaxationTime,                  // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   heightChannel,                          // charPhysLength: reference length of simulation geometry
    (T)   charPhysViscosity * Re / heightChannel, // physVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   charPhysViscosity,                      // charPhysViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   physDensity                             // physDensity: physical density in __kg / m^3__
  );
  converter.print();

  const T radiusEllipsoidY = controls[0];
  const T radiusEllipsoidZ = controls[1];
  const T volumeEllipsoid  = 0.001;
  const T radiusEllipsoidX = 0.75 * volumeEllipsoid / (M_PI*radiusEllipsoidY*radiusEllipsoidZ);

  Vector<T,3> center( centerEllipsoidX, centerChannelYZ, centerChannelYZ ); // Fixed center of ellipsoid
  Vector<T,3> radius( radiusEllipsoidX, radiusEllipsoidY, radiusEllipsoidZ );
  clout << "Radius X: " << radiusEllipsoidX  << std::endl;
  clout << "Radius Y: " << radiusEllipsoidY  << std::endl;
  clout << "Radius Z: " << radiusEllipsoidZ  << std::endl;
  Vector<T,3> radius_plus_eps = radius;
  for (int i=0; i < 3; ++i) {
    radius_plus_eps[i] += epsilon * converter.getConversionFactorLength();
  }

  // Indicator ellipsoid
  std::shared_ptr<IndicatorF3D<T>> ellipsoid = std::make_shared<IndicatorEllipsoid3D<T>>( center, radius_plus_eps );

  // Smooth ellipsoid to set porosity field
  std::shared_ptr<AnalyticalF3D<T,T>> smoothEllipsoid = std::make_shared<SmoothIndicatorSigmoidEllipsoid3D<T,T,false>>( center, radius, epsilon * converter.getConversionFactorLength() );

  // === Initial 2nd Step: Prepare Geometry ===
  // Construct channel
  Vector<T,3> origin( 0., 0., 0. );
  Vector<T,3> extent( lengthChannel, heightChannel, heightChannel );
  IndicatorCuboid3D<T> channel( extent, origin );

  // Instantiation of a cuboidGeometry with weights
  CuboidDecomposition3D<T> cuboidGeometry( channel, converter.getConversionFactorLength(), singleton::mpi().getSize() );

  // Periodic boundaries
  cuboidGeometry.setPeriodicity({false, true, true});

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,3> superGeometry( cuboidGeometry, loadBalancer );
  prepareGeometry( converter, channel, superGeometry);

  // === Initial 3rd Step: Prepare Lattice ===
  SuperLattice<T,DESCRIPTOR> sLattice( converter, superGeometry );
  prepareLattice( sLattice, converter, superGeometry, smoothEllipsoid);

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation" << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  auto maxiT = converter.getLatticeTime( maxPhysT );
  timer.start();

  for (std::size_t iT = 0; iT < maxiT; ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions === can be skipped for forced flow type
    setBoundaryValues( sLattice, converter, superGeometry, iT );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer );
  }
  timer.stop();
  timer.printSummary();
  return util::pow(objective_details(sLattice, converter, superGeometry), 2.);
}
