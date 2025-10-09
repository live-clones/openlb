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

/* ellipsoidShape3d.cpp:
 * This example examines a steady flow past an ellipsoid placed inside a rectangular channel.
 * The half-axes are alligned along the x, y and z coordinate directions and periodic boundary treatment
 * is applied in y and z directions.
 *
 * Contains:
 * - Standard channel flow simulation
 * - Flow simulation while computing derivatives with respect to ellipsoid radii in y and z direction using ADf
 * - Shape optimization using LBFGS:
 *   radii in y and z direction are optimized to minimize L2-Norm of dissipation over the entire channel
 *   radius in x direction is adjusted to enforce fixed ellipsoid volume
 */

// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include "olb.h"

using namespace olb;
using namespace olb::opti;

using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q19<>
>;
using MyADfCase = Case<
  NavierStokes, Lattice<util::ADf<double,2>, descriptors::D3Q19<>
>;
using MyOptiCase = OptiCaseADf<
  Controlled, MyCase,
  Derivatives, MyADfCase
>;

namespace olb::parameters{

struct CONTROLS : public descriptors::FIELD_BASE<2>();

}

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

template <typename PARAMETERS>
auto createMesh(PARAMETERS& parameters) {
  using T = PARAMETERS::value_t;
  const Vector extent = parameters.template get<parameters::DOMAIN_EXTENT>();
  const Vector origin{0., 0., 0.};
  IndicatorCuboid3D<T> channel( extent, origin );

  Mesh<T,3> mesh(channel,
                 extent[0] / parameters.template get<parameters::RESOLUTION>(),
		 singleton::mpi().getSize());
  mesh.setOverlap(parameters.template get<parameters::OVERLAP>());
  return mesh;
}

template<typename CASE>
void prepareGeometry(CASE& myCase)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;
  using T = CASE::value_t;
  auto& parameters = myCase.getParmeters();
  auto& superGeometry = myCase.getGeometry();
  const T delta = parameters.get<parameters::DOMAIN_EXTENT>()[0] / parameters.get<parameters::RESOLUTION>() / 2.0;

  const Vector extent = parameters.template get<parameters::DOMAIN_EXTENT>();
  const Vector origin{0., 0., 0.};
  IndicatorCuboid3D<T> channel( extent, origin );
  superGeometry.rename( 0, 1, channel );

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

template<typename CASE>
void prepareLattice( SuperLattice<T,DESCRIPTOR>& sLattice,
                     SuperGeometry<T,3>& superGeometry,
                     std::shared_ptr<AnalyticalF3D<T,T>> smoothEllipsoid)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  sLattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    parameters.get<parameters::RESOLUTION>()          // resolution: number of voxels per charPhysL
    parameters.get<parameters::RELAXATION_TIME>(),    // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    parameters.get<parameters::DOMAIN_EXTENT>()[0],   // charPhysLength: reference length of simulation geometry
    parameters.get<parameters::PHYS_CHAR_VELOCITY>(), // physVelocity: maximal/highest expected velocity during simulation in __m / s__
    parameters.get<parameters::PHYS_CHAR_VISCOSITY>(),// charPhysViscosity: physical kinematic viscosity in __m^2 / s__
    parameters.get<parameters::PHYS_CHAR_DENSITY>()   // physDensity: physical density in __kg / m^3__
  );
  auto& converter = sLattice.getUnitConverter();
  converter.print();
  const T omega = converter.getLatticeRelaxationFrequency();

  // Bulk dynamics for all materials
  auto bulkIndicator = superGeometry.getMaterialIndicator( { 1,3,4,5 } );
  sLattice.template defineDynamics<PorousBGKdynamics<T,DESCRIPTOR>>(bulkIndicator);

  // boundary conditions
  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);

  // Set required porosities
  Vector radius = parameters.get<parameters::ELLIPSOID_RADII>();
  Vector radius_plus_eps = radius;
  const T eps = parameters.get<parameters::SMOOTH_LAYER_THICKNESS>() *
	        parameters.get<parameters::EXTENT_DOMAIN>()[0] /
	        parameters.get<parameters::RESOLUTION>(); 
  for (int i=0; i < 3; ++i) {
    radius_plus_eps += eps;
  }
  std::shared_ptr<IndicatorF3D<T>> ellipsoid = std::make_shared<IndicatorEllipsoid3D<T>>( parameters.get<parameters::ELLIPSOID_POS>(), radius_plus_eps );

  // Smooth ellipsoid to set porosity field
  std::shared_ptr<AnalyticalF3D<T,T>> smoothEllipsoid = std::make_shared<SmoothIndicatorSigmoidEllipsoid3D<T,T,false>>( parameters.get<parameters::ELLIPSOID_RADII>(), radius, eps);
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

template<typename CASE>
void setBoundaryValues( SuperLattice<T,DESCRIPTOR>& sLattice,
                        SuperGeometry<T,3>& superGeometry, int iT )
{
  auto& converter = sLattice.getUnitConverter();
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

template<typename CASE>
T objective_details( SuperLattice<T,DESCRIPTOR>& sLattice,
                     SuperGeometry<T,3>& superGeometry )
{
  OstreamManager clout( std::cout,"objective_details" );
  auto& converter = sLattice.getUnitConverter();
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
template<typename CASE>
void getResults( SuperLattice<T, DESCRIPTOR>& sLattice,
                 int iT,
                 SuperGeometry<T,3>& superGeometry,
                 util::Timer<T>& timer)
{
  OstreamManager clout( std::cout,"getResults" );
  auto& converter = sLattice.getUnitConverter();

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

int main( int argc, char* argv[] )
{
  initialize( &argc, &argv );
  using T = S;
  OstreamManager clout( std::cout,"simulateEllipsoid" );
  MyCase::ParametersD myCaseParametersD;
  {
    using namespace parameters;
    myCaseParametersD.set<DOMAIN_EXTENT>({0.4, 1.2, 1.2});
    myCaseParametersD.set<RESOLUTION>(20);
    myCaseParametersD.set<RELAXATION_TIME>(0.55);
    myCaseParametersD.set<PHYS_CHAR_VISCOSITY>(0.001);
    myCaseParametersD.set<REYNOLDS>(50);
    myCaseParametersD.set<PHYS_CHAR_VELOCITY>([&]() {
      return myCaseParameters.get<PHYS_CHAR_VISCOSITY>() *
             myCaseParameters.get<REYNOLDS>() /
	     myCaseParameters.get<EXTENT_DOMAIN>()[0];
    });
    myCaseParametersD.set<PHYS_CHAR_DENSITY>(1.0);
    myCaseParametersD.set<CONTROLS>({0.05, 0.05});
    myCaseParametersD.set<ELLIPSOID_VOLUME>(0.001);
    myCaseParametersD.set<SMOOTH_LAYER_THICKNESS>(5);
    myCaseParametersD.set<ELLIPSOID_RADII>([&]() {
      const T radiusZ = 0.75*myCaseParametersD.get<ELLIPSOID_VOLUME>() /
	      (M_PI*myCaseParametersD.get<CONTROLS>()[0]*myCaseParametersD.get<CONTROLS>()[1]);
      return {myCaseParametersD.get<CONTROLS>()[0],
              myCaseParametersD.get<CONTROLS>()[1],
	      radiusZ};
    });
    myCaseParametersD.set<ELLIPSOID_POS>([&]() {
      return {0.5,
              myCaseParametersD.get<DOMAIN_EXTENT>()[1] / 2.0,
	      myCaseParametersD.get<DOMAIN_EXTENT>()[2] / 2.0};
    });
  }

  auto mesh = createMesh(myCaseParametersD);

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

  // Instantiation of a superGeometry
  prepareGeometry( channel, superGeometry);

  // === Initial 3rd Step: Prepare Lattice ===
  SuperLattice<T,DESCRIPTOR> sLattice( superGeometry );
  prepareLattice( sLattice, superGeometry, smoothEllipsoid);

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation" << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  auto maxiT = converter.getLatticeTime( maxPhysT );
  timer.start();

  for (std::size_t iT = 0; iT < maxiT; ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions === can be skipped for forced flow type
    setBoundaryValues( sLattice, superGeometry, iT );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, iT, superGeometry, timer );
  }
  timer.stop();
  timer.printSummary();
  clout << util::pow(objective_details(sLattice, superGeometry), 2.) << std::endl;

  /*if constexpr (false) {
    std::vector<U> radiusEllipsoidYZ({.05, .05});
    radiusEllipsoidYZ[0].setDiffVariable(0);
    radiusEllipsoidYZ[1].setDiffVariable(1);
    clout << simulateEllipsoid3D<U>( radiusEllipsoidYZ ) << std::endl;
  }

  if constexpr (true){
    solver::OptiCaseAD<S,2,std::vector> optiCase(
      objective<S>,
      objective<U>);
    OptimizerLBFGS<S,std::vector<S>> optimizer(
    2, 1.e-16, 10, .01, 10, "StrongWolfe", 20, 1.e-4, true, "", "log",
    true, 0.19, true, 0.01, false, 0., true,
    {OptimizerLogType::value, OptimizerLogType::control, OptimizerLogType::derivative});

    std::vector<S> startValue(0.08, 0.08);
    optimizer.setControl(startValue);
    optimizer.optimize(optiCase);
  }*/
}
