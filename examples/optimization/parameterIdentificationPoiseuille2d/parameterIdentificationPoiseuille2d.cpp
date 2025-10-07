/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2007, 2012, 2022 Jonas Latt, Mathias J. Krause
 *  Vojtech Cvrcek, Peter Weisbrod, Julius Jeßberger
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

/** \file A simple two-dimensional fluid flow optimization problem is solved.
 * The setup is a planar channel flow, similar to the example poiseuille2d.
 * For a given pressure drop, the steady flow is simulated and the mass flow
 * rate is computed. In the optimization problem, the inlet pressure is
 * determined s.t. a pre-defined mass flow rate is achieved.
 */

// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::opti;

using S = double;
using U = util::ADf<S,1>;
using DESCRIPTOR = descriptors::D2Q9<>;

const int N = 50;             // resolution
const S lx  = 2.;             // length of the channel
const S ly  = 1.;             // height of the channel
const S Re = 10.;             // Reynolds number
const S maxPhysT = 30.;       // max. simulation time in s, SI unit
const S physInterval = 0.25;  // interval for the convergence check in s
const S residuum = 1e-9;      // residuum for the convergence check
const S wantedMassFlow = 0.00026159;
S inletPressure = 0.000659;

using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D2Q9<>>
>;
using MyOptiCase = OptiCaseFDQ<
  Controlled, MyCase
>;

template <typename PARAMETERS>
auto createMesh(PARAMETERS& parameters) {
  using T = PARAMETERS::value_t;
  const Vector<T,2> extend( lx, ly );
  const Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid( extend, origin );

  Mesh<T,2> mesh(cuboid, ly / N, singleton::mpi().getSize());
  mesh.setOverlap(parameters.template get<parameters::OVERLAP>());
  return mesh; 
}

template<typename CASE>
void prepareGeometry(CASE& myCase)
{
  using T = CASE::value_t;
  auto& superGeometry = myCase.getGeometry();

  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,{1,1} );

  const T physSpacing = ly / N;
  const Vector<T,2> extend    {physSpacing / T(2), ly};
  const Vector<T,2> originIn  {-physSpacing / T(4), 0};
  const Vector<T,2> originOut {lx-physSpacing / T(4), 0};

  IndicatorCuboid2D<T> inflow( extend, originIn );
  superGeometry.rename( 2,3,1,inflow );

  IndicatorCuboid2D<T> outflow( extend, originOut );
  superGeometry.rename( 2,4,1,outflow );

  superGeometry.clean(false);
  superGeometry.innerClean(false);
  superGeometry.checkForErrors(false);
}

template<typename CASE>
void prepareLattice(CASE& myCase)
{
  using T = CASE::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  sLattice.template setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    int {N},     // resolution: number of voxels per charPhysL
    (T)   0.8,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   1,     // charPhysLength: reference length of simulation geometry
    (T)   1,     // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1./Re, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0    // physDensity: physical density in __kg / m^3__
  );
  auto& converter = sLattice.getUnitConverter();
  const T omega = converter.getLatticeRelaxationFrequency();

  sLattice.template defineDynamics<BGKdynamics<T,DESCRIPTOR>>(superGeometry, 1);
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);

  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);

  sLattice.template setParameter<descriptors::OMEGA>(omega);
}

template <typename CASE>
void setInitialValues(CASE& myCase) {
  using T = CASE::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  auto& converter = sLattice.getUnitConverter();
 
  // Initial conditions
  AnalyticalLinear2D<T,T> rho(
    - inletPressure / lx * descriptors::invCs2<T,DESCRIPTOR>(),
    0,
    inletPressure * descriptors::invCs2<T,DESCRIPTOR>() + 1 );

  const T Lx = converter.getLatticeLength( lx ) - 1;
  const T Ly = converter.getLatticeLength( ly ) - 1;
  const T maxVelocity = inletPressure * Ly * Ly
    / (8.0 * converter.getLatticeViscosity() * Lx);
  const T radius = T(0.5) * (ly - converter.getPhysDeltaX());
  const std::vector<T> axisPoint { lx/T(2), ly/T(2) };
  const std::vector<T> axisDirection { 1, 0 };
  Poiseuille2D<T> u( axisPoint, axisDirection, maxVelocity, radius );

  const std::vector<T> zero(2, T());
  AnalyticalConst2D<T, T> u0(zero);

  sLattice.defineRhoU(superGeometry, 0, rho, u0);
  sLattice.iniEquilibrium(superGeometry, 0, rho, u0);

  const auto domain = superGeometry.getMaterialIndicator({1,2,3,4});
  sLattice.defineRhoU( domain, rho, u );
  sLattice.iniEquilibrium( domain, rho, u );

  sLattice.initialize();
}

template<typename CASE>
void getResults(CASE& myCase, std::size_t iT, util::Timer<typename CASE::value_t>& timer, bool hasConverged) {
  using T = CASE::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& converter = sLattice.getUnitConverter();

  SuperVTMwriter2D<T> vtmWriter( "poiseuille2d" );
  const bool lastTimeStep
   = ( hasConverged || (iT + 1 == converter.getLatticeTime( maxPhysT )) );
  const int vtmIter  = converter.getLatticeTime( maxPhysT/20. );

  SuperLatticePhysVelocity2D<T,DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure( sLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  if (iT == 0) {
      sLattice.communicate();
      SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice );
      SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice );
      vtmWriter.write( cuboid );
      vtmWriter.write( rank );
      vtmWriter.createMasterFile();
  }

  // Output on the console
  timer.update( iT );
  timer.printStep();
  sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

  if ( iT%vtmIter==0 || lastTimeStep )
  {
    sLattice.communicate();
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write( iT );
  }
}

template <typename CASE>
void simulate(CASE& myCase) {
  using T = CASE::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  auto& converter = sLattice.getUnitConverter();

  std::cout << "starting simulation with pressure = " << inletPressure << "..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ),
    superGeometry.getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( converter.getLatticeTime( physInterval ), residuum );
  timer.start();
  for (std::size_t iT=0 ; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {
    if ( converge.hasConverged() ) {
      std::cout << "Simulation converged." << std::endl;
      getResults(myCase, iT, timer, converge.hasConverged());
      break;
    }
    sLattice.collideAndStream();
    converge.takeValue(sLattice.getStatistics().getMaxU(), false );
  } 
  
  timer.stop();
  timer.printSummary();
}

template <typename CASE>
void setInitialControls(MyOptiCase& optiCase) {
  using T = CASE::value_t;
  std::vector<T> control({0.0001});
  optiCase.getController().set(control);
}

template <typename CASE>
void applyControls(MyOptiCase& optiCase) {
  inletPressure = optiCase.getController()[0];
}

/// Compute squared error between simulated and wanted mass flow rate
template <typename CASE>
CASE::value_t massFlowError(MyOptiCase& optiCase)
{
  using T = CASE::value_t;
  auto& myCase = optiCase.getCase(Controlled{});
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  auto& converter = sLattice.getUnitConverter();
  sLattice.setProcessingContext(ProcessingContext::Evaluation);

  SuperLatticeVelocity2D velocity(sLattice);
  SuperLatticeDensity2D density(sLattice);
  SuperPlaneIntegralFluxMass2D<T> massFlowRate(
    velocity, density, superGeometry, converter.getConversionFactorMass(),
    converter.getConversionFactorTime(), Vector<T,2>({T(0.5)*lx, T(0.5)*ly}),
    Vector<T,2>({0, 1}), BlockDataReductionMode::Analytical
  );
  const int input[3] = {0};
  T mFlow[4] = {0.};
  massFlowRate(mFlow, input);
  std::cout << "Mass flow rate = " << mFlow[0] << std::endl;
  const T res = mFlow[0];
  const T wantedRes {wantedMassFlow};
  return 0.5 * (res - wantedRes) * (res - wantedRes);
}

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  initialize( &argc, &argv );

  // === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParametersD;
  
  // === Step 3: Create Mesh ===
  auto mesh = createMesh(myCaseParametersD);

  // === Step 4: Create Case ===
  MyCase myCase(myCaseParametersD, mesh);

  // Create OptiCase
  MyOptiCase optiCase;

  // Set Case
  optiCase.setCase<Controlled>(myCase);

  // Define initial values for the control
  setInitialControls<MyCase>(optiCase);

  // Define objective computation routine
  optiCase.setObjective([&](MyOptiCase& optiCase) {
    myCase.resetLattices();
    applyControls<MyCase>(optiCase);
    prepareGeometry(myCase);
    prepareLattice(myCase);
    setInitialValues(myCase);
    simulate(myCase);
    return massFlowError<MyCase>(optiCase);
  });

  OptimizerLBFGS<S,std::vector<S>> optimizer(
    1, 1.e-7, 20, 1., 10, "StrongWolfe", 20, 1.e-4, true, "", "log",
    true, 0.01, true, 0., false, 0., true,
    {OptimizerLogType::value, OptimizerLogType::control, OptimizerLogType::derivative});
  optimizer.optimize(optiCase);
}
