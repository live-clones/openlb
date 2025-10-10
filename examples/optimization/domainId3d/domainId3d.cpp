/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Shota Ito
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

#include <olb.h>
#include "../helper.h"  // Will be removed once SuperLatticeFieldReductionO enables indicator support

using namespace olb;

using T = double;
using DESCRIPTOR = descriptors::DualPorousD3Q19Descriptor;
using DYNAMICS = PorousBGKdynamics<T,DESCRIPTOR>;
using DUAL_DYNAMICS = DualPorousBGKDynamics<T,DESCRIPTOR>;     // Dynamics of dual problem
using OBJECTIVE = functors::L2DistanceF<functors::VelocityF>;  // Objective functional
using CONTROLS = descriptors::POROSITY;                        // Controlled field


const int resolution      = 20;
const T relaxationTime    = 0.8;
const T physLength        = 1.0;
const T physVelocity      = 1.0;
const T physViscosity     = 0.1;
const T physDensity       = 1.0;
const T physMaxT          = 6.0;
const T physStartT        = 1.0;
const T physBoundaryT     = 0.02;

IndicatorCuboid3D<T> designDomain({0.4, 0.4, 0.4}, {0.3, 0.3, 0.3});
IndicatorCuboid3D<T> referenceObject({0.2, 0.2, 0.2}, {0.4, 0.4, 0.4});

void prepareGeometry(const UnitConverter<T,DESCRIPTOR>& converter,
                     SuperGeometry<T,3>& superGeometry)
{
  superGeometry.rename(0, 2);
  superGeometry.rename(2, 1, {1,1,1});

  Vector<T,3> origin{0.,0.,0.};
  origin[0] += converter.getPhysDeltaX()/2.;
  origin[1] -= converter.getPhysDeltaX()/2.;
  origin[2] += converter.getPhysDeltaX()/2.;

  Vector<T,3> extend{physLength, 0., physLength};
  extend[0] -= 2.*converter.getPhysDeltaX()/2.;
  extend[1] +=    converter.getPhysDeltaX()/2.;
  extend[2] -= 2.*converter.getPhysDeltaX()/2.;

  IndicatorCuboid3D<T> inflow(extend, origin);
  superGeometry.rename(2, 3, inflow);

  origin = {0., physLength, 0.};
  origin[0] += converter.getPhysDeltaX()/2.;
  origin[1] -= converter.getPhysDeltaX()/2.;
  origin[2] += converter.getPhysDeltaX()/2.;

  extend = {physLength, physLength, physLength};
  extend[0] -= 2.*converter.getPhysDeltaX()/2.;
  extend[1] +=    converter.getPhysDeltaX()/2.;
  extend[2] -= 2.*converter.getPhysDeltaX()/2.;

  IndicatorCuboid3D<T> outflow(extend, origin);
  superGeometry.rename(2, 4, outflow);

  // Indicators for material 6 (designDomain)
  superGeometry.rename(1, 6, designDomain);

  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
}

void prepareLattice(const UnitConverter<T,DESCRIPTOR>& converter,
                    SuperLattice<T, DESCRIPTOR>& sLattice,
                    SuperGeometry<T,3>& superGeometry)
{
  auto bulkIndicator = superGeometry.getMaterialIndicator({1,3,4,6});
  sLattice.template defineDynamics<PorousBGKdynamics>(bulkIndicator);

  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);
  boundary::set<boundary::LocalVelocity>(sLattice, superGeometry.getMaterialIndicator({3}));
  boundary::set<boundary::LocalPressure>(sLattice, superGeometry.getMaterialIndicator({4}));

  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  AnalyticalConst3D<T,T> rhoF(1);
  Vector<T,3> velocity;
  AnalyticalConst3D<T,T> uF(velocity);

  sLattice.defineRhoU(bulkIndicator, rhoF, uF);
  sLattice.iniEquilibrium(bulkIndicator, rhoF, uF);

  AnalyticalConst3D<T,T> one(1.);
  AnalyticalConst3D<T,T> zero(0.);

  sLattice.template defineField<descriptors::POROSITY>(bulkIndicator, one);
  sLattice.template defineField<descriptors::POROSITY>(superGeometry, referenceObject, zero);

  sLattice.initialize();
}

void setBoundaryValues(const UnitConverter<T,DESCRIPTOR>& converter,
                       SuperLattice<T, DESCRIPTOR>& sLattice,
                       std::size_t iT, SuperGeometry<T,3>& superGeometry)
{
  const std::size_t itStart = converter.getLatticeTime(physStartT);
  const std::size_t itUpdate = converter.getLatticeTime(physBoundaryT);
  if (iT <= itStart && iT % itUpdate == 0) {
    PolynomialStartScale<T,std::size_t> StartScale( itStart, T( 1 ) );
    std::size_t iTvec[1] = {iT};
    T frac[1] = {};
    StartScale( frac,iTvec );

    AnalyticalConst3D<T,T> uF(0., frac[0] * converter.getCharLatticeVelocity(), 0.);
    sLattice.defineU(superGeometry, 3, uF);
    sLattice.template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }
}

void getResults(const UnitConverter<T,DESCRIPTOR>& converter,
                SuperLattice<T, DESCRIPTOR>& sLattice,
                LatticeResults<T, DESCRIPTOR>& results,
                std::size_t iT, util::Timer<T> timer)
{
  const std::size_t iTvtk = converter.getLatticeTime(physMaxT/5);
  const std::size_t iTlog = converter.getLatticeTime(physMaxT/5);

  // Writes the VTK files
  if (iT%iTvtk == 0) {
    results.write(converter, sLattice, iT);
  }

  // Get statistics
  if (iT%iTlog == 0) {
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    timer.print(iT);
  }
}

auto build(std::string name) {
  // === 1st Step: Initialization ===
  OstreamManager clout(std::cout, name);
  clout << "Building simulation setup..." << std::endl;
  auto lData = std::make_unique<LatticeData<T,DESCRIPTOR>>(name);

  // Provide the unit converter the characteristic entities
  auto& converter = lData->create<UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    resolution,        // resolution
    relaxationTime,    // relaxationTime
    physLength,        // charPhysLength: reference length of simulation geometry in [m]
    physVelocity,      // charPhysVelocity: highest expected velocity during simulation in [m/s]
    physViscosity,     // physViscosity: physical kinematic viscosity in [m^2/s]
    physDensity        // physDensity: physical density [kg/m^3]
  );
  converter.print();

  // === 2nd Step: Prepare Geometry ===
  Vector extend{physLength, physLength, physLength};
  Vector origin{0, 0, 0};
  IndicatorCuboid3D<T> cuboid(extend, origin);
  auto& cuboidDecomposition = lData->create<CuboidDecomposition3D<T>>(cuboid, converter.getPhysDeltaX(),singleton::mpi().getSize());

  auto& loadBalancer = lData->create<HeuristicLoadBalancer<T>>(cuboidDecomposition);

  auto& superGeometry = lData->create<SuperGeometry<T,3>>(cuboidDecomposition, loadBalancer);
  prepareGeometry(converter, superGeometry);

  lData->create<SuperLattice<T,DESCRIPTOR>>(converter, superGeometry);
  clout << "Done." << std::endl;
  return lData;
}

void simulate(std::unique_ptr<LatticeData<T,DESCRIPTOR>>& lData) {
  auto& converter = lData->getUnitConverter();
  auto& superGeometry = lData->getSuperGeometry();
  auto& sLattice = lData->getSuperLattice();
  auto& results = lData->getLatticeResults();

  OstreamManager clout(std::cout, lData->getName());
  clout << "Starting simulation..." << std::endl;

  // === 4th Step: Main Loop with Timer ===
  const std::size_t iTmax = converter.getLatticeTime(physMaxT);
  util::Timer<T> timer(iTmax, superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(converter, sLattice, iT, superGeometry);
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults(converter, sLattice, results, iT, timer);
  }

  timer.stop();
  timer.printSummary();
  clout << "Done." << std::endl;
}

void prepareReference(auto& reference) {
  prepareLattice(reference->getUnitConverter(), reference->getSuperLattice(), reference->getSuperGeometry());
  reference->getLatticeResults().template add<SuperLatticePhysPressure3D<T,DESCRIPTOR>,
                                              SuperLatticePhysVelocity3D<T,DESCRIPTOR>,
                                              SuperLatticePorosity3D<T,DESCRIPTOR>>();
}

void preparePrimal(auto& primal, const std::vector<T>& controls) {
  auto& converter = primal->getUnitConverter();
  auto& superGeometry = primal->getSuperGeometry();
  auto& sLattice = primal->getSuperLattice();
  auto bulkIndicator = superGeometry.getMaterialIndicator({1,2,3,4,6});

  // Define primal physics
  sLattice.template defineDynamics<DYNAMICS>(superGeometry.getMaterialIndicator({1,6}));
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);
  boundary::set<boundary::LocalVelocity>(sLattice, superGeometry.getMaterialIndicator({3}));
  boundary::set<boundary::LocalPressure>(sLattice, superGeometry.getMaterialIndicator({4}));
  sLattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  // Initialize primal problem
  AnalyticalConst3D<T,T> rhoF(1);
  Vector<T,3> velocity;
  AnalyticalConst3D<T,T> uF(velocity);
  sLattice.defineRhoU(bulkIndicator, rhoF, uF);
  sLattice.iniEquilibrium(bulkIndicator, rhoF, uF);
  AnalyticalConst3D<T,T> one(1.);
  sLattice.template defineField<CONTROLS>(bulkIndicator, one);

  // Convert control vector to porosity field
  std::vector<T> porosity = opti::applyProjection<opti::projection::Sigmoid<T>>(controls);
  setFieldFromSerialized<CONTROLS>(porosity, sLattice, superGeometry.getMaterialIndicator({6}));
  sLattice.initialize();

  primal->getLatticeResults().template add<SuperLatticePhysPressure3D<T,DESCRIPTOR>,
                                           SuperLatticePhysVelocity3D<T,DESCRIPTOR>,
                                           SuperLatticePorosity3D<T,DESCRIPTOR>>();
}

T computeObjective(auto& reference, auto& primal) {
  const auto& converter = reference->getUnitConverter();
  auto& superGeometry = reference->getSuperGeometry();
  auto& refLattice = reference->getSuperLattice();
  auto& primalLattice = primal->getSuperLattice();
  auto objectiveDomain = superGeometry.getMaterialIndicator({1,6});

  // Evaluate functor for objective computation
  auto objectiveO = makeWriteFunctorO<OBJECTIVE,opti::J>(primalLattice);
  objectiveO->restrictTo(objectiveDomain);
  writePhysFunctorTo<functors::VelocityF,OBJECTIVE::Reference>(refLattice,
                                                               objectiveDomain,
                                                               converter.getConversionFactorVelocity());
  // Get solution from the reference simulation for the inverse problem
  copyFields<OBJECTIVE::Reference,OBJECTIVE::Reference>(refLattice, primalLattice);
  objectiveO->template setParameter<descriptors::CONVERSION>(converter.getConversionFactorVelocity());
  objectiveO->template setParameter<descriptors::NORMALIZE>(norm(refLattice, converter, objectiveDomain));
  objectiveO->apply();

  // Compute source term for the dual simulation
  auto objectiveDerivativeO = makeWriteFunctorO<functors::DerivativeF<OBJECTIVE,descriptors::POPULATION,DYNAMICS>,
                                                opti::DJDF>(primalLattice);
  objectiveDerivativeO->restrictTo(objectiveDomain);
  objectiveDerivativeO->template setParameter<descriptors::CONVERSION>(converter.getConversionFactorVelocity());
  objectiveDerivativeO->template setParameter<descriptors::NORMALIZE>(norm(refLattice, converter, objectiveDomain));
  objectiveDerivativeO->template setParameter<descriptors::DX>(converter.getPhysDeltaX());
  objectiveDerivativeO->apply();

  primalLattice.setProcessingContext(ProcessingContext::Evaluation);
  return integrate<opti::J>(primalLattice, objectiveDomain)[0];
}

void prepareDual(auto& dual, auto& primal, auto& reference) {
  auto& converter = dual->getUnitConverter();
  auto& superGeometry = dual->getSuperGeometry();
  auto& sLattice = dual->getSuperLattice();
  auto& primalLattice = primal->getSuperLattice();
  auto bulkIndicator = superGeometry.getMaterialIndicator({1,2,3,4,6});

  // Define dual physics
  sLattice.template defineDynamics<DUAL_DYNAMICS>(superGeometry.getMaterialIndicator({1,6}));
  boundary::set<boundary::BounceBack>(sLattice, superGeometry.getMaterialIndicator({2,3,4}));
  sLattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  // Initialize dual problem
  AnalyticalConst3D<T,T> rhoF(1);
  Vector<T,3> velocity;
  AnalyticalConst3D<T,T> uF(velocity);
  sLattice.defineRhoU(bulkIndicator, rhoF, uF);
  sLattice.iniEquilibrium(bulkIndicator, rhoF, uF);

  // This needs to be before copying the fields as otherwise the the copied fields will be overwritten
  sLattice.initialize();

  // Provide fields required by the dual collision operator
  copyFields<CONTROLS,CONTROLS>(primalLattice, sLattice);
  writeFunctorTo<functors::PopulationF,opti::F>(primalLattice, bulkIndicator);
  copyFields<opti::F,opti::F>(primalLattice, sLattice);
  copyFields<opti::DJDF,opti::DJDF>(primalLattice, sLattice);

  dual->getLatticeResults().template add<SuperLatticePhysPressure3D<T,DESCRIPTOR>,
                                         SuperLatticePhysVelocity3D<T,DESCRIPTOR>,
                                         SuperLatticePorosity3D<T,DESCRIPTOR>,
                                         SuperLatticeField3D<T,DESCRIPTOR,opti::DJDF>>();
}

std::vector<T> computeSensitivity(auto& dual, auto& primal, const std::vector<T>& controls) {
  const auto& converter = primal->getUnitConverter();
  auto& superGeometry = primal->getSuperGeometry();
  auto& primalLattice = primal->getSuperLattice();
  auto& dualLattice = dual->getSuperLattice();

  // Evaluate optimality condition
  auto optimalityO = makeWriteFunctorO<functors::OptimalityF<DYNAMICS,CONTROLS>,
                                       opti::SENSITIVITY<CONTROLS>>(dualLattice);
  // Derivative of applied projection
  std::vector<T> projectionD = opti::applyProjection<opti::projection::SigmoidD<T>>(controls);
  setFieldFromSerialized<opti::DPROJECTIONDALPHA<CONTROLS>>(projectionD, dualLattice, superGeometry.getMaterialIndicator({6}));
  dualLattice.template setProcessingContext<Array<opti::DPROJECTIONDALPHA<CONTROLS>>>(ProcessingContext::Simulation);
  // Compute jacobian of collision operator regarding control variable
  auto dCDalphaO = makeWriteFunctorO<functors::DerivativeF<functors::CollisionF<DYNAMICS>,CONTROLS,DYNAMICS>,
                                     opti::DCDALPHA<CONTROLS>>(primalLattice);
  dCDalphaO->template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  dCDalphaO->template setParameter<descriptors::DX>(1.0);
  dCDalphaO->apply();
  // Jacobian is computed on primal lattice as jacobian is evaluated for primal populations
  copyFields<opti::DCDALPHA<CONTROLS>,opti::DCDALPHA<CONTROLS>>(primalLattice, dualLattice);
  optimalityO->template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  optimalityO->apply();

  // Return serial vector containing total derivatives of objective regarding controls
  dualLattice.setProcessingContext(ProcessingContext::Evaluation);
  return getSerializedFromField<opti::SENSITIVITY<CONTROLS>>(dualLattice, superGeometry.getMaterialIndicator({6}));
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);

  // Create Simulation instances
  auto refSimulation = build(" REFERENCE ");
  auto primalSimulation = build(" PRIMAL ");
  auto dualSimulation = build(" DUAL ");

  // Create instances for IO for each simulation
  refSimulation->template create<LatticeResults<T,DESCRIPTOR>>("refSolution");
  primalSimulation->template create<LatticeResults<T,DESCRIPTOR>>("primalSolution");
  dualSimulation->template create<LatticeResults<T,DESCRIPTOR>>("dualSolution");

  // Simulate reference solution
  prepareReference(refSimulation);
  simulate(refSimulation);

  // Define routine for objective and derivative computation
  opti::OptiCaseAnalytical<T,std::vector<T>> optiCase;
  optiCase.setObjective([&](const std::vector<T>& controls) -> T {
    primalSimulation->resetLattice();
    preparePrimal(primalSimulation, controls);
    simulate(primalSimulation);
    return computeObjective(refSimulation, primalSimulation);
  });
  optiCase.setDerivative([&](const std::vector<T>& controls, std::vector<T>& derivatives) {
    dualSimulation->resetLattice();
    prepareDual(dualSimulation, primalSimulation, refSimulation);
    simulate(dualSimulation);
    derivatives = computeSensitivity(dualSimulation, primalSimulation, controls);
  });

  // Optimizer setup
  const std::size_t dimCtrl = getSerializedFieldSize<CONTROLS>(primalSimulation->getSuperLattice(),
                                                               primalSimulation->getSuperGeometry().getMaterialIndicator({6}));
  opti::OptimizerLBFGS<T,std::vector<T>> optimizer(
      dimCtrl, 1.e-10, 10, 1., 20, "Wolfe", 20, 1.e-4);
  optimizer.setStartValue(opti::projection::getInitialControl(1e-2,
    opti::projection::Sigmoid<T>(),
    refSimulation->getUnitConverter(),
    opti::Permeability));
  optimizer.optimize(optiCase);
}
