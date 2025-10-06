/* Lattice Boltzmann sample, written in C++, using the OpenLB
 * library

 * Copyright (C) 2019-2022 Mathias J. Krause, Julius Jeßberger,
 * E-mail contact: info@openlb.net
 * The most recent release of OpenLB can be downloaded at
 * <http:  //www.openlb.net/>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;
using namespace olb::opti;
using namespace olb::parameters;

using T = double;
using NSDESCRIPTOR = D3Q19<>;
using ADDESCRIPTOR = D3Q7<>;

using MyCase = Case<
  NavierStokes,       Lattice<T,NSDESCRIPTOR>,
  Concentration0, Lattice<T,ADDESCRIPTOR>>;

using MyOptiCase = OptiCaseFDQ<Controlled, MyCase>;

util::TimeIntegrator<T> density1(0.0, 1.0, 0.1);
util::TimeIntegrator<T> density2(0.0, 1.0, 0.1);
util::TimeIntegrator<T> variance(0.0, 1.0, 0.1);
util::TimeIntegrator<T> varianceNorm(0.0, 1.0, 0.1);

namespace olb::parameters {

struct CHAR_PHYS_LENGTH : public descriptors::FIELD_BASE<1> { };
struct CHAR_PHYS_U : public descriptors::FIELD_BASE<1> { };
struct DIFFUSION : public descriptors::FIELD_BASE<1> { };
struct LATTICE_U : public descriptors::FIELD_BASE<1> { };
struct PHYS_START_PERIOD : public descriptors::FIELD_BASE<1> { };
struct PHYS_START_T : public descriptors::FIELD_BASE<1> { };

// control variables
struct PHYS_PERIOD  : public descriptors::FIELD_BASE<1> { };
struct AMPLITUDE_PHYS_PRESSURE : public descriptors::FIELD_BASE<1> { };
struct DIFFERENCE_PERIOD : public descriptors::FIELD_BASE<1> { };
}

Mesh<T,MyCase::d> createMesh(MyCase::ParametersD& params) {
  const int nC = util::max(16, 4 * singleton::mpi().getSize());
  const T dx = params.get<DX>();
  STLreader<T> stlReader("microMixer3d.stl", dx, T{1});
  IndicatorLayer3D<T> extendedDomain(stlReader, dx);
  Mesh<T,MyCase::d> mesh(extendedDomain, dx, nC);
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();
  const T dx = params.get<DX>();

  STLreader<T> stlReader("microMixer3d.stl", dx, T{1.});
  geometry.rename(0, 2, stlReader);
  geometry.rename(2, 1, {1, 1, 1});

  // Returns the minimum phys position in each direction for material 2
  Vector<T,3> minR = geometry.getStatistics().getMinPhysR(2);
  Vector<T,3> maxR = geometry.getStatistics().getMaxPhysR(2);
  Vector<T,3> centerR = geometry.getStatistics().getCenterPhysR(2);

  // sets circle of both inflows and the outflow, with direction and radius
  IndicatorCircle3D<T> inflow1(minR[0], minR[1], centerR[2], 0., -1., 0.,
                               (maxR[0] - minR[0]) / 3.);
  IndicatorCircle3D<T> inflow2(maxR[0], minR[1], centerR[2], 0., -1., 0.,
                               (maxR[0] - minR[0]) / 3.);
  IndicatorCircle3D<T> outflow((maxR[0] + minR[0]) / 2., maxR[1], centerR[2],
                               0., 1., 0., (maxR[0] - minR[0]) / 3.);

  // sets cylinder on that in-/out-flow circles with length
  IndicatorCylinder3D<T> layerInflow1(inflow1, dx);
  IndicatorCylinder3D<T> layerInflow2(inflow2, dx);
  IndicatorCylinder3D<T> layerOutflow(outflow, dx);

  geometry.rename(2, 3, 1, layerInflow1); // layer of inflow1 gets mat = 3
  geometry.rename(2, 4, 1, layerInflow2); // layer of inflow2 gets mat = 4
  geometry.rename(2, 5, 1, layerOutflow); // layer of outflow gets mat = 5

  geometry.clean(false);
  geometry.innerClean(false);
  geometry.checkForErrors(false);

  geometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;
  auto& params = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& sLatticeAD = myCase.getLattice(Concentration0{});

  sLattice.setUnitConverter<UnitConverterFromResolutionAndLatticeVelocity<T,NSDESCRIPTOR>>(
    (int)   params.get<RESOLUTION>(),        //resolution
    ( T )   params.get<LATTICE_U>(),         //charLatticeVelocity
    ( T )   params.get<CHAR_PHYS_LENGTH>(),  //charPhysLength
    ( T )   params.get<CHAR_PHYS_U>(),       //charPhysVelocity
    ( T )   params.get<VISCOSITY>(),         //physViscosity
    ( T )   params.get<DENSITY>()            //physDensity
  );
  sLatticeAD.setUnitConverter<UnitConverterFromResolutionAndLatticeVelocity<T,ADDESCRIPTOR>>(
    (int)   params.get<RESOLUTION>(),        //resolution
    ( T )   params.get<LATTICE_U>(),         //charLatticeVelocity
    ( T )   params.get<CHAR_PHYS_LENGTH>(),  //charPhysLength
    ( T )   params.get<CHAR_PHYS_U>(),       //charPhysVelocity
    ( T )   params.get<VISCOSITY>(),         //physViscosity
    ( T )   params.get<DENSITY>()            //physDensity
  );

  auto& converter = sLattice.getUnitConverter();
  converter.print();

  auto bulkIndicator = geometry.getMaterialIndicator({1,3,4,5});
  // dynamics for fluid
  sLattice.template defineDynamics<BGKdynamics>(bulkIndicator);

  // boundary conditions for fluid
  boundary::set<boundary::BounceBack>(sLattice, geometry, 2);
  boundary::set<boundary::InterpolatedPressure>(sLattice, geometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, geometry, 4);
  boundary::set<boundary::InterpolatedVelocity>(sLattice, geometry, 5);


  // dynamics for adsorptive
  auto bulkIndicatorAD = geometry.getMaterialIndicator({1,3,5});
  sLatticeAD.template defineDynamics<ParticleAdvectionDiffusionBGKdynamics>(bulkIndicatorAD);

  // boundary for adsorptive
  boundary::set<boundary::BounceBack>(sLatticeAD, geometry, 2);
  boundary::set<boundary::BounceBack>(sLatticeAD, geometry, 4);
  boundary::set<boundary::AdvectionDiffusionDirichlet>(sLatticeAD, geometry, 3);
  setZeroGradientBoundary<T,ADDESCRIPTOR>(sLatticeAD, geometry.getMaterialIndicator({5}));

  // set parameters
  sLattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  sLatticeAD.template setParameter<descriptors::OMEGA>(
    converter.getLatticeRelaxationFrequencyFromDiffusivity<ADDESCRIPTOR>(
      params.get<DIFFUSION>()));

  // define lattice coupling operator
  auto& coupling = myCase.setCouplingOperator(
    "NavierStokesAdvectionDiffusionCoupling",
    NavierStokesAdvectionDiffusionVelocityCoupling{},
    names::NavierStokes{},   sLattice,
    names::Concentration0{}, sLatticeAD
  );
  coupling.restrictTo(geometry.getMaterialIndicator({1}));

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  auto& geometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& sLatticeAD = myCase.getLattice(Concentration0{});

  // initialisation for fluid
  AnalyticalConst3D<T, T> rho1(1.);
  AnalyticalConst3D<T, T> u0(0., 0., 0.);

  auto initIndicator = geometry.getMaterialIndicator({1,2,3,4,5});
  sLattice.defineRhoU(initIndicator, rho1, u0);
  sLattice.iniEquilibrium(initIndicator, rho1, u0);

  // initialisation for adsorptive
  auto initIndicatorAD = geometry.getMaterialIndicator({1,2,4,5});
  AnalyticalConst3D<T, T> rhoSmall(1.e-8);
  sLatticeAD.defineRhoU(initIndicatorAD, rhoSmall, u0);
  sLatticeAD.defineRhoU(geometry, 3, rho1, u0);
  sLatticeAD.iniEquilibrium(initIndicatorAD, rhoSmall, u0);
  sLatticeAD.iniEquilibrium(geometry, 3, rho1, u0);

  // Make the lattice ready for simulation
  sLattice.initialize();
  sLatticeAD.initialize();
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT) {
  OstreamManager clout(std::cout, "setBoundaryValues");
  auto& params = myCase.getParameters();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& converter = sLattice.getUnitConverter();

  std::vector<T> maxVelocity(3, T());
  const T distanceToBoundary = converter.getConversionFactorLength() / 2.;
  const T latticeVelNS = converter.getLatticeVelocity(converter.getCharPhysVelocity());
  const std::size_t itStartTime = converter.getLatticeTime(params.get<PHYS_START_T>());

  if (iT <= itStartTime && iT % 50 == 0) {
    SinusStartScale<T,int> startScale(itStartTime, T(1));
    int help[1] = {(int) iT};
    T frac[3] = {T()};
    startScale(frac, help);

    // set lattice velocity on boundary
    maxVelocity[1] = latticeVelNS * frac[0];
    RectanglePoiseuille3D<T> u5(geometry, 5, maxVelocity,
                                distanceToBoundary, distanceToBoundary, distanceToBoundary);
    sLattice.defineU(geometry, 5, u5);
  }

  const int itStartPeriodTime = converter.getLatticeTime(params.get<PHYS_START_PERIOD>());
  const T amplitude = params.get<AMPLITUDE_PHYS_PRESSURE>() / converter.getConversionFactorPressure();

  T rho = 1.;
  if ((iT <= itStartTime + 0.5 * itStartPeriodTime) && (iT > itStartTime)) {
    Cosinus<T,T> cos(params.get<PHYS_START_PERIOD>(), T(0.5) * amplitude);
    T help[1] = {converter.getPhysTime(iT - itStartTime)};
    T frac[1] = {T()};
    cos(frac, help);

    rho = util::densityFromPressure<T,NSDESCRIPTOR>(T(-0.5) * amplitude + frac[0]);
    AnalyticalConst3D<T,T> rhovar(rho);
    sLattice.defineRho(geometry, 4, rhovar);
  }

  if (iT > itStartTime + 0.5 * itStartPeriodTime)
  {
    CosinusComposite<T,T> cosComp(params.get<PHYS_PERIOD>(), amplitude, params.get<DIFFERENCE_PERIOD>());
    T help[1] = { converter.getPhysTime(iT - itStartTime - 0.5 * itStartPeriodTime) };
    T frac[1] = { T() };
    cosComp(frac, help);

    rho = util::densityFromPressure<T,NSDESCRIPTOR>(-frac[0]);
    AnalyticalConst3D<T, T> rhovar( rho );
    sLattice.defineRho(geometry, 4, rhovar);
  }
}

void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT) {
  auto& params = myCase.getParameters();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& sLatticeAD = myCase.getLattice(Concentration0{});
  auto& geometry = myCase.getGeometry();
  auto& converter = sLattice.getUnitConverter();
  const unsigned vtkSaveT = converter.getLatticeTime(0.5);
  const unsigned statIter = converter.getLatticeTime(0.5);

  SuperLatticePhysVelocity3D<T,NSDESCRIPTOR> velocityNS(sLattice, converter);
  SuperLatticePhysPressure3D<T,NSDESCRIPTOR> pressureNS(sLattice, converter);

  SuperLatticeDensity3D<T,ADDESCRIPTOR> adsorptive(sLatticeAD);

  SuperVTMwriter3D<T> vtmWriter("microMixer3d");
  vtmWriter.addFunctor(velocityNS);
  vtmWriter.addFunctor(adsorptive);

  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  if (iT % vtkSaveT == 0) {
    vtmWriter.write(iT);
  }

  if (iT % statIter == 0) {
    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
  }

  const T time = converter.getPhysTime(iT);
  const T dt = converter.getPhysDeltaT();
  const T physMaxTime = params.get<MAX_PHYS_T>();
  const T physPeriod = params.get<PHYS_PERIOD>();

  if (iT == 0) {
    density1.reset(physMaxTime - T(2) * physPeriod, physMaxTime - physPeriod, dt);
    density2.reset(physMaxTime - physPeriod, physMaxTime, dt);
    variance.reset(physMaxTime - physPeriod, physMaxTime, dt);
    varianceNorm.reset(physMaxTime - physPeriod, physMaxTime, dt);
  }

  if (time >= physMaxTime - 2.0 * physPeriod - dt) {

    int input[1] = { };
    T output[adsorptive.getTargetDim()+1];  // for concentration

    // average concentration at the outlet
    SuperAverage3D<T>(adsorptive, geometry, 5).operator()(output, input);
    density1.takeValue(iT, output[0] / physPeriod);
    density2.takeValue(iT, output[0] / physPeriod);
  }

  if (time >= physMaxTime - 1.0 * physPeriod - dt) {

    int input[1] = { };
    T output1[1+1];  // for variance

    // variance of the density at the outlet
    const T mu = density1.getResult();
    SuperConst3D<T> expectedValue(geometry, mu);
    SuperAverage3D<T>((adsorptive - expectedValue) * (adsorptive - expectedValue), geometry, 5).operator()(output1, input);
    variance.takeValue(iT, output1[0] / physPeriod);
    // variance of the density at the outlet (with normalization)
    SuperConst3D<T> factor(geometry, 0.25 / (mu * mu));
    SuperAverage3D<T>(factor * (adsorptive - expectedValue) * (adsorptive - expectedValue), geometry, 5).operator()(output1, input);
    varianceNorm.takeValue(iT, output1[0] / physPeriod);
  }
}

void simulate(MyCase& myCase) {
  OstreamManager clout(std::cout, "simulate");
  auto& params = myCase.getParameters();

  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(
    params.get<MAX_PHYS_T>());
  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (size_t iT = 0; iT < iTmax; ++iT) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    myCase.getLattice(NavierStokes{}).collideAndStream();
    myCase.getOperator("NavierStokesAdvectionDiffusionCoupling").apply();
    myCase.getLattice(Concentration0{}).collideAndStream();

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, timer, iT);
  }

  timer.stop();
  timer.printSummary();
}


/*
void setInitialControl(MyOptiCase& optiCase) {
  optiCase.getController().set({0., 4., 0.5});
}

void applyControl(MyOptiCase& optiCase) {
  auto& control = optiCase.getController();
  physPeriod = control[0];
  amplitudePhysPressure = control[1];
  differencePeriod = control[2];
}

/// Return time-averaged std. deviation of the adsorptive density
T objectiveF(MyOptiCase& optiCase) {
  OstreamManager clout (std::cout, "objectiveF");
  clout << "Optimize by segregation intensity." << std::endl;

  // average over time
  // number of time steps = length of period / step length
  const T conversionTime = optiCase.getCase(Controlled{}).getLattice(NavierStokes{}).getUnitConverter().getConversionFactorTime();
  const T timeStepsPeriod(util::floor(optiCase.getController().get()[0] / conversionTime));
  clout << "Time steps per period            = " << timeStepsPeriod << std::endl;

  T refVariance = density1.getResult() * (T(1) - density1.getResult());
  T segIntensity = variance.getResult() / refVariance;
  clout << "Average density over first int.  = " << std::setprecision (12) << T(density1.getResult()) << std::endl;
  clout << "Av. density over second interval = " << std::setprecision (12) << T(density2.getResult()) << std::endl;
  clout << "Average variance                 = " << std::setprecision (12) << T(variance.getResult()) << std::endl;
  clout << "Norm. av. variance               = " << std::setprecision (12) << T(varianceNorm.getResult()) << std::endl;
  clout << "Danckwerts reference variance    = " << std::setprecision (12) << refVariance << std::endl;
  clout << "Danckwerts segregation intensity = " << std::setprecision (12)<< segIntensity << std::endl << std::endl;

  return segIntensity;
}

T computeObjective(MyOptiCase& optiCase) {
  auto& controlledCase = optiCase.getCase(Controlled{});

  // Set control variables
  controlledCase.resetLattices();
  applyControl(optiCase);

  // Prepare Case
  prepareGeometry(controlledCase);
  prepareLattice(controlledCase);
  setInitialValues(controlledCase);
  simulate(controlledCase);

  // Evaluate objective functor to compute objective value
  return objectiveF(optiCase);
}
*/
int main(int argc, char* argv[])
{
  // Step 2: Initialization
  initialize(&argc, &argv);

  // Step 2.1: Set parameters
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<LATTICE_U          >(0.1);
    myCaseParameters.set<CHAR_PHYS_U        >(0.0332);
    myCaseParameters.set<VISCOSITY          >(1.e-6);
    myCaseParameters.set<CHAR_PHYS_LENGTH   >(0.00133);
    myCaseParameters.set<MAX_PHYS_T         >(5.);   // time for fluid simulation
    myCaseParameters.set<PHYS_START_T       >(0.6);  // time to start fluid pulsation
    myCaseParameters.set<PHYS_START_PERIOD  >(0.4);  // time to start fluid pulsation
    myCaseParameters.set<DENSITY            >(1000);
    myCaseParameters.set<DIFFUSION          >(1.e-6);  // 1.e-9
    myCaseParameters.set<RESOLUTION         >(7);    // resolution of the hydraulic diameter  // 36
    myCaseParameters.set<DX>(myCaseParameters.get<CHAR_PHYS_LENGTH>() / myCaseParameters.get<RESOLUTION>());

    // to be removed for optimization?
    myCaseParameters.set<PHYS_PERIOD        >(0.6);
    myCaseParameters.set<AMPLITUDE_PHYS_PRESSURE>(3.);
    myCaseParameters.set<DIFFERENCE_PERIOD  >(0.5);
  }
  myCaseParameters.fromCLI(argc, argv);

  // Step 3: Create Mesh
  Mesh mesh = createMesh(myCaseParameters);

  // Step 4: Create Case
  MyCase myCase(myCaseParameters, mesh);

  MyOptiCase optiCase;
  optiCase.setCase<Controlled>(myCase);




  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase);

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);
/*
  // Initialize control parameters
  setInitialControl(optiCase);

  // set objective
  optiCase.setObjective(computeObjective);

  /// @li Step E: Create an Optimizer
  OptimizerLBFGS<T,std::vector<T>> optimizer(
    optiCase.getController().size(), 1.e-7, 20, 0.5, 10, "StrongWolfe", 20, 1.e-6, true, "", "log",
    true, 2.0, true, 0.1, false, 0., true,
    {OptimizerLogType::value, OptimizerLogType::control, OptimizerLogType::derivative});

  /// @li Step F: Optimize
  optimizer.optimize(optiCase);
*/
}
