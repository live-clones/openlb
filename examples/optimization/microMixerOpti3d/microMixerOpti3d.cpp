/* Lattice Boltzmann sample, written in C++, using the OpenLB
 * library
 *
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

// TODO: Case is not running yet!

//#define WriteVTK
//#define PULSATING
//#define ADD_INFOS
//#define VALIDATION_SENSITIVITIES
//#define ADT

#include <olb.h>

using namespace olb;
using namespace names;
using namespace opti;

using T = FLOATING_POINT_TYPE;
using NSDESCRIPTOR = descriptors::D3Q19<>;
using TDESCRIPTOR = descriptors::D3Q7<>;
using MyCase = Case<NavierStokes, Lattice<T,NSDESCRIPTOR>,
                    AdvectionDiffusion, Lattice<T,TDESCRIPTOR>>;
using MyOptiCase = OptiCaseFDQ<Controlled, MyCase>;

const T latticeU = 0.1;
const T charU = 0.0332;   // charU
const T charNu = 1.e-6;   // charNu
const T charL = 0.00133;  // hydraulic diameter = 4 * surface / perimeter
const T charRho = 1000;   // charRhoFluid
const int N = 36;         // resolution of the hydraulic diameter
const T physDeltaX = charL/N;

const T diffusion = 1.e-9;
const T physMaxTime = 3;   // time for fluid simulation
const T physStartTime = .6;    // time to start fluid pulsation

// Control variables
T physPeriod;
T amplitudePhysPressure;
T differencePeriod;

util::TimeIntegrator<T> density1(0.0, 1.0, 0.1);
util::TimeIntegrator<T> density2(0.0, 1.0, 0.1);
util::TimeIntegrator<T> variance(0.0, 1.0, 0.1);
util::TimeIntegrator<T> varianceNorm(0.0, 1.0, 0.1);
#ifdef ADD_INFOS
  util::TimeIntegrator<T> velocityIn1(0.0, 1.0, 0.1);
  util::TimeIntegrator<T> velocityIn2(0.0, 1.0, 0.1);
  util::TimeIntegrator<T> velocityOut(0.0, 1.0, 0.1);
  util::TimeIntegrator<T> pressureIn1(0.0, 1.0, 0.1);
  util::TimeIntegrator<T> pressureIn2(0.0, 1.0, 0.1);
  util::TimeIntegrator<T> pressurePlain(0.0, 1.0, 0.1);
  util::TimeIntegrator<T> densityPlain(0.0, 1.0, 0.1);
  util::TimeIntegrator<T> constant1(0.0, 1.0, 0.1);
  util::TimeIntegrator<T> constant2(0.0, 1.0, 0.1);
  util::TimeIntegrator<T> density1_in1(0.0, 1.0, 0.1);
  util::TimeIntegrator<T> density2_in1(0.0, 1.0, 0.1);
  util::TimeIntegrator<T> density1_in2(0.0, 1.0, 0.1);
  util::TimeIntegrator<T> density2_in2(0.0, 1.0, 0.1);
#endif
T minInletVelocity {0};

// for stability
const T physStartPeriod {0.4};
const T maxVelocityFactor {1.0};

Mesh<T,MyCase::d> createMesh() {
  const int nC = util::max(16, 4 * singleton::mpi().getSize());
  STLreader<T> stlReader("microMixer3d.stl", physDeltaX, T{1.});
  IndicatorLayer3D<T> extendedDomain(stlReader, physDeltaX);
  return Mesh<T,MyCase::d>(extendedDomain, physDeltaX, nC);
}

void prepareGeometry(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;
  auto& superGeometry = myCase.getGeometry();

  STLreader<T> stlReader("microMixer3d.stl", physDeltaX, T{1.});
  superGeometry.rename(0, 2, stlReader);
  superGeometry.rename(2, 1, {1, 1, 1});

// Returns the minimum phys position in each direction for material 2
  Vector<T,3> minR = superGeometry.getStatistics().getMinPhysR(2);
  Vector<T,3> maxR = superGeometry.getStatistics().getMaxPhysR(2);
  Vector<T,3> centerR = superGeometry.getStatistics().getCenterPhysR(2);

  // sets circle of both inflows and the outflow, with direction and radius
  IndicatorCircle3D<T> inflow1(minR[0], minR[1], centerR[2], 0., -1., 0.,
                               (maxR[0] - minR[0]) / 3.);
  IndicatorCircle3D<T> inflow2(maxR[0], minR[1], centerR[2], 0., -1., 0.,
                               (maxR[0] - minR[0]) / 3.);
  IndicatorCircle3D<T> outflow((maxR[0] + minR[0]) / 2., maxR[1], centerR[2],
                               0., 1., 0., (maxR[0] - minR[0]) / 3.);

  // sets cylinder on that in-/out-flow circles with length
  IndicatorCylinder3D<T> layerInflow1(inflow1, physDeltaX);
  IndicatorCylinder3D<T> layerInflow2(inflow2, physDeltaX);
  IndicatorCylinder3D<T> layerOutflow(outflow, physDeltaX);

  superGeometry.rename(2, 3, 1, layerInflow1); // layer of inflow1 gets mat = 3
  superGeometry.rename(2, 4, 1, layerInflow2); // layer of inflow2 gets mat = 4
  superGeometry.rename(2, 5, 1, layerOutflow); // layer of outflow gets mat = 5

  superGeometry.clean(false);
  superGeometry.innerClean(false);
  superGeometry.checkForErrors(false);

  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;
  auto& superGeometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& sLatticeAD = myCase.getLattice(AdvectionDiffusion{});

  sLattice.setUnitConverter<UnitConverterFromResolutionAndLatticeVelocity<T,NSDESCRIPTOR>>(
    (int)   N,              //resolution
    ( T )   latticeU,       //charLatticeVelocity
    ( T )   charL,          //charPhysLength
    ( T )   charU,          //charPhysVelocity
    ( T )   charNu,         //physViscosity
    ( T )   charRho         //physDensity
  );
  clout << "navier-stokes converter " << std::endl;

  sLatticeAD.setUnitConverter<UnitConverterFromResolutionAndLatticeVelocity<T,TDESCRIPTOR>>(
    (int)   N,              //resolution
    ( T )   latticeU,       //charLatticeVelocity
    ( T )   charL,          //charPhysLength
    ( T )   charU,          //charPhysVelocity
    ( T )   diffusion,      //physViscosity
    ( T )   charRho         //physDensity
  );
  clout << "advection-diffusion converter " << std::endl;
  auto& converterNS = sLattice.getUnitConverter();
  auto& converterAD = sLatticeAD.getUnitConverter();
  converterNS.print();
  converterAD.print();

  auto bulkIndicator = superGeometry.getMaterialIndicator({1,3,4,5});
  // dynamics for fluid
  sLattice.template defineDynamics<BGKdynamics<T,NSDESCRIPTOR>>(bulkIndicator);

  // boundary conditions for fluid
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);
  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 5);


  // dynamics for adsorptive
  auto bulkIndicatorAD = superGeometry.getMaterialIndicator({1,3,5});
  sLatticeAD.template defineDynamics<ParticleAdvectionDiffusionBGKdynamics<T,TDESCRIPTOR>>(bulkIndicatorAD);

  // boundary for adsorptive
  boundary::set<boundary::BounceBack>(sLatticeAD, superGeometry, 2);
  boundary::set<boundary::BounceBack>(sLatticeAD, superGeometry, 4);
  boundary::set<boundary::AdvectionDiffusionDirichlet>(sLatticeAD, superGeometry, 3);
  setZeroGradientBoundary<T,TDESCRIPTOR>(sLatticeAD, superGeometry.getMaterialIndicator({5}));

  // set parameters
  sLattice.template setParameter<descriptors::OMEGA>(converterNS.getLatticeRelaxationFrequency());
  sLatticeAD.template setParameter<descriptors::OMEGA>(converterAD.getLatticeRelaxationFrequency());

  // define lattice coupling operator
  auto& coupling = myCase.setCouplingOperator(
    "NavierStokesAdvectionDiffusionCoupling",
    NavierStokesAdvectionDiffusionCoupling{},
    names::NavierStokes{}, sLattice,
    names::Temperature{}, sLatticeAD
  );
  coupling.restrictTo(superGeometry.getMaterialIndicator({1}));

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  auto& superGeometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& sLatticeAD = myCase.getLattice(AdvectionDiffusion{});

  // initialisation for fluid
  AnalyticalConst3D<T, T> rho1(1.);
  AnalyticalConst3D<T, T> u0(0., 0., 0.);

  auto initIndicator = superGeometry.getMaterialIndicator({1,2,3,4,5});
  sLattice.defineRhoU(initIndicator, rho1, u0);
  sLattice.iniEquilibrium(initIndicator, rho1, u0);

  // initialisation for adsorptive
  auto initIndicatorAD = superGeometry.getMaterialIndicator({1,2,4,5});
  AnalyticalConst3D<T, T> rhoSmall(1.e-8);
  sLatticeAD.defineRhoU(initIndicatorAD, rhoSmall, u0);
  sLatticeAD.defineRhoU(superGeometry, 3, rho1, u0);
  sLatticeAD.iniEquilibrium(initIndicatorAD, rhoSmall, u0);
  sLatticeAD.iniEquilibrium(superGeometry, 3, rho1, u0);

  // Make the lattice ready for simulation
  sLattice.initialize();
  sLatticeAD.initialize();
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT) {
  OstreamManager clout(std::cout, "setBoundaryValues");
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  auto& converterNS = sLattice.getUnitConverter();

  std::vector <T> maxVelocity(3, T());
  const T distanceToBoundary = converterNS.getConversionFactorLength() / 2.;
  const T latticeVelNS = converterNS.getLatticeVelocity(converterNS.getCharPhysVelocity());
  const std::size_t itStartTime = converterNS.getLatticeTime(physStartTime);

  if (iT <= itStartTime && iT % 50 == 0) {
    SinusStartScale<T,int> startScale(itStartTime, T(1));
    int help[1] = {(int) iT};
    T frac[3] = {T()};
    startScale(frac, help);

    // set lattice velocity on boundary
    maxVelocity[1] = latticeVelNS * frac[0];
    RectanglePoiseuille3D<T> u5(superGeometry, 5, maxVelocity,
                                distanceToBoundary, distanceToBoundary, distanceToBoundary);
    sLattice.defineU(superGeometry, 5, u5);
  }

  if (iT == 0) {
#ifdef ADD_INFOS
    const S dt (converterNS.getPhysDeltaT());
    pressurePlain.reset(physMaxTime - physPeriod, physMaxTime, dt);
    densityPlain.reset(physMaxTime - physPeriod, physMaxTime, dt);
#endif
  }

#ifdef PULSATING
  const int itStartPeriodTime = converterNS.getLatticeTime(physStartPeriod);
  const T amplitude = amplitudePhysPressure / converterNS.getConversionFactorPressure();

  T rho = 1.;
  if ((iT <= itStartTime + 0.5 * itStartPeriodTime) && (iT > itStartTime)) {
    Cosinus<T,T> cos(physStartPeriod, T(0.5) * amplitude);
    T help[1] = {converterNS.getPhysTime(iT - itStartTime)};
    T frac[1] = {T()};
    cos(frac, help);

    rho = util::densityFromPressure<T,NSDESCRIPTOR>(T(-0.5) * amplitude + frac[0]);
    AnalyticalConst3D<T,T> rhovar(rho);
    sLattice.defineRho(superGeometry, 4, rhovar);
  }

  if (iT > itStartTime + 0.5 * itStartPeriodTime)
  {
    CosinusComposite<T,T> cosComp(physPeriod, amplitude, differencePeriod);
    T help[1] = { converterNS.getPhysTime(iT - itStartTime - 0.5 * itStartPeriodTime) };
    T frac[1] = { T() };
    cosComp(frac, help);

    rho = util::densityFromPressure<T,NSDESCRIPTOR>(-frac[0]);
    AnalyticalConst3D<T, T> rhovar( rho );
    sLattice.defineRho(superGeometry, 4, rhovar);
#ifdef ADD_INFOS
    pressurePlain.takeValue(iT, converterNS.getPhysPressure(frac[0]) / physPeriod);
    densityPlain.takeValue(iT, converterNS.getPhysDensity(util::densityFromPressure<T,NSDESCRIPTOR>(frac[0])) / physPeriod);
#endif
  }
#endif

#ifdef ADD_INFOS
  if (iT % 25 == 0) {
    if (singleton::mpi().getRank() == 0) {
      std::ofstream file;
      file.open("pulse.txt", std::ofstream::out | std::ofstream::app);
      if ( file.is_open() ) {
        file << time << "\t" << velocity << "\t" << rho << "\n";
        file.close();
      }
    }
  }
#endif
}

void getResults(MyCase& myCase, std::size_t iT) {
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& sLatticeAD = myCase.getLattice(AdvectionDiffusion{});
  auto& superGeometry = myCase.getGeometry();
  auto& converterNS = sLattice.getUnitConverter();

  SuperLatticePhysVelocity3D<T,NSDESCRIPTOR> velocityNS(sLattice, converterNS);
  SuperLatticePhysPressure3D<T,NSDESCRIPTOR> pressureNS(sLattice, converterNS);

  SuperLatticeDensity3D<T,TDESCRIPTOR> adsorptive(sLatticeAD);
#ifdef VALIDATION_SENSITIVITIES
  Vector<T,3> ll{0.005,0.0197,-0.0003};
  Vector<T,3> ur{0.009,0.0205,0.0011};
  IndicatorCuboid3D<T> centralCube_indi (ur-ll, ll);
  AnalyticalFfromIndicatorF3D<T,T> centralCube_analytical (centralCube_indi);
  SuperLatticeFfromAnalyticalF3D<T,NSDESCRIPTOR> centralCube (centralCube_analytical, sLattice);
  SuperIndicatorFfromIndicatorF3D<T> centralCube_sIndi (centralCube_indi, superGeometry);
#endif

#ifdef WriteVTK
  SuperVTMwriter3D<T> vtmWriter("microMixer3d");
  vtmWriter.addFunctor(velocityNS);
  vtmWriter.addFunctor(adsorptive);
#ifdef VALIDATION_SENSITIVITIES
  vtmWriter.addFunctor(centralCube);
#endif

  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  if (iT % converterNS.getLatticeTime(vtkSaveT) == 0) {
    vtmWriter.write(iT);
  }
#endif

  const T time = converterNS.getPhysTime(iT);
  const T dt = converterNS.getPhysDeltaT();

  if (iT == 0) {
    density1.reset(physMaxTime - T(2) * physPeriod, physMaxTime - physPeriod, dt);
    density2.reset(physMaxTime - physPeriod, physMaxTime, dt);
    variance.reset(physMaxTime - physPeriod, physMaxTime, dt);
    varianceNorm.reset(physMaxTime - physPeriod, physMaxTime, dt);
    #ifdef ADD_INFOS
      velocityIn1.reset(physMaxTime - physPeriod, physMaxTime, dt);
      velocityIn2.reset(physMaxTime - physPeriod, physMaxTime, dt);
      velocityOut.reset(physMaxTime - physPeriod, physMaxTime, dt);
      pressureIn1.reset(physMaxTime - physPeriod, physMaxTime, dt);
      pressureIn2.reset(physMaxTime - physPeriod, physMaxTime, dt);
      constant1.reset(physMaxTime - T(2) * physPeriod, physMaxTime - physPeriod, dt);
      constant2.reset(physMaxTime - physPeriod, physMaxTime, dt);
      density1_in1.reset(physMaxTime - T(2) * physPeriod, physMaxTime - physPeriod, dt);
      density2_in1.reset(physMaxTime - physPeriod, physMaxTime, dt);
      density1_in2.reset(physMaxTime - T(2) * physPeriod, physMaxTime - physPeriod, dt);
      density2_in2.reset(physMaxTime - physPeriod, physMaxTime, dt);
    #endif
    minInletVelocity = 0;
  }

  if (time >= physMaxTime - 2.0 * physPeriod - dt) {

    int input[1] = { };
    T output[adsorptive.getTargetDim()+1];  // for concentration

    // average concentration at the outlet
    SuperAverage3D<T>(adsorptive, superGeometry, 5).operator()(output, input);
    density1.takeValue(iT, output[0] / physPeriod);
    density2.takeValue(iT, output[0] / physPeriod);
    #ifdef ADD_INFOS
      constant1.takeValue(iT, T(1) / physPeriod);
      constant2.takeValue(iT, T(1) / physPeriod);
      SuperAverage3D<T>(adsorptive, superGeometry, 3).operator()(output, input);
      density1_in1.takeValue(iT, output[0] / physPeriod);
      density2_in1.takeValue(iT, output[0] / physPeriod);
      SuperAverage3D<T>(adsorptive, superGeometry, 4).operator()(output, input);
      density1_in2.takeValue(iT, output[0] / physPeriod);
      density2_in2.takeValue(iT, output[0] / physPeriod);
    #endif
  }

  if (time >= physMaxTime - 1.0 * physPeriod - dt) {

    int input[1] = { };
    T output1[1+1];  // for variance
    #ifdef ADD_INFOS
      T output4[velocityNS.getTargetDim()+1];  // for velocity
      T output6[pressureNS.getTargetDim()+1];  // for inlet pressure

      // average velocity at the inlet
      SuperAverage3D<T>(velocityNS, superGeometry, 3).operator()(output4, input);
      velocityIn1.takeValue(iT, output4[1] / physPeriod);
      SuperAverage3D<T>(velocityNS, superGeometry, 4).operator()(output4, input);
      velocityIn2.takeValue(iT, output4[1] / physPeriod);
      // average velocity at the outlet
      SuperAverage3D<T>(velocityNS, superGeometry, 5).operator()(output4, input);
      velocityOut.takeValue(iT, output4[1] / physPeriod);
      // av. pressure at the inlet
      SuperAverage3D<T>(pressureNS, superGeometry, 3).operator()(output6, input);
      pressureIn1.takeValue(iT, output6[0] / physPeriod);
      SuperAverage3D<T>(pressureNS, superGeometry, 4).operator()(output6, input);
      pressureIn2.takeValue(iT, output6[0] / physPeriod);
    #endif

    // variance of the density at the outlet
    const T mu = density1.getResult();
    SuperConst3D<T> expectedValue(superGeometry, mu);
    SuperAverage3D<T>((adsorptive - expectedValue) * (adsorptive - expectedValue), superGeometry, 5).operator()(output1, input);
    variance.takeValue(iT, output1[0] / physPeriod);
    // variance of the density at the outlet (with normalization)
    SuperConst3D<T> factor(superGeometry, 0.25 / (mu * mu));
    SuperAverage3D<T>(factor * (adsorptive - expectedValue) * (adsorptive - expectedValue), superGeometry, 5).operator()(output1, input);
    varianceNorm.takeValue(iT, output1[0] / physPeriod);
  }

#ifdef ADD_INFOS
  int input[1] = { };
  T output7[3+1];
  SuperMin3D<T>(velocityNS, superGeometry, 3).operator()(output7, input);
  minInletVelocity = util::min(minInletVelocity, output7[1]);
#endif

#ifdef VALIDATION_SENSITIVITIES
  if (time >= physMaxTime - dt) {
    int input2[1] = { };
    T output8[3+1];
    T output9[3+1];
    T output10[1+1];
    T output11[1+1];

    SuperL2Norm3D<T> (velocityNS, superGeometry, 3).operator()(output8, input2);
    SuperL2Norm3D<T> (velocityNS, centralCube_sIndi).operator()(output9, input2);
    SuperL2Norm3D<T> (adsorptive, superGeometry, 5).operator()(output10, input2);
    SuperL2Norm3D<T> (adsorptive, centralCube_sIndi).operator()(output11, input2);


    if (singleton::mpi().getRank() == 0) {
      std::ofstream file;
      file.open("sensitivities.txt", std::ofstream::out | std::ofstream::app);
      if ( file.is_open() ) {
        // period length, av. vel. inlet1, av. vel. box, av. conc. outlet, av. conc. box
        file << std::setprecision(12) << optiCase.getController().getControl(0) << ' '
          << output8[0] << ' '
          << output9[0] << ' '
          << output10[0] << ' '
          << output11[0] << '\n' << std::endl;
        file.close();
      }
    }
  }
#endif
}

void simulate(MyCase& myCase) {
  OstreamManager clout(std::cout, "simulate");

  /// === 4th Step: Fluid Main Loop with Timer ===
  util::Timer<T> timer(myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(physMaxTime),
                       myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (size_t iT = 0; iT < myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(physMaxTime); ++iT) {
    setTemporalValues(myCase, iT);

    myCase.getLattice(NavierStokes{}).collideAndStream();
    myCase.getOperator("NavierStokesAdvectionDiffusionCoupling").apply();
    myCase.getLattice(AdvectionDiffusion{}).collideAndStream();

    getResults(myCase, iT);
  }

  timer.stop();
  timer.printSummary();
}

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

  #ifdef ADD_INFOS
    clout << "Av. int. over 1/p, first inverval= " << std::setprecision (12) << T(constant1.getResult()) << std::endl;
    clout << "Av. int. over 1/p, sec. inverval = " << std::setprecision (12) << T(constant2.getResult()) << std::endl;
    clout << "Av. velocity at inlet 1          = " << std::setprecision (12)<< T(velocityIn1.getResult()) << std::endl;
    clout << "Av. velocity at inlet 2          = " << std::setprecision (12)<< T(velocityIn2.getResult()) << std::endl;
    clout << "Av. velocity at outlet           = " << std::setprecision (12)<< T(velocityOut.getResult()) << std::endl;
    clout << "Av. pressure at inlet 1          = " << std::setprecision (12)<< T(pressureIn1.getResult()) << std::endl;
    clout << "Av. pressure at inlet 2          = " << std::setprecision (12)<< T(pressureIn2.getResult()) << std::endl;
    clout << "Av. pressure at inlet 2 plain    = " << std::setprecision (12)<< T(pressurePlain.getResult()) << std::endl;
    clout << "Av. density at inlet 2 plain     = " << std::setprecision (12)<< T(densityPlain.getResult()) << std::endl;
    clout << "Av. conc. at inlet 1, first int. = " << std::setprecision (12)<< T(density1_in1.getResult()) << std::endl;
    clout << "Av. conc. at inlet 1, sec. int.  = " << std::setprecision (12)<< T(density2_in1.getResult()) << std::endl;
    clout << "Av. conc. at inlet 2, first int. = " << std::setprecision (12)<< T(density1_in2.getResult()) << std::endl;
    clout << "Av. conc. at inlet 2, sec. int.  = " << std::setprecision (12)<< T(density2_in2.getResult()) << std::endl;
    clout << "Min inlet velocity = " << T(minInletVelocity) << std::endl;
  #endif

  // clout << "Optimize by norm. av. variance." << std::endl;
  // output[0] = varianceNorm.getResult();

  if (singleton::mpi().getRank() == 0) {
    std::ofstream file;
    file.open("periods.txt", std::ofstream::out | std::ofstream::app);
    if ( file.is_open() ) {
      file << optiCase.getController().get()[0] << " " << std::setprecision (12)
        << T(density1.getResult()) << " " << T(density2.getResult()) << std::endl;
      file.close();
    }
  }

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

int main(int argc, char* argv[])
{
  // Step 2: Initialization
  initialize(&argc, &argv);

  // Step 2.1: Set parameters
  MyCase::ParametersD myCaseParametersD;
  {
    using namespace parameters;
    // list here parameters
  }

  // Step 3: Create Mesh
  Mesh mesh = createMesh();

  // Step 4: Create Case
  MyCase myCase(myCaseParametersD, mesh);

  MyOptiCase optiCase;
  optiCase.setCase<Controlled>(myCase);

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
}
