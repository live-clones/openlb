/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Florian Kaiser
 *                2008 Orestis Malaspinas
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

/* Header related to squareCavity2dTurbulent.cpp
 * The reference is the paper in "Gaedtke, M., Wachter, S., Raedle, M., Nirschl, H., & Krause, M. J. (2018).
 * Application of a lattice Boltzmann method combined with a Smagorinsky turbulence model to spatially resolved heat flux inside a refrigerated vehicle.
 * Computers & Mathematics with Applications, 76(10), 2315-2329."
 */

// natural convection of air in a square cavity in 2D

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::names;

namespace olb::parameters {
  struct SIM_VALUES  : public descriptors::FIELD_BASE<5> { };
}

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D2Q9<descriptors::FORCE>>,
  Temperature,  Lattice<double, descriptors::D2Q5<descriptors::VELOCITY>>
>;

/// @brief Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters){
    using T = MyCase::value_t;

    const T dx = parameters.get<parameters::PHYS_DELTA_X>();
    const Vector domainExtend = parameters.get<parameters::DOMAIN_EXTENT>();
    const Vector extend = {domainExtend[0] + dx, domainExtend[1] + dx};
    Vector origin{-dx / 2., -dx / 2.};
    IndicatorCuboid2D<T> cuboid(extend, origin);

    Mesh<T,MyCase::d> mesh(cuboid, dx, singleton::mpi().getSize());
    mesh.setOverlap(parameters.get<parameters::OVERLAP>());
    mesh.getCuboidDecomposition().setPeriodicity({false,false});

    return mesh;
}

void prepareGeometry(MyCase& myCase){
    OstreamManager clout(std::cout, "preprareGeometry");

    clout << "Prepare Geometry ..." << std::endl;

    using T = MyCase::value_t;
    auto& geometry = myCase.getGeometry();
    auto& parameters = myCase.getParameters();

    const T dx                = parameters.get<parameters::PHYS_DELTA_X>();
    const Vector domainExtend = parameters.get<parameters::DOMAIN_EXTENT>();
    const T lx                = domainExtend[0];

    geometry.rename(0, 4);

    Vector origin {0, 0};
    Vector extend {lx, lx};
    IndicatorCuboid2D<T> cuboid2(extend, origin);

    geometry.rename(4, 1, cuboid2);

    Vector extendwallleft{dx, lx + dx};
    Vector originwallleft{-dx / 2., -dx / 2.};
    IndicatorCuboid2D<T> wallleft(extendwallleft, originwallleft);

    Vector extendwallright{dx, lx + dx};
    Vector originwallright{lx, -dx / 2.};
    IndicatorCuboid2D<T> wallright(extendwallright, originwallright);

    geometry.rename(4,2,1,wallleft);
    geometry.rename(4,3,1,wallright);

    geometry.clean();
    geometry.innerClean();
    geometry.checkForErrors();

    geometry.print();
    clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase){
    OstreamManager clout(std::cout,"prepareLattice");
    clout << "Prepare Lattice ..." << std::endl;

    using T = MyCase::value_t;
    using NSEDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
    using ADEDESCRIPTOR = MyCase::descriptor_t_of<Temperature>;

    auto& geometry = myCase.getGeometry();
    auto& parameters = myCase.getParameters();

    auto& NSElattice = myCase.getLattice(NavierStokes{});
    auto& ADElattice = myCase.getLattice(Temperature{});

    const T physCharLength          = parameters.get<parameters::PHYS_CHAR_LENGTH>();
    const int N             = parameters.get<parameters::RESOLUTION>();
    const T physViscosity           = parameters.get<parameters::PHYS_KINEMATIC_VISCOSITY>();
    const T physDeltaX              = parameters.get<parameters::PHYS_DELTA_X>();
    const T physCharVelocity        = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
    const T physDeltaT              = 2. * 0.056 / physCharVelocity * physCharLength / N;
    const T physDensity             = parameters.get<parameters::PHYS_CHAR_DENSITY>();
    const T physThermalExpansion    = parameters.get<parameters::PHYS_THERMAL_EXPANSION>();
    const T physThermalConductivity = parameters.get<parameters::PHYS_THERMAL_CONDUCTIVITY>();
    const T physHeatCapacity        = parameters.get<parameters::PHYS_HEAT_CAPACITY>();
    const T g                       = parameters.get<parameters::GRAVITATIONAL_ACC>();
    const T smagoConst              = parameters.get<parameters::SMAGORINSKY_CONST>();
    const T prTurb                  = parameters.get<parameters::PRANDTL_TURB>();
    const T Tcold                   = parameters.get<parameters::T_COLD>();
    const T Thot                    = parameters.get<parameters::T_HOT>();

    NSElattice.setUnitConverter<ThermalUnitConverter<T,NSEDESCRIPTOR,ADEDESCRIPTOR>>(
        (T) physDeltaX,
        (T) physDeltaT,
        (T) physCharLength,
        (T) physCharVelocity,
        (T) physViscosity,
        (T) physDensity,
        (T) physThermalConductivity,
        (T) physHeatCapacity,
        (T) physThermalExpansion,
        (T) Tcold,
        (T) Thot
    );
    const auto& converter = NSElattice.getUnitConverter();
    converter.print();

    ADElattice.setUnitConverter(converter);

    NSElattice.defineDynamics<ExternalTauEffLESForcedBGKdynamics<T,NSEDESCRIPTOR,momenta::AdvectionDiffusionBulkTuple>>(geometry.getMaterialIndicator({1, 2, 3}));
    ADElattice.defineDynamics<ExternalTauEffLESBGKadvectionDiffusionDynamics>(geometry.getMaterialIndicator({1, 2, 3}));

    boundary::set<boundary::BounceBack>(ADElattice, geometry, 4);
    boundary::set<boundary::BounceBack>(NSElattice, geometry, 4);

    /// sets boundary
    boundary::set<boundary::AdvectionDiffusionDirichlet>(ADElattice, geometry.getMaterialIndicator({2, 3}));
    boundary::set<boundary::LocalVelocity>(NSElattice, geometry.getMaterialIndicator({2, 3}));

    T NSEomega  =  converter.getLatticeRelaxationFrequency();
    T ADEomega  =  converter.getLatticeThermalRelaxationFrequency();

    AnalyticalConst2D<T,T> tauNSE(1. / NSEomega);
    AnalyticalConst2D<T,T> tauADE(1. / ADEomega);

    NSElattice.defineField<descriptors::TAU_EFF>( geometry.getMaterialIndicator({1, 2, 3}), tauNSE );
    ADElattice.defineField<descriptors::TAU_EFF>( geometry.getMaterialIndicator({1, 2, 3}), tauADE );

    T boussinesqForcePrefactor = g / converter.getConversionFactorVelocity() * converter.getConversionFactorTime() *
                               converter.getCharPhysTemperatureDifference() * converter.getPhysThermalExpansionCoefficient();

    const T preFactor = smagoConst*smagoConst
                    * descriptors::invCs2<T,NSEDESCRIPTOR>()*descriptors::invCs2<T,NSEDESCRIPTOR>()
                    * 2*util::sqrt(2);

    auto& coupling = myCase.setCouplingOperator(
        "Boussinesq",
        SmagorinskyBoussinesqCoupling{},
        names::NavierStokes{}, NSElattice,
        names::Temperature{},  ADElattice
    );
    coupling.setParameter<SmagorinskyBoussinesqCoupling::T0>(
      converter.getLatticeTemperature(Tcold));
    coupling.setParameter<SmagorinskyBoussinesqCoupling::FORCE_PREFACTOR>(
      boussinesqForcePrefactor * Vector<T,2>{0.0,1.0});
    coupling.setParameter<SmagorinskyBoussinesqCoupling::SMAGORINSKY_PREFACTOR>(preFactor);
    coupling.setParameter<SmagorinskyBoussinesqCoupling::PR_TURB>(prTurb);
    coupling.setParameter<SmagorinskyBoussinesqCoupling::OMEGA_NSE>(
      NSEomega);
    coupling.setParameter<SmagorinskyBoussinesqCoupling::OMEGA_ADE>(
      ADEomega);

    clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase){
    OstreamManager clout(std::cout,"setInitialValues");
    clout << "Set initial values ..." << std::endl;

    using T = MyCase::value_t;

    auto& geometry = myCase.getGeometry();
    auto& NSElattice = myCase.getLattice(NavierStokes{});
    auto& ADElattice = myCase.getLattice(Temperature{});

    const auto& converter = NSElattice.getUnitConverter();
    T NSEomega = converter.getLatticeRelaxationFrequency();
    T ADEomega = converter.getLatticeThermalRelaxationFrequency();

    T Tcold = converter.getCharPhysLowTemperature();
    T Thot  = converter.getCharPhysHighTemperature();
    T Tmean = (Thot + Tcold) / 2.;

    /// define initial conditions
    AnalyticalConst2D<T,T> rho(1.);
    AnalyticalConst2D<T,T> u0(0.0, 0.0);
    AnalyticalConst2D<T,T> T_cold(converter.getLatticeTemperature(Tcold));
    AnalyticalConst2D<T,T> T_hot(converter.getLatticeTemperature(Thot));
    AnalyticalConst2D<T,T> T_mean(converter.getLatticeTemperature(Tmean));

    /// for each material set Rho, U and the Equilibrium
    NSElattice.defineRhoU(geometry.getMaterialIndicator({1, 2, 3}), rho, u0);
    NSElattice.iniEquilibrium(geometry.getMaterialIndicator({1, 2, 3}), rho, u0);

    ADElattice.defineRho(geometry, 1, T_mean);
    ADElattice.iniEquilibrium(geometry, 1, T_mean, u0);
    ADElattice.defineRho(geometry, 2, T_hot);
    ADElattice.iniEquilibrium(geometry, 2, T_hot, u0);
    ADElattice.defineRho(geometry, 3, T_cold);
    ADElattice.iniEquilibrium(geometry, 3, T_cold, u0);

    NSElattice.setParameter<descriptors::OMEGA>(NSEomega);
    ADElattice.setParameter<descriptors::OMEGA>(ADEomega);

    /// Make the lattice ready for simulation
    NSElattice.initialize();
    ADElattice.initialize();

    clout << "Set initial values ... OK" << std::endl;
}

void computeNusselt(MyCase& myCase){
  OstreamManager clout(std::cout, "computeNusselt");

  using T = MyCase::value_t_of<NavierStokes>;
  auto& geometry      = myCase.getGeometry();
  auto& NSElattice    = myCase.getLattice(NavierStokes{});
  auto& ADElattice    = myCase.getLattice(Temperature{});
  auto& converter     = NSElattice.getUnitConverter();
  auto& parameters    = myCase.getParameters();

  const int N = converter.getResolution();
  int material = 0;
  T T_x = 0, T_xplus1 = 0, T_xplus2 = 0, q = 0, voxel = 0;

  for (int iC = 0; iC < NSElattice.getLoadBalancer().size(); iC++) {
    int ny = NSElattice.getBlock(iC).getNy();
    int iX = 0;
    for (int iY = 0; iY < ny; ++iY) {
      material = geometry.getBlockGeometry(iC).getMaterial(iX,iY);

      T_x = ADElattice.getBlock(iC).get(iX,iY).computeRho();
      T_xplus1 = ADElattice.getBlock(iC).get(iX+1,iY).computeRho();
      T_xplus2 = ADElattice.getBlock(iC).get(iX+2,iY).computeRho();

      if ( material == 2 ) {
        q += (3.0*T_x - 4.0*T_xplus1 + 1.0*T_xplus2)/2.0*N;
        voxel++;
      }
    }
  }

  #ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(q, MPI_SUM);
    singleton::mpi().reduceAndBcast(voxel, MPI_SUM);
  #endif
  parameters.set<parameters::NUSSELT>(q / (T)voxel);
}

void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                int iT)
{
  OstreamManager clout(std::cout, "getResults");
  using NSEDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ADEDESCRIPTOR = MyCase::descriptor_t_of<Temperature>;
  using T = MyCase::value_t_of<NavierStokes>;
  auto& NSElattice        = myCase.getLattice(NavierStokes{});
  auto& ADElattice        = myCase.getLattice(Temperature{});
  const auto& converter   = NSElattice.getUnitConverter();
  auto& parameters        = myCase.getParameters();
  const int statIter      = converter.getLatticeTime(parameters.get<parameters::STAT_ITER>());
  const int vtkIter       = converter.getLatticeTime(parameters.get<parameters::VTK_ITER>());
  const bool converged    = parameters.get<parameters::CONVERGED>();

  const T Thot            = parameters.get<parameters::T_HOT>();
  const T Tcold           = parameters.get<parameters::T_COLD>();
  const T lx              = parameters.get<parameters::DOMAIN_EXTENT>()[0];

  SuperVTMwriter2D<T> vtkWriter("squareCavity2dTurbulent");

  SuperLatticePhysVelocity2D<T, NSEDESCRIPTOR> velocity(NSElattice, converter);
  SuperLatticePhysPressure2D<T, NSEDESCRIPTOR> pressure(NSElattice, converter);
  SuperLatticePhysTemperature2D<T, NSEDESCRIPTOR, ADEDESCRIPTOR> temperature(ADElattice, converter);
  AnalyticalFfromSuperF2D<T> interpolation(velocity, true);

  if (iT == 0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, NSEDESCRIPTOR> cuboid(NSElattice);
    SuperLatticeRank2D<T, NSEDESCRIPTOR> rank(NSElattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.createMasterFile();
    vtkWriter.write(iT);
  }

  if ((iT % vtkIter == 0 && iT > 0) || converged) {
    vtkWriter.addFunctor(pressure);
    vtkWriter.addFunctor(velocity);
    vtkWriter.addFunctor(temperature);
    vtkWriter.write(iT);
  }

  if (iT % statIter == 0 || converged) {
    NSElattice.setProcessingContext(ProcessingContext::Evaluation);
    ADElattice.setProcessingContext(ProcessingContext::Evaluation);

    timer.update(iT);
    timer.printStep();

    /// NSLattice statistics console output
    NSElattice.getStatistics().print(iT,converter.getPhysTime(iT));
    /// ADElattice statistics console output
    ADElattice.getStatistics().print(iT,converter.getPhysTime(iT));


    BlockReduction2D2D<T> planeReduction(temperature, 600, BlockDataSyncMode::ReduceOnly);
    BlockGifWriter<T> gifWriter;
    gifWriter.write(planeReduction, Tcold-1, Thot+1, iT, "temperature");

    SuperEuklidNorm2D<T, NSEDESCRIPTOR> normVel( velocity );
    BlockReduction2D2D<T> planeReduction2(normVel, 600, BlockDataSyncMode::ReduceOnly);
    BlockGifWriter<T> gifWriter2;
    gifWriter2.write( planeReduction2, iT, "velocity" );
  }

  if ( converged ) {
    computeNusselt(myCase);

    /// Initialize vectors for data output
    T xVelocity[2] = { T() };
    T outputVelX[2] = { T() };
    T yVelocity[2] = { T() };
    T outputVelY[2] = { T() };
    const int outputSize = 512;
    Vector<T, outputSize> velX;
    Vector<T, outputSize> posX;
    Vector<T, outputSize> velY;
    Vector<T, outputSize> posY;

    /// loop for the resolution of the cavity at x = lx/2 in yDirection and vice versa
    for (int n = 0; n < outputSize; ++n) {
      T yPosition[2] = { lx / 2, lx * n / (T) outputSize };
      T xPosition[2] = { lx * n / (T) outputSize, lx / 2 };

      /// Interpolate xVelocity at x = lx/2 for each yPosition
      interpolation(xVelocity, yPosition);
      interpolation(yVelocity, xPosition);
      /// Store the interpolated values to compare them among each other in order to detect the maximum
      velX[n] = xVelocity[0];
      posY[n] = yPosition[1];
      velY[n] = yVelocity[1];
      posX[n] = xPosition[0];

      /// Initialize output with the corresponding velocities and positions at the origin
      if (n == 0) {
        outputVelX[0] = velX[0];
        outputVelX[1] = posY[0];
        outputVelY[0] = velY[0];
        outputVelY[1] = posX[0];
      }
      /// look for the maximum velocity in xDirection and the corresponding position in yDirection
      if (n > 0 && velX[n] > outputVelX[0]) {
        outputVelX[0] = velX[n];
        outputVelX[1] = posY[n];
      }
      /// look for the maximum velocity in yDirection and the corresponding position in xDirection
      if (n > 0 && velY[n] > outputVelY[0]) {
        outputVelY[0] = velY[n];
        outputVelY[1] = posX[n];
      }
    }

    parameters.set<parameters::SIM_VALUES>(
      {outputVelX[0],
       outputVelY[0],
       outputVelX[1],
       outputVelY[1],
       parameters.get<parameters::NUSSELT>()}
    );
  }
}

void simulate(MyCase& myCase){
    OstreamManager clout(std::cout,"Simulation");
    clout << "Starting Simulation ..." << std::endl;

    using T = MyCase::value_t;
    auto& parameters      = myCase.getParameters();
    auto& NSElattice      = myCase.getLattice(NavierStokes{});
    auto& ADElattice      = myCase.getLattice(Temperature{});
    const auto& converter = NSElattice.getUnitConverter();

    const int iTmax       = converter.getLatticeTime(parameters.get<parameters::MAX_PHYS_T>());

    util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());

    timer.start();

    const int convIter = parameters.get<parameters::CONV_ITER>();
    util::ValueTracer<T> converge(6, parameters.get<parameters::CONVERGENCE_PRECISION>());

    for (int iT = 0; iT < iTmax; ++iT) {

      if (converge.hasConverged() && !parameters.get<parameters::CONVERGED>()) {
        parameters.set<parameters::CONVERGED>(true);
        clout << "Simulation converged." << std::endl;
        clout << "Time " << iT << "." << std::endl;

        getResults(myCase, timer, iT);
        break;
      }

      NSElattice.collideAndStream();
      ADElattice.collideAndStream();
      myCase.getOperator("Boussinesq").apply();
      getResults(myCase, timer, iT);
      if(!parameters.get<parameters::CONVERGED>() && iT % convIter == 0){
        ADElattice.setProcessingContext(ProcessingContext::Evaluation);
        computeNusselt(myCase);
        converge.takeValue(parameters.get<parameters::NUSSELT>(), true);
      }
    }

    timer.stop();
    timer.printSummary();
}
