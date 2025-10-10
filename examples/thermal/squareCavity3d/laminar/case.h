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

/* Header related to squareCavity3dLaminar.cpp
 * The reference is the paper in "Gaedtke, M., Wachter, S., Raedle, M., Nirschl, H., & Krause, M. J. (2018).
 * Application of a lattice Boltzmann method combined with a Smagorinsky turbulence model to spatially resolved heat flux inside a refrigerated vehicle.
 * Computers & Mathematics with Applications, 76(10), 2315-2329."
 */

// natural convection of air in a square cavity in 3D

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::names;

namespace olb::parameters {
  struct SIM_VALUES  : public descriptors::FIELD_BASE<5> { };
  struct N_CELLS_Z   : public descriptors::FIELD_BASE<1> { };
}

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q19<descriptors::FORCE>>,
  Temperature,  Lattice<double, descriptors::D3Q7<descriptors::VELOCITY>>
>;

/// @brief Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters){
    using T = MyCase::value_t;

    const T dx                = parameters.get<parameters::PHYS_DELTA_X>();
    const Vector domainExtend = parameters.get<parameters::DOMAIN_EXTENT>();
    const T cellsZ            = parameters.get<parameters::N_CELLS_Z>();
    const Vector extend{domainExtend[0] + dx, domainExtend[1] + dx, cellsZ * dx };
    const Vector origin{-dx / 2., -dx / 2., -dx / 2.};
    IndicatorCuboid3D<T> cuboid(extend, origin);

    Mesh<T,MyCase::d> mesh(cuboid, dx, singleton::mpi().getSize());
    mesh.setOverlap(parameters.get<parameters::OVERLAP>());
    mesh.getCuboidDecomposition().setPeriodicity({false, false, true});

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
    const T cellsZ            = parameters.get<parameters::N_CELLS_Z>();
    const T lx                = domainExtend[0];

    geometry.rename(0, 4);

    Vector origin {0, 0, 0};
    Vector extend {lx, lx, cellsZ * dx};
    IndicatorCuboid3D<T> cuboid2(extend, origin);

    geometry.rename(4, 1, cuboid2);

    Vector extendwallleft{dx, lx + dx, cellsZ * dx};
    Vector originwallleft{-dx / 2., -dx / 2., -dx / 2.};
    IndicatorCuboid3D<T> wallleft(extendwallleft, originwallleft);

    Vector extendwallright{dx, lx + dx, cellsZ * dx};
    Vector originwallright{lx - dx / 2., -dx / 2., -dx / 2.};
    IndicatorCuboid3D<T> wallright(extendwallright, originwallright);

    geometry.rename(4, 2, 1, wallleft  );
    geometry.rename(4, 3, 1, wallright );

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
    const T tau                     = parameters.get<parameters::LATTICE_RELAXATION_TIME>();
    const T physViscosity           = parameters.get<parameters::PHYS_KINEMATIC_VISCOSITY>();
    const T physDeltaX              = parameters.get<parameters::PHYS_DELTA_X>();
    const T physDeltaT              = (tau - 0.5) / descriptors::invCs2<T,NSEDESCRIPTOR>() * physDeltaX * physDeltaX / physViscosity;
    const T physCharVelocity        = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
    const T physDensity             = parameters.get<parameters::PHYS_CHAR_DENSITY>();
    const T physThermalExpansion    = parameters.get<parameters::PHYS_THERMAL_EXPANSION>();
    const T physThermalConductivity = parameters.get<parameters::PHYS_THERMAL_CONDUCTIVITY>();
    const T physHeatCapacity        = parameters.get<parameters::PHYS_HEAT_CAPACITY>();
    const T g                       = parameters.get<parameters::GRAVITATIONAL_ACC>();
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

    NSElattice.defineDynamics<ForcedBGKdynamics>(geometry.getMaterialIndicator({1, 2, 3}));
    ADElattice.defineDynamics<AdvectionDiffusionBGKdynamics>(geometry.getMaterialIndicator({1, 2, 3}));

    boundary::set<boundary::BounceBack>(ADElattice, geometry, 4);
    boundary::set<boundary::BounceBack>(NSElattice, geometry, 4);

    boundary::set<boundary::AdvectionDiffusionDirichlet>(ADElattice, geometry.getMaterialIndicator({2, 3}));
    boundary::set<boundary::LocalVelocity>(NSElattice, geometry.getMaterialIndicator({2, 3}));

    T boussinesqForcePrefactor = g / converter.getConversionFactorVelocity() * converter.getConversionFactorTime() *
                               converter.getCharPhysTemperatureDifference() * converter.getPhysThermalExpansionCoefficient();

    auto& coupling = myCase.setCouplingOperator(
      "Boussinesq",
      NavierStokesAdvectionDiffusionCoupling{},
      names::NavierStokes{}, NSElattice,
      names::Temperature{},  ADElattice
    );
    coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::T0>(
      converter.getLatticeTemperature(Tcold));
    coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::FORCE_PREFACTOR>(
      boussinesqForcePrefactor * Vector{0.0,1.0, 0.0}
    );

    clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase){
    OstreamManager clout(std::cout,"setInitialValues");
    clout << "Set initial values ..." << std::endl;

    using T               = MyCase::value_t;

    auto& geometry        = myCase.getGeometry();
    auto& NSElattice      = myCase.getLattice(NavierStokes{});
    auto& ADElattice      = myCase.getLattice(Temperature{});
    const auto& converter = NSElattice.getUnitConverter();

    T NSEomega = converter.getLatticeRelaxationFrequency();
    T ADEomega = converter.getLatticeThermalRelaxationFrequency();

    T Tcold = converter.getCharPhysLowTemperature();
    T Thot  = converter.getCharPhysHighTemperature();
    T Tmean = (Thot + Tcold) / 2.;

    /// define initial conditions
    AnalyticalConst3D<T,T> rho(1.);
    AnalyticalConst3D<T,T> u0(0.0, 0.0, 0.0);
    AnalyticalConst3D<T,T> T_cold(converter.getLatticeTemperature(Tcold));
    AnalyticalConst3D<T,T> T_hot(converter.getLatticeTemperature(Thot));
    AnalyticalConst3D<T,T> T_mean(converter.getLatticeTemperature(Tmean));

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
  int material = 0, voxel = 0;
  T T_x = 0, T_xplus1 = 0, T_xplus2 = 0, q = 0;

  for (int iC = 0; iC < NSElattice.getLoadBalancer().size(); iC++) {
    int ny = NSElattice.getBlock(iC).getNy();

    int iX = 0;
    int iZ = 1;

    for (int iY = 0; iY < ny; ++iY) {
      material = geometry.getBlockGeometry(iC).getMaterial(iX,iY,iZ);

      T_x = ADElattice.getBlock(iC).get(iX,iY,iZ).computeRho();
      T_xplus1 = ADElattice.getBlock(iC).get(iX+1,iY,iZ).computeRho();
      T_xplus2 = ADElattice.getBlock(iC).get(iX+2,iY,iZ).computeRho();

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
  const int statIter      = converter.getLatticeTime(parameters.get<parameters::PHYS_STAT_ITER_T>());
  const int vtkIter       = converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T>());
  const bool converged    = parameters.get<parameters::CONVERGED>();

  const T Thot            = parameters.get<parameters::T_HOT>();
  const T Tcold           = parameters.get<parameters::T_COLD>();
  const T lx              = parameters.get<parameters::DOMAIN_EXTENT>()[0];

  SuperVTMwriter3D<T> vtkWriter("squareCavity3dLaminar");

  SuperLatticePhysVelocity3D<T, NSEDESCRIPTOR> velocity(NSElattice, converter);
  SuperLatticePhysPressure3D<T, NSEDESCRIPTOR> pressure(NSElattice, converter);
  SuperLatticePhysTemperature3D<T, NSEDESCRIPTOR, ADEDESCRIPTOR> temperature(ADElattice, converter);
  AnalyticalFfromSuperF3D<T> interpolation(velocity, true);

  if (iT == 0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, NSEDESCRIPTOR> cuboid(NSElattice);
    SuperLatticeRank3D<T, NSEDESCRIPTOR> rank(NSElattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.createMasterFile();
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

    /// NSElattice statistics console output
    NSElattice.getStatistics().print(iT,converter.getPhysTime(iT));
    /// ADElattice statistics console output
    ADElattice.getStatistics().print(iT,converter.getPhysTime(iT));

    const double a[3] = {0, 0, 1.};
    BlockReduction3D2D<T> planeReduction(temperature, a);
    BlockGifWriter<T> gifWriter;
    gifWriter.write(planeReduction, Tcold*0.98, Thot*1.02, iT, "temperature");

    SuperEuklidNorm3D<T> normVel( velocity );
    BlockReduction3D2D<T> planeReduction2(normVel, {0, 0, 1});
    BlockGifWriter<T> gifWriter2;
    gifWriter2.write( planeReduction2, iT, "velocity" );
  }

  if ( converged ) {
    computeNusselt(myCase);

    T xVelocity[3] = { T() };
    T outputVelX[3] = { T() };
    T yVelocity[3] = { T() };
    T outputVelY[3] = { T() };

    const int outputSize = 512;
    Vector<T, outputSize> velX;
    Vector<T, outputSize> posX;
    Vector<T, outputSize> velY;
    Vector<T, outputSize> posY;

    const T N           = parameters.get<parameters::RESOLUTION>();

    /// loop for the resolution of the cavity at x = lx/2 in yDirection and vice versa
    for (int n = 0; n < outputSize; ++n) {
      T yPosition[3] = { lx / 2, lx * n / (T) outputSize, lx / N * 2 / 2 };
      T xPosition[3] = { lx * n / (T) outputSize, lx / 2, lx / N * 2 / 2 };

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
    auto& parameters = myCase.getParameters();
    auto& NSElattice = myCase.getLattice(NavierStokes{});
    auto& ADElattice = myCase.getLattice(Temperature{});
    const auto& converter = NSElattice.getUnitConverter();

    const int iTmax = converter.getLatticeTime(parameters.get<parameters::MAX_PHYS_T>());

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
