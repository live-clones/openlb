/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2008 Orestis Malaspinas
 *                2025 Adrian Kummerlaender
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

#include <olb.h>

#include "fringeZone.h"

const std::string xmlPath = "decomposition/n9.xml"; 

using namespace olb;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q27<>>
>;

enum class BulkModel {

  RLB = 0,
  SmagorinskyBGK = 1,
  LocalSmagorinskyBGK = 2,
  ShearSmagorinskyBGK = 3,
  ConsistentStrainSmagorinskyBGK = 4,
  KrauseBGK = 5,

};

namespace olb::parameters {

  struct BULK_MODEL              : public descriptors::TYPED_FIELD_BASE<int, 1> { };
  struct INLET_CYLINDER_SIZE     : public descriptors::FIELD_BASE<2> { }; // 0 -> length; 1 -> radius and char length

  struct SMAGORINSKY             : public descriptors::FIELD_BASE<1> { };
  struct LATTICE_RELAXATION_TIME : public descriptors::FIELD_BASE<1> { };

  struct TURBULENCE_INTENSITY    : public descriptors::FIELD_BASE<1> { };
  struct TURBULENCE_N_SEEDS      : public descriptors::FIELD_BASE<1> { };
  struct TURBULENCE_N_TIME       : public descriptors::FIELD_BASE<1> { };
  struct TURBULENCE_SIGMA        : public descriptors::FIELD_BASE<1> { };

  struct STAT_ITER               : public descriptors::TYPED_FIELD_BASE<int, 1> { };
  struct VTK_SAVE_ITER           : public descriptors::TYPED_FIELD_BASE<int, 1> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  std::vector<T> origin(3,T());
  IndicatorCuboid3D<T> cuboid(extent, origin);

  const T physDeltaX = params.get<parameters::PHYS_DELTA_X>(); 
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  return mesh;
}

/// Constructs indicator of inlet tube geometry
std::shared_ptr<IndicatorF3D<MyCase::value_t>> makeInletI(MyCase::ParametersD& params)
{
  using T = MyCase::value_t;
  Vector<T,3> origin(
    T(),
    params.get<parameters::DOMAIN_EXTENT>()[1] / 2.0 + params.get<parameters::PHYS_DELTA_X>(),
    params.get<parameters::DOMAIN_EXTENT>()[2] / 2.0 + params.get<parameters::PHYS_DELTA_X>()
  );
  Vector<T,3> extend(
    params.get<parameters::INLET_CYLINDER_SIZE>()[0] + 5 * params.get<parameters::PHYS_DELTA_X>(),
    params.get<parameters::DOMAIN_EXTENT>()[1] / 2.0 + params.get<parameters::PHYS_DELTA_X>(),
    params.get<parameters::DOMAIN_EXTENT>()[2] / 2.0 + params.get<parameters::PHYS_DELTA_X>()
  );
  return std::shared_ptr<IndicatorF3D<T>>(
    new IndicatorCylinder3D<T>(extend, origin, params.get<parameters::INLET_CYLINDER_SIZE>()[1])
  );
}

/// Constructs indicator of injection tube geometry
std::shared_ptr<IndicatorF3D<MyCase::value_t>> makeInjectionTubeI(MyCase::ParametersD& params)
{
  using T = MyCase::value_t;
  Vector<T,3> origin(
    params.get<parameters::INLET_CYLINDER_SIZE>()[0],
    params.get<parameters::DOMAIN_EXTENT>()[1] / 2.0 + params.get<parameters::PHYS_DELTA_X>(),
    params.get<parameters::DOMAIN_EXTENT>()[2] / 2.0 + params.get<parameters::PHYS_DELTA_X>()
  );
  Vector<T,3> extend(
    params.get<parameters::DOMAIN_EXTENT>()[0],
    params.get<parameters::DOMAIN_EXTENT>()[1] / 2.0 + params.get<parameters::PHYS_DELTA_X>(),
    params.get<parameters::DOMAIN_EXTENT>()[2] / 2.0 + params.get<parameters::PHYS_DELTA_X>()
  );
  return std::shared_ptr<IndicatorF3D<T>>(
    new IndicatorCylinder3D<T>(
      extend, 
      origin, 
      std::min(params.get<parameters::DOMAIN_EXTENT>()[1], params.get<parameters::DOMAIN_EXTENT>()[2]))
  );
}

std::shared_ptr<IndicatorF3D<MyCase::value_t>> makeNozzleI(MyCase::ParametersD& params)
{
  auto inletCylinder = makeInletI(params);
  auto injectionTube = makeInjectionTubeI(params);
  return inletCylinder + injectionTube;
}

void prepareGeometry(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  geometry.rename(0, 2);

  auto inletCylinder = makeInletI(myCase.getParameters());
  geometry.rename(2, 1, *inletCylinder);

  auto injectionTube = makeInjectionTubeI(myCase.getParameters());
  geometry.rename(2, 1, *injectionTube);

  {
    Vector<T,3> origin(T(),
      params.get<parameters::DOMAIN_EXTENT>()[1] / 2.0  + params.get<parameters::PHYS_DELTA_X>(),
      params.get<parameters::DOMAIN_EXTENT>()[2] / 2.0  + params.get<parameters::PHYS_DELTA_X>()
    );
    Vector<T,3> extend(params.get<parameters::PHYS_DELTA_X>(),
      params.get<parameters::DOMAIN_EXTENT>()[1] / 2.0  + params.get<parameters::PHYS_DELTA_X>(),
      params.get<parameters::DOMAIN_EXTENT>()[2] / 2.0  + params.get<parameters::PHYS_DELTA_X>()
    );

    IndicatorCylinder3D<T> cylinderIN(extend, origin, params.get<parameters::INLET_CYLINDER_SIZE>()[1]);
    geometry.rename(1,3, cylinderIN);
  }

  {
    Vector<T,3> origin(params.get<parameters::INLET_CYLINDER_SIZE>()[0] - params.get<parameters::PHYS_DELTA_X>(),
      params.get<parameters::DOMAIN_EXTENT>()[1] / 2.0  + params.get<parameters::PHYS_DELTA_X>(),
      params.get<parameters::DOMAIN_EXTENT>()[2] / 2.0  + params.get<parameters::PHYS_DELTA_X>()
    );
    Vector<T,3> extend(params.get<parameters::INLET_CYLINDER_SIZE>()[0],
      params.get<parameters::DOMAIN_EXTENT>()[1] / 2.0  + params.get<parameters::PHYS_DELTA_X>(),
      params.get<parameters::DOMAIN_EXTENT>()[2] / 2.0  + params.get<parameters::PHYS_DELTA_X>()
    );

    IndicatorCylinder3D<T> cylinderOUT(extend, origin, params.get<parameters::DOMAIN_EXTENT>()[2]);
    geometry.rename(1,4, cylinderOUT);
  }

  geometry.clean();
  geometry.innerClean();
  geometry.checkForErrors();

  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase) {
  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});

  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

  const T resolution = params.get<parameters::RESOLUTION>();
  const T latticeRelaxationTime = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T physCharLength = params.get<parameters::INLET_CYLINDER_SIZE>()[1];
  const T physCharVelocity = params.get<parameters::PHYS_CHAR_VELOCITY>();
  const T physViscosity = params.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T physDensity = params.get<parameters::PHYS_CHAR_DENSITY>();

  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    resolution,              // resolution: number of voxels per charPhysL
    latticeRelaxationTime,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    physCharLength,          // charPhysLength: reference length of simulation geometry
    physCharVelocity,        // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    physViscosity,           // physViscosity: physical kinematic viscosity in __m^2 / s__
    physDensity              // physDensity: physical density in __kg / m^3__
  );

  const auto& converter = lattice.getUnitConverter();
  converter.print();

  // Material=0 -->do nothing
  lattice.defineDynamics<NoDynamics>(geometry, 0);

  // Material=1 -->bulk dynamics
  // Material=3 -->bulk dynamics (inflow)
  // Material=4 -->bulk dynamics (outflow)
  auto bulkIndicator = geometry.getMaterialIndicator({1, 3, 4});

  switch (params.get<parameters::BULK_MODEL>()) {
    case static_cast<int>(BulkModel::RLB):
      lattice.defineDynamics<RLBdynamics>(bulkIndicator);
      break;
    case static_cast<int>(BulkModel::ShearSmagorinskyBGK):
      lattice.defineDynamics<ShearSmagorinskyBGKdynamics>(bulkIndicator);
      lattice.setParameter<collision::LES::SMAGORINSKY>(params.get<parameters::SMAGORINSKY>());
      break;
    case static_cast<int>(BulkModel::KrauseBGK):
      lattice.defineDynamics<KrauseBGKdynamics>(bulkIndicator);
      lattice.setParameter<collision::LES::SMAGORINSKY>(params.get<parameters::SMAGORINSKY>());
    case static_cast<int>(BulkModel::ConsistentStrainSmagorinskyBGK):
      lattice.defineDynamics<ConStrainSmagorinskyBGKdynamics>(bulkIndicator);
      lattice.setParameter<collision::LES::SMAGORINSKY>(params.get<parameters::SMAGORINSKY>());
      break;
    case static_cast<int>(BulkModel::SmagorinskyBGK):
      lattice.defineDynamics<SmagorinskyBGKdynamics>(bulkIndicator);
      lattice.setParameter<collision::LES::SMAGORINSKY>(params.get<parameters::SMAGORINSKY>());
      break;
    case static_cast<int>(BulkModel::LocalSmagorinskyBGK):
    default:
      lattice.defineDynamics<LocalSmagorinskyBGKdynamics>(bulkIndicator);
      FringeZoneSmagorinskyConstant smagorinskyFringe(converter, T{0.15});
      lattice.defineField<collision::LES::SMAGORINSKY>(bulkIndicator, smagorinskyFringe);
      break;
  }

  // Material=2 -->bounce back
  lattice.defineDynamics<BounceBack>(geometry, 2);

  const T omega = converter.getLatticeRelaxationFrequency();
  boundary::set<boundary::LocalVelocity>(lattice, geometry, 3);
  boundary::set<boundary::InterpolatedConvection>(lattice, geometry, 4);

  lattice.setParameter<descriptors::OMEGA>(omega);

  {
    auto& communicator = lattice.getCommunicator(stage::PostCollide());
    communicator.clearRequestedCells();
    communicator.requestOverlap(1, geometry.getMaterialIndicator({1,2,3,4}));
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(
  MyCase& myCase
) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});

  AnalyticalConst3D<T,T> rhoF(1);
  Vector<T,3> velocity{};
  AnalyticalConst3D<T,T> uF(velocity);

  lattice.defineRhoU(geometry.getMaterialIndicator({1, 2, 3, 4}), rhoF, uF);
  lattice.iniEquilibrium(geometry.getMaterialIndicator({1, 2, 3, 4}), rhoF, uF);
  lattice.initialize();
  geometry.updateStatistics();

}

void setTemporalValues(
  MyCase& myCase,
  std::size_t iT,
  VortexMethodTurbulentVelocityBoundary<MyCase::value_t, MyCase::descriptor_t_of<NavierStokes>>& vortex
) { 
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});
  const auto& converter = lattice.getUnitConverter();

  if (params.get<parameters::BULK_MODEL>() == static_cast<int>(BulkModel::ShearSmagorinskyBGK)) {
    lattice.setParameter<descriptors::LATTICE_TIME>(iT);
  }

  const T startUpFactor = 0.001 * 1.0 / converter.getCharLatticeVelocity();
  const T maxStartTPhys = (60*converter.getCharPhysLength()/converter.getCharPhysVelocity()) * startUpFactor;
  const auto maxStartT = converter.getLatticeTime(maxStartTPhys);

  auto uSol = std::shared_ptr<AnalyticalF3D<T,T>>(new CirclePowerLaw3D<T>(geometry, 3, converter.getCharLatticeVelocity(), 8, T()));

  if (iT <= maxStartT) {
    PolynomialStartScale<T,std::size_t> scale(maxStartT, 1);
    T frac{};
    scale(&frac, &iT);

    vortex.setVelocityProfile(converter.getConversionFactorVelocity() * frac * uSol);
    lattice.template setProcessingContext<Array<U_PROFILE>>(ProcessingContext::Simulation);
    vortex.apply(iT);
  } else {
    vortex.apply(iT);
  }
}

void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& params = myCase.getParameters();
  const auto& converter = lattice.getUnitConverter();

  SuperVTMwriter3D<T> vtkWriter("nozzle3d");

  if (iT==0) {
    SuperLatticeCuboid3D cuboid(lattice);
    SuperLatticeRank3D rank(lattice);
    SuperLatticePlatform platform(lattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.write(platform);
    vtkWriter.createMasterFile();
  }

  if (iT % params.get<parameters::STAT_ITER>() == 0) {
    timer.update(iT);
    timer.printStep();
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
    if (std::isnan(lattice.getStatistics().getAverageRho())) {
      std::exit(-1);
    }
  }

  // Writes the vtk files
  if (iT % params.get<parameters::VTK_SAVE_ITER>() == 0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    lattice.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtkWriter("nozzle3d");
      SuperLatticePhysVelocity3D velocity(lattice, converter);
      SuperLatticePhysPressure3D pressure(lattice, converter);
      vtkWriter.addFunctor(velocity);
      vtkWriter.addFunctor(pressure);
      task(vtkWriter, iT);
    });
  }
}

void simulate(
  MyCase& myCase
) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& params = myCase.getParameters();
  const auto& converter = lattice.getUnitConverter();

  /// === Vortex Method Setup (Turbulent inlet condition)
  Vector<T,3> originI(
    params.get<parameters::PHYS_DELTA_X>(), 
    params.get<parameters::DOMAIN_EXTENT>()[1] / 2.0 + params.get<parameters::PHYS_DELTA_X>(), 
    params.get<parameters::DOMAIN_EXTENT>()[2] / 2.0 + params.get<parameters::PHYS_DELTA_X>()
  );
  Vector<T,3> extendI(
    T(), 
    params.get<parameters::DOMAIN_EXTENT>()[1] / 2.0 + params.get<parameters::PHYS_DELTA_X>(), 
    params.get<parameters::DOMAIN_EXTENT>()[2] / 2.0 + params.get<parameters::PHYS_DELTA_X>()
  );
  Vector<T,3> inflowAxis{1, 0, 0};
  IndicatorCylinder3D<T> cylinderIN(extendI, originI, 1.);
  VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR> vortex(
    geometry.getMaterialIndicator(3),
    cylinderIN,
    converter,
    lattice,
    params.get<parameters::TURBULENCE_N_SEEDS>(),                                
    params.get<parameters::TURBULENCE_N_TIME>(),                               
    params.get<parameters::TURBULENCE_SIGMA>(), 
    inflowAxis
  );
  
  auto intensity = std::shared_ptr<AnalyticalF3D<T,T>>(
    new AnalyticalConst3D<T,T>(params.get<parameters::TURBULENCE_INTENSITY>())
  );
  vortex.setIntensityProfile(intensity);

  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(
    params.get<parameters::MAX_PHYS_T>()
  );

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT, vortex);

    /// === Step 8.2: Collide and Stream Execution ===
    myCase.getLattice(NavierStokes{}).collideAndStream();

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, timer, iT);
  }

  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[]) {
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<BULK_MODEL>(static_cast<int>(BulkModel::LocalSmagorinskyBGK));

    myCaseParameters.set<DOMAIN_EXTENT>({60, 11, 11});
    myCaseParameters.set<INLET_CYLINDER_SIZE>({4, 1});

    myCaseParameters.set<RESOLUTION   >(5);
    myCaseParameters.set<PHYS_DELTA_X >(
      myCaseParameters.get<INLET_CYLINDER_SIZE>()[1] / myCaseParameters.get<RESOLUTION>()
    );
    myCaseParameters.set<PHYS_CHAR_VELOCITY>(1);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(0.0002);
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1);
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.500018);
    myCaseParameters.set<SMAGORINSKY>(0.15);
    myCaseParameters.set<MAX_PHYS_T   >(150);

    myCaseParameters.set<TURBULENCE_INTENSITY>(0.05);
    myCaseParameters.set<TURBULENCE_N_SEEDS>(50);
    myCaseParameters.set<TURBULENCE_N_TIME>(0.1);
    myCaseParameters.set<TURBULENCE_SIGMA>(
      0.1 * myCaseParameters.get<INLET_CYLINDER_SIZE>()[1]
    );

    myCaseParameters.set<STAT_ITER>(200);
    myCaseParameters.set<VTK_SAVE_ITER>(500);
  }
  myCaseParameters.fromCLI(argc, argv);

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase);

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);
}
