/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2019 Mathias J. Krause, Adrian Kummerlaender
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net
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

#include "olb.h"

using namespace olb;
using namespace olb::names;
//ys using DESCRIPTOR   = descriptors::D3Q19<>;
//ys using T            = float;
//ys using BulkDynamics = BGKdynamics<T, DESCRIPTOR>;

// === Step 1: Declarations ===
using MyCase = Case<NavierStokes, Lattice<double, descriptors::D3Q19<>>>;

namespace olb::parameters {
struct NO_EXPORT_RESULTS : public descriptors::TYPED_FIELD_BASE<bool, 1> {};
struct CUBOIDS_PER_PROCESS : public descriptors::TYPED_FIELD_BASE<std::size_t, 1> {};

} // namespace olb::parameters

// Undefine to test a minimal bounce back cavity
#define LID_DRIVEN
//// Use bounce back (velocity) boundaries instead of local velocity
//#define LID_DRIVEN_BOUNCE_BACK
//// Use single fused collision kernel instead of individual dispatch on GPUs
//#define GPU_USE_FUSED_COLLISION

/// @brief Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters)
{
  using T                 = MyCase::value_t;
  const T      physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();
  const T      length = parameters.get<parameters::PHYS_CHAR_LENGTH>(); // length of the cavity in x- and z-direction
  Vector<T, 3> origin((T)0);
  Vector<T, 3> extend(length + 0.5 * physDeltaX);
  IndicatorCuboid3D<T> cuboid(extend, origin);

  Mesh<T, MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

/// @brief Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
/// @note The material numbers will be used to assign physics to lattice nodes
void prepareGeometry(MyCase& myCase)
{

  OstreamManager clout(std::cout, "prepareGeometry");
  //clout << "Prepare Geometry ..." << std::endl;
  using T          = MyCase::value_t;
  auto& sGeometry  = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();
  const T length     = parameters.get<parameters::PHYS_CHAR_LENGTH>(); // length of the cavity in x- and z-direction

  Vector<T, 3>         origin((T)0);
  Vector<T, 3>         extend(length + 0.25 * physDeltaX);
  IndicatorCuboid3D<T> cuboid(extend, origin);

  sGeometry.rename(0, 2, cuboid);
  sGeometry.rename(2, 1, {1, 1, 1});

  T                    eps = physDeltaX;
  Vector<T, 3>         libOrigin(-eps, length - eps, -eps);
  Vector<T, 3>         libExtend(length + 2 * eps, 2 * eps, length + 2 * eps);
  IndicatorCuboid3D<T> lid(libExtend, libOrigin);

  sGeometry.rename(2, 3, 1, lid);

  const bool verbose = false;
  sGeometry.clean(verbose);

  sGeometry.innerClean(verbose);
  sGeometry.checkForErrors(verbose);

  //  sGeometry.print();

  //clout << "Prepare Geometry ... OK" << std::endl;
  return;
}

/// @brief Set lattice dynamics
/// @param myCase The Case instance which keeps the simulation data
void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  using T          = MyCase::value_t;
  auto& sGeometry  = myCase.getGeometry();
  auto& sLattice   = myCase.getLattice(NavierStokes {});
  auto& parameters = myCase.getParameters();
  {
    using namespace olb::parameters;
    sLattice.setUnitConverter<UnitConverterFromResolutionAndLatticeVelocity<T, MyCase::descriptor_t_of<NavierStokes>>>(
        parameters.get<RESOLUTION>(), // resolution: number of voxels per charPhysL
        0.01,                         //lattice_char_velocity : characteristic non-dimensional(lattice) valocity
                                      //ys         parameters
        //ys .get<LATTICE_CHAR_VELOCITY>(),  //lattice_char_velocity : characteristic non-dimensional(lattice) valocity
        parameters.get<PHYS_CHAR_LENGTH>(), // charPhysLength: reference length of simulation geometry
        1.0,        // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
        1. / 1000., // physViscosity: physical kinematic viscosity in __m^2 / s__
        1.0         // physDensity: physical density in __kg / m^3__
    );
  }

  /// Material=1 -->bulk dynamics
  sLattice.defineDynamics<BGKdynamics>(sGeometry, 1);

#ifdef LID_DRIVEN
#ifdef LID_DRIVEN_BOUNCE_BACK
  boundary::set<boundary::BounceBack>(sLattice, sGeometry, 2);
  sLattice.defineDynamics<BounceBackVelocity>(sGeometry, 3);
#else // Local velocity boundaries
  boundary::set<boundary::LocalVelocity>(sLattice, sGeometry, 2);
  boundary::set<boundary::LocalVelocity>(sLattice, sGeometry, 3);
#endif
#else
  boundary::set<boundary::BounceBack>(sLattice, sGeometry, 2);
  boundary::set<boundary::BounceBack>(sLattice, sGeometry, 3);
#endif

  // Alternative GPU-specific performance tuning option
#if defined(PLATFORM_GPU_CUDA) && defined(GPU_USE_FUSED_COLLISION)
#ifdef LID_DRIVEN_BOUNCE_BACK
  // Enable non-virtual dispatching of common collision operators (optional, improves performance)
  sLattice.forBlocksOnPlatform<Platform::GPU_CUDA>([](auto& block) {
    block.setCollisionO(gpu::cuda::getFusedCollisionO<T, DESCRIPTOR, BulkDynamics, BounceBack<T, DESCRIPTOR>,
                                                      BounceBackVelocity<T, DESCRIPTOR>>());
  });
#else // Local velocity boundaries
  sLattice.forBlocksOnPlatform<Platform::GPU_CUDA>([](auto& block) {
    block.setCollisionO(
        gpu::cuda::getFusedCollisionO<
            T, DESCRIPTOR, BulkDynamics,
            CombinedRLBdynamics<T, DESCRIPTOR, BulkDynamics, momenta::RegularizedVelocityBoundaryTuple<0, -1>>,
            CombinedRLBdynamics<T, DESCRIPTOR, BulkDynamics, momenta::RegularizedVelocityBoundaryTuple<0, 1>>,
            CombinedRLBdynamics<T, DESCRIPTOR, BulkDynamics, momenta::RegularizedVelocityBoundaryTuple<1, -1>>,
            CombinedRLBdynamics<T, DESCRIPTOR, BulkDynamics, momenta::RegularizedVelocityBoundaryTuple<1, 1>>,
            CombinedRLBdynamics<T, DESCRIPTOR, BulkDynamics, momenta::RegularizedVelocityBoundaryTuple<2, -1>>,
            CombinedRLBdynamics<T, DESCRIPTOR, BulkDynamics, momenta::RegularizedVelocityBoundaryTuple<2, 1>>>());
  });
#endif
#endif // PLATFORM_GPU_CUDA
}

/// Set initial condition for primal variables (velocity and density)
/// @param myCase The Case instance which keeps the simulation data
/// @note Be careful: initial values have to be set using lattice units
void setInitialValues(MyCase& myCase)
{
  OstreamManager clout(std::cout, "Initialization");
  //clout << "lattice initialization ..." << std::endl;

  using T = MyCase::value_t;

  auto& sLattice = myCase.getLattice(NavierStokes {});

  auto& sGeometry = myCase.getGeometry();

  AnalyticalConst3D<T, T> rhoF((T)1);
  AnalyticalConst3D<T, T> uF((T)0, (T)0, (T)0);

  auto bulkIndicator = sGeometry.getMaterialIndicator({1, 2, 3});
  sLattice.iniEquilibrium(bulkIndicator, rhoF, uF);
  sLattice.defineRhoU(bulkIndicator, rhoF, uF);

#ifdef LID_DRIVEN
  AnalyticalConst3D<T, T> uTop(sLattice.getUnitConverter().getCharLatticeVelocity(), (T)0, (T)0);
  sLattice.defineU(sGeometry, 3, uTop);
#endif

  const T omega = sLattice.getUnitConverter().getLatticeRelaxationFrequency();

  sLattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  sLattice.initialize();
  //clout << "Initialization ... OK" << std::endl;
  return;
}

void getResults(MyCase& myCase)
{
  using T        = MyCase::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes {});

  SuperVTMwriter3D<T> vtmWriter("cavity3d");
  const auto&         converter = sLattice.getUnitConverter();

  SuperLatticeCuboid3D       cuboidF(sLattice);
  SuperLatticeRank3D         rankF(sLattice);
  SuperLatticePhysVelocity3D velocityF(sLattice, converter);
  SuperLatticePhysPressure3D pressureF(sLattice, converter);

  vtmWriter.write(cuboidF);
  vtmWriter.write(rankF);
  vtmWriter.write(velocityF);
  vtmWriter.write(pressureF);
}

/// @brief Execute simulation: set initial values and run time loop
/// @param myCase The Case instance which keeps the simulation data
void simulate(MyCase& myCase)
{
  OstreamManager clout(std::cout, "Time marching");
  using T          = MyCase::value_t;
  auto& sLattice   = myCase.getLattice(NavierStokes {});
  auto& parameters = myCase.getParameters();
  sLattice.setStatisticsOff();
  if (!parameters.get<parameters::NO_EXPORT_RESULTS>()) {
    sLattice.writeSummary();
    getResults(myCase);
  }

  //this time marching to make MLUPs idling.
  for (std::size_t iT = 0; iT < 10; ++iT) {
    sLattice.collideAndStream();
  }

#ifdef PLATFORM_GPU_CUDA
  gpu::cuda::device::synchronize();
#endif

  util::Timer<T> timer(parameters.get<parameters::TIME_STEPS>(), myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT < parameters.get<parameters::TIME_STEPS>(); ++iT) {
    sLattice.collideAndStream();
  }

#ifdef PLATFORM_GPU_CUDA
  gpu::cuda::device::synchronize();
#endif

  timer.stop();
  timer.update(parameters.get<parameters::TIME_STEPS>());

  sLattice.setProcessingContext(ProcessingContext::Evaluation);

  if (singleton::mpi().isMainProcessor()) {
    std::cout << parameters.get<parameters::RESOLUTION>() << ", " << parameters.get<parameters::TIME_STEPS>() << ", "
              << parameters.get<parameters::CUBOIDS_PER_PROCESS>() << ", " << singleton::mpi().getSize() << ", "
#ifdef PARALLEL_MODE_OMP
              << singleton::omp().getSize() << ", "
#else
              << 1 << ", "
#endif
              << timer.getTotalMLUPs() << std::endl;
  }

  if (!parameters.get<parameters::NO_EXPORT_RESULTS>()) {
    getResults(myCase);
  }
}

int main(int argc, char** argv)
{
  initialize(&argc, &argv, false, false);
  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION>(100);
    myCaseParameters.set<TIME_STEPS>(100);
    myCaseParameters.set<PHYS_CHAR_LENGTH>(1.0);
    myCaseParameters.set<NO_EXPORT_RESULTS>(true);
    myCaseParameters.set<CUBOIDS_PER_PROCESS>(1);
    myCaseParameters.set<PHYS_DELTA_X>([&] {
      return myCaseParameters.get<PHYS_CHAR_LENGTH>() / myCaseParameters.get<RESOLUTION>();
    });
  }
  myCaseParameters.fromCLI(argc, argv);

  if (!myCaseParameters.get<parameters::NO_EXPORT_RESULTS>()) {
    singleton::directories().setOutputDir("./tmp/");
  }

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

  return 0;
}
