/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Jan Marquardt
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

/*
 * TODO:
 * - Test setup with different configurations
 * - Update formatting of output as needed
 * - Use enumeration instead of preprocessor directives
 * - Improve description below
 */

/*
 * This example is based on 10.1016/j.cpc.2024.109321.
 * The limestone shapes are from the online particle database PARROT
 * (https://parrot.tu-freiberg.de/).
 * The number of triangles has been reduced.
 * For the limestones, please reduce relaxation time tau to at least 0.55
 */

#include "olb.h"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <memory>
#include <vector>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::particles;
using namespace olb::particles::dynamics;
using namespace olb::particles::contact;
using namespace olb::particles::communication;
using namespace olb::particles::access;
using namespace olb::util;

#define WriteVTK
#define WithContact

typedef double T;

typedef enum {
    LIMESTONES,
    CUBES,
    SPHERES
} ParticleType;

const ParticleType particleType = SPHERES;
T  wantedParticleVolumeFraction = 0.15;

// Define lattice type
typedef D3Q19<POROSITY, VELOCITY_NUMERATOR, VELOCITY_DENOMINATOR,
#ifdef WithContact
              CONTACT_DETECTION,
#endif
              FORCE>
    DESCRIPTOR;

// Define particle type
typedef PARTICLE_DESCRIPTOR<
    DESCRIPTOR::d, GENERAL_EXTENDABLE<DESCRIPTOR::d>,
    MOBILITY_VERLET<DESCRIPTOR::d>, SURFACE_RESOLVED_PARALLEL<DESCRIPTOR::d>,
    FORCING_RESOLVED<DESCRIPTOR::d>, PHYSPROPERTIES_RESOLVED<DESCRIPTOR::d>,
#ifdef WithContact
    MECHPROPERTIES_COLLISION<DESCRIPTOR::d>,
    NUMERICPROPERTIES_RESOLVED_CONTACT<DESCRIPTOR::d>,
#endif
    PARALLELIZATION_RESOLVED>
    PARTICLETYPE;

// Define particle-particle contact type
typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
    PARTICLECONTACTTYPE;

// Define particle-wall contact type
typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
    WALLCONTACTTYPE;

constexpr T nonDimensionalTime = 200.0;
constexpr T gravit             = 9.81;

T maxPhysT;

// Discretization Settings
int         res = 9;
constexpr T tau = 0.6;


// Time Settings
int iTpurge;
int iTwrite;

// Fluid Settings
constexpr T fluidDensity = 1000;
T           dynamicViscosity;
T           kinematicViscosity;

// Particle Settings
unsigned int particleNumber = 0;
T            ArchimedesNumber(1000.);
T            densityRatio(3.3);
constexpr T  radius             = 1.5e-3;
constexpr T  diameter           = 2 * radius;
T            equivalentDiameter = diameter;
T            particleDensity;
T            particleVolumeFraction;

// Domain Settings
constexpr T        extentToDiameter = 15.0;
Vector<T, 3>       extent(extentToDiameter* diameter);
const Vector<T, 3> origin(T {0});

// Characteristic Quantities
constexpr T charPhysLength = diameter;

// Contact Settings
constexpr T coefficientOfRestitution   = 0.926;
constexpr T coefficientKineticFriction = 0.16;
constexpr T coefficientStaticFriction  = T {2} * coefficientKineticFriction;
constexpr T staticKineticTransitionVelocity         = 0.001;
constexpr T YoungsModulusParticle                   = 5.0e3;
constexpr T PoissonRationParticle                   = 0.245;
T           particleEnlargementForContact           = 0;
constexpr unsigned contactBoxResolutionPerDirection = 8;


// STLs with scaling dimensions to almost reach equivalent volume
std::vector<std::pair<std::string, T>> limestoneStlFiles = {
    std::make_pair("./limestone/312_reduced.stl", 2.518e-4),
    std::make_pair("./limestone/529_reduced.stl", 2.625e-4),
    std::make_pair("./limestone/1076_reduced.stl", 1.012e-4),
    std::make_pair("./limestone/1270_reduced.stl", 1.214e-4),
    std::make_pair("./limestone/1810_reduced.stl", 0.862e-4)
  };


#ifdef PARALLEL_MODE_MPI
MPI_Comm
    averageParticleVelocityComm; /// Communicator for calculation of average particle velocity
MPI_Comm
    numberParticleCellsComm; /// Communicator for calculation of current number of particle cells
#endif // PARALLEL_MODE_MPI

std::string getParticleIdentifier(const std::size_t& pID)
{
  return std::to_string(pID);
}

T evalTerminalVelocitySingleParticleStokes()
{
  return (2. / 9.) * (particleDensity - fluidDensity) * gravit * radius *
         radius / dynamicViscosity;
};


T evalAverageSettlingVelocity(XParticleSystem<T, PARTICLETYPE>& xParticleSystem)
{
  T averageSettlingVelocity {0};

  communication::forParticlesInSuperParticleSystem<
      T, PARTICLETYPE, conditions::valid_particle_centres>(
      xParticleSystem,
      [&](Particle<T, PARTICLETYPE>&       particle,
          ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
        averageSettlingVelocity += getVelocity(particle)[2];
      });

#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(averageSettlingVelocity, MPI_SUM,
                                  singleton::mpi().bossId(),
                                  averageParticleVelocityComm);
#endif // PARALLEL_MODE_MPI

  return averageSettlingVelocity / particleNumber;
}

void updateBodyForce(SuperLattice<T, DESCRIPTOR>&        sLattice,
                     XParticleSystem<T, PARTICLETYPE>&   xParticleSystem,
                     SuperGeometry<T, 3>&                superGeometry,
                     UnitConverter<T, DESCRIPTOR> const& converter)
{
  SuperLatticePhysExternalPorosity3D<T, DESCRIPTOR> porosity(sLattice,
                                                             converter);
  SuperAverage3D<T> avPorosityF(porosity, superGeometry, 1);
  int               input[1]{};
  T                 fluidVolumeFraction[1];
  avPorosityF(fluidVolumeFraction, input);
  const T volumeRatio = particleVolumeFraction / fluidVolumeFraction[0];

  // Apply equal to submerged weight of the particles to the fluid
  std::vector<T> balancingAcceleration(3, T(0.));
  balancingAcceleration[2] =
      gravit * volumeRatio * (particleDensity / fluidDensity - 1);

  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
  SuperAverage3D<T>                         avgVel(velocity, superGeometry, 1);
  T                                         vel[3]{};
  avgVel(vel, input);
  // Avoid numerical drift
  balancingAcceleration[2] -=
      vel[2] * converter.getCharPhysVelocity() / converter.getCharPhysLength();

  const T conversionFactor = converter.getConversionFactorTime() *
                             converter.getConversionFactorTime() /
                             converter.getConversionFactorLength();
  balancingAcceleration[2] *= conversionFactor;

  AnalyticalConst3D<T, T> acc(balancingAcceleration);
  sLattice.defineField<descriptors::FORCE>(superGeometry, 1, acc);
}

template <typename T>
T calculateCubeEdgeLengthFromSphereRadius()
{
  // Calculate the volume of the sphere
  const T sphereVolume = (4.0 / 3.0) * M_PI * util::pow(radius, 3);

  // Calculate the edge length of the cube
  const T cubeEdgeLength = util::pow(sphereVolume, 1.0 / 3.0);

  return cubeEdgeLength;
}

// Prepare geometry
void prepareGeometry(SuperGeometry<T, 3>&                superGeometry,
                     UnitConverter<T, DESCRIPTOR> const& converter)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0, 1);

  superGeometry.clean();
  superGeometry.innerClean();

  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(SuperLattice<T, DESCRIPTOR>&        sLattice,
                    UnitConverter<T, DESCRIPTOR> const& converter,
                    SuperGeometry<T, 3>&                superGeometry)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  /// Material=0 -->do nothing
  sLattice.defineDynamics<
      PorousParticleKupershtokhForcedBGKdynamics<T, DESCRIPTOR>>(superGeometry,
                                                                 1);
  sLattice.defineDynamics<BounceBack>(superGeometry, 2);

  sLattice.setParameter<descriptors::OMEGA>(
      converter.getLatticeRelaxationFrequency());

  {
    auto& communicator = sLattice.getCommunicator(stage::PostPostProcess());
    communicator
        .requestFields<POROSITY, VELOCITY_NUMERATOR, VELOCITY_DENOMINATOR>();
    communicator.requestOverlap(sLattice.getOverlap());
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Set Boundary Values
void setBoundaryValues(SuperLattice<T, DESCRIPTOR>&        sLattice,
                       UnitConverter<T, DESCRIPTOR> const& converter, int iT,
                       SuperGeometry<T, 3>& superGeometry)
{
  OstreamManager clout(std::cout, "setBoundaryValues");

  if (iT == 0) {
    AnalyticalConst3D<T, T> zero(0.);
    AnalyticalConst3D<T, T> one(1.);
    sLattice.defineField<descriptors::POROSITY>(
        superGeometry.getMaterialIndicator(1), one);
    // Set initial condition
    AnalyticalConst3D<T, T>    ux(0.);
    AnalyticalConst3D<T, T>    uy(0.);
    AnalyticalConst3D<T, T>    uz(0.);
    AnalyticalComposed3D<T, T> u(ux, uy, uz);

    AnalyticalConst3D<T, T> rho(1.);

    // Initialize all values of distribution functions to their local equilibrium
    sLattice.defineRhoU(superGeometry, 1, rho, u);
    sLattice.iniEquilibrium(superGeometry, 1, rho, u);

    // Make the lattice ready for simulation
    sLattice.initialize();
  }
}

void getResults(SuperLattice<T, DESCRIPTOR>&        sLattice,
                UnitConverter<T, DESCRIPTOR> const& converter, int iT,
                SuperGeometry<T, 3>& superGeometry, Timer<double>& timer,
                XParticleSystem<T, PARTICLETYPE>& xParticleSystem)
{
  OstreamManager clout(std::cout, "getResults");

  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
#ifdef WriteVTK
  SuperVTMwriter3D<T>                       vtkWriter("sedimentation");
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
  SuperLatticePhysExternalPorosity3D<T, DESCRIPTOR> externalPor(sLattice,
                                                                converter);

  vtkWriter.addFunctor(pressure);
  vtkWriter.addFunctor(velocity);
  vtkWriter.addFunctor(externalPor);

  if (iT == 0) {
    /// Writes the converter log file
    SuperLatticeCuboid3D<T, DESCRIPTOR>   cuboid(sLattice);
    SuperLatticeRank3D<T, DESCRIPTOR>     rank(sLattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.createMasterFile();
  }

  if (iT % iTwrite == 0) {
    vtkWriter.write(iT);
  }
#endif // WriteVTK

  /// Writes output on the console
  if (iT % iTwrite == 0) {
    // TODO: Update formatting as needed
    clout << "Average settling velocity: "
          << evalAverageSettlingVelocity(xParticleSystem) << " in m/s"
          << std::endl;

    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
  }
}

int main(int argc, char* argv[])
{
  /// === 1st Step: Initialization ===
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");
  #ifdef PARALLEL_MODE_MPI //Check if MPI is available, otherwise throw error
  std::vector<std::string> cmdInput;

  if (argc > 1) {
    cmdInput.assign(argv + 1, argv + argc);
    res = std::stoi(cmdInput[0]);
  }
  if (argc > 2) {
    wantedParticleVolumeFraction = std::stod(cmdInput[1]);
  }
  if (argc > 3) {
    densityRatio = std::stod(cmdInput[2]);
  }
  if (argc > 4) {
    ArchimedesNumber = std::stod(cmdInput[3]);
  }

  const std::string positionsfilename =[]{
    if constexpr (particleType==LIMESTONES){
      return std::string("limestone/particlepositions_limestone_15_0.150021");
    }
    else return std::string("particlepositions_sphere_1.5mm_15_0.30");
  }();

  const T domainVolume = extent[0] * extent[1] * extent[2];
  particleDensity      = densityRatio * fluidDensity;
  kinematicViscosity =
      util::sqrt(gravit * (densityRatio - 1) *
                 util::pow(equivalentDiameter, 3) / ArchimedesNumber);

  const Vector<T, 3> externalAcceleration = {
      .0, .0, -gravit * (1. - fluidDensity / particleDensity)};

  // Estimation maximal velocity using Stoke's law
  dynamicViscosity         = kinematicViscosity * fluidDensity;
  const T charPhysVelocity = evalTerminalVelocitySingleParticleStokes();
  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
      (int)res,              // resolution
      (T)tau,                // latticeRelaxationTime
      (T)charPhysLength,     // charPhysLength
      (T)charPhysVelocity,   // charPhysVelocity
      (T)kinematicViscosity, // physViscosity
      (T)fluidDensity        // fluidDensity
  );
  converter.write();
  converter.print();


  if (MPI_Comm_dup(MPI_COMM_WORLD, &averageParticleVelocityComm) !=
      MPI_SUCCESS) {
    throw std::runtime_error("Unable to duplicate MPI communicator");
  }
  if (MPI_Comm_dup(MPI_COMM_WORLD, &numberParticleCellsComm) != MPI_SUCCESS) {
    throw std::runtime_error("Unable to duplicate MPI communicator");
  }


  std::size_t iT = 0;
  maxPhysT       = nonDimensionalTime * equivalentDiameter / charPhysVelocity;
  particleEnlargementForContact = converter.getPhysDeltaX() / T {5};

  /// === 2rd Step: Prepare Geometry ===
  /// Instantiation of a cuboidGeometry with weights
  IndicatorCuboid3D<T> cuboid(extent, origin);
  constexpr auto       getPeriodicity = []() {
    return Vector<bool, 3>(true, true, true);
  };

  const unsigned numberProcesses = singleton::mpi().getSize();

CuboidDecomposition3D<T> cuboidGeometry(
      cuboid, converter.getConversionFactorLength(), numberProcesses);

cuboidGeometry.setPeriodicity({ true, false, false });

  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);
  SuperGeometry<T, 3>      superGeometry(cuboidGeometry, loadBalancer, 2);
  prepareGeometry(superGeometry, converter);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice(converter, superGeometry);

  prepareLattice(sLattice, converter, superGeometry);

  std::vector<SolidBoundary<T, DESCRIPTOR::d>> solidBoundaries;
  const T epsilon     = T {0.5} * converter.getConversionFactorLength();
  const T halfEpsilon = T {0.5} * epsilon;

  constexpr T overlapSecurityFactor = 1.0;
  T maxCircumRadius = T {0};

  // Create smooth indicators for limestones
  std::vector<std::shared_ptr<STLreader<T>>> limestoneSTLreaders;
  std::vector<std::unique_ptr<olb::SmoothIndicatorCustom3D<T, T, true>>>
      limestoneIndicators;
  const T latticeSpacingDiscreteParticle =
    T {0.2} * converter.getConversionFactorLength();

  if constexpr (particleType==LIMESTONES){
    clout << "Initializing limestone ..." << std::endl;

    {
      unsigned i = 0;
      for (auto& STLfile : limestoneStlFiles) {
        limestoneSTLreaders.push_back(std::make_shared<STLreader<T>>(
            STLfile.first, converter.getConversionFactorLength(),
            STLfile.second));
        limestoneIndicators.push_back(
            std::make_unique<olb::SmoothIndicatorCustom3D<T, T, true>>(
                latticeSpacingDiscreteParticle, limestoneSTLreaders[i],
                olb::Vector<T, 3>(T {}), epsilon, olb::Vector<T, 3>(T {})));
        limestoneIndicators.back()->calcMofiAndMass(particleDensity);
        ++i;
      }
    }

    clout << "Initializing limestone ... OK" << std::endl;

    // Create Particle Dynamics
    // Create ParticleSystem

    for (auto& STLsurface : limestoneIndicators) {
      maxCircumRadius = util::max(maxCircumRadius, STLsurface->getCircumRadius());
    }
  }
  if constexpr (particleType==SPHERES){
    maxCircumRadius = radius + halfEpsilon;
  }


  if constexpr (access::providesContactMaterial<PARTICLETYPE>()) {
    const T detectionDistance =
        T {0.5} * util::sqrt(PARTICLETYPE::d) * converter.getPhysDeltaX();
    maxCircumRadius = maxCircumRadius - halfEpsilon +
                      util::max(halfEpsilon, detectionDistance);
  }
  maxCircumRadius *= overlapSecurityFactor;

  // ensure parallel mode is enabled

    SuperParticleSystem<T, PARTICLETYPE> xParticleSystem(superGeometry,
                                                       maxCircumRadius);
    particles::communication::checkCuboidSizes(xParticleSystem);


  // Create ParticleContact container
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  // Container containg contact properties
  ContactProperties<T, 1> contactProperties;
  contactProperties.set(
      0, 0,
      evalEffectiveYoungModulus(YoungsModulusParticle, YoungsModulusParticle,
                                PoissonRationParticle, PoissonRationParticle),
      coefficientOfRestitution, coefficientKineticFriction,
      coefficientStaticFriction, staticKineticTransitionVelocity);

  // Create and assign resolved particle dynamics
  xParticleSystem.defineDynamics<VerletParticleDynamics<T, PARTICLETYPE>>();

  //Create particle manager handling coupling, gravity and particle dynamics
  ParticleManager<T, DESCRIPTOR, PARTICLETYPE> particleManager(
      xParticleSystem, superGeometry, sLattice, converter, externalAcceleration, getPeriodicity());

  // Create Communicators

  const communication::ParticleCommunicator& particleCommunicator =
      particleManager.getParticleCommunicator();


  const std::function<T(const std::size_t&)> getCircumRadius =
      [&](const std::size_t& pID) {
        if constexpr (particleType==LIMESTONES){
          return limestoneIndicators[pID % limestoneStlFiles.size()]
              ->getCircumRadius();
        }
        if constexpr (particleType==SPHERES){
          return radius + T {0.5} * epsilon;
        }
        if constexpr (particleType==CUBES){

          return T {0.5} *
                (calculateCubeEdgeLengthFromSphereRadius<T>() * util::sqrt(3) +
                  epsilon);
        }

      };
  const std::function<T(const std::size_t&)> getParticleVolume =
      [&](const std::size_t& pID) {
        if constexpr (particleType==LIMESTONES){
          return limestoneIndicators[pID % limestoneStlFiles.size()]->getVolume();
        }
        if constexpr (particleType==SPHERES){
          return T {4} / T {3} * M_PI * radius * radius * radius;
        }
        if constexpr (particleType==CUBES){
          return util::pow(calculateCubeEdgeLengthFromSphereRadius<T>(), 3);
        }

      };

  const std::function<void(
      const particles::creators::SpawnData<T, DESCRIPTOR::d>&,
      const std::string&)>
      createParticleFromString =
          [&](const particles::creators::SpawnData<T, DESCRIPTOR::d>& data,
              const std::string&                                      pID) {
            const PhysR<T, 3>  physPosition   = data.position;
            const Vector<T, 3> angleInDegrees = data.angleInDegree;

          if constexpr (particleType==LIMESTONES){
            if (particleNumber < limestoneStlFiles.size()) {
              std::shared_ptr <STLreader<T>> limestoneIndicator =
                  std::make_shared<STLreader<T>>(
                      limestoneStlFiles[particleNumber].first,
                      converter.getConversionFactorLength(),
                      limestoneStlFiles[particleNumber].second);
              creators::addResolvedArbitraryShape3D<T, PARTICLETYPE>(
                  xParticleSystem,  physPosition,
                  latticeSpacingDiscreteParticle,
                  limestoneIndicator,
                  epsilon, particleDensity );
            }
            else {
              creators::addResolvedObject<T, PARTICLETYPE>(
                  xParticleSystem,
                  particleNumber % limestoneStlFiles.size(), physPosition,
                  particleDensity, angleInDegrees);
            }
          }
          if constexpr (particleType==SPHERES){
            creators::addResolvedSphere3D(xParticleSystem,
                                          physPosition, radius, epsilon,
                                          particleDensity);
          }
          if constexpr (particleType==CUBES){
            creators::addResolvedCuboid3D(
                xParticleSystem, physPosition,
                Vector<T, 3>(calculateCubeEdgeLengthFromSphereRadius<T>()),
                epsilon, particleDensity);
          }
            ++particleNumber;
          };

  if (positionsfilename != "") {
    clout << "Spawning particles from " << positionsfilename << " ..."
          << std::endl;
    if (std::filesystem::exists(positionsfilename)) {
      std::vector<particles::creators::SpawnData<T, DESCRIPTOR::d>>
          tmpSpawnData;
          tmpSpawnData = particles::creators::setParticles<T, 3>(
          positionsfilename, wantedParticleVolumeFraction, cuboid, domainVolume,
          getParticleVolume, createParticleFromString);

      T tmpVolume = T {0};
      for (unsigned pID = 0; pID < tmpSpawnData.size(); ++pID) {
        tmpVolume += getParticleVolume(pID);
      }
      particleVolumeFraction = tmpVolume / domainVolume;
    }
    else {
      OLB_ASSERT(false, positionsfilename + " does not exist.");
    }
  }
  else {
    OLB_ASSERT(false, "No particle positions file given.");
  }
 if constexpr (particleType==LIMESTONES){
    limestoneIndicators.clear();
  }

  clout << "Spawning particles from " << positionsfilename << " ... OK"
        << std::endl;

  maxCircumRadius = 0.;
  forParticlesInSuperParticleSystem(
      xParticleSystem,
      [&](Particle<T, PARTICLETYPE>&       particle,
          ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
#ifdef WithContact
        particle.setField<MECHPROPERTIES, MATERIAL>(0);
        particle.setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(
            particleEnlargementForContact);
#endif // WithContact

        const T currCircumRadius = access::getRadius(particle);
        maxCircumRadius          = util::max(currCircumRadius, maxCircumRadius);
      });

  singleton::mpi().reduceAndBcast(maxCircumRadius, MPI_MAX,
                                  singleton::mpi().bossId(),
                                  particleCommunicator.contactTreatmentComm);

  xParticleSystem.updateOffsetFromCircumRadius(overlapSecurityFactor *
                                               maxCircumRadius);

  /// === 5th Step: Definition of Initial and Boundary Conditions ===
  setBoundaryValues(sLattice, converter, 0, superGeometry);

  /// === 4th Step: Main Loop with Timer ===
  clout << "MaxIT: " << converter.getLatticeTime(maxPhysT) << std::endl;
  iTwrite = util::max(0.02 * converter.getLatticeTime(maxPhysT), 1);
  iTpurge = util::max(util::ceil(0.06 * converter.getLatticeTime(maxPhysT)), 1);

  Timer<double> timer(converter.getLatticeTime(maxPhysT),
                      superGeometry.getStatistics().getNvoxel());
  timer.start();
  for (iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {
#ifndef WithContact
    // Execute particle manager
    particleManager.execute<
        couple_lattice_to_parallel_particles<T, DESCRIPTOR, PARTICLETYPE>,
        communicate_parallel_surface_force<T, PARTICLETYPE>,
        apply_gravity<T, PARTICLETYPE>,
        process_dynamics_parallel<T, PARTICLETYPE>,
        update_particle_core_distribution<T, PARTICLETYPE>>();

    particles::dynamics::coupleResolvedParticlesToLattice<
        T, DESCRIPTOR, PARTICLETYPE, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
        xParticleSystem, contactContainer, superGeometry, sLattice, converter,
        solidBoundaries, getPeriodicity);

#else  // WithContact
    // Couple lattice to particles
    couple_lattice_to_parallel_particles<T, DESCRIPTOR, PARTICLETYPE>::execute(
        xParticleSystem, superGeometry, sLattice, converter, getPeriodicity());

    communicate_parallel_surface_force<T, PARTICLETYPE>::execute(
        xParticleSystem, particleCommunicator);

    // Apply contacts
    particles::contact::processContacts<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                        WALLCONTACTTYPE,
                                        ContactProperties<T, 1>>(
        xParticleSystem, solidBoundaries, contactContainer, contactProperties,
        superGeometry,
        particleCommunicator.contactTreatmentComm,
        contactBoxResolutionPerDirection, T {4. / (3 * util::sqrt(M_PI))},
        getPeriodicity);

    // Apply gravity
    communication::forParticlesInSuperParticleSystem<
        T, PARTICLETYPE,
        conditions::valid_particle_centres //only consider center for resolved
        >(xParticleSystem,
          [&](Particle<T, PARTICLETYPE>&       particle,
              ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
            apply_gravity<T, PARTICLETYPE>::execute(xParticleSystem, particle,
                                                    externalAcceleration,
                                                    converter.getPhysDeltaT());
          });


    // Process particles (Verlet algorithm)
    communication::forParticlesInSuperParticleSystem<
        T, PARTICLETYPE,
        conditions::valid_particle_centres //only consider center for resolved
        >(xParticleSystem,
          [&](Particle<T, PARTICLETYPE>&       particle,
              ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
            particle.process(converter.getPhysDeltaT());
          });

    communicatePostContactTreatmentContacts(
        contactContainer, xParticleSystem, converter.getPhysDeltaX(),
        particleCommunicator.particleContactDetectionComm,
        particleCommunicator.wallContactDetectionComm,
        getPeriodicity());

    update_particle_core_distribution<T, PARTICLETYPE>::execute(
        xParticleSystem, converter.getPhysDeltaX(), particleCommunicator,
        getPeriodicity());

    // Couple particles to lattice
    particles::dynamics::coupleResolvedParticlesToLattice<
        T, DESCRIPTOR, PARTICLETYPE, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
        xParticleSystem, contactContainer, superGeometry, sLattice, converter,
        solidBoundaries, getPeriodicity);

    // Communicate found contacts
    communicateContactsParallel(
        contactContainer, xParticleSystem, converter.getPhysDeltaX(),
        particleCommunicator.particleContactDetectionComm,
        particleCommunicator.wallContactDetectionComm,
        getPeriodicity());

    if constexpr (isPeriodic(getPeriodicity())) {
      accountForPeriodicParticleBoundary(xParticleSystem, contactContainer,
                                         superGeometry, getPeriodicity);
    }
#endif // WithContact

    if (iT % iTpurge == 0) {
      purgeInvalidParticles<T, PARTICLETYPE>(xParticleSystem);
    }

    // Get Results
    getResults(sLattice, converter, iT, superGeometry, timer, xParticleSystem);

    updateBodyForce(sLattice, xParticleSystem, superGeometry, converter);

    // Collide and stream
    sLattice.collideAndStream();
  }

  timer.update(iT);
  timer.stop();
  timer.printSummary();
  #else
    throw std::runtime_error(
      "This example is not designed to run in serial mode.");
  #endif // PARALLEL_MODE_MPI
}
