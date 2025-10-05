/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006-2025 Mathias J. Krause, Fabian Klemens,
 *  Julius Jessberger, Shota Ito
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
#include "analyticalSolutionTestFlow3D.h"
#include "../helper.h"   // Will be removed once SuperLatticeFieldReductionO enables indicator support

using namespace olb;
using namespace olb::names;
using namespace olb::opti;

/// @brief Step 1: Declare simulation structure.
/// Model name and lattice type are collected in a Case class
using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = descriptors::D3Q19<descriptors::FORCE>;

/// State the type(s) of the used simulation cases and give a name to each case
template<typename T>
using MyCase = Case<NavierStokes, Lattice<T, DESCRIPTOR>>;
using MyOptiCase = OptiCaseADfFromDim<3,
  Controlled, MyCase<T>,
  Reference, MyCase<T>
>;

using OBJECTIVE = functors::L2DistanceF<functors::VelocityF>;
using CONTROLS = descriptors::FORCE;

// Physical parameters
const T physLength        = 1;      ///< length of the squared cube [m]
const T physVelocity      = 1;      ///< characteristic velocity [m/s]
const T physViscosity     = 0.1;    ///< kinetic viscosity of fluid [m*m/s]
const T physDensity       = 1;      ///< fluid density [kg/(m*m*m)]
const T physMaxT          = 6.0;    ///< maximal simulation time [s]
const T physStartT        = 4.0;    ///< duration of smooth start up of simulation [s]

// Discretization parameters
const int resolution      = 7;     ///< number of cells at the side of the cube
const T physDeltaX = physLength / resolution;      ///< spatial mesh width [m]
const T latticeVelocity   = 0.07;   ///< quotient of grid spacing to temporal spacing


/// @brief Step 3: Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
template<typename T>
Mesh<T,MyCase<T>::d> createMesh()
{
  /// @par Contents
  /// @li Define the spatial simulation domain
  const Vector<T,3> origin{0, 0, 0};
  const Vector<int,3> extend{resolution+1, resolution+1, resolution+1};

  /// @li Create the mesh decomposed into `singleton::mpi().getSize()` cuboids
  return Mesh<T,MyCase<T>::d>(origin, physDeltaX, extend, singleton::mpi().getSize());
}

/// @brief Step 5: Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
/// @note The material numbers are used to assign physics to lattice nodes
template<typename T>
void prepareGeometry(MyCase<T>& myCase) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  /// @par Contents
  /// @li Store a reference for simplified access
  auto& geometry = myCase.getGeometry();

  /// @li Set material numbers
  geometry.rename(0, 2);
  geometry.rename(2, 1, {1,1,1});

  /// @li Print some information on geometry
  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

/// @brief Step 6: Set lattice dynamics and boundary conditions
/// @param myCase The Case instance which keeps the simulation data
template<typename T>
void prepareLattice(MyCase<T>& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  /// @par Contents
  /// @li Store references for simplified access
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});

  /// @li Set up a unit converter with the characteristic physical units
  lattice.template setUnitConverter<UnitConverterFromResolutionAndLatticeVelocity<T,DESCRIPTOR>>(
    resolution,        // resolution
    latticeVelocity,   // charLatticeVelocity
    physLength,        // charPhysLength: reference length of simulation geometry in [m]
    physVelocity,      // charPhysVelocity: highest expected velocity during simulation in [m/s]
    physViscosity,     // physViscosity: physical kinematic viscosity in [m^2/s]
    physDensity        // physDensity: physical density [kg/m^3]
  );
  auto& converter = lattice.getUnitConverter();
  converter.print();

  /// @li Material=1 --> bulk dynamics
  lattice.template defineDynamics<ForcedBGKdynamics>(geometry.getMaterialIndicator({1}));
  /// @li Material=2,3 --> velocity boundary
  boundary::set<boundary::LocalVelocity>(lattice, geometry.getMaterialIndicator({2}));
  /// @li Set lattice relaxation frequency
  lattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
}

// Step 7: Set initial values for primal variables (e.g. velocity, density) and additional fields
/// @param myCase The Case instance which keeps the simulation data
/// @note Initial values have to be set using lattice units
template<typename T>
void setInitialValues(MyCase<T>& myCase)
{
  /// @par Contents
  /// @li Store references for simplified access
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& converter = lattice.getUnitConverter();

  /// @li Initialize density to one everywhere
  /// @li Initialize velocity to be tangential at the lid and zero in the bulk and at the walls
  AnalyticalConst3D<T,T> rhoF(1);
  const Vector<T,3> u{0,0,0};
  AnalyticalConst3D<T,T> uF(u);

  /// @li Initialize populations to equilibrium state
  lattice.defineRhoU(geometry.getMaterialIndicator({1,2}), rhoF, uF);
  lattice.iniEquilibrium(geometry.getMaterialIndicator({1,2}), rhoF, uF);

  /// @li Set force field
  ForceTestFlow3D<T,T,DESCRIPTOR> forceF(converter);
  const T latticeScaling(converter.getConversionFactorMass() / converter.getConversionFactorForce());
  AnalyticalScaled3D<T,T> scaledForceF(forceF, latticeScaling);  // conversion to lattice units
  lattice.template defineField<descriptors::FORCE>(geometry.getMaterialIndicator({1}), scaledForceF);

  /// @li Initialize lattice
  lattice.initialize();
}

/// Step 8.1: Update boundary values at times (and additional fields, if needed)
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
/// In order to start the simulation slowly, we start with zero boundary velocity and gradually increase it
/// @note Boundary values have to be set using lattice units
template<typename T>
void setTemporalValues(MyCase<T>& myCase,
                       std::size_t iT)
{
  /// @par Contents
  /// @li Store references for simplified access
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& converter = lattice.getUnitConverter();

  const std::size_t itStart = converter.getLatticeTime(physStartT);
  if (iT <= itStart) {
    /// @li Compute scaling factor
    PolynomialStartScale<T,T> startScaleF(itStart, T(1));
    const T iTvec[1] = {T(iT)};
    T frac[1] = {};
    startScaleF(frac, iTvec);

    /// @li Take analytical velocity solution, scale it to lattice units, set the boundary data
    VelocityTestFlow3D<T,T,DESCRIPTOR> velocityF(converter);
    AnalyticalScaled3D<T,T> uBoundaryStartF(velocityF, frac[0] / converter.getConversionFactorVelocity());
    lattice.defineU(geometry, 2, uBoundaryStartF);

    /// @li Communicate the new boundary velocity to GPU (if needed)
    lattice.template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }
}

/// Step 8.3: Compute simulation results at times
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
template<typename T>
void getResults(MyCase<T>& myCase,
                util::Timer<T>& timer,
                std::size_t iT)
{
  /// @par Contents
  /// @li Store references for simplified access
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();

  /// @li Write vtk plots from time to time
  const std::size_t iTvtk = converter.getLatticeTime(physMaxT/5);
  const std::size_t iTlog = converter.getLatticeTime(physMaxT/5);

  SuperVTMwriter3D<T> vtmWriter("testFlow3d");
  SuperLatticePhysVelocity3D velocityF(lattice, converter);
  SuperLatticePhysPressure3D pressureF(lattice, converter);
  vtmWriter.addFunctor(velocityF);
  vtmWriter.addFunctor(pressureF);

  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  // Writes the VTK files
  if (iT%iTvtk == 0 && iT > 0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }

  /// @li Print some (numerical and computational) statistics
  if (iT%iTlog == 0) {
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
    timer.print(iT);
  }
}

/// Step 8: Execute simulation
/// @param myCase The Case instance which keeps the simulation data
/// Run time loop
template<typename T>
void simulate(MyCase<T>& myCase)
{
  /// @par Contents
  /// @li Setup timer
  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(physMaxT);
  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT)
  {
    /// @li Step 8.1: Update the Boundary Values and Fields at Times
    setTemporalValues<T>(myCase, iT);

    /// @li Step 8.2: Collide and Stream Execution
    myCase.getLattice(NavierStokes{}).collideAndStream();

    /// @li Stripe off density offset due to Dirichlet boundary conditions
    myCase.getLattice(NavierStokes{}).stripeOffDensityOffset(
      myCase.getLattice(NavierStokes{}).getStatistics().getAverageRho() - T{1});

    /// @li Step 8.3: Computation and Output of the Results
    getResults<T>(myCase, timer, iT);
  }

  /// @li Evaluate timer
  timer.stop();
  timer.printSummary();
}

void setInitialControl(MyOptiCase& optiCase) {
  // Intialize the three components of control
  optiCase.getController().set({0, 0, 0});
}

template<typename T>
void applyControl(MyOptiCase& optiCase) {
  // decide whether we solve for value or derivatives
  auto& controlledCase = optiCase.getCaseByType<T>();
  auto& lattice = controlledCase.getLattice(NavierStokes{});
  auto designDomain = controlledCase.getGeometry().getMaterialIndicator({1});
  auto& converter = lattice.getUnitConverter();

  /// Scale the force field component-wise by the control
  std::vector<T> control = optiCase.getController<T>().get();
  std::shared_ptr<AnalyticalF<3,T,T>> controlF
   = std::make_shared<AnalyticalConst3D<T,T>>(control);
  ForceTestFlow3D<T,T,DESCRIPTOR> forceF(converter);
  const T latticeScaling(converter.getConversionFactorMass() / converter.getConversionFactorForce());
  std::shared_ptr<AnalyticalF<3,T,T>> latticeForceF
   = std::make_shared<AnalyticalScaled3D<T,T>>(forceF, latticeScaling);
  latticeForceF = latticeForceF * controlF;
  lattice.template defineField<descriptors::FORCE>(designDomain, *latticeForceF);
  lattice.setProcessingContext(ProcessingContext::Simulation);
}

template<typename T>
T objectiveF(MyOptiCase& optiCase) {
  // decide whether we solve for value or derivatives
  auto& controlledCase = optiCase.getCaseByType<T>();
  auto& controlledLattice = controlledCase.getLattice(NavierStokes{});
  auto& converter = controlledLattice.getUnitConverter();
  auto objectiveDomain = controlledCase.getGeometry().getMaterialIndicator({1});

  auto& referenceLattice = optiCase.getCase(Reference{}).getLattice(NavierStokes{});
  auto& refConverter = referenceLattice.getUnitConverter();
  auto refObjectiveDomain = optiCase.getCase(Reference{}).getGeometry().getMaterialIndicator({1});

  // Evaluate functor for objective computation
  auto objectiveO = makeWriteFunctorO<OBJECTIVE,opti::J>(controlledLattice);
  objectiveO->restrictTo(objectiveDomain);

  // Write velocity field in the reference case
  writePhysFunctorTo<functors::VelocityF,OBJECTIVE::Reference>(referenceLattice,
                                                               refObjectiveDomain,
                                                               refConverter.getConversionFactorVelocity());

  // Get solution from the reference simulation for the inverse problem
  copyFields<OBJECTIVE::Reference,OBJECTIVE::Reference>(referenceLattice, controlledLattice);
  objectiveO->template setParameter<descriptors::CONVERSION>(converter.getConversionFactorVelocity());
  const T normalize = norm(referenceLattice, refConverter, refObjectiveDomain);
  objectiveO->template setParameter<descriptors::NORMALIZE>(normalize);
  objectiveO->apply();

  controlledLattice.setProcessingContext(ProcessingContext::Evaluation);
  return integrate<opti::J>(controlledLattice, objectiveDomain)[0];
}

template<typename T>
T computeObjective(MyOptiCase& optiCase) {
  // decide whether we solve for value or derivatives
  auto& controlledCase = optiCase.getCaseByType<T>();
  // Reset prior simulation lattice
  controlledCase.resetLattices();

  // Prepare new lattices
  prepareLattice<T>(controlledCase);
  setInitialValues<T>(controlledCase);

  // Set updated controls in the simulation
  applyControl<T>(optiCase);

  // Execute simulation
  simulate<T>(controlledCase);

  // Compute Objective value from simulation results
  return objectiveF<T>(optiCase);
}

/// Steps 2-8: Set up and run a simulation
int main(int argc, char* argv[])
{
  /// @par Contents
  /// @li Step 2: Initialization
  initialize(&argc, &argv);
  using U = util::ADf<T,3>;
  MyCase<T>::ParametersD myCaseParametersD;
  MyCase<U>::ParametersD myAdfCaseParametersD;

  /// @li Step 3: Create Mesh
  Mesh mesh = createMesh<T>();
  Mesh adfMesh = createMesh<U>();

  // ==== BELOW HERE OPTIMIZATION SPECIFIC ====
  /// @li Step B: Create Cases and prepare
  // Prepare controlled case
  MyCase<T> myCase(myCaseParametersD, mesh);
  prepareGeometry<T>(myCase);

  // Create ADf-typed case for gradient computation
  MyCase<U> adfCase(myAdfCaseParametersD, adfMesh);
  prepareGeometry<U>(adfCase);

  // Compute solution for the objective functional
  MyCase<T> referenceCase(myCaseParametersD, mesh);
  prepareGeometry<T>(referenceCase);
  prepareLattice<T>(referenceCase);
  setInitialValues<T>(referenceCase);
  simulate<T>(referenceCase);

  /// @li Step C: Create OptiCase and set cases
  MyOptiCase optiCase;
  optiCase.setCase<Controlled>(myCase);
  optiCase.setCase<Derivatives>(adfCase);
  optiCase.setCase<Reference>(referenceCase);

  /// @li Step D: Set initial control
  setInitialControl(optiCase);

  /// @li Step E: Define objective routine
  optiCase.setObjective(computeObjective<T>, computeObjective<U>);

  /// @li Step G: Create an Optimizer
  OptimizerLBFGS<T,std::vector<T>> optimizer(
    optiCase.getController().size(), 1.e-10, 10, 1., 20, "Wolfe", 20, 1.e-4);

  /// @li Step H: Optimize
  optimizer.optimize(optiCase);
}
