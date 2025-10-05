/*  Lattice Boltzmann sample, written in C++, using the OpenLBlibrary
 *
 *  Copyright (C) 2025 Mathias J. Krause, David Heidenthal, Adrian Kummerländer, Michael Grinschewski, Fedor Bukreev
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

/* airfoil2d.cpp:
 * This example examines a steady flow past a NACA airfoil placed in a channel.
 * At the inlet, a random Turbulent profile is imposed on the velocity,
 * whereas the outlet implements a Dirichlet pressure condition set by p = 0.
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;

using T          = FLOATING_POINT_TYPE;
using DESCRIPTOR = D2Q9<POROSITY,VELOCITY>;

// Parameters for the simulation setup
const int N              = 250;             // resolution of the model
const T   Re             = 1400000;         // Reynolds number
const T   physViscosity  = 1.e-5;           // physical kinematic viscosity
const T   intensity      = 0.02;            // turbulence intensity for inlet
const T   maxPhysT       = 16.;             // max. simulation time in s, SI unit
const T   chordLength    = 1.;              // Length of airfoil from leading edge to trailing edge
const T   L              = chordLength / N; // latticeL
const T   lengthX        = 6.;
const T   lengthY        = 2. + L;
const T   centerAirfoilX = lengthX / 4.;
const T   centerAirfoilY = lengthY / 2.;
const T   angleOfAttack  = 5.*std::numbers::pi/90.;

// The four digits in the index are the digits of the following numbers so NACA 1410 would be
// camber = 0.01, camberPos = 0.4, thicknessPercentage = 0.10
const T camber              = 0.01;
const T camberPos           = 0.4;
const T thicknessPercentage = 0.1; // max thickness to chord length ratio

/// Wallfunction parameters
//  Used method for density reconstruction
//  0: use local density
//  1: extrapolation (Guo)
//  2: constant (rho = 1.)
const int rhoMethod = 0;
//  Used wall profile
//  0: power law profile
//  1: Spalding profile
const int wallFunctionProfile = 1;
// check if descriptor with body force is used
const bool bodyForce = false;
// interpolate sampling velocity along given normal between lattice voxels
const bool interpolateSampleVelocity = true;
// use van Driest damping function for turbulent viscosity in boundary cell
const bool useVanDriest = false;
//  distance from cell to real wall in lattice units if no geometry indicator is given as input
const T latticeWallDistance = 0.5;
//  distance from cell to velocity sampling point in lattice units
const T samplingCellDistance = 3.5;
const bool movingWall = false;
const bool averageVelocity = false;

using FringeDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::BulkTuple,
  equilibria::ThirdOrder,
  collision::ParameterFromCell<collision::LES::SMAGORINSKY,
                               collision::SmagorinskyEffectiveOmega<collision::RLBThirdOrder>>
>;

struct SmoothInflowUpdateO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct VELOCITY : public descriptors::FIELD_BASE<0,1> { };

  using parameters = meta::list<VELOCITY>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    auto& cell = cells.template get<names::NavierStokes>();
    auto u = parameters.template get<VELOCITY>();
    cell.defineU(u.data());
  }
};

// Stores geometry information in form of material numbers
void prepareGeometry(SuperGeometry<T, 2>&                superGeometry,
                     IndicatorF2D<T>&                    airfoil)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  Vector<T, 2> extend(lengthX, lengthY);
  Vector<T, 2> origin;

  superGeometry.rename(0, 2);

  superGeometry.rename(2, 1, {1, 1});

  // Set material number for inflow
  extend[0] = 2. * L;
  origin[0] = -L;
  IndicatorCuboid2D<T> inflow(extend, origin);
  superGeometry.rename(2, 3, 1, inflow);

  {
    IndicatorCuboid2D<T> cornerI(L, Vector{-L/2, lengthY});
    superGeometry.rename(2, 6, cornerI);
  }
  {
    IndicatorCuboid2D<T> cornerI(L, Vector{-L/2, -L/2});
    superGeometry.rename(2, 6, cornerI);
  }

  {
    IndicatorCuboid2D<T> fringeI(Vector{20*L, lengthY}, Vector{lengthX-20*L, 0});
    superGeometry.rename(1, 7, fringeI);
  }

  // Set material number for outflow
  origin[0] = lengthX - L;
  IndicatorCuboid2D<T> outflow(extend, origin);
  superGeometry.rename(2, 4, outflow);
  // Set material number for airfoil
  superGeometry.rename(1, 5, airfoil);

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean(false, {1,7});
  superGeometry.checkForErrors();
  superGeometry.communicate();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(SuperLattice<T, DESCRIPTOR>&        sLattice,
                    SuperGeometry<T, 2>&                superGeometry,
                    IndicatorF2D<T>&                    airfoil,
                    WallModelParameters<T>&             wallModelParameters)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const auto& converter = sLattice.getUnitConverter();
  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  // auto bulkIndicator = superGeometry.getMaterialIndicator({1});
  setTurbulentWallModelDynamics(sLattice, superGeometry, 1,
                                wallModelParameters);

  // Material=2 -->Full Slip
  boundary::set<boundary::FullSlip>(sLattice, superGeometry, 2); // <-

  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 6);

  // Setting of the boundary conditions

  //if boundary conditions are chosen to be local

  //if boundary conditions are chosen to be interpolated
  boundary::set<boundary::LocalVelocity>(sLattice, superGeometry, 3);
  boundary::set<T,DESCRIPTOR,boundary::InterpolatedPressure<T,DESCRIPTOR,FringeDynamics>>(sLattice,
                                                                                          superGeometry.getMaterialIndicator(4),
                                                                                          superGeometry.getMaterialIndicator(7),
                                                                                          superGeometry.getMaterialIndicator(0));

  setBouzidiBoundary(sLattice, superGeometry, 5, airfoil);
  setTurbulentWallModel(sLattice, superGeometry, 5, wallModelParameters);

  // Initial conditions
  AnalyticalConst2D<T, T> rhoF(1);
  AnalyticalConst2D<T, T> rho0(0);
  std::vector<T>          velocity(2, T(0));
  AnalyticalConst2D<T, T> uF(velocity);

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU(superGeometry.getMaterialIndicator({1,2,3,4,5,6,7}), rhoF, uF);
  sLattice.iniEquilibrium(superGeometry.getMaterialIndicator({1,2,3,4,5,6,7}), rhoF, uF);
  sLattice.defineField<VELOCITY>(superGeometry.getMaterialIndicator({1,2,3,4,5,6,7}), uF);
  sLattice.defineField<POROSITY>(superGeometry.getMaterialIndicator({1,2,3,4,6,7}), rhoF);
  sLattice.defineField<POROSITY>(superGeometry, 5, rho0);

  sLattice.defineDynamics<FringeDynamics>(superGeometry, 7);
  {
    Vector<T,2> origin(lengthX-20*converter.getPhysDeltaX(),
                       -converter.getPhysDeltaX());
    Vector<T,2> extend(20*converter.getPhysDeltaX(),
                       lengthY + 2*converter.getPhysDeltaX());
    IndicatorCuboid2D<T> out(extend, origin);
    SuperIndicatorFfromIndicatorF2D<T> outI(out, superGeometry);

    AnalyticalFfromCallableF<2,T,T> smagorinskyF([&](Vector<T,2> physR) -> Vector<T,1> {
      return 0.3 + ((physR[0] - origin[0]) / (20*converter.getPhysDeltaX())) * (2-0.3);
    });

    sLattice.defineField<collision::LES::SMAGORINSKY>(outI, smagorinskyF);
  }

  sLattice.setParameter<descriptors::OMEGA>(omega);
  sLattice.setParameter<collision::LES::SMAGORINSKY>(0.3);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Computes the pressure drop between the voxels before and after the cylinder
void getResults(SuperLattice<T, DESCRIPTOR>& sLattice,
                std::size_t iT, SuperGeometry<T, 2>& sGeometry,
                util::Timer<T>& timer, IndicatorF2D<T>& airfoil)
{
  OstreamManager clout(std::cout, "getResults");

  const int vtkIter  = sLattice.getUnitConverter().getLatticeTime(.3);
  const int statIter = sLattice.getUnitConverter().getLatticeTime(.1);

  if (iT == 0) {
    SuperVTMwriter2D<T> vtmWriter("airfoil2d");
    vtmWriter.createMasterFile();
  }

  if (iT % statIter == 0) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    T point[2] = {};
    point[0]   = centerAirfoilX + chordLength;
    point[1]   = centerAirfoilY;
    SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure(sLattice, sLattice.getUnitConverter());
    AnalyticalFfromSuperF2D<T> intpolateP(pressure, true);
    T                          p;
    intpolateP(&p, point);

    // Timer console output
    timer.update(iT);
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print(iT, sLattice.getUnitConverter().getPhysTime(iT));

    // Drag, lift, pressure drop
    AnalyticalFfromSuperF2D<T>            intpolatePressure(pressure, true);
    SuperLatticePhysDrag2D<T, DESCRIPTOR> drag(sLattice, sGeometry, 5,
                                               sLattice.getUnitConverter());

    T point1[2] = {};
    T point2[2] = {};

    point1[0] = centerAirfoilX - chordLength / 2;
    point1[1] = centerAirfoilY;

    point2[0] = centerAirfoilX + chordLength / 2;
    point2[1] = centerAirfoilY;

    T p1, p2;
    intpolatePressure(&p1, point1);
    intpolatePressure(&p2, point2);

    clout << "pressure1=" << p1;
    clout << "; pressure2=" << p2;

    T pressureDrop = p1 - p2;
    clout << "; pressureDrop=" << pressureDrop;

    int input[3] = {};
    T   _drag[drag.getTargetDim()];
    drag(_drag, input);
    clout << "; drag=" << _drag[0] << "; lift=" << _drag[1] << std::endl;
  }

  if (iT % vtkIter == 0) {
    SuperVTMwriter2D<T>                       vtmWriter("airfoil2d");
    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(sLattice, sLattice.getUnitConverter());
    SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure(sLattice, sLattice.getUnitConverter());
    SuperLatticeField2D<T, DESCRIPTOR, descriptors::Y1> y1(sLattice);
    SuperLatticePhysField2D<T, DESCRIPTOR, descriptors::WMVELOCITY> wmvelocity(sLattice, sLattice.getUnitConverter().getConversionFactorVelocity());
    wmvelocity.getName() = "wmvelocity";
    SuperLatticePhysWallShearStress2D<T, DESCRIPTOR> wss(sLattice, sGeometry, 5, sLattice.getUnitConverter(), airfoil);

    vtmWriter.addFunctor(velocity);
    vtmWriter.addFunctor(pressure);
    vtmWriter.addFunctor(y1);
    vtmWriter.addFunctor(wmvelocity);
    vtmWriter.addFunctor(wss);

    vtmWriter.write(iT);
  }
}

int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");

  const UnitConverterFromResolutionAndLatticeVelocity<T, DESCRIPTOR> converter(
      int {N},                           // resolution: number of voxels per charPhysL
      (T) 0.03,                           // Max cell speed?
      (T) chordLength,                    // charPhysLength: reference length of simulation geometry
      (T) Re * physViscosity/chordLength, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
      (T) physViscosity,                  // physViscosity: physical kinematic viscosity in __m^2 / s__
      (T) 1.0                             // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("airfoil2d");

  // === 2rd Step: Prepare Geometry ===
  Vector<T, 2>         extend(lengthX, lengthY);
  Vector<T, 2>         origin;
  IndicatorCuboid2D<T> cuboid(extend, origin);

  #ifdef PARALLEL_MODE_MPI
    const int noOfCuboids = singleton::mpi().getSize();
  #else
    const int noOfCuboids = 1;
  #endif

  // Instantiation of a cuboidGeometry with weights
  CuboidDecomposition2D<T> cuboidGeometry(cuboid, L, noOfCuboids);

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  // Instantiation of a superGeometry
  SuperGeometry<T, 2> superGeometry(cuboidGeometry, loadBalancer);

  Vector<T, 2>          center(centerAirfoilX, centerAirfoilY);
  IndicatorAirfoil2D<T> airfoilI = IndicatorAirfoil2D<T>(
      center, chordLength, camber, camberPos, thicknessPercentage);
  IndicatorRotate<T,2> airfoil(center, angleOfAttack, airfoilI);
  prepareGeometry(superGeometry, airfoil);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T,DESCRIPTOR> sLattice(converter, superGeometry);
  WallModelParameters<T>      wallModelParameters;
  wallModelParameters.bodyForce                 = bodyForce;
  wallModelParameters.rhoMethod                 = rhoMethod;
  wallModelParameters.samplingCellDistance      = samplingCellDistance;
  wallModelParameters.interpolateSampleVelocity = interpolateSampleVelocity;
  wallModelParameters.useVanDriest              = useVanDriest;
  wallModelParameters.wallFunctionProfile       = wallFunctionProfile;
  wallModelParameters.latticeWallDistance       = latticeWallDistance;
  wallModelParameters.movingWall                = movingWall;
  wallModelParameters.averageVelocity           = averageVelocity;

  prepareLattice(sLattice, superGeometry, airfoil,
                 wallModelParameters);

  SuperLatticeCoupling smoothInflowUpdateO(
    SmoothInflowUpdateO{},
    names::NavierStokes{}, sLattice);
  smoothInflowUpdateO.restrictTo(superGeometry.getMaterialIndicator(3));

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer(converter.getLatticeTime(maxPhysT),
                       superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {
    {
      const std::size_t iTmaxStart = converter.getLatticeTime(0.1 * maxPhysT);
      const std::size_t iTupdate   = converter.getLatticeTime(0.001);
      if (iT % iTupdate == 0 && iT <= iTmaxStart) {
        T frac[1] = {};
        PolynomialStartScale<T,T> scale(iTmaxStart, T(1));
        T iTvec[1] = {T( iT )};
        scale(frac, iTvec);
        smoothInflowUpdateO.setParameter<SmoothInflowUpdateO::VELOCITY>({frac[0] * converter.getCharLatticeVelocity(), 0});
        smoothInflowUpdateO.apply();
      }
    }

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults(sLattice, iT, superGeometry, timer, airfoil);
  }

  timer.stop();
  timer.printSummary();
}
