/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2025 Christoph Gaul, Mathias J. Krause
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

/* aneurysm
 * This example represents flow in an aneurysm.
 * The geometry is read from a STL file. Stress on the aneurysm wall can be visualized with the Mesh.pvd file.
 */

#include <olb.h>

using namespace olb;

using T          = FLOATING_POINT_TYPE;
using DESCRIPTOR = descriptors::D3Q19<>;

// Parameters for the simulation setup
const int N        = 120;      // Number of cells in the length of the aneurysm
const T   Re       = 708.;     // Reynolds number
const T   maxPhysT = 16.;      // max. simulation time in s, SI unit
const T   L        = 0.05 / N; // latticeL
const T   CFL      = 0.02;     // CFL number

// Stores data from stl file in geometry in form of material numbers
void prepareGeometry(UnitConverter<T, DESCRIPTOR> const& converter, IndicatorF3D<T>& indicator, STLreader<T>& aneurysm,
                     SuperGeometry<T, 3>& superGeometry)
{
  superGeometry.rename(0, 2, indicator);
  superGeometry.rename(2, 1, aneurysm);
  superGeometry.clean();

  std::vector<Vector<T, 3>> normals(5);
  normals[0] = {-0.026646, 0.963492, -0.266408};
  normals[1] = {0.906283, -0.171557, 0.386290};
  normals[2] = {-0.159757, 0.968927, -0.18883};
  normals[3] = {0.994170, -0.101600, -0.036091};
  normals[4] = {0.02343, 0.98639, -0.1627};

  std::vector<Vector<T, 3>> centerpoints(5);
  centerpoints[0] = {0.005072, 0.034537, -0.020140};
  centerpoints[1] = {0.032770, 0.018303, 0.014303};
  centerpoints[2] = {-0.0224, 0.0185, -0.02975};
  centerpoints[3] = {0.03653, 0.0152, -0.005};
  centerpoints[4] = {-0.017263, 0.0355, 0.0193};

  std::vector<T> radii(5);
  radii[0] = 0.001111;
  radii[1] = 0.001421;
  radii[2] = 0.002200;
  radii[3] = 0.001394;
  radii[4] = 0.002518;
  IndicatorCircle3D<T> inletCircle(centerpoints[0], normals[0], radii[0]);
  IndicatorCircle3D<T> inletCircle2(centerpoints[1], normals[1], radii[1]);
  IndicatorCircle3D<T> inletCircle3(centerpoints[2], normals[2], radii[2]);
  IndicatorCircle3D<T> outletCircle(centerpoints[3], normals[3], radii[3]);
  IndicatorCircle3D<T> outletCircle2(centerpoints[4], normals[4], radii[4]);

  IndicatorCylinder3D<T> inletCylinder(inletCircle, 2. * converter.getPhysDeltaX());
  IndicatorCylinder3D<T> inletCylinder2(inletCircle2, 2. * converter.getPhysDeltaX());
  IndicatorCylinder3D<T> inletCylinder3(inletCircle3, 2. * converter.getPhysDeltaX());
  IndicatorCylinder3D<T> outletCylinder(outletCircle, 2. * converter.getPhysDeltaX());
  IndicatorCylinder3D<T> outletCylinder2(outletCircle2, 2. * converter.getPhysDeltaX());

  superGeometry.rename(2, 4, 1, inletCylinder);
  superGeometry.rename(2, 4, 1, inletCylinder2);
  superGeometry.rename(2, 3, 1, inletCylinder3);
  superGeometry.rename(2, 4, 1, outletCylinder);
  superGeometry.rename(2, 4, 1, outletCylinder2);

  superGeometry.clean();
  superGeometry.checkForErrors();
  superGeometry.print();
}

// Set up the geometry of the simulation
void prepareLattice(SuperLattice<T, DESCRIPTOR>& sLattice, UnitConverter<T, DESCRIPTOR> const& converter,
                    SuperGeometry<T, 3>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  // Material=1 -->bulk dynamics
  sLattice.defineDynamics<BGKdynamics>(superGeometry, 1);

  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);

  // Material=3 -->fixed velocity
  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 3);

  // Material=4 -->fixed pressure
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);

  // Material=5 -->bounce back
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 5);

  // Initial conditions
  AnalyticalConst3D<T, T> rhoF(1);
  AnalyticalConst3D<T, T> uF(T(0), T(0), T(0));

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU(superGeometry, 1, rhoF, uF);
  sLattice.iniEquilibrium(superGeometry, 1, rhoF, uF);

  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues(SuperLattice<T, DESCRIPTOR>& sLattice, UnitConverter<T, DESCRIPTOR> const& converter, int iT,
                       SuperGeometry<T, 3>& superGeometry)
{
  OstreamManager clout(std::cout, "setBoundaryValues");

  // No of time steps for smooth start-up
  int iTmaxStart = converter.getLatticeTime(maxPhysT * 0.4);
  int iTupdate   = 30;

  if (iT % iTupdate == 0 && iT <= iTmaxStart) {
    PolynomialStartScale<T, int> StartScale(iTmaxStart, T(1));

    // Creates and sets the Poiseuille inflow profile using functors
    int iTvec[1] = {iT};
    T   frac[1]  = {};
    StartScale(frac, iTvec);
    std::vector<T> maxVelocity(3, 0);
    maxVelocity[0]                      = 2.25 * frac[0] * converter.getCharLatticeVelocity();
    T                     distance2Wall = converter.getPhysDeltaX() / 2.;
    CirclePoiseuille3D<T> poiseuilleU(superGeometry, 3, maxVelocity[0], distance2Wall);
    sLattice.defineU(superGeometry, 3, poiseuilleU);

    // Update velocity on GPU
    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
        ProcessingContext::Simulation);
  }
}

// Computes the pressure drop between the voxels before and after the cylinder
void getResults(SuperLattice<T, DESCRIPTOR>& sLattice, UnitConverter<T, DESCRIPTOR> const& converter, int iT,
                SuperGeometry<T, 3>& superGeometry, util::Timer<T>& timer, STLreader<T>& stlReader,
                VTUsurfaceWriter<T>& vtuWriter)
{

  OstreamManager clout(std::cout, "getResults");

  const std::size_t vtkIter  = converter.getLatticeTime(.4);
  const std::size_t statIter = converter.getLatticeTime(.4);

  SuperVTMwriter3D<T> vtmWriter("aneurysm");

  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
  SuperLatticeDensity3D                     densityF(sLattice);
  AnalyticalFfromSuperF3D                   smoothDensityF(densityF);

  SuperLatticeStress3D    stressF(sLattice);
  AnalyticalFfromSuperF3D smoothStressF(stressF);

  PhysWallShearStressOnSurface3D<T, DESCRIPTOR> interpolatedWssF(converter, smoothDensityF, smoothStressF, stlReader);
  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(pressure);
  interpolatedWssF.getName() = "WallShearStress";

  if (iT == 0) {
    vtmWriter.createMasterFile();
    vtuWriter.createMasterFile();
    vtuWriter.addFunctor(interpolatedWssF);
    vtuWriter.addFunctor(velocity);
    vtuWriter.addFunctor(pressure);
  }

  // Writes the vtu files
  if (iT % vtkIter == 0) {
    // Send values from GPU to CPU for evaluation
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    vtuWriter.write(iT);
  }

  // Writes the vtk files
  if ( iT%vtkIter==0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    sLattice.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriter("aneurysm");
      SuperLatticePhysVelocity3D velocity(sLattice, sLattice.getUnitConverter());
      SuperLatticePhysPressure3D pressure(sLattice, sLattice.getUnitConverter());
      vtmWriter.addFunctor(velocity);
      vtmWriter.addFunctor(pressure);
      task(vtmWriter, iT);
  });
}

  // Writes output on the console
  if (iT % statIter == 0) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    // Timer console output
    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
  }
}

int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");

  UnitConverter<T, DESCRIPTOR> converter(
      (T)L,               // physDeltaX: spacing between two lattice cells in [m]
      (T)CFL * L / 0.2,   // physDeltaT: time step in [s]
      (T)0.007,           // charPhysLength: reference length of simulation geometry in [m]
      (T)0.2,             // charPhysVelocity: highest expected velocity during simulation in [m/s]
      (T)0.2 * 0.05 / Re, // physViscosity: physical kinematic viscosity in [m^2/s] 3.7e-6
      (T)1080.0           // physDensity: physical density in [kg/m^3]
  );
  converter.print();

  // === 2nd Step: Prepare Geometry ===
  // Instantiation of the STLreader class
  // file name, voxel size in meter, stl unit in meter, outer voxel no., inner voxel no.
  // Reading the pipe STL
  STLreader<T> aneurysm("aneurysm3d.stl", converter.getPhysDeltaX(), 1);
  // Extending the geometry with 1 cell for outer boundaries
  IndicatorLayer3D<T> extendedDomain(aneurysm, 2. * converter.getPhysDeltaX());

  // Instantiation of a cuboidDecomposition with weights
  CuboidDecomposition3D<T> cuboidDecomposition(extendedDomain, converter.getPhysDeltaX(), singleton::mpi().getSize());

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  // Instantiation of a superGeometry
  SuperGeometry<T, 3> superGeometry(cuboidDecomposition, loadBalancer);

  prepareGeometry(converter, extendedDomain, aneurysm, superGeometry);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice(converter, superGeometry);

  //prepareLattice and set boundaryCondition
  prepareLattice(sLattice, converter, superGeometry);

  VTUsurfaceWriter<T> vtuWriter("Mesh", sLattice.getCuboidDecomposition(), sLattice.getLoadBalancer());
  vtuWriter.addSTL(aneurysm);
  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(sLattice, converter, iT, superGeometry);

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults(sLattice, converter, iT, superGeometry, timer, aneurysm, vtuWriter);
  }

  timer.stop();
  timer.printSummary();
}
