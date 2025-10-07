/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C)
 *  2023      Nicolas Hafen, Frantisek Prinz, Mathias J. Krause
 *  2011-2016 Thomas Henn, Mathias J. Krause, Marie-Luise Maier
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

/* bifurcation3d.cpp:
 * This example examines a steady particulate flow past a bifurcation. At the inlet,
 * an inflow condition with grad_n u = 0 and rho = 1 is implemented.
 * At both outlets, a Poiseuille profile is imposed on the velocity.
 * After a start time, particles are put into the bifurcation at the
 * inlet and experience a stokes drag force.
 *
 * A publication using the same geometry can be found here:
 * http://link.springer.com/chapter/10.1007/978-3-642-36961-2_5
 *  *
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace olb::names;
using namespace particles;
using namespace particles::subgrid;
using namespace particles::communication;
using namespace particles::dynamics;
using namespace particles::creators;
using namespace particles::io;

using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q19<>>
>;

// === Step 1: Declarations ===
namespace olb::parameters {

struct PART_RADIUS : public descriptors::FIELD_BASE<1> {};
struct PART_RHO : public descriptors::FIELD_BASE<1> {};
struct IT_PERIOD : public descriptors::FIELD_BASE<1> {};
struct INLET_RADIUS : public descriptors::FIELD_BASE<1> {};
struct INLET_CENTER : public descriptors::FIELD_BASE<0, 1> {};
struct INLET_NORMAL : public descriptors::FIELD_BASE<0, 1> {};
struct OUTLET_RADIUS0 : public descriptors::FIELD_BASE<1> {};
struct OUTLET_CENTER0 : public descriptors::FIELD_BASE<0, 1> {};
struct OUTLET_NORMAL0 : public descriptors::FIELD_BASE<0, 1> {};
struct OUTLET_RADIUS1 : public descriptors::FIELD_BASE<1> {};
struct OUTLET_CENTER1 : public descriptors::FIELD_BASE<0, 1> {};
struct OUTLET_NORMAL1 : public descriptors::FIELD_BASE<0, 1> {};

} // namespace olb::parameters

Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters)
{
  using T                 = MyCase::value_t;
  const T      physDeltaX = 2 * parameters.get<parameters::INLET_RADIUS>() / parameters.get<parameters::RESOLUTION>();
  STLreader<T> stlReader("../bifurcation3d.stl", physDeltaX);
  IndicatorLayer3D<T> extendedDomain(stlReader, physDeltaX);
  std::size_t         noOfCuboids = util::max(20, singleton::mpi().getSize());
  Mesh<T, MyCase::d>  mesh(extendedDomain, physDeltaX, noOfCuboids);
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

using T = FLOATING_POINT_TYPE;
typedef D3Q19<> DESCRIPTOR;
typedef SubgridParticle3D PARTICLETYPE;

#define BOUZIDI

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const T Re = 50;                    // Reynolds number
int N = 19;                         // resolution of the model
const T radius = 1.5e-4;            // particles radius
const T partRho = 998.2;            // particles density

const T fluidMaxPhysT = T( 5 );     // max. fluid simulation time in s, SI unit
const T particleMaxPhysT = T( 20 ); // max. particle simulation time in s, SI unit

std::size_t noOfParticles = 1000;   // total number of inserted particles

//Set capture method:
// materialCapture: based on material number
// wallCapture:     based on more accurate stl description
typedef enum {materialCapture, wallCapture} ParticleDynamicsSetup;
const ParticleDynamicsSetup particleDynamicsSetup = wallCapture;

// center of inflow and outflow regions [m]
Vector<T, 3> inletCenter( T(), T(), 0.0786395 );
Vector<T, 3> outletCenter0( -0.0235929682287551, -0.000052820468762797,
                            -0.021445708949909 );
Vector<T, 3> outletCenter1( 0.0233643529416147, 0.00000212439067050152,
                            -0.0211994104877918 );

// radii of inflow and outflow regions [m]
T inletRadius = 0.00999839;
T outletRadius0 = 0.007927;
T outletRadius1 = 0.00787134;

// normals of inflow and outflow regions
Vector<T, 3> inletNormal( T(), T(), T( -1 ) );
Vector<T, 3> outletNormal0( 0.505126, -0.04177, 0.862034 );
Vector<T, 3> outletNormal1( -0.483331, -0.0102764, 0.875377 );

//Ensure that parallel mode is used
#ifdef PARALLEL_MODE_MPI

void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T          = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();

  const T      physDeltaX = 2 * parameters.get<parameters::INLET_RADIUS>() / parameters.get<parameters::RESOLUTION>();
  STLreader<T> stlReader("../bifurcation3d.stl", physDeltaX);
  IndicatorLayer3D<T> extendedDomain(stlReader, physDeltaX);
  geometry.rename(0, 2, extendedDomain);
  geometry.rename(2, 1, stlReader);
  geometry.clean();

  // rename the material at the inlet
  IndicatorCircle3D<T>   inletCircle(parameters.get<parameters::INLET_CENTER>(),
                                     parameters.get<parameters::INLET_NORMAL>(),
                                     parameters.get<parameters::INLET_RADIUS>());
  IndicatorCylinder3D<T> inlet(inletCircle, 2 * physDeltaX);

  geometry.rename(2, 3, 1, inlet);

  // rename the material at the outlet0
  IndicatorCircle3D<T>   outletCircle0(parameters.get<parameters::OUTLET_CENTER0>(),
                                       parameters.get<parameters::OUTLET_NORMAL0>(),
                                       0.95 * parameters.get<parameters::OUTLET_RADIUS0>());
  IndicatorCylinder3D<T> outlet0(outletCircle0, 4 * physDeltaX);
  geometry.rename(2, 4, outlet0);

  // rename the material at the outlet1
  IndicatorCircle3D<T>   outletCircle1(parameters.get<parameters::OUTLET_CENTER1>(),
                                       parameters.get<parameters::OUTLET_NORMAL1>(),
                                       0.95 * parameters.get<parameters::OUTLET_RADIUS1>());
  IndicatorCylinder3D<T> outlet1(outletCircle1, 4 * physDeltaX);
  geometry.rename(2, 5, outlet1);

  IndicatorCircle3D<T>   inletCircleExtended(parameters.get<parameters::INLET_CENTER>(),
                                             parameters.get<parameters::INLET_NORMAL>(),
                                             parameters.get<parameters::INLET_RADIUS>() + 2 * physDeltaX);
  IndicatorCylinder3D<T> inletExtended(inletCircleExtended, 2 * physDeltaX);
  geometry.rename(2, 6, inletExtended);

  geometry.clean();
  geometry.innerClean(true);
  geometry.checkForErrors();

  clout << "Prepare Geometry ... OK" << std::endl;
  return;
}


void prepareLattice( MyCase& myCase )
{

  OstreamManager clout( std::cout, "prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;
  using T            = MyCase::value_t;
  auto& parameters   = myCase.getParameters();
  auto& geometry     = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(NavierStokes {});
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    int {N}, // resolution: number of voxels per charPhysL
    (T)   0.557646,                    // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   inletRadius*2.,              // charPhysLength: reference length of simulation geometry
    (T)   Re*1.5e-5/( inletRadius*2 ), // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1.5e-5,                      // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.225                        // physDensity: physical density in __kg / m^3__
  );
  auto& converter = lattice.getUnitConverter();
  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  lattice.defineDynamics<BGKdynamics>(superGeometry, 1);

  #ifdef BOUZIDI
  setBouzidiBoundary(lattice, superGeometry, 2, stlReader);
  #else
  boundary::set<boundary::BounceBack>(lattice, superGeometry, 2);
  #endif

  // Material=3 -->bulk dynamics (inflow)
  lattice.defineDynamics<BGKdynamics>(superGeometry, 3);

  // Material=4 -->bulk dynamics (outflow)
  lattice.defineDynamics<BGKdynamics>(superGeometry, 4);
  lattice.defineDynamics<BGKdynamics>(superGeometry, 5);

  // Setting of the boundary conditions
  boundary::set<boundary::InterpolatedPressure>(lattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedVelocity>(lattice, superGeometry, 4);
  boundary::set<boundary::InterpolatedVelocity>(lattice, superGeometry, 5);

  lattice.setParameter<descriptors::OMEGA>(omega);

  clout << "Prepare Lattice ... OK" << std::endl;
  return;
}


// Generates a slowly increasing sinusoidal inflow for the first iTMax timesteps
void setBoundaryValues( MyCase& myCase, std::size_t iT )
{
  OstreamManager clout( std::cout, "setBoundaryValues" );
  using T            = MyCase::value_t;
  auto& parameters   = myCase.getParameters();
  auto& geometry     = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(NavierStokes {});
  auto& converter = lattice.getUnitConverter();
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  // No of time steps for smooth start-up
  std::size_t iTmaxStart = converter.getLatticeTime( 0.8 * parameters::get<parameters::FLUID_MAX_PHYS_T>() );
  std::size_t iTperiod = parameters::get<parameters::IT_PERIOD>()/converter.getPhysDeltaT();  // amount of timesteps when new boundary conditions are reset

  if ( iT == 0 ) {

    AnalyticalConst3D<T, T> rhoF( 1 );
    std::vector<T> velocity( 3, T() );
    AnalyticalConst3D<T, T> uF( velocity );

    lattice.iniEquilibrium( geometry, 1, rhoF, uF );
    lattice.iniEquilibrium( geometry, 2, rhoF, uF );
    lattice.iniEquilibrium( geometry, 3, rhoF, uF );
    lattice.iniEquilibrium( geometry, 4, rhoF, uF );
    lattice.iniEquilibrium( geometry, 5, rhoF, uF );

    lattice.defineRhoU( geometry, 1, rhoF, uF );
    lattice.defineRhoU( geometry, 2, rhoF, uF );
    lattice.defineRhoU( geometry, 3, rhoF, uF );
    lattice.defineRhoU( geometry, 4, rhoF, uF );
    lattice.defineRhoU( geometry, 5, rhoF, uF );

    // Make the lattice ready for simulation
    lattice.initialize();
  }

  else if ( iT <= iTmaxStart && iT % iTperiod == 0 ) {
    SinusStartScale<T, int> startScale( iTmaxStart, T( 1 ) );
    std::size_t iTvec[1] = { iT };
    T frac[1] = { T( 0 ) };
    startScale( frac, iTvec );
    T maxVelocity = frac[0] * converter.getCharLatticeVelocity() * 3. / 4.
                    * util::pow( parameters::get<parameters::INLET_RADIUS>(), 2 ) / util::pow( parameters::get<parameters::OUTLET_RADIUS_0>(), 2 );

    Vector<T,3> outletNormal0 = parameters::get<parameters::OUTLET_NORMAL0>();
    Vector<T,3> outletNormal1 = parameters::get<parameters::OUTLET_NORMAL1>();
    Vector<T,3> outletCenter0 = parameters::get<parameters::OUTLET_CENTER0>();
    Vector<T,3> outletCenter1 = parameters::get<parameters::OUTLET_CENTER1>();
    T outletRadius0 = parameters::get<parameters::OUTLET_RADIUS0>();
    T outletRadius1 = parameters::get<parameters::OUTLET_RADIUS1>();
    CirclePoiseuille3D<T> poiseuilleU4( outletCenter0[0], outletCenter0[1],
                                        outletCenter0[2], outletNormal0[0],
                                        outletNormal0[1], outletNormal0[2],
                                        outletRadius0 * 0.95, -maxVelocity );

    CirclePoiseuille3D<T> poiseuilleU5( outletCenter1[0], outletCenter1[1],
                                        outletCenter1[2], outletNormal1[0],
                                        outletNormal1[1], outletNormal1[2],
                                        outletRadius1 * 0.95, -maxVelocity );

    superLattice.defineU( geometry, 4, poiseuilleU4 );
    superLattice.defineU( geometry, 5, poiseuilleU5 );
  }
}

//Prepare particles
void prepareParticles( MyCase& myCase,
  SuperParticleSystem<T,PARTICLETYPE>& superParticleSystem,
  SolidBoundary<T,3>& wall,
  ParticleDynamicsSetup particleDynamicsSetup )
{
  OstreamManager clout( std::cout, "prepareParticles" );
  clout << "Prepare Particles ..." << std::endl;

  //Create material indicators for particle boundaries
  std::vector<int> wallMaterials {2};
  std::shared_ptr<SuperIndicatorMaterial<T,3>> wallMaterialIndicator =
    std::make_shared<SuperIndicatorMaterial<T,3>>(superGeometry, wallMaterials);
  std::vector<int> outletMaterials {4,5};
  std::shared_ptr<SuperIndicatorMaterial<T,3>> outletMaterialIndicator =
    std::make_shared<SuperIndicatorMaterial<T,3>>(superGeometry, outletMaterials);

  //Add selected particle dynamics
  if (particleDynamicsSetup==wallCapture){
    //Create verlet dynamics with material aware wall capture
    superParticleSystem.defineDynamics<
      VerletParticleDynamicsMaterialAwareWallCaptureAndEscape<T,PARTICLETYPE>>(
        wall, wallMaterialIndicator, outletMaterialIndicator);
  } else {
    //Create verlet dynamics with material capture
    superParticleSystem.defineDynamics<
      VerletParticleDynamicsMaterialCaptureAndEscape<T,PARTICLETYPE>>(
        wallMaterialIndicator, outletMaterialIndicator);
  }

  // particles generation at inlet3
  Vector<T, 3> center( parameters::get<parameters::INLET_CENTER>() );
  center[2] = 0.074;
  IndicatorCircle3D<T> inflowCircle( center, parameters::get<parameters::INLET_NORMAL>(), parameters::get<parameters::INLET_RADIUS>() -
                                     converter.getPhysDeltaX() * 2.5 );
  IndicatorCylinder3D<T> inletCylinder( inflowCircle, 0.01 *
                                        converter.getPhysDeltaX() );

  //Add particles
  Randomizer<T> randomizer;
  addParticles( superParticleSystem, inletCylinder, parameters::get<parameters::PARTICLE_DENSITY>(),
                parameters::get<parameters::PARTICLE_RADIUS>(), parameters::get<parameters::NO_OF_PARTICLES>(), randomizer );

  //Print super particle system summary
  superParticleSystem.print();
  clout << "Prepare Particles ... OK" << std::endl;
}

bool getResults( MyCase& myCase, std::size_t iT,
                 Timer<double>& fluidTimer, STLreader<T>& stlReader, bool fluidExists,
                 SuperParticleSystem<T,PARTICLETYPE>& superParticleSystem,
                 Timer<double>& particleTimer )
{
  OstreamManager clout( std::cout, "getResults" );
  SuperVTMwriter3D<T> vtmWriter( "bifurcation3d" );
  SuperVTMwriter3D<T> vtmWriterStartTime( "startingTimeBifurcation3d" );

  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( superLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( superLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  vtmWriterStartTime.addFunctor( velocity );
  vtmWriterStartTime.addFunctor( pressure );

  VTUwriter<T,PARTICLETYPE,true> superParticleWriter("particles_master",false,false);

  SuperParticleGroupedFieldF<T,PARTICLETYPE,GENERAL,POSITION>
    particlePosition( superParticleSystem );
  SuperParticleGroupedFieldF<T,PARTICLETYPE,PHYSPROPERTIES,MASS>
    particleMass( superParticleSystem );
  SuperParticleGroupedFieldF<T,PARTICLETYPE,PHYSPROPERTIES,RADIUS>
    particleRadius( superParticleSystem );
  SuperParticleGroupedFieldF<T,PARTICLETYPE,MOBILITY,VELOCITY>
    particleVelocity( superParticleSystem );
  SuperParticleGroupedFieldF<T,PARTICLETYPE,FORCING,FORCE>
    particleForce( superParticleSystem );
  SuperParticleGroupedFieldF<T,PARTICLETYPE,DYNBEHAVIOUR,ACTIVE>
    particleActivity( superParticleSystem );
  superParticleWriter.addFunctor(particlePosition,"Position");
  superParticleWriter.addFunctor(particleMass,"Mass");
  superParticleWriter.addFunctor(particleRadius,"Radius");
  superParticleWriter.addFunctor(particleVelocity,"Velocity");
  superParticleWriter.addFunctor(particleForce,"Force");
  superParticleWriter.addFunctor(particleActivity,"Active");

  std::size_t fluidMaxT = converter.getLatticeTime( parameters::get<parameters::FLUID_MAX_PHYS_T>() );

  if ( iT == 0 ) {
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( superLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( superLattice );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
    superParticleWriter.createMasterFile();
    vtmWriterStartTime.createMasterFile();

    clout << "N=" << parameters::get<parameters::RESOUTION>() <<"; maxTimeSteps(fluid)="
          << converter.getLatticeTime( fluidMaxPhysT ) << "; noOfCuboid="
          << geometry.getCuboidDecomposition().size() << "; Re=" << parameters::get<parameters::REYNOLDS>()
          <<  "; noOfParticles=" << parameters::get<parameters::NO_OF_PARTICLES>() << "; maxTimeSteps(particle)="
          << converter.getLatticeTime( parameters::get<parameters::PARTICLE_MAX_PHYS_T>() )
          << "; St=" << ( 2.* parameters::get<parameters::PARTICLE_DENSITY>()
              * parameters::get<parameters::PARTICLE_RADIUS>() *parameters::get<parameters::PARTICLE_RADIUS>()  * converter.getCharPhysVelocity() )
              / ( 9.*converter.getPhysViscosity()*converter.getPhysDensity()*converter.getCharPhysLength() ) << std::endl;
  }

  // Writes the .vtk and .gif files
  if ( iT % iTperiod == 0 ) {
    if ( !fluidExists && iT <= fluidMaxT ) {
      vtmWriterStartTime.write(iT);
      SuperEuklidNorm3D<T> normVel( velocity );
      BlockReduction3D2D<T> planeReduction( normVel, {0, -1, 0}, 600, BlockDataSyncMode::ReduceOnly );
      // write output as JPEG
      heatmap::write(planeReduction, iT);
    }
    if (iT > fluidMaxT) {
      // only write .vtk-files after the fluid calculation is finished
      vtmWriter.write(iT - fluidMaxT);
    }
  }

  // Writes output on the console for the fluid phase
  if (iT < converter.getLatticeTime( fluidMaxPhysT ) && iT%iTperiod == 0 ) {

    // Timer statics
    fluidTimer.update( iT );
    fluidTimer.printStep();

    // Lattice statistics
    superLattice.getStatistics().print( iT, converter.getPhysTime( iT ) );

    // Flux at the inlet and outlet regions
    const std::vector<int> materials = { 1, 3, 4, 5 };

    IndicatorCircle3D<T> inlet(
      inletCenter + 2. * converter.getPhysDeltaX() * inletNormal,
      inletNormal, inletRadius + 2. * converter.getPhysDeltaX() );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxInflow( superLattice, converter, superGeometry, inlet, materials );
    vFluxInflow.print( "inflow", "ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxInflow( superLattice, converter, superGeometry, inlet, materials );
    pFluxInflow.print( "inflow", "N", "Pa" );

    IndicatorCircle3D<T> outlet0(
      outletCenter0 + 2. * converter.getPhysDeltaX() * outletNormal0,
      outletNormal0, outletRadius0 + 2. * converter.getPhysDeltaX() );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOutflow0( superLattice, converter, superGeometry, outlet0, materials );
    vFluxOutflow0.print( "outflow0", "ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxOutflow0( superLattice, converter, superGeometry, outlet0, materials );
    pFluxOutflow0.print( "outflow0", "N", "Pa" );

    IndicatorCircle3D<T> outlet1(
      outletCenter1 + 2. * converter.getPhysDeltaX() * outletNormal1,
      outletNormal1, outletRadius1 + 2. * converter.getPhysDeltaX() );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOutflow1( superLattice, converter, superGeometry, outlet1, materials );
    vFluxOutflow1.print( "outflow1", "ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxOutflow1( superLattice, converter, superGeometry, outlet1, materials );
    pFluxOutflow1.print( "outflow1", "N", "Pa" );
  }

  // Writes output on the console for the fluid phase
  if ( iT >= converter.getLatticeTime( fluidMaxPhysT ) &&
       (iT%iTperiod == 0 || iT == converter.getLatticeTime( fluidMaxPhysT )) ) {
    // advance particle timer
    particleTimer.print( iT - fluidMaxT );

    //delete invalidated particles
    purgeInvalidParticles<T,PARTICLETYPE>( superParticleSystem );

    //Define materials for capture rate
    std::vector<int> materialsOutout {4,5};
    SuperIndicatorMaterial<T,3> materialIndicatorOutput (superGeometry, materialsOutout);

    std::size_t noActive;
    captureStatistics( superParticleSystem, materialIndicatorOutput, noActive );


    superParticleWriter.write( iT - fluidMaxT );

    // true as long as certain amount of active particles
    if ( noActive < 0.001 * parameters::get<parameters::NO_OF_PARTICLES>()
         && iT > 0.9*converter.getLatticeTime( parameters::get<parameters::MAX_PHYS_T>() + particleMaxPhysT ) ) {
      return false;
    }
    //Additional criterion added 02.02.23
    if ( noActive==0 ) {
      return false;
    }
  }
  return true;
}


#endif //PARALLEL_MODE_MPI





int main( int argc, char* argv[] )
{

  initialize( &argc, &argv );

  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout, "main" );

  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<REYNOLDS>(50);
    myCaseParameters.set<RESOLUTION>(19);
    myCaseParameters.set<IT_PERIOD>(0.1419);   // time interval in s after which boundary conditions are updated
    myCaseParameters.set<PART_RADIUS>(1.5e-4); // particles radius
    myCaseParameters.set<PART_RHO>(998.2);     // particles density
    myCaseParameters.set<FLUID_MAX_PHYS_T>(5);      // max. simulation time in s, SI unit
    myCaseParameters.set<PARTICLE_MAX_PHYS_T>(20); // max. particle simulation time in s, SI unit
    myCaseParameters.set<NO_OF_PARTICLES>(1000); // total number of inserted particles
    myCaseParameters.set<INLET_CENTER>({0., 0., 0.0786395});
    myCaseParameters.set<OUTLET_CENTER0>({-0.0235929682287551, -0.000052820468762797, -0.021445708949909});
    myCaseParameters.set<OUTLET_CENTER1>({0.0233643529416147, 0.00000212439067050152, -0.0211994104877918});
    myCaseParameters.set<INLET_RADIUS>(0.00999839);
    myCaseParameters.set<OUTLET_RADIUS0>(0.007927);
    myCaseParameters.set<OUTLET_RADIUS1>(0.00787134);
    myCaseParameters.set<INLET_NORMAL>({0., 0., -1.});
    myCaseParameters.set<OUTLET_NORMAL0>({0.505126, -0.04177, 0.862034});
    myCaseParameters.set<OUTLET_NORMAL1>({-0.483331, -0.0102764, 0.875377});
  }
  myCaseParameters.fromCLI(argc, argv);

  Mesh mesh = createMesh(myCaseParameters);
  MyCase myCase = MyCase(myCaseParameters, mesh);



#ifdef PARALLEL_MODE_MPI
  STLreader<T> stlReader( "../bifurcation3d.stl", converter.getPhysDeltaX() );
  IndicatorLayer3D<T> extendedDomain( stlReader,
                                      converter.getPhysDeltaX() );

  const unsigned latticeMaterial = 2; //Material number of wall
  const unsigned contactMaterial = 0; //Material identifier (only relevant for contact model)
  SolidBoundary<T,3> wall( std::make_unique<IndicInverse<T,3>>(stlReader),
                           latticeMaterial, contactMaterial );

  prepareGeometry( myCase );

  prepareLattice(  myCase );

  SuperParticleSystem<T,PARTICLETYPE> superParticleSystem(superGeometry);

  ParticleManager<T,DESCRIPTOR,PARTICLETYPE> particleManager(
    superParticleSystem, superGeometry, superLattice, converter);

  prepareParticles( myCase, superParticleSystem, wall,
    particleDynamicsSetup, randomizer );



  Timer<double> fluidTimer( converter.getLatticeTime( fluidMaxPhysT ),
                            superGeometry.getStatistics().getNvoxel() );

  Timer<double> particleTimer( converter.getLatticeTime( particleMaxPhysT ),
                               noOfParticles );
  fluidTimer.start();

  std::size_t iT = 0;
  // amount of timesteps when getResults rewrites data
  int iTperiod = converter.getLatticeTime( .2 );

  bool fluidExists = true;

  // checks whether there is already data of the fluid from an earlier calculation
  if ( !( superLattice.load( "fluidSolution_N"+std::to_string(N) ) ) ) {

    fluidExists = false;

    // if there is no data available, it is generated
    for ( ; iT <= converter.getLatticeTime( fluidMaxPhysT ); ++iT ) {

      // during run up time boundary values are set, collide and stream step,
      // results of fluid, afterwards only particles are simulated
      setBoundaryValues( myCase, iT );
      superLattice.collideAndStream();

      getResults( myCase, iT, iTperiod,
                  fluidTimer, stlReader, fluidExists,
                  superParticleSystem, particleTimer );
    }

    fluidTimer.stop();
    fluidTimer.printSummary();

    superLattice.communicate();
    // calculated results are written in a file
    superLattice.save( "fluidSolution_N"+std::to_string(N) );
  }

  // if there exists already data of the fluid from an earlier calculation, this is used
  else {

    iT = converter.getLatticeTime( fluidMaxPhysT );
    getResults( superLattice, converter, iT,
                iTperiod, superGeometry, fluidTimer, stlReader, fluidExists,
                superParticleSystem, particleTimer);

  }

  // initialize particle velocity
  initializeParticleVelocity( superLattice, superGeometry, converter, superParticleSystem );

  // after the fluid calculation, particle simulation starts
  particleTimer.start();

  for ( ; iT <= converter.getLatticeTime( fluidMaxPhysT + particleMaxPhysT ); ++iT ) {
    // particles simulation starts after run up time is over
    particleManager.execute<
      couple_lattice_to_particles<T,DESCRIPTOR,PARTICLETYPE>,
      process_dynamics<T,PARTICLETYPE>,
      update_particle_core_distribution<T,PARTICLETYPE>
    >();

    if ( !getResults( superLattice, converter, iT,
                      iTperiod, superGeometry, fluidTimer,
                      stlReader, fluidExists,
                      superParticleSystem,
                      particleTimer) ){
      break;
    }
  }
  particleTimer.stop();
  particleTimer.printSummary();
#else
  std::cerr << std::endl
            << "ERROR: Subgrid particles can only be used with MPI!"
            << std::endl << std::endl;
#endif //PARALLEL_MODE_MPI
}
