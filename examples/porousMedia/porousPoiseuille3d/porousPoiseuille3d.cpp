/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2019 Fabian Klemens, Davide Dapelo, Mathias J. Krause
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

/* porousPoiseuille3d.cpp:
 * This example examines a 3D Poseuille flow with porous media.
 * Two porous media LB methods can be used here:
 * Spaid and Phelan (doi:10.1063/1.869392), or
 * Guo and Zhao (doi:10.1103/PhysRevE.66.036304)
 * 
 * Case-specific arguments:
 * BOUNDARY_TYPE: 0=bounceBack, 1=local, 2=interpolated
 * POROSITY_TYPE: 0=BGK, 1=SpaidPhelan, 2=GuoZhao
 * PERMEABILIY: default 1e-2
 * Default: Resolution=50, Permeability=1e-2, BoundaryTpe=interpolated, PorosityType=GuoZhao
 */

#include <olb.h>
using namespace olb;
using namespace olb::names;
enum class BoundaryType: int {
  bounceBack    = 0,
  local         = 1,
  interpolated  = 2
};
enum class PorosityType: int {
  BGK             = 0,
  SPAID_PHELAN    = 1,
  GUO_ZHAO        = 2
};
namespace olb::parameters {
  struct BOUNDARY_TYPE  : public descriptors::TYPED_FIELD_BASE<BoundaryType,1> { };
  struct POROSITY_TYPE  : public descriptors::TYPED_FIELD_BASE<PorosityType,1> { };
  struct CONVERGED      : public descriptors::TYPED_FIELD_BASE<bool,1> { };
  struct PERMEABILITY   : public descriptors::FIELD_BASE<1> { };
  struct EPSILON        : public descriptors::FIELD_BASE<1> { };
  struct INITIAL_PRESSURE_L : public descriptors::FIELD_BASE<1> { };
  struct PRESSURE_GRADIENT  : public descriptors::FIELD_BASE<1> { };
  struct VISCOSITY       : public descriptors::FIELD_BASE<1> { };
  struct CONVERGENCE_CHECK_T        : public descriptors::FIELD_BASE<1> { };
  struct CONVERGENCE_CHECK_RESIDUUM : public descriptors::FIELD_BASE<1> { };
}

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D2Q9<>>
>;

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector origin{0, 0};
  IndicatorCuboid2D<T> cuboid(extent, origin);

  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

template <typename T, typename DESCRIPTOR>
class PhysicalToLatticeVelocityF3D: public AnalyticalF3D<T,T> {
protected:
  AnalyticalF3D<T,T>* f;
  UnitConverter<T,DESCRIPTOR> converter;

public:
  PhysicalToLatticeVelocityF3D(AnalyticalF3D<T,T>* f_, UnitConverter<T,DESCRIPTOR> const& converter_)
    : AnalyticalF3D<T,T>(3), f(f_), converter(converter_) {};

  bool operator()(T output[], const T x[]) override
  {
    (*f)(output, x);
    for (int i=0; i<3; ++i) {
      output[i] = converter.getLatticeVelocity( output[i] );
    }
    return true;
  };
};


// Parameters for the simulation setup
const T length  = 2.;         // length of the pie
const T diameter  = 1.;       // diameter of the pipe
int N = 21;                   // resolution of the model
const T physU = 1.;           // physical velocity
const T Re = 1.;              // Reynolds number
const T physRho = 1.;         // physical density
const T tau = 0.8;            // lattice relaxation time
const T maxPhysT = 20.;       // max. simulation time in s, SI unit
const T residuum = 1e-5;      // residuum for the convergence check

// Scaled Parameters
const T radius  = diameter/2.;            // radius of the pipe
const T physInterval = 0.0125*maxPhysT;   // interval for the convergence check in s

// Approximation of the modified Bessel function (doi:10.1088/1742-6596/1043/1/012003)
template <typename T>
T besselApprox( T x )
{
  return util::cosh(x) / util::pow( 1 + 0.25*util::pow(x,2), 0.25 ) * ( 1 + 0.24273*util::pow(x,2) )/( 1 + 0.43023*util::pow(x,2) );
}

/// Functional to calculate velocity profile on pipe with porous media.
template <typename T>
class PorousPoiseuille3D : public AnalyticalF3D<T,T> {
protected:
  T Kin, mu, dp, epsilon, radius;
  bool trunc;

public:
  PorousPoiseuille3D( MyCase& myCase, T radius_ )
    : AnalyticalF3D<T,T>(3), radius(radius_) 
  {
    auto& parameters = myCase.getParameters();
    Kin       = parameters.get<parameters::PERMEABILITY>();
    dp        = parameters.get<parameters::PRESSURE_GRADIENT>();
    mu        = parameters.get<parameters::VISCOSITY>();
    epsilon   = parameters.get<parameters::EPSILON>();
  }

  bool operator()(T output[], const T x[]) override
  {
    T r = util::sqrt(epsilon/Kin);
    T dist = util::sqrt( util::pow(x[1]-radius, 2.) + util::pow(x[2]-radius, 2.) );
    output[0] = Kin / mu * dp * ( 1. - besselApprox(r*dist) / besselApprox(r*radius)  );
    output[1] = 0.;
    output[2] = 0.;

    if ( dist > radius ) {
      output[0] = 0.;
    }

    return true;
  };
};

void prepareGeometry( MyCase& myCase ) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry    = myCase.getGeometry();
  auto& parameters  = myCase.getParameters();
  T phyDeltaX       = parameters.get<parameters::PHYS_DELTA_X>();
  T physDeltaX2     = parameters.get<parameters::PHYS_DELTA_X>() / 2.;
  T radius          = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  T length          = parameters.get<parameters::DOMAIN_EXTENT>()[1];

  Vector<T, 3> center0(-phyDeltaX * 0.2, radius, radius);
  Vector<T, 3> center1(length, radius, radius);
  IndicatorCylinder3D<T> pipe(center0, center1, radius);

  geometry.rename(0, 2);
  geometry.rename(2, 1, pipe);

  Vector<T, 3> origin(0, radius, radius);
  Vector<T, 3> extend = origin;

  // Set material number for inflow
  origin[0] = -physDeltaX2;
  extend[0] = physDeltaX2 * 2;
  IndicatorCylinder3D<T> inflow(origin, extend, radius);
  geometry.rename(2, 3, 1, inflow);

  // Set material number for outflow
  origin[0] = length - physDeltaX2;
  extend[0] = physDeltaX;
  IndicatorCylinder3D<T> outflow(extend, origin, radius);
  geometry.rename(2, 4, 1, outflow);

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  // Removes all not needed boundary voxels inside the surface
  geometry.innerClean();
  geometry.checkForErrors();

  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( MyCase& myCase ) {
  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});
  size_t  N   = parameters.get<parameters::RESOLUTION>();
  T       Re  = parameters.get<parameters::REYNOLDS>();

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  #ifdef GUO_ZHAO
  Dynamics<T,DESCRIPTOR>* dynamics = new DYNAMICS<T,DESCRIPTOR>(omega);
  sLattice.defineDynamics(geometry, 1, dynamics);
  #else
  sLattice.defineDynamics<DYNAMICS>(geometry, 1);
  #endif

  Vector<T, 3> center0(0, radius, radius);
  Vector<T, 3> center1(length, radius, radius);

  AnalyticalConst3D<T,T> zero(0.);
  AnalyticalConst3D<T,T> one(1.);

  T nu = (tau-0.5)/3.;
  T h = converter.getPhysDeltaX();
#ifdef SPAID_PHELAN
  T d = 1. - (h*h*nu*tau/Kin);
  clout << "Lattice Porosity: " << d << std::endl;
  clout << "Kmin: " << h*h*nu*tau << std::endl;
  if (Kin < h*h*nu*tau) {
    clout << "WARNING: Chosen K is too small!" << std::endl;
    exit(1);
  }
  AnalyticalConst3D<T,T> porosity( d );
  for (int i: {
         0,1,2,3,4
       }) {
    sLattice.defineField<POROSITY>(geometry, i, porosity);
  }
#elif defined GUO_ZHAO
  AnalyticalConst3D<T,T> Nu(nu);
  AnalyticalConst3D<T,T> k(Kin/(h*h));
  AnalyticalConst3D<T,T> eps(epsilon);
  for (int i: {
         0,1,2,3,4
       }) {
    sLattice.defineField<EPSILON>(geometry, i, eps);
    sLattice.defineField<NU>(geometry, i, Nu);
    sLattice.defineField<K>(geometry, i, k);
  }
#endif

  // Bouzidi
  center0[0] -= 0.5*converter.getPhysDeltaX();
  center1[0] += 0.5*converter.getPhysDeltaX();
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  setBouzidiBoundary(sLattice, geometry, 2, pipe);
  // Interp
  //sLattice.defineDynamics<DYNAMICS>(geometry, 2);
  boundary::set<boundary::InterpolatedVelocity>(sLattice, geometry, 2);

  // Material=3 --> bulk dynamics
  #ifdef GUO_ZHAO
  sLattice.defineDynamics(geometry, 3, dynamics);
  #else
  sLattice.defineDynamics<DYNAMICS>(geometry, 3);
  #endif
  boundary::set<boundary::InterpolatedVelocity>(sLattice, geometry, 3);

  // Material=4 --> bulk dynamics
  #ifdef GUO_ZHAO
  sLattice.defineDynamics(geometry, 4, dynamics);
  #else
  sLattice.defineDynamics<DYNAMICS>(geometry, 4);
  #endif
  boundary::set<boundary::InterpolatedPressure>(sLattice, geometry, 4);

  // Initial conditions
  // Pressure for Poiseuille flow with maximum velocity of charU at K->infty
  p0 = 4. * converter.getPhysViscosity() * converter.getCharPhysVelocity() * length / (radius * radius);
  T p0L = converter.getLatticePressure(p0);
  AnalyticalLinear3D<T, T> rho(-p0L / length * invCs2<T,DESCRIPTOR>(), 0, 0, p0L * invCs2<T,DESCRIPTOR>() + 1);

  dp = p0/length;
  mu = converter.getPhysViscosity()*converter.getPhysDensity();

  //CirclePoiseuille3D<T> uSol( {0., radius, radius}, {1, 0, 0}, converter.getCharPhysVelocity(), radius );
  PorousPoiseuille3D<T> uSol( Kin, mu, dp, radius );
  PhysicalToLatticeVelocityF3D<T> u(&uSol, converter);

  // Initialize all values of distribution functions to their local equilibrium
  for (int i: {
         0,1,2,3,4
       }) {
    sLattice.defineRhoU(geometry, i, rho, u);
    sLattice.iniEquilibrium(geometry, i, rho, u);
  }

  sLattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

/// Compute error norms
void error( SuperGeometry<T,3>& geometry,
            SuperLattice<T, DESCRIPTOR>& sLattice,
            UnitConverter<T,DESCRIPTOR> const& converter,
            AnalyticalF3D<T,T>& uSol)
{

  OstreamManager clout( std::cout,"error" );


  int tmp[]= { };
  T result[2]= { };
  auto indicatorF = geometry.getMaterialIndicator(1);
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> u( sLattice, converter );

  // Velocity error
  SuperAbsoluteErrorL1Norm3D<T> absVelocityErrorNormL1(u, uSol, indicatorF);
  absVelocityErrorNormL1(result, tmp);
  clout << "velocity-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm3D<T> relVelocityErrorNormL1(u, uSol, indicatorF);
  relVelocityErrorNormL1(result, tmp);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm3D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, tmp);
  clout << "velocity-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm3D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, tmp);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm3D<T> absVelocityErrorNormLinf(u, uSol, indicatorF);
  absVelocityErrorNormLinf(result, tmp);
  clout << "velocity-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm3D<T> relVelocityErrorNormLinf(u, uSol, indicatorF);
  relVelocityErrorNormLinf(result, tmp);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;

  // Pressure error
  T p0L = converter.getLatticePressure(p0);
  AnalyticalLinear3D<T,T> pressureSol( -converter.getPhysPressure( p0L )/length, 0, 0, converter.getPhysPressure( p0L ) );
  SuperLatticePhysPressure3D<T,DESCRIPTOR> pressure( sLattice,converter );

  SuperAbsoluteErrorL1Norm3D<T> absPressureErrorNormL1(pressure, pressureSol, indicatorF);
  absPressureErrorNormL1(result, tmp);
  clout << "pressure-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm3D<T> relPressureErrorNormL1(pressure, pressureSol, indicatorF);
  relPressureErrorNormL1(result, tmp);
  clout << "; pressure-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm3D<T> absPressureErrorNormL2(pressure, pressureSol, indicatorF);
  absPressureErrorNormL2(result, tmp);
  clout << "pressure-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm3D<T> relPressureErrorNormL2(pressure, pressureSol, indicatorF);
  relPressureErrorNormL2(result, tmp);
  clout << "; pressure-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm3D<T> absPressureErrorNormLinf(pressure, pressureSol, indicatorF);
  absPressureErrorNormLinf(result, tmp);
  clout << "pressure-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm3D<T> relPressureErrorNormLinf(pressure, pressureSol, indicatorF);
  relPressureErrorNormLinf(result, tmp);
  clout << "; pressure-Linf-error(rel)=" << result[0] << std::endl;
}


// Output to console and files
void getResults( SuperLattice<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry<T,3>& geometry, util::Timer<T>& timer, bool hasConverged )
{

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter3D<T> vtmWriter( "porousPoiseuille3d" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  //CirclePoiseuille3D<T> uSol( {0., radius, radius}, {1, 0, 0}, converter.getCharPhysVelocity(), radius );
  PorousPoiseuille3D<T> uSol( Kin, mu, dp, radius );
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> uSolF( uSol, sLattice );
  vtmWriter.addFunctor( uSolF );

  const int vtmIter  = converter.getLatticeTime( maxPhysT/20. );
  const int statIter = converter.getLatticeTime( maxPhysT/20. );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }

  // Writes the vtm files and profile text file
  if ( iT%vtmIter==0 || hasConverged ) {
    vtmWriter.write( iT );
  }


  // Writes output on the console
  if ( iT%statIter==0 || hasConverged ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Calculate inflow and outflow flux
    std::vector<int> materials = { 1, 3, 4 };
    Vector<T,3> normal( 1, 0, 0 );
    auto mode = BlockDataReductionMode::Discrete;
    Vector<T,3> posInflow = geometry.getStatistics().getMinPhysR( 1 );
    Vector<T,3> posOutflow = geometry.getStatistics().getMaxPhysR( 1 );

    SuperPlaneIntegralFluxVelocity3D<T> vFluxIn( sLattice, converter,
        geometry, posInflow, normal, materials, mode );
    SuperPlaneIntegralFluxPressure3D<T> pFluxIn( sLattice, converter,
        geometry, posInflow, normal, materials, mode );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOut( sLattice, converter,
        geometry, posOutflow, normal, materials, mode );
    SuperPlaneIntegralFluxPressure3D<T> pFluxOut( sLattice, converter,
        geometry, posOutflow, normal, materials, mode );

    vFluxIn.print( "Inflow" );
    pFluxIn.print( "Inflow" );
    vFluxOut.print( "Outflow" );
    pFluxOut.print( "Outflow" );

    // Error norms
    AnalyticalFfromSuperF3D<T> intpolatePressure( pressure, true );

    T point1[3] = {0, radius, radius};
    T point2[3] = {0, radius, radius};

    point1[0] = length*0.5 - length*0.01;
    point2[0] = length*0.5 + length*0.01;

    T p1, p2;
    intpolatePressure( &p1,point1 );
    intpolatePressure( &p2,point2 );

    clout << "pressure1=" << p1;
    clout << "; pressure2=" << p2;

    T pressureDrop = p1-p2;
    clout << "; pressureDrop=" << pressureDrop << std::endl;

    error(geometry, sLattice, converter, uSol);

    // Gnuplot
    Gnuplot<T> gplot( "velocityProfile" );
    T uAnalytical[3] = {};
    T uNumerical[3] = {};
    AnalyticalFfromSuperF3D<T> intpolateVelocity( velocity, true );
    for (int i=0; i<101; i++) {
      T yInput = diameter*T(i)/T(100);
      T input[3] = {length*T(0.5), yInput, radius};
      uSol(uAnalytical, input);
      intpolateVelocity(uNumerical, input);
      gplot.setData( yInput, {uAnalytical[0], uNumerical[0]}, {"analytical","numerical"} );
    }

    // Create PNG file
    gplot.writePNG();
  }
}

int main( int argc, char* argv[] )
{

  // === 1st Step: Initialization ===
  initialize( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  if (argc > 1) {
    if (argv[1][0]=='-'&&argv[1][1]=='h') {
      OstreamManager clout( std::cout,"help" );
      clout<<"Usage: program [Resolution] [Permeability]" <<std::endl;
      clout<<"Default: Resolution=21, Permeability=1e-2" <<std::endl;
      return 0;
    }
  }

  if (argc > 1) {
    N = atoi(argv[1]);
    if (N < 1) {
      std::cerr << "Fluid domain is too small" << std::endl;
      return 1;
    }
  }

  if (argc > 2) {
    Kin = atof(argv[2]);
    if (Kin < 0) {
      std::cerr << "Permeabilty must be greater than 0" << std::endl;
      return 2;
    }
  }

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {N},                  // resolution: number of voxels per charPhysL
    (T)   tau,                // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   diameter,           // charPhysLength: reference length of simulation geometry
    (T)   physU,              // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   diameter*physU/Re,  // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   physRho             // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  //converter.print();
  // Writes the converter log in a file
  converter.write("porousPoiseuille3d");


  // === 2nd Step: Prepare Geometry ===

  Vector<T, 3> center0(0, radius, radius);
  Vector<T, 3> center1(length + 0.5 * converter.getPhysDeltaX(), radius, radius);
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  IndicatorLayer3D<T> extendedDomain(pipe, converter.getPhysDeltaX());

  // Instantiation of a cuboidDecomposition with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = 2*singleton::mpi().getSize();
#else // ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = 6;
#endif // ifdef PARALLEL_MODE_MPI
  CuboidDecomposition3D<T> cuboidDecomposition(extendedDomain, converter.getPhysDeltaX(), noOfCuboids);

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  // Instantiation of a geometry
  SuperGeometry<T,3> geometry(cuboidDecomposition, loadBalancer);

  prepareGeometry(converter, geometry);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( converter, geometry );

  //prepareLattice and setBoundaryConditions
  prepareLattice(sLattice, converter, geometry);

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), geometry.getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( converter.getLatticeTime( physInterval ), residuum );
  timer.start();

  for ( std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << std::endl;
      getResults( sLattice, converter, iT, geometry, timer, converge.hasConverged() );

      break;
    }

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    // in this application no boundary conditions have to be adjusted

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, geometry, timer, converge.hasConverged()  );
    converge.takeValue( sLattice.getStatistics().getAverageEnergy(), true );
  }

  timer.stop();
  timer.printSummary();
}
