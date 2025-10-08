/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2018 Marc Haußmann, Mathias J. Krause, Jonathan Jeppener-Haltenhoff
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

/* poiseuille3d.cpp:
 * This example examines a 3D Poseuille flow
 * It illustrates the usage of different flow setups, boundary conditions
 * and computation of error norms.
 * The following simulation options can be combined freely:
 * - if the compiler flag ENABLE_MRT is set, mrt collision operators are used
 * - forced/ nonForced flow
 * - different boundary conditions
 * - simulation only or eoc convergence analysis
 */

// the main code of the simulation is in poiseuille3d.h as it is also used by the
// example ../../pdeSolverEoc/poiseuille3dEoc

// set flag in order to use mrt collision operators instead of bgk
//#define ENABLE_MRT

/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2018 Marc Haußmann, Mathias J. Krause, Jonathan Jeppener-Haltenhoff
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


using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::names;

using DESCRIPTOR = D3Q19<FORCE>;

//typedef enum {forced, nonForced} FlowType;//TODO porpoer enums
//typedef enum {bounceBack, local, interpolated, bouzidi, freeSlip, partialSlip} BoundaryType;
  
enum FlowType : int {
  FORCED = 0,
  NON_FORCED = 1
};
enum BoundaryType : int {
  BOUNCE_BACK = 0,
  LOCAL = 1,
  INTERPOLATED = 2,
  BOUZIDI= 3 ,
  FREE_SLIP = 4,
  PARTIAL_SLIP = 5
};

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q19<FORCE>>
>;

namespace olb::parameters {  
  struct FLOW_TYPE : public TYPED_FIELD_BASE<FlowType,1> { };
  struct BOUNDARY_TYPE : public TYPED_FIELD_BASE<BoundaryType,1> { };
  
  struct RE : public descriptors::FIELD_BASE<1> { }; // reynolds number
  //struct TAU : public descriptors::FIELD_BASE<1> { }; //lattice relaxation time
  struct DIAMETER : public descriptors::FIELD_BASE<1> { }; // diameter of the pipe
  struct LENGTH : public descriptors::FIELD_BASE<1> { }; // length of the pipe
  struct INTERVAL : public descriptors::FIELD_BASE<1> { }; // interval for the convergence check in s
  struct RESIDUUM : public descriptors::FIELD_BASE<1> { }; //residuum for the convergence check
  //struct TUNE : public descriptors::FIELD_BASE<1> { }; // for partialSlip only: 0->bounceBack, 1->freeSlip
}

/// @brief Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const T radius = parameters.get<parameters::DIAMETER>()/2;
  const T length = parameters.get<parameters::LENGTH>();
  Vector<T, 3> center0(0, radius, radius);
  Vector<T, 3> center1(length + 0.5 * parameters.get<parameters::PHYS_DELTA_X>(), radius, radius);
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  IndicatorLayer3D<T> extendedDomain(pipe, parameters.get<parameters::PHYS_DELTA_X>());
  Mesh<T, MyCase::d> mesh(extendedDomain, parameters.get<parameters::PHYS_DELTA_X>(), singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

// Stores geometry information in form of material numbers
void prepareGeometry(MyCase& myCase)
{
  using T = MyCase::value_t;
  OstreamManager clout(std::cout, "prepareGeometry");

  clout << "Prepare Geometry ..." << std::endl;
  
  auto&parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  
  const T radius = parameters.get<parameters::DIAMETER>()/2.0;
  const T length = parameters.get<parameters::LENGTH>();
  FlowType flowType = parameters.get<parameters::FLOW_TYPE>();

  const T dx = parameters.get<parameters::PHYS_DELTA_X>();

  Vector<T, 3> center0(-dx * 0.2, radius, radius);
  Vector<T, 3> center1(length, radius, radius);
  if (flowType == FORCED) {
    center0[0] -= 3.*dx;
    center1[0] += 3.*dx;
  }
  IndicatorCylinder3D<T> pipe(center0, center1, radius);

  geometry.rename(0, 2);

  geometry.rename(2, 1, pipe);

  if (flowType == NON_FORCED) {
    geometry.clean();
    Vector<T, 3> origin(0, radius, radius);
    Vector<T, 3> extend = origin;

    // Set material number for inflow
    origin[0] = -dx * 2;
    extend[0] = dx * 2;
    IndicatorCylinder3D<T> inflow(origin, extend, radius);
    geometry.rename(2, 3, 1, inflow);

    // Set material number for outflow
    origin[0] = length - 2 * dx;
    extend[0] = length + 2 * dx;
    IndicatorCylinder3D<T> outflow(extend, origin, radius);
    geometry.rename(2, 4, 1, outflow);
  }

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  // Removes all not needed boundary voxels inside the surface
  geometry.innerClean();
  geometry.checkForErrors();

  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(MyCase& myCase)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;
  using T = MyCase::value_t;
  using BulkDynamics       = BGKdynamics<T,DESCRIPTOR>;
  using ForcedBulkDynamics = ForcedBGKdynamics<T,DESCRIPTOR>;
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});
  
  const T radius = parameters.get<parameters::DIAMETER>()/2.0;
  const T length = parameters.get<parameters::LENGTH>();
  
  const T diameter = parameters.get<parameters::DIAMETER>();
  const T physU = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();
  
  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    int {parameters.get<parameters::RESOLUTION>()},                  // resolution: number of voxels per charPhysL
    (T)   parameters.get<TAU>(), // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   diameter,           // charPhysLength: reference length of simulation geometry
    (T)   physU,// charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   diameter*physU/parameters.get<parameters::RE>(),  // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   parameters.get<parameters::PHYS_CHAR_DENSITY>() // physDensity: physical density in __kg / m^3__
  );

  

  const T omega = lattice.getUnitConverter().getLatticeRelaxationFrequency();
  FlowType flowType = parameters.get<parameters::FLOW_TYPE>();
  BoundaryType boundaryType = parameters.get<parameters::BOUNDARY_TYPE>();

  if (flowType == NON_FORCED) {
    lattice.defineDynamics<BulkDynamics>(myCase.getGeometry(), 1);
  } else {
    lattice.defineDynamics<ForcedBulkDynamics>(myCase.getGeometry(), 1);
  }

  Vector<T, 3> center0(0, radius, radius);
  Vector<T, 3> center1(length, radius, radius);

  std::vector<T> origin = { length, radius, radius};
  std::vector<T> axis = { 1, 0, 0 };

  CirclePoiseuille3D<T> poiseuilleU(origin, axis, lattice.getUnitConverter().getCharLatticeVelocity(), radius);

  if (boundaryType == BOUNCE_BACK) {
    boundary::set<boundary::BounceBack>(lattice, myCase.getGeometry(), 2);
  }
  else if (boundaryType == FREE_SLIP) {
    boundary::set<boundary::FullSlip>(lattice, myCase.getGeometry(), 2);
  }
  else if (boundaryType == PARTIAL_SLIP) {
    boundary::set<boundary::PartialSlip>(lattice, myCase.getGeometry(), 2);
    lattice.template setParameter<descriptors::TUNER>(parameters.get<TUNER>());
  }
  else if (boundaryType == BOUZIDI) {
    center0[0] -= 0.5*physDeltaX;
    center1[0] += 0.5*physDeltaX;
    if (flowType == FORCED) {
      center0[0] -= 3.*physDeltaX;
      center1[0] += 3.*physDeltaX;
    }
    IndicatorCylinder3D<T> pipe(center0, center1, radius);
    setBouzidiBoundary<T, DESCRIPTOR, BouzidiPostProcessor>(lattice, myCase.getGeometry(), 2, pipe);
  }
  else {
    if (boundaryType == LOCAL) {
      boundary::set<boundary::LocalVelocity>(lattice, myCase.getGeometry(), 2);
    }
    else {
      boundary::set<boundary::InterpolatedVelocity>(lattice, myCase.getGeometry(), 2);
    }
  }

  if (flowType == NON_FORCED) {
    if (boundaryType == BOUZIDI) {
      IndicatorCylinder3D<T> pipe(center0, center1, radius);
      setBouzidiBoundary<T, DESCRIPTOR, BouzidiVelocityPostProcessor>(lattice, myCase.getGeometry(), 3, pipe);
    }
    else {
      // Material=3 -->bulk dynamics
      if (boundaryType == LOCAL) {
        boundary::set<boundary::LocalVelocity>(lattice, myCase.getGeometry(), 3);
      }
      else {
        boundary::set<boundary::InterpolatedVelocity>(lattice, myCase.getGeometry(), 3);
      }
    }
    // Material=4 -->bulk dynamics
    if (boundaryType == LOCAL) {
      boundary::set<boundary::LocalPressure>(lattice, myCase.getGeometry(), 4);
    }
    else {
      boundary::set<boundary::LocalPressure>(lattice, myCase.getGeometry(), 4);
    }
  }

  if (flowType == FORCED) {
    // Initial conditions
    T D = lattice.getUnitConverter().getLatticeLength(diameter);

    std::vector<T> poiseuilleForce(3, T());
    poiseuilleForce[0] = 4. * lattice.getUnitConverter().getLatticeViscosity() * lattice.getUnitConverter().getCharLatticeVelocity() / (D * D / 4. );
    AnalyticalConst3D<T,T> force( poiseuilleForce );

    // Initialize force
    lattice.defineField<FORCE>(myCase.getGeometry(), 1, force);
    lattice.defineField<FORCE>(myCase.getGeometry(), 2, force );


    AnalyticalConst3D<T, T> rhoF(1);

    lattice.defineRhoU(myCase.getGeometry(), 1, rhoF, poiseuilleU);
    lattice.iniEquilibrium(myCase.getGeometry(), 1, rhoF, poiseuilleU);
    lattice.defineRhoU(myCase.getGeometry(), 2, rhoF, poiseuilleU);
    lattice.iniEquilibrium(myCase.getGeometry(), 2, rhoF, poiseuilleU);
  }
  else {
    // Initial conditions
    T p0 = 4. * lattice.getUnitConverter().getPhysViscosity() * physU * length / (radius * radius);

    p0 = lattice.getUnitConverter().getLatticePressure(p0);
    AnalyticalLinear3D<T, T> rho(-p0 / length * descriptors::invCs2<T,DESCRIPTOR>(), 0, 0, p0 * descriptors::invCs2<T,DESCRIPTOR>() + 1);

    std::vector<T> velocity(3, T());
    AnalyticalConst3D<T, T> uF(velocity);

    // Initialize all values of distribution functions to their local equilibrium
    lattice.defineRhoU(myCase.getGeometry(), 0, rho, uF);
    lattice.iniEquilibrium(myCase.getGeometry(), 0, rho, uF);
    lattice.defineRhoU(myCase.getGeometry(), 1, rho, poiseuilleU);
    lattice.iniEquilibrium(myCase.getGeometry(), 1, rho, poiseuilleU);
    lattice.defineRhoU(myCase.getGeometry(), 2, rho, poiseuilleU);
    lattice.iniEquilibrium(myCase.getGeometry(), 2, rho, poiseuilleU);
    lattice.defineRhoU(myCase.getGeometry(), 3, rho, poiseuilleU);
    lattice.iniEquilibrium(myCase.getGeometry(), 3, rho, poiseuilleU);
    if(boundaryType == BOUZIDI) {
      setBouzidiVelocity(lattice, myCase.getGeometry(), 3, poiseuilleU);
    }
    lattice.defineRhoU(myCase.getGeometry(), 4, rho, poiseuilleU);
    lattice.iniEquilibrium(myCase.getGeometry(), 4, rho, poiseuilleU);
  }

  // Set relaxation time to omega in all dynamics
  lattice.setParameter<descriptors::OMEGA>(omega);
  // Make the lattice ready for simulation
  lattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Compute error norms
void error_calc(MyCase& myCase)
{
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& superGeometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  OstreamManager clout( std::cout,"error" );
  const UnitConverter<T,DESCRIPTOR>& converter = sLattice.getUnitConverter();
  sLattice.setProcessingContext(ProcessingContext::Evaluation);
  
  const T radius = parameters.get<parameters::DIAMETER>()/2.0;
  const T length = parameters.get<parameters::LENGTH>();
  
  // set up size-increased indicator and instantiate wall shear stress functor (wss)
  Vector<T, 3> center0Extended(-converter.getPhysDeltaX() * 0.2, radius, radius);
  Vector<T, 3> center1Extended(length, radius, radius);
  if (parameters.get<parameters::FLOW_TYPE>() == FORCED) {
    center0Extended[0] -= 4.*converter.getPhysDeltaX();
    center1Extended[0] += 4.*converter.getPhysDeltaX();
  }
  IndicatorCylinder3D<T> pipeExtended(center0Extended, center1Extended, radius);
  IndicatorLayer3D<T> indicatorExtended (pipeExtended, 0.9*converter.getPhysDeltaX()*parameters.get<parameters::RESOLUTION>()/11.);
  SuperLatticePhysWallShearStress3D<T,DESCRIPTOR> wss(sLattice, superGeometry, 2, converter, indicatorExtended);

  int tmp[]= { };
  T result[2]= { };

  // velocity error
  const T maxVelocity = converter.getCharPhysVelocity();
  std::vector<T> axisPoint = {length, radius, radius};
  std::vector<T> axisDirection = { 1, 0, 0 };
  CirclePoiseuille3D<T> uSol(axisPoint, axisDirection, maxVelocity, radius);
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> u( sLattice,converter );
  auto indicatorF = superGeometry.getMaterialIndicator(1);

  SuperAbsoluteErrorL1Norm3D<T> absVelocityErrorNormL1(u, uSol, indicatorF);
  absVelocityErrorNormL1(result, tmp);
  clout << "velocity-L1-error(abs)=" << result[0];
  T velocityL1AbsError = result[0];
  SuperRelativeErrorL1Norm3D<T> relVelocityErrorNormL1(u, uSol, indicatorF);
  relVelocityErrorNormL1(result, tmp);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm3D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, tmp);
  clout << "velocity-L2-error(abs)=" << result[0];
  T velocityL2AbsError = result[0];
  SuperRelativeErrorL2Norm3D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, tmp);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm3D<T> absVelocityErrorNormLinf(u, uSol, indicatorF);
  absVelocityErrorNormLinf(result, tmp);
  clout << "velocity-Linf-error(abs)=" << result[0];
  T velocityLinfAbsError = result[0];
  SuperRelativeErrorLinfNorm3D<T> relVelocityErrorNormLinf(u, uSol, indicatorF);
  relVelocityErrorNormLinf(result, tmp);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;

  // strainRate error
  CirclePoiseuilleStrainRate3D<T, DESCRIPTOR> sSol( converter, radius );
  SuperLatticePhysStrainRate3D<T,DESCRIPTOR> s( sLattice,converter );

  SuperAbsoluteErrorL1Norm3D<T> absStrainRateErrorNormL1(s, sSol, indicatorF);
  absStrainRateErrorNormL1(result, tmp);
  clout << "strainRate-L1-error(abs)=" << result[0];
  T strainRateL1AbsError = result[0];
  SuperRelativeErrorL1Norm3D<T> relStrainRateErrorNormL1(s, sSol, indicatorF);
  relStrainRateErrorNormL1(result, tmp);
  clout << "; strainRate-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm3D<T> absStrainRateErrorNormL2(s, sSol, indicatorF);
  absStrainRateErrorNormL2(result, tmp);
  clout << "strainRate-L2-error(abs)=" << result[0];
  T strainRateL2AbsError = result[0];
  SuperRelativeErrorL2Norm3D<T> relStrainRateErrorNormL2(s, sSol, indicatorF);
  relStrainRateErrorNormL2(result, tmp);
  clout << "; strainRate-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm3D<T> absStrainRateErrorNormLinf(s, sSol, indicatorF);
  absStrainRateErrorNormLinf(result, tmp);
  clout << "strainRate-Linf-error(abs)=" << result[0];
  T strainRateLinfAbsError = result[0];
  SuperRelativeErrorLinfNorm3D<T> relStrainRateErrorNormLinf(s, sSol, indicatorF);
  relStrainRateErrorNormLinf(result, tmp);
  clout << "; strainRate-Linf-error(rel)=" << result[0] << std::endl;

  // wallShearStress error
  AnalyticalConst3D<T,T> wssSol(4. * converter.getPhysViscosity() * converter.getPhysDensity() * maxVelocity / parameters.get<parameters::DIAMETER>());
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> wssSolLattice (wssSol, sLattice);

  auto indicatorB = superGeometry.getMaterialIndicator(2);

  SuperAbsoluteErrorL1Norm3D<T> absWallShearStressErrorNormL1(wss, wssSol, indicatorB);
  absWallShearStressErrorNormL1(result, tmp);
  clout << "wss-L1-error(abs)=" << result[0];
  T wssL1AbsError = result[0];
  SuperRelativeErrorL1Norm3D<T> relWallShearStressErrorNormL1(wss, wssSol, indicatorB);
  relWallShearStressErrorNormL1(result, tmp);
  clout << "; wss-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm3D<T> absWallShearStressErrorNormL2(wss, wssSol, indicatorB);
  absWallShearStressErrorNormL2(result, tmp);
  clout << "wss-L2-error(abs)=" << result[0];
  T wssL2AbsError = result[0];
  SuperRelativeErrorL2Norm3D<T> relWallShearStressErrorNormL2(wss, wssSol, indicatorB);
  relWallShearStressErrorNormL2(result, tmp);
  clout << "; wss-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm3D<T> absWallShearStressErrorNormLinf(wss, wssSol, indicatorB);
  absWallShearStressErrorNormLinf(result, tmp);
  clout << "wss-Linf-error(abs)=" << result[0];
  T wssLinfAbsError = result[0];
  SuperRelativeErrorLinfNorm3D<T> relWallShearStressErrorNormLinf(wss, wssSol, indicatorB);
  relWallShearStressErrorNormLinf(result, tmp);
  clout << "; wss-Linf-error(rel)=" << result[0] << std::endl;

  if (parameters.get<parameters::FLOW_TYPE>() == NON_FORCED) {
    // pressure error
    T p0 = 4. * converter.getPhysViscosity() * maxVelocity * length / (radius * radius);
    AnalyticalLinear3D<T, T> pressureSol(-p0 / length, 0, 0, p0);
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);

    SuperAbsoluteErrorL1Norm3D<T> absPressureErrorNormL1(pressure, pressureSol, indicatorF);
    absPressureErrorNormL1(result, tmp);
    clout << "pressure-L1-error(abs)=" << result[0];
    T pressureL1AbsError = result[0];
    SuperRelativeErrorL1Norm3D<T> relPressureErrorNormL1(pressure, pressureSol, indicatorF);
    relPressureErrorNormL1(result, tmp);
    clout << "; pressure-L1-error(rel)=" << result[0] << std::endl;

    SuperAbsoluteErrorL2Norm3D<T> absPressureErrorNormL2(pressure, pressureSol, indicatorF);
    absPressureErrorNormL2(result, tmp);
    clout << "pressure-L2-error(abs)=" << result[0];
    T pressureL2AbsError = result[0];
    SuperRelativeErrorL2Norm3D<T> relPressureErrorNormL2(pressure, pressureSol, indicatorF);
    relPressureErrorNormL2(result, tmp);
    clout << "; pressure-L2-error(rel)=" << result[0] << std::endl;

    SuperAbsoluteErrorLinfNorm3D<T> absPressureErrorNormLinf(pressure, pressureSol, indicatorF);
    absPressureErrorNormLinf(result, tmp);
    clout << "pressure-Linf-error(abs)=" << result[0];
    T pressureLinfAbsError = result[0];
    SuperRelativeErrorLinfNorm3D<T> relPressureErrorNormLinf(pressure, pressureSol, indicatorF);
    relPressureErrorNormLinf(result, tmp);
    clout << "; pressure-Linf-error(rel)=" << result[0] << std::endl;
  }
}

// Output to console and files
void getResults( std::size_t iT,
                 util::Timer<MyCase::value_t>& timer, bool hasConverged, MyCase& myCase)
{
  using T = MyCase::value_t;
  //Initialize Gnuplot
  static Gnuplot<T> gplot("centerVelocity");

  auto& parameters = myCase.getParameters();
  auto& superGeometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  OstreamManager clout( std::cout,"getResults" );
  const UnitConverter<T,DESCRIPTOR>& converter = sLattice.getUnitConverter();
  const bool lastTimeStep = ( hasConverged || (iT + 1 == converter.getLatticeTime( parameters.get<parameters::MAX_PHYS_T>() )) );
  const bool noslipBoundary = ((parameters.get<parameters::BOUNDARY_TYPE>() != FREE_SLIP) && (parameters.get<parameters::BOUNDARY_TYPE>() != PARTIAL_SLIP));
  const int statIter = converter.getLatticeTime( parameters.get<parameters::MAX_PHYS_T>()/20. );
  
  const T radius = parameters.get<parameters::DIAMETER>()/2;
  const T length = parameters.get<parameters::LENGTH>();
  
  // set up size-increased indicator and instantiate wall shear stress functor (wss)
  Vector<T, 3> center0Extended(-converter.getPhysDeltaX() * 0.2, radius, radius);
  Vector<T, 3> center1Extended(length, radius, radius);
  if (parameters.get<parameters::FLOW_TYPE>() == FORCED) {
    center0Extended[0] -= 4.*converter.getPhysDeltaX();
    center1Extended[0] += 4.*converter.getPhysDeltaX();
  }
  IndicatorCylinder3D<T> pipeExtended(center0Extended, center1Extended, radius);
  IndicatorLayer3D<T> indicatorExtended (pipeExtended, 0.9*converter.getPhysDeltaX()*parameters.get<parameters::RESOLUTION>()/11.);
  SuperLatticePhysWallShearStress3D<T,DESCRIPTOR> wss(sLattice, superGeometry, 2, converter, indicatorExtended);
  
  SuperVTMwriter3D<T> vtmWriter( "poiseuille3d" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );
  vtmWriter.addFunctor( wss );
  const T maxVelocity = converter.getCharPhysVelocity();
  std::vector<T> axisPoint = {length, radius, radius};
  std::vector<T> axisDirection = { 1, 0, 0 };
  CirclePoiseuille3D<T> uSol(axisPoint, axisDirection, maxVelocity, radius);
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> analyticalVelocityLattice(uSol, sLattice);
  analyticalVelocityLattice.getName() = "analytical solution";
  vtmWriter.addFunctor(analyticalVelocityLattice);

  const int vtmIter  = converter.getLatticeTime( parameters.get<parameters::MAX_PHYS_T>()/20. );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    SuperLatticeDiscreteNormal3D<T, DESCRIPTOR> discreteNormal( sLattice, superGeometry, superGeometry.getMaterialIndicator({2, 3}) );
    SuperLatticeDiscreteNormalType3D<T, DESCRIPTOR> discreteNormalType( sLattice, superGeometry, superGeometry.getMaterialIndicator({2, 3, 4, 5}) );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.write( discreteNormal );
    vtmWriter.write( discreteNormalType );

    vtmWriter.createMasterFile();
  }

  // Writes the vtm files and profile text file
  if ( iT%vtmIter==0 || lastTimeStep ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    vtmWriter.write( iT );

    SuperEuklidNorm3D<T> normVel( velocity );
    BlockReduction3D2D<T> planeReduction( normVel, Vector<T,3>({0,0,1}), 600, BlockDataSyncMode::ReduceOnly );
    // write output as JPEG
    heatmap::write(planeReduction, iT);

  }
  

  // Writes output on the console
  if ( iT%statIter==0 || lastTimeStep ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Error norms
    if (noslipBoundary) {
      if ( lastTimeStep ) {
        error_calc(myCase);
      }
    }
  }

  // Gnuplot output
  if ((noslipBoundary) && (lastTimeStep)) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    // plot velocity magnitude over line through the center of the simulation domain
    const T maxVelocity = converter.getCharPhysVelocity();
    T D = converter.getLatticeLength( parameters.get<parameters::DIAMETER>() );
    T dx = 1. / T(converter.getResolution());
    T point[3] { };
    point[0] = length/2.;
    point[2] = ( T )radius;
    std::vector<T> axisPoint {length, radius, radius};
    std::vector<T> axisDirection { 1, 0, 0 };
    CirclePoiseuille3D<T> uSol(axisPoint, axisDirection, maxVelocity, radius);
    T analytical[3] { };
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
    AnalyticalFfromSuperF3D<T> intpolateVelocity( velocity, true, 1 );
    T numerical[3] { };
    for ( int iY=0; iY<=D; ++iY ) {
      point[1] = ( T )converter.getPhysLength(iY);
      uSol( analytical,point );
      intpolateVelocity( numerical,point );
      gplot.setData( iY*dx, {analytical[0],numerical[0]}, {"analytical","numerical"} );
    }
    // Create PNG file
    gplot.writePNG();
  }
}

void simulate(MyCase& myCase) {
  OstreamManager clout( std::cout,"simulate" );
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  auto& converter = sLattice.getUnitConverter();
  
  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(
    parameters.get<parameters::MAX_PHYS_T>());
  
  

  const T radius = parameters.get<parameters::DIAMETER>()/2.0;
  const T length = parameters.get<parameters::LENGTH>();

  // set up size-increased indicator and instantiate wall shear stress functor (wss)
  Vector<T, 3> center0Extended(-converter.getPhysDeltaX() * 0.2, radius, radius);
  Vector<T, 3> center1Extended(length, radius, radius);
  if (parameters.get<parameters::FLOW_TYPE>() == FORCED) {
    center0Extended[0] -= 4.*converter.getPhysDeltaX();
    center1Extended[0] += 4.*converter.getPhysDeltaX();
  }
  IndicatorCylinder3D<T> pipeExtended(center0Extended, center1Extended, radius);
  IndicatorLayer3D<T> indicatorExtended (pipeExtended, 0.9*converter.getPhysDeltaX()*parameters.get<parameters::RESOLUTION>()/11.);
  SuperLatticePhysWallShearStress3D<T,DESCRIPTOR> wss(sLattice, superGeometry, 2, converter, indicatorExtended);
  
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer(iTmax, superGeometry.getStatistics().getNvoxel());
    util::ValueTracer<T> converge( converter.getLatticeTime( parameters.get<parameters::INTERVAL>() ), parameters.get<parameters::RESIDUUM>() );
  timer.start();
  
  for ( std::size_t iT = 0; iT < converter.getLatticeTime( parameters.get<parameters::MAX_PHYS_T>() ); ++iT ) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << std::endl;
      getResults( iT, timer,
        converge.hasConverged(), myCase);
      break;
    }

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    // in this application no boundary conditions have to be adjusted

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( iT, timer,
      converge.hasConverged(), myCase );
    converge.takeValue( sLattice.getStatistics().getAverageEnergy(), true );
  }

  timer.stop();
  timer.printSummary();
  
}


int main( int argc, char* argv[] )
{
  //int N = 21;
  //bool eoc = false;
  
  // === 1st Step: Initialization ===
  initialize( &argc, &argv );
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION         >(           21);//TODO check N?
    //myCaseParameters.set<PHYS_DELTA_X       >(       0.1/21);//TODO hier schon defininieren?
    //myCaseParameters.set<PHYS_DELTA_T       >(0.00078125);
    myCaseParameters.set<PHYS_CHAR_VELOCITY >(          1.0);
    //myCaseParameters.set<PHYS_CHAR_VISCOSITY>(     0.001);
    myCaseParameters.set<PHYS_CHAR_DENSITY  >(          1.0);
    //myCaseParameters.set<DOMAIN_EXTENT      >({1.0, 1.0});
    myCaseParameters.set<MAX_PHYS_T         >(          20.);
    myCaseParameters.set<LENGTH             >(           2.);
    myCaseParameters.set<DIAMETER           >(           1.);
    myCaseParameters.set<RE                 >(          10.);
    myCaseParameters.set<TAU                >(          0.8);
    myCaseParameters.set<RESIDUUM           >(         1e-5);
    myCaseParameters.set<TUNER              >(           0.);
    myCaseParameters.set<FLOW_TYPE          >(    NON_FORCED);
    myCaseParameters.set<BOUNDARY           >(  INTERPOLATED);
  }
  myCaseParameters.fromCLI(argc, argv);
  {
    using namespace olb::parameters;
    myCaseParameters.set<INTERVAL           >(myCaseParameters.get<MAX_PHYS_T>()*0.0125); //TODO CHECK
    myCaseParameters.set<PHYS_DELTA_X       >(myCaseParameters.get<DIAMETER>()/myCaseParameters.get<RESOLUTION>()*0.0125);//TODO CHECk defininieren?
  }
  
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  
  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);
  
  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);
  
  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);
  
  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase);
  
  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  //setInitialValues(myCase); TODO sections in prepareLattice?
  
  /// === Step 8: Simulate ===
  simulate(myCase);

  // Instantiation of a superGeometry
  //const int overlap = (flowType == FORCED) ? 2 : 3;
  //SuperGeometry<T,3> superGeometry(cuboidDecomposition, loadBalancer, overlap);

}
