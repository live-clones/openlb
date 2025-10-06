/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Florian Kaiser, Stephan Simonis
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
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierCauchy,
  Lattice<double,
          descriptors::D2Q8<descriptors::tag::MRT,
                            descriptors::FORCE,
                            descriptors::MOMENT_VECTOR,
                            descriptors::DISP_SOLID,
                            descriptors::SIGMA_SOLID>>
>;

namespace olb::parameters {
  struct KAPPA  : public descriptors::FIELD_BASE<1> { };
}

template <typename T, typename DESCRIPTOR>
class ForceField2D : public AnalyticalF2D<T, T> {

protected:
  MyCase _myCase;

public:
  ForceField2D(MyCase& myCase) : AnalyticalF2D<T, T>(2),
  _myCase(myCase) {};

  bool operator()(T output[], const T input[]) override
  {
    auto& lattice = _myCase.getLattice(NavierCauchy{});
    T pi = std::numbers::pi_v<double>;

    const auto& converter = lattice.getUnitConverter();
    T dx = converter.getPhysDeltaX();
    T dt = converter.getPhysDeltaT();
    T kappa = converter.getDampingFactor();
    T lambda = converter.getLatticeLambda();
    T mu = converter.getLatticeShearModulus();
    T charU = converter.getCharPhysDisplacement();

    T x = input[0];
    T y = input[1];

    T x_sx = -.3;
    T x_sy =  .4;
    T y_sx =  .1;
    T y_sy = -.7;

    T omega_11 = 1. / (mu / descriptors::invCs2<T,DESCRIPTOR>() + 0.5);
    T omega_d  = 1. / (2 * mu / (1 - descriptors::invCs2<T,DESCRIPTOR>()) + 0.5);
    T omega_s  = 1. / (2 * (mu + lambda) / (1 + descriptors::invCs2<T,DESCRIPTOR>()) + 0.5);

    T tau_11 = 1. / omega_11 - 0.5;
    T tau_s = 1. / omega_d - 0.5;
    T tau_p = 1. / omega_s - 0.5;
    T tau_12 = (1. / tau_11 + 4. * tau_11) / 8.;

    // These are correct
    T d1 = -1. / 4. - 1. / 2. * descriptors::invCs2<T,DESCRIPTOR>() * tau_12 * tau_s + pow(tau_s, 2) / 2.
         -  1. / 2. * descriptors::invCs2<T,DESCRIPTOR>() * pow(tau_s, 2) + 1. / 2. * descriptors::invCs2<T,DESCRIPTOR>() * tau_12 * tau_p
         +  pow(tau_p, 2) / 2. + 1. / 2. * descriptors::invCs2<T,DESCRIPTOR>() * pow(tau_p, 2);
    T d2 = -descriptors::invCs2<T,DESCRIPTOR>() / 2. + descriptors::invCs2<T,DESCRIPTOR>() * pow(tau_11, 2) + descriptors::invCs2<T,DESCRIPTOR>() * tau_11 * tau_12
         + 1. / 2. * descriptors::invCs2<T,DESCRIPTOR>() * tau_12 * tau_s - pow(tau_s, 2) / 2.
         + 1. / 2. * descriptors::invCs2<T,DESCRIPTOR>() * pow(tau_s, 2) + 1. / 2. * descriptors::invCs2<T,DESCRIPTOR>() * tau_12 * tau_p
         + pow(tau_p, 2) / 2. + 1. / 2. * descriptors::invCs2<T,DESCRIPTOR>() * pow(tau_p, 2);
    T d3 = -descriptors::invCs2<T,DESCRIPTOR>() / 4.
         +  descriptors::invCs2<T,DESCRIPTOR>() * pow(tau_11, 2)
         +  descriptors::invCs2<T,DESCRIPTOR>() * tau_11 * tau_12;

    T mu_phys     =     mu * (dx * dx * kappa) / dt;
    T lambda_phys = lambda * (dx * dx * kappa) / dt;

    output[0] = bx(0, x, y, x_sx, x_sy, y_sx, y_sy, dx, mu_phys, lambda_phys, charU, pi) * dt / kappa
              + (d1 * d2bx_dx2(0, x, y, x_sx, x_sy, y_sx, y_sy, dx, mu_phys, lambda_phys, charU, pi)
              + d2 * d2by_dxdy(0, x, y, x_sx, x_sy, y_sx, y_sy, dx, mu_phys, lambda_phys, charU, pi)
              + d3 * d2bx_dy2(0, x, y, x_sx, x_sy, y_sx, y_sy, dx, mu_phys, lambda_phys, charU, pi)) * dx * dx * dt / kappa;

    output[1] = by(0, x, y, x_sx, x_sy, y_sx, y_sy, dx, mu_phys, lambda_phys, charU, pi) * dt / kappa
              + (d1 * d2by_dy2(0, x, y, x_sx, x_sy, y_sx, y_sy, dx, mu_phys, lambda_phys, charU, pi)
              + d2 * d2bx_dxdy(0, x, y, x_sx, x_sy, y_sx, y_sy, dx, mu_phys, lambda_phys, charU, pi)
              + d3 * d2by_dx2(0, x, y, x_sx, x_sy, y_sx, y_sy, dx, mu_phys, lambda_phys, charU, pi)) * dx * dx * dt / kappa;

    return true;
  };

  // Output 0
  T bx(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T dx, T mu, T lambda,
       T charU, T pi)
  {
    return 0.0036*mu*pi*pi*sin(2*pi*(y + y_sx))*cos(2*pi*(x + x_sx)) - 0.0028*pi*pi*(lambda + mu)*sin(2*pi*(x + x_sy))*sin(2*pi*(y + y_sy)) + 0.0036*pi*pi*(lambda + 2*mu)*sin(2*pi*(y + y_sx))*cos(2*pi*(x + x_sx));
  }

  T d2bx_dx2(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T dx, T mu,
             T lambda, T charU, T pi)
  {
    return pi*pi*pi*pi*(-0.0144*mu*sin(2*pi*(y + y_sx))*cos(2*pi*(x + x_sx)) + 0.0112*(lambda + mu)*sin(2*pi*(x + x_sy))*sin(2*pi*(y + y_sy)) - 0.0144*(lambda + 2*mu)*sin(2*pi*(y + y_sx))*cos(2*pi*(x + x_sx)));
  }

  T d2by_dxdy(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T dx, T mu,
              T lambda, T charU, T pi)
  {
    return pi*pi*pi*pi*(0.0112*mu*sin(2*pi*(x + x_sy))*sin(2*pi*(y + y_sy)) - 0.0144*(lambda + mu)*sin(2*pi*(y + y_sx))*cos(2*pi*(x + x_sx)) + 0.0112*(lambda + 2*mu)*sin(2*pi*(x + x_sy))*sin(2*pi*(y + y_sy)));

  }

  T d2bx_dy2(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T dx, T mu,
             T lambda, T charU, T pi)
  {
    return pi*pi*pi*pi*(-0.0144*mu*sin(2*pi*(y + y_sx))*cos(2*pi*(x + x_sx)) + 0.0112*(lambda + mu)*sin(2*pi*(x + x_sy))*sin(2*pi*(y + y_sy)) - 0.0144*(lambda + 2*mu)*sin(2*pi*(y + y_sx))*cos(2*pi*(x + x_sx)));
  }

  // Output 1
  T by(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T dx, T mu, T lambda,
       T charU, T pi)
  {
    return 0.0028*mu*pi*pi*cos(2*pi*(x + x_sy))*cos(2*pi*(y + y_sy)) + 0.0036*pi*pi*(lambda + mu)*sin(2*pi*(x + x_sx))*cos(2*pi*(y + y_sx)) + 0.0028*pi*pi*(lambda + 2*mu)*cos(2*pi*(x + x_sy))*cos(2*pi*(y + y_sy));
  }

  T d2by_dy2(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T dx, T mu,
             T lambda, T charU, T pi)
  {
    return -pi*pi*pi*pi*(0.0112*mu*cos(2*pi*(x + x_sy))*cos(2*pi*(y + y_sy)) + 0.0144*(lambda + mu)*sin(2*pi*(x + x_sx))*cos(2*pi*(y + y_sx)) + 0.0112*(lambda + 2*mu)*cos(2*pi*(x + x_sy))*cos(2*pi*(y + y_sy)));
  }

  T d2bx_dxdy(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T dx, T mu,
              T lambda, T charU, T pi)
  {
    return -pi*pi*pi*pi*(0.0144*mu*sin(2*pi*(x + x_sx))*cos(2*pi*(y + y_sx)) + 0.0112*(lambda + mu)*cos(2*pi*(x + x_sy))*cos(2*pi*(y + y_sy)) + 0.0144*(lambda + 2*mu)*sin(2*pi*(x + x_sx))*cos(2*pi*(y + y_sx)));
  }

  T d2by_dx2(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T dx, T mu,
             T lambda, T charU, T pi)
  {
    return -pi*pi*pi*pi*(0.0112*mu*cos(2*pi*(x + x_sy))*cos(2*pi*(y + y_sy)) + 0.0144*(lambda + mu)*sin(2*pi*(x + x_sx))*cos(2*pi*(y + y_sx)) + 0.0112*(lambda + 2*mu)*cos(2*pi*(x + x_sy))*cos(2*pi*(y + y_sy)));
  }
};

template <typename T, typename DESCRIPTOR>
class ManufacturedSolutionU2D : public AnalyticalF2D<T, T> {

protected:

public:
  ManufacturedSolutionU2D(MyCase& myCase) : AnalyticalF2D<T, T>(2) {};

  bool operator()(T output[], const T input[]) override
  {
    T pi = std::numbers::pi_v<double>;

    T x = input[0];
    T y = input[1];

    T x_sx = -.3;
    T x_sy =  .4;
    T y_sx =  .1;
    T y_sy = -.7;

    // lattice units
    output[0] = 9. / 10000. * util::cos(2. * pi * (x + x_sx)) * util::sin(2. * pi * (y + y_sx));
    output[1] = 7. / 10000. * util::cos(2. * pi * (x + x_sy)) * util::cos(2. * pi * (y + y_sy));

    return true;
  };
};

template <typename T, typename DESCRIPTOR>
class ManufacturedSolutionStress2D : public AnalyticalF2D<T, T> {

protected:
  MyCase _myCase;

public:
  ManufacturedSolutionStress2D(MyCase& myCase) : AnalyticalF2D<T, T>(4),
  _myCase(myCase) {};

  bool operator()(T output[], const T input[]) override
  {
    T pi = std::numbers::pi_v<double>;

    const auto& converter = _myCase.getLattice(NavierCauchy{}).getUnitConverter();
    T dx = converter.getPhysDeltaX();
    T dt = converter.getPhysDeltaT();
    T kappa = converter.getDampingFactor();
    T lambda = converter.getLatticeLambda();
    T mu = converter.getLatticeShearModulus();

    T x = input[0];
    T y = input[1];

    T x_sx = -.3;
    T x_sy =  .4;
    T y_sx =  .1;
    T y_sy = -.7;

    T latticeFactor = dt / (kappa * dx);

    // xx
    output[0] = latticeFactor * ((2. * mu + lambda) * dux_dx(0, x, y, x_sx, x_sy, y_sx, y_sy, pi)
              + lambda * duy_dy(0, x, y, x_sx, x_sy, y_sx, y_sy, pi));
    // xy
    output[1] = latticeFactor * (mu * (dux_dy(0, x, y, x_sx, x_sy, y_sx, y_sy, pi)
              + duy_dx(0, x, y, x_sx, x_sy, y_sx, y_sy, pi)));

    // yx
    output[2] = latticeFactor * (mu * (dux_dy(0, x, y, x_sx, x_sy, y_sx, y_sy, pi)
              + duy_dx(0, x, y, x_sx, x_sy, y_sx, y_sy, pi)));

    // yy
    output[3] = latticeFactor * ((2. * mu + lambda) * duy_dy(0, x, y, x_sx, x_sy, y_sx, y_sy, pi)
              + lambda * dux_dx(0, x, y, x_sx, x_sy, y_sx, y_sy, pi));

    return true;
  };

  T dux_dx(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T pi) {
    return -(9. * pi * sin(2. * pi * (x + x_sx)) * sin(2. * pi * (y + y_sx))) / 5000.;
  }

  T dux_dy(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T pi) {
    return (9. * pi * cos(2. * pi * (x + x_sx)) * cos(2. * pi * (y + y_sx))) / 5000.;
  }

  T duy_dx(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T pi) {
    return -(7. * pi * cos(2. * pi * (y + y_sy)) * sin(2. * pi * (x + x_sy))) / 5000.;
  }

  T duy_dy(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T pi) {
    return -(7. * pi * cos(2. * pi * (x + x_sy)) * sin(2. * pi * (y + y_sy))) / 5000.;
  }
};

/// @brief Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector origin{0, 0};
  IndicatorCuboid2D<T> cuboid(extent, origin);

  const T physDeltaX = extent[0] / parameters.get<parameters::RESOLUTION>();
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({true,true});
  return mesh;
}

/// @brief Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
/// @note The material numbers will be used to assign physics to lattice nodes
void prepareGeometry(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();

  /// Set material numbers
  geometry.rename(0, 1);

  geometry.getStatistics().print();
}

/*
void prepareOmegas(LinElaUnitConverter<T, LDESCRIPTOR> converter,
                   T descriptors::invCs2<T,DESCRIPTOR>(),
                   std::vector<T>& allOmegas)
{
  T omega_11 = 1. / (converter.getLatticeShearModulus() / descriptors::invCs2<T,DESCRIPTOR>() + 0.5);
  T omega_d = 1. / (2. * converter.getLatticeShearModulus()/(1. - descriptors::invCs2<T,DESCRIPTOR>()) + 0.5);
  T omega_s = 1. / (2. * (converter.getLatticeShearModulus() + converter.getLatticeLambda()) / (1. + descriptors::invCs2<T,DESCRIPTOR>()) + 0.5);

  T tau_11 = 1. / omega_11 - 0.5;
  T tau_s = 1. / omega_d - 0.5;
  T tau_p = 1. / omega_s - 0.5;
  // T tau_12 = 0.5;
  T tau_12 = (1. / tau_11 + 4. * tau_11) / 8.;
  T tau_21 = tau_12;
  // T tau_22 = 0.5;
  T tau_22 = -(12. * pow(tau_11, 3) + 28. * pow(tau_11, 2) * tau_p + tau_11 * (-5. + 32. * pow(tau_p, 2)) + tau_p * (-9. + 16. * tau_12 * tau_p + 64. * pow(tau_p, 2))) / (3. * (1. + 4. * pow(tau_11, 2) - 8. * tau_12 * tau_p));
  T omega_12 = 1. / (tau_12 + 0.5);
  T omega_21 = 1. / (tau_21 + 0.5);
  T omega_22 = 1. / (tau_22 + 0.5);

  allOmegas = {omega_11, omega_d, omega_s, omega_12, omega_21, omega_22};
  return;
}
*/

/*
void prepareLattice(SuperLattice<T, LDESCRIPTOR>& lattice,
SuperGeometry<T, 2>& superGeometry,
LinElaUnitConverter<T, LDESCRIPTOR> converter,
std::vector<T>& allOmegas,
ForceField2D<T, T>& forceSol,
T descriptors::invCs2<T,DESCRIPTOR>())
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  lattice.defineDynamics<BoolakeeLinearElasticity>(superGeometry.getMaterialIndicator({1}));

  auto bulkIndicator = superGeometry.getMaterialIndicator({ 1 });

  std::vector<T>          iniPop = {0., 0., 0., 0., 0., 0., 0., 0.};
  std::vector<T>          iniMom = {0., 0.};
  std::vector<T>          iniStress = {0., 0., 0.};

  AnalyticalConst2D<T, T> initialPopulationF(iniPop);
  AnalyticalConst2D<T, T> initialMomentsF(iniMom);
  AnalyticalConst2D<T, T> initialStressF(iniStress);

  // dx, dt, descriptors::invCs2<T,DESCRIPTOR>(), mü, lambda, kappa, epsilon
  T magic[8] = {converter.getConversionFactorLength(),
  converter.getConversionFactorTime(),
  descriptors::invCs2<T,DESCRIPTOR>(),
  converter.getLatticeShearModulus(),
  converter.getLatticeLambda(),
  converter.getDampingFactor(),
  converter.getCharPhysVelocity(),
  converter.getEpsilon()};

  lattice.setParameter<descriptors::MAGIC_SOLID>(magic);

  lattice.setParameter<descriptors::OMEGA_SOLID>(allOmegas);
  lattice.defineField<FORCE>(        superGeometry, 1, forceSol);
  lattice.defineField<POPULATION>(   superGeometry, 1, initialPopulationF);
  lattice.defineField<MOMENT_VECTOR>(superGeometry, 1, initialMomentsF);
  lattice.defineField<SIGMA_SOLID>(  superGeometry, 1, initialStressF);

  /// Make the lattice ready for simulation
  lattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}
*/

/// @brief Set lattice dynamics
/// @param myCase The Case instance which keeps the simulation data
void prepareLattice(MyCase& myCase) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierCauchy>;
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierCauchy{});
  auto& geometry = myCase.getGeometry();

  /// Material=1 -->bulk dynamics
  lattice.defineDynamics<BoolakeeLinearElasticity>(geometry, 1);

  const T physDeltaX            = parameters.get<parameters::PHYS_DELTA_X>();
  const T physDeltaT            = parameters.get<parameters::PHYS_DELTA_T>();
  const T physCharLength        = parameters.get<parameters::PHYS_CHAR_LENGTH>();
  const T physCharDisplacement  = parameters.get<parameters::PHYS_CHAR_DISPLACEMENT>();
  const T physYoungsModulus     = parameters.get<parameters::YOUNGS_MODULUS>();
  const T physPoissonRatio      = parameters.get<parameters::POISSON_RATIO>();
  const T kappa                 = parameters.get<parameters::KAPPA>();

  // Set up a unit converter with the characteristic physical units
  lattice.setUnitConverter<LinElaUnitConverter<T, DESCRIPTOR>>(
    physDeltaX,
    physDeltaT,
    physCharLength,
    physCharDisplacement,
    physYoungsModulus,
    physPoissonRatio,
    kappa
  );
  const auto& converter = lattice.getUnitConverter();
  converter.print();

  T omega_11 = 1. / (converter.getLatticeShearModulus() / descriptors::invCs2<T,DESCRIPTOR>() + 0.5);
  T omega_d = 1. / (2. * converter.getLatticeShearModulus()/(1. - descriptors::invCs2<T,DESCRIPTOR>()) + 0.5);
  T omega_s = 1. / (2. * (converter.getLatticeShearModulus() + converter.getLatticeLambda()) / (1. + descriptors::invCs2<T,DESCRIPTOR>()) + 0.5);

  T tau_11 = 1. / omega_11 - 0.5;
  T tau_p = 1. / omega_s - 0.5;
  T tau_12 = (1. / tau_11 + 4. * tau_11) / 8.;
  T tau_21 = tau_12;
  T tau_22 = -(12. * pow(tau_11, 3) + 28. * pow(tau_11, 2) * tau_p + tau_11 * (-5. + 32. * pow(tau_p, 2)) + tau_p * (-9. + 16. * tau_12 * tau_p + 64. * pow(tau_p, 2))) / (3. * (1. + 4. * pow(tau_11, 2) - 8. * tau_12 * tau_p));
  T omega_12 = 1. / (tau_12 + 0.5);
  T omega_21 = 1. / (tau_21 + 0.5);
  T omega_22 = 1. / (tau_22 + 0.5);

  lattice.setParameter<descriptors::OMEGA_SOLID>({omega_11, omega_d, omega_s, omega_12, omega_21, omega_22});
  T magic[8] = {
    converter.getPhysDeltaX(),
    converter.getPhysDeltaT(),
    descriptors::invCs2<T,DESCRIPTOR>(),
    converter.getLatticeShearModulus(),
    converter.getLatticeLambda(),
    converter.getDampingFactor(),
    converter.getCharPhysDisplacement(),
    converter.getEpsilon()
  };

  lattice.setParameter<descriptors::MAGIC_SOLID>(magic);

  ForceField2D<T, DESCRIPTOR> force(myCase);
  lattice.defineField<descriptors::FORCE>(geometry, 1, force);
}

/// Set initial condition for primal variables (velocity and density)
/// @param myCase The Case instance which keeps the simulation data
/// @note Be careful: initial values have to be set using lattice units
void setInitialValues(MyCase& myCase)
{
  // Nothing to do here, because simulation is initialized with 0
}

/// Update boundary values at times (and external fields, if they exist)
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
/// @note Be careful: boundary values have to be set using lattice units
void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{
  // Nothing to do here, because simulation does not depend on time
}

/*
std::vector<T> getResults( SuperLattice<T, LDESCRIPTOR>& lattice,
                            LinElaUnitConverter<T, LDESCRIPTOR> converter,
                            SuperGeometry<T, 2>& superGeometry,
                            int iT,
                            int maxIt)
{
  OstreamManager clout(std::cout, "getResults");

  SuperVTMwriter2D<T> vtmWriter("periodicPlate");

  CSV<T> csvWriter("TEMP_CSV");

  SuperGeometryF<T,2> materials(superGeometry);
  vtmWriter.addFunctor( materials );

  if (iT == 0) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, LDESCRIPTOR>   cuboid(lattice);
    SuperLatticeRank2D<T, LDESCRIPTOR>     rank(lattice);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();
  }

  SuperLatticePhysVelocity2D<T, LDESCRIPTOR>          velocityPhys(lattice, converter);
  SuperLatticeVelocity2D<T, LDESCRIPTOR>              velocity(lattice);
  SuperLatticeDensity2D<T, LDESCRIPTOR>               density(lattice);
  SuperLatticeFpop2D<T, LDESCRIPTOR>                  population(lattice);
  SuperLatticeField2D<T, LDESCRIPTOR, FORCE>          forceField(lattice);

  ManufacturedSolutionU2D<T, LDESCRIPTOR>             dispSol(converter);
  SuperLatticeFfromAnalyticalF2D<T, LDESCRIPTOR>      dispSolLattice(dispSol, lattice);

  ManufacturedSolutionStress2D<T, LDESCRIPTOR>        stressSol(converter);
  SuperLatticeFfromAnalyticalF2D<T, LDESCRIPTOR>      stressSolLattice(stressSol, lattice);

  // Fields for error calc
  SuperLatticeField2D<T, LDESCRIPTOR, DISP_SOLID>     moments(lattice);
  SuperLatticeField2D<T, LDESCRIPTOR, SIGMA_SOLID>    stress(lattice);

  auto indicatorF = superGeometry.getMaterialIndicator({1});

  vtmWriter.addFunctor(population,       "population");
  vtmWriter.addFunctor(moments,          "numerical disp");
  vtmWriter.addFunctor(forceField,       "force");
  vtmWriter.addFunctor(dispSolLattice,   "analytical disp");
  vtmWriter.addFunctor(stressSolLattice, "analytical stress");
  vtmWriter.addFunctor(stress,           "numerical stress");

  vtmWriter.write(iT);

  T   l2UResult[2]        = {T(), T()};
  T   lInfUResult[2]      = {T(), T()};
  T   l2StressResult[2]   = {T(), T()};
  T   lInfStressResult[2] = {T(), T()};
  int tmp[]               = {int()};

  SuperRelativeErrorL2Norm2D<T>   relUErrorL2Norm(lattice, moments, dispSol, indicatorF);
  SuperRelativeErrorLinfNorm2D<T> relUErrorLinfNorm(lattice, moments, dispSol, indicatorF);
  SuperRelativeErrorL2Norm2D<T>   relStressErrorL2Norm( lattice, stress, stressSol, indicatorF );
  SuperRelativeErrorLinfNorm2D<T> relStressErrorLinfNorm( lattice, stress, stressSol, indicatorF );

  relUErrorL2Norm(l2UResult, tmp);
  relUErrorLinfNorm(lInfUResult, tmp);
  relStressErrorL2Norm(l2StressResult, tmp);
  relStressErrorLinfNorm(lInfStressResult, tmp);

  csvWriter.writeDataFile(iT, l2UResult[0],  "l2UErr");
  csvWriter.writeDataFile(iT, lInfUResult[0], "lInfUErr");
  csvWriter.writeDataFile(iT, l2StressResult[0],  "l2StressErr");
  csvWriter.writeDataFile(iT, lInfStressResult[0], "lInfStressErr");

  std::vector<T> returnVec = {l2UResult[0], lInfUResult[0], l2StressResult[0], lInfStressResult[0]};
  clout << "N\t" << "L2 U Error\t" << "LInf U Error\t" << "L2 Stress Error\t" << "LInf Stress Error" << std::endl;
  clout << converter.getResolution() << "\t" << returnVec[0] << "\t" << returnVec[1] << "\t" << returnVec[2] << "\t\t" << returnVec[3] << std::endl;

  return returnVec;

}
*/

/// Compute simulation results at times
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  OstreamManager clout(std::cout, "getResults");
  /// Write vtk plots every 0.3 seconds (of phys. simulation time)
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierCauchy>;
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierCauchy{});
  auto& converter = lattice.getUnitConverter();
  auto& geometry = myCase.getGeometry();

  const std::size_t iTlog = parameters.get<parameters::IT_LOG>();
  const std::size_t iTvtk = parameters.get<parameters::IT_VTK>();
  const std::size_t iTMax = converter.getLatticeTime(parameters.get<parameters::MAX_PHYS_T>());

  SuperVTMwriter2D<T> vtmWriter("periodicPlate");

  ManufacturedSolutionU2D<T, DESCRIPTOR>             dispSol(myCase);
  SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR>      dispSolLattice(dispSol, lattice);

  ManufacturedSolutionStress2D<T, DESCRIPTOR>        stressSol(myCase);
  SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR>      stressSolLattice(stressSol, lattice);

  // Fields for error calc
  SuperLatticeField2D<T, DESCRIPTOR, olb::descriptors::DISP_SOLID>     disp(lattice);
  SuperLatticeField2D<T, DESCRIPTOR, olb::descriptors::SIGMA_SOLID>    stress(lattice);

  if (iT == 0) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, DESCRIPTOR>   cuboid(lattice);
    SuperLatticeRank2D<T, DESCRIPTOR>     rank(lattice);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();
  }

  // Writes the VTK files
  if (iT%iTvtk == 0 && iT > 0) {
    vtmWriter.addFunctor(disp,             "numerical disp");
    vtmWriter.addFunctor(dispSolLattice,   "analytical disp");
    vtmWriter.addFunctor(stress,           "numerical stress");
    vtmWriter.addFunctor(stressSolLattice, "analytical stress");
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }

  /// Print some (numerical and computational) statistics
  if (iT%iTlog == 0 || iT % iTMax == 0) {
    T   l2UResult[2]        = {T(), T()};
    T   lInfUResult[2]      = {T(), T()};
    T   l2StressResult[2]   = {T(), T()};
    T   lInfStressResult[2] = {T(), T()};
    int tmp[]               = {int()};

    auto indicatorF = geometry.getMaterialIndicator({ 1 });
    SuperRelativeErrorL2Norm2D<T>   relUErrorL2Norm(lattice, disp, dispSol, indicatorF);
    SuperRelativeErrorLinfNorm2D<T> relUErrorLinfNorm(lattice, disp, dispSol, indicatorF);
    SuperRelativeErrorL2Norm2D<T>   relStressErrorL2Norm( lattice, stress, stressSol, indicatorF );
    SuperRelativeErrorLinfNorm2D<T> relStressErrorLinfNorm( lattice, stress, stressSol, indicatorF );

    relUErrorL2Norm(l2UResult, tmp);
    relUErrorLinfNorm(lInfUResult, tmp);
    relStressErrorL2Norm(l2StressResult, tmp);
    relStressErrorLinfNorm(lInfStressResult, tmp);

    clout << "N\t" << "L2 U Error\t" << "LInf U Error\t" << "L2 Stress Error\t" << "LInf Stress Error" << std::endl;
    clout << converter.getResolution() << "\t" << l2UResult[0] << "\t" << lInfUResult[0] << "\t" << l2StressResult[0] << "\t\t" << lInfStressResult[0] << std::endl;
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
    timer.print(iT);
  }
}

/// @brief Execute simulation: set initial values and run time loop
/// @param myCase The Case instance which keeps the simulation data
void simulate(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  const T physMaxT = parameters.get<parameters::MAX_PHYS_T>();

  const std::size_t iTmax = myCase.getLattice(NavierCauchy{}).getUnitConverter().getLatticeTime(physMaxT);
  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    myCase.getLattice(NavierCauchy{}).collideAndStream();

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, timer, iT);
  }

  timer.stop();
  timer.printSummary();
}

/*
int main(int argc, char **argv)
{
  // === 1st Step: Initialization ===
  OstreamManager clout(std::cout, "main");
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

  CLIreader args(argc, argv);
  const std::size_t N  = args.getValueOrFallback<std::size_t>("--res", 40);

  T simTime = 3.2;

  // Characteristical physicial values
  T dx = 1. / N;
  T dt = dx * dx;
  T charL = 1.;
  T charU = 1.;

  // Poisson Ratio
  T nu = 0.8;

  // Damping Factor
  T kappa = 1.;

  // Young Modulus
  T ELattice = 0.11;
  T E = ELattice * (dx * dx * kappa) / dt; // phys

  // Length of square plate
  T length = 1.0;

  T descriptors::invCs2<T,DESCRIPTOR>() = T(1) / descriptors::invCs2<T,LDESCRIPTOR>();

  LinElaUnitConverter<T, LDESCRIPTOR> converter(
    dx, // deltaX
    dt, // deltaT
    charL, // charL
    charU, // charU
    E, // phys
    nu, // phys
    kappa // phys
  );

  converter.print();

  const int noOfCuboids = singleton::mpi().getSize();

  std::vector<T> extend(2,T());
  std::vector<T> origin(2,T());

  extend[0] = length - dx;
  extend[1] = length - dx;

  origin[0] = dx / 2.;
  origin[1] = dx / 2.;

  IndicatorCuboid2D<T> cuboid(extend, origin);

  CuboidDecomposition<T, 2> cuboidDecomposition(cuboid, converter.getConversionFactorLength(), noOfCuboids);
  // cuboidDecomposition2D<T> cuboidDecomposition(0.0, 0.0, converter.getConversionFactorLength(), N, N, noOfCuboids);
  cuboidDecomposition.setPeriodicity({true, true});

  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  SuperGeometry<T, 2> superGeometry(cuboidDecomposition, loadBalancer);
  prepareGeometry(superGeometry);

  SuperLattice<T, LDESCRIPTOR> lattice(converter, superGeometry);

  ForceField2D<T, T> forceSol(descriptors::invCs2<T,DESCRIPTOR>(), converter);

  std::vector<T> allOmegas = {0., 0., 0., 0., 0., 0.};
  prepareOmegas(converter, descriptors::invCs2<T,DESCRIPTOR>(), allOmegas);

  prepareLattice(lattice, superGeometry, converter, allOmegas, forceSol, descriptors::invCs2<T,DESCRIPTOR>());

  converter.print();

  int maxIt = simTime / converter.getConversionFactorTime();

  clout << "Awaiting " << maxIt << " Time steps. Starting simulation..." << std::endl;
  util::Timer<T> timer(converter.getLatticeTime(simTime),
                       superGeometry.getStatistics().getNvoxel());
  timer.start();

  std::vector<T> errors = {};
  std::ofstream fout;
  std::string dataFile = singleton::directories().getLogOutDir() + std::to_string(N) +  "_Err.dat";
  fout.open(dataFile.c_str(), std::ios::trunc);
  fout << "N;it;l2UErr;lInfUErr;l2StressErr;lInfStressErr" << std::endl;

  int numDataPoints = 400.;

  for (int iT = 0; iT <= maxIt; ++iT) {

    timer.update(iT);

    if (iT % (maxIt / numDataPoints) == 0) {
      errors = getResults(lattice, converter, superGeometry, iT, maxIt);
      timer.printStep();
      if (singleton::mpi().isMainProcessor()) {
        fout << N << ";" << converter.getPhysTime(iT) << ";" << errors[0] << ";" << errors[1] << ";" <<errors[2] << ";" << errors[3] << std::endl;
      }
    }

    lattice.collideAndStream();
  }
  fout.close();
  timer.stop();
}
*/

/// Setup and run a simulation
int main(int argc, char* argv[]) {
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<KAPPA                  >(  1.          );
    myCaseParameters.set<RESOLUTION             >( 40           );
    myCaseParameters.set<YOUNGS_MODULUS         >(  0.11        );
    myCaseParameters.set<POISSON_RATIO          >(  0.7         );
    myCaseParameters.set<PHYS_CHAR_DISPLACEMENT >(  1.0         );
    myCaseParameters.set<PHYS_CHAR_LENGTH       >(  1.0         );
    myCaseParameters.set<DOMAIN_EXTENT          >( { 1.0, 1.0 } );
    myCaseParameters.set<MAX_PHYS_T             >(  3.2         );
    myCaseParameters.set<IT_LOG                 >(100           );
    myCaseParameters.set<IT_VTK                 >(100           );
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