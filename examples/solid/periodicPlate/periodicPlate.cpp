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

#include "olb.h"
#include <cmath>

using namespace olb;
using namespace olb::descriptors;

using T = FLOATING_POINT_TYPE;

using LDESCRIPTOR = D2Q8<tag::MRT,FORCE,MOMENT_VECTOR,DISP_SOLID,SIGMA_SOLID>;

constexpr T pi = M_PI;

template <typename T, typename _DESCRIPTOR>
class ForceField2D : public AnalyticalF2D<T, T> {

protected:
  T dt;
  T dx;
  T charU;
  T lambda;
  T kappa;
  T mu;
  T theta;

public:
  ForceField2D(T inTheta, LinElaUnitConverter<T, LDESCRIPTOR> converter)
      : AnalyticalF2D<T, T>(2),
      dt(converter.getConversionFactorTime()),
      dx(converter.getConversionFactorLength()),
      charU(converter.getCharPhysVelocity()),
      lambda(converter.getLatticeLambda()),
      kappa(converter.getDampingFactor()),
      mu(converter.getLatticeShearModulus()),
      theta(inTheta) {};

  bool operator()(T output[], const T input[]) override
  {

    T x = input[0];
    T y = input[1];

    T x_sx = -.3;
    T x_sy =  .4;
    T y_sx =  .1;
    T y_sy = -.7;

    T omega_11 = 1. / (mu / theta + 0.5);
    T omega_d  = 1. / (2 * mu / (1 - theta) + 0.5);
    T omega_s  = 1. / (2 * (mu + lambda) / (1 + theta) + 0.5);

    T tau_11 = 1. / omega_11 - 0.5;
    T tau_s = 1. / omega_d - 0.5;
    T tau_p = 1. / omega_s - 0.5;
    T tau_12 = (1. / tau_11 + 4. * tau_11) / 8.;
    T tau_21 = tau_12;
    T tau_22 = -(12. * pow(tau_11, 3)
             + 28. * pow(tau_11, 2) * tau_p + tau_11 * (-5. + 32. * pow(tau_p, 2))
             + tau_p * (-9. + 16. * tau_12 * tau_p + 64. * pow(tau_p, 2))) / (3. * (1. + 4. * pow(tau_11, 2) - 8. * tau_12 * tau_p));
    T omega_12 = 1. / (tau_12 + 0.5);
    T omega_21 = 1. / (tau_21 + 0.5);
    T omega_22 = 1. / (tau_22 + 0.5);

    // These are correct
    T d1 = -1. / 4. - 1. / 2. * theta * tau_12 * tau_s + pow(tau_s, 2) / 2.
         -  1. / 2. * theta * pow(tau_s, 2) + 1. / 2. * theta * tau_12 * tau_p
         +  pow(tau_p, 2) / 2. + 1. / 2. * theta * pow(tau_p, 2);
    T d2 = -theta / 2. + theta * pow(tau_11, 2) + theta * tau_11 * tau_12
         + 1. / 2. * theta * tau_12 * tau_s - pow(tau_s, 2) / 2.
         + 1. / 2. * theta * pow(tau_s, 2) + 1. / 2. * theta * tau_12 * tau_p
         + pow(tau_p, 2) / 2. + 1. / 2. * theta * pow(tau_p, 2);
    T d3 = -theta / 4.
         +  theta * pow(tau_11, 2)
         +  theta * tau_11 * tau_12;

    T mu_phys     =     mu * (dx * dx * kappa) / dt;
    T lambda_phys = lambda * (dx * dx * kappa) / dt;

    output[0] = bx(0, x, y, x_sx, x_sy, y_sx, y_sy, dx, mu_phys, lambda_phys, charU) * dt / kappa
              + (d1 * d2bx_dx2(0, x, y, x_sx, x_sy, y_sx, y_sy, dx, mu_phys, lambda_phys, charU)
              + d2 * d2by_dxdy(0, x, y, x_sx, x_sy, y_sx, y_sy, dx, mu_phys, lambda_phys, charU)
              + d3 * d2bx_dy2(0, x, y, x_sx, x_sy, y_sx, y_sy, dx, mu_phys, lambda_phys, charU)) * dx * dx * dt / kappa;

    output[1] = by(0, x, y, x_sx, x_sy, y_sx, y_sy, dx, mu_phys, lambda_phys, charU) * dt / kappa
              + (d1 * d2by_dy2(0, x, y, x_sx, x_sy, y_sx, y_sy, dx, mu_phys, lambda_phys, charU)
              + d2 * d2bx_dxdy(0, x, y, x_sx, x_sy, y_sx, y_sy, dx, mu_phys, lambda_phys, charU)
              + d3 * d2by_dx2(0, x, y, x_sx, x_sy, y_sx, y_sy, dx, mu_phys, lambda_phys, charU)) * dx * dx * dt / kappa;

    return true;
  };

  // Output 0
  T bx(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T dx, T mu, T lambda,
       T charU)
  {
    return 0.0036*mu*pi*pi*sin(2*pi*(y + y_sx))*cos(2*pi*(x + x_sx)) - 0.0028*pi*pi*(lambda + mu)*sin(2*pi*(x + x_sy))*sin(2*pi*(y + y_sy)) + 0.0036*pi*pi*(lambda + 2*mu)*sin(2*pi*(y + y_sx))*cos(2*pi*(x + x_sx));
  }

  T d2bx_dx2(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T dx, T mu,
             T lambda, T charU)
  {
    return pi*pi*pi*pi*(-0.0144*mu*sin(2*pi*(y + y_sx))*cos(2*pi*(x + x_sx)) + 0.0112*(lambda + mu)*sin(2*pi*(x + x_sy))*sin(2*pi*(y + y_sy)) - 0.0144*(lambda + 2*mu)*sin(2*pi*(y + y_sx))*cos(2*pi*(x + x_sx)));
  }

  T d2by_dxdy(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T dx, T mu,
              T lambda, T charU)
  {
    return pi*pi*pi*pi*(0.0112*mu*sin(2*pi*(x + x_sy))*sin(2*pi*(y + y_sy)) - 0.0144*(lambda + mu)*sin(2*pi*(y + y_sx))*cos(2*pi*(x + x_sx)) + 0.0112*(lambda + 2*mu)*sin(2*pi*(x + x_sy))*sin(2*pi*(y + y_sy)));

  }

  T d2bx_dy2(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T dx, T mu,
             T lambda, T charU)
  {
    return pi*pi*pi*pi*(-0.0144*mu*sin(2*pi*(y + y_sx))*cos(2*pi*(x + x_sx)) + 0.0112*(lambda + mu)*sin(2*pi*(x + x_sy))*sin(2*pi*(y + y_sy)) - 0.0144*(lambda + 2*mu)*sin(2*pi*(y + y_sx))*cos(2*pi*(x + x_sx)));
  }

  // Output 1
  T by(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T dx, T mu, T lambda,
       T charU)
  {
    return 0.0028*mu*pi*pi*cos(2*pi*(x + x_sy))*cos(2*pi*(y + y_sy)) + 0.0036*pi*pi*(lambda + mu)*sin(2*pi*(x + x_sx))*cos(2*pi*(y + y_sx)) + 0.0028*pi*pi*(lambda + 2*mu)*cos(2*pi*(x + x_sy))*cos(2*pi*(y + y_sy));
  }

  T d2by_dy2(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T dx, T mu,
             T lambda, T charU)
  {
    return -pi*pi*pi*pi*(0.0112*mu*cos(2*pi*(x + x_sy))*cos(2*pi*(y + y_sy)) + 0.0144*(lambda + mu)*sin(2*pi*(x + x_sx))*cos(2*pi*(y + y_sx)) + 0.0112*(lambda + 2*mu)*cos(2*pi*(x + x_sy))*cos(2*pi*(y + y_sy)));
  }

  T d2bx_dxdy(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T dx, T mu,
              T lambda, T charU)
  {
    return -pi*pi*pi*pi*(0.0144*mu*sin(2*pi*(x + x_sx))*cos(2*pi*(y + y_sx)) + 0.0112*(lambda + mu)*cos(2*pi*(x + x_sy))*cos(2*pi*(y + y_sy)) + 0.0144*(lambda + 2*mu)*sin(2*pi*(x + x_sx))*cos(2*pi*(y + y_sx)));
  }

  T d2by_dx2(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy, T dx, T mu,
             T lambda, T charU)
  {
    return -pi*pi*pi*pi*(0.0112*mu*cos(2*pi*(x + x_sy))*cos(2*pi*(y + y_sy)) + 0.0144*(lambda + mu)*sin(2*pi*(x + x_sx))*cos(2*pi*(y + y_sx)) + 0.0112*(lambda + 2*mu)*cos(2*pi*(x + x_sy))*cos(2*pi*(y + y_sy)));
  }
};

template <typename T, typename _DESCRIPTOR>
class ManufacturedSolutionU2D : public AnalyticalF2D<T, T> {

protected:
  T charU;

public:
  ManufacturedSolutionU2D(LinElaUnitConverter<T, LDESCRIPTOR> converter) : AnalyticalF2D<T, T>(2),
      charU(converter.getCharPhysVelocity())
      {};

  bool operator()(T output[], const T input[]) override
  {
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

template <typename T, typename _DESCRIPTOR>
class ManufacturedSolutionStress2D : public AnalyticalF2D<T, T> {

protected:
  T dt;
  T dx;
  T lambda;
  T mu;
  T kappa;

public:
  ManufacturedSolutionStress2D(LinElaUnitConverter<T, LDESCRIPTOR> converter) : AnalyticalF2D<T, T>(4),
      dt(converter.getConversionFactorTime()),
      dx(converter.getConversionFactorLength()),
      mu(converter.getPhysShearModulus()),
      lambda(converter.getPhysLambda()),
      kappa(converter.getDampingFactor()) {};

  bool operator()(T output[], const T input[]) override
  {
    T x = input[0];
    T y = input[1];

    T x_sx = -.3;
    T x_sy =  .4;
    T y_sx =  .1;
    T y_sy = -.7;

    T latticeFactor = dt / (kappa * dx);

    // xx
    output[0] = latticeFactor * ((2. * mu + lambda) * dux_dx(0, x, y, x_sx, x_sy, y_sx, y_sy)
              + lambda * duy_dy(0, x, y, x_sx, x_sy, y_sx, y_sy));
    // xy
    output[1] = latticeFactor * (mu * (dux_dy(0, x, y, x_sx, x_sy, y_sx, y_sy)
              + duy_dx(0, x, y, x_sx, x_sy, y_sx, y_sy)));

    // yx
    output[2] = latticeFactor * (mu * (dux_dy(0, x, y, x_sx, x_sy, y_sx, y_sy)
              + duy_dx(0, x, y, x_sx, x_sy, y_sx, y_sy)));

    // yy
    output[3] = latticeFactor * ((2. * mu + lambda) * duy_dy(0, x, y, x_sx, x_sy, y_sx, y_sy)
              + lambda * dux_dx(0, x, y, x_sx, x_sy, y_sx, y_sy));

    return true;
  };

  T dux_dx(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy) {
    return -(9. * pi * sin(2. * pi * (x + x_sx)) * sin(2. * pi * (y + y_sx))) / 5000.;
  }

  T dux_dy(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy) {
    return (9. * pi * cos(2. * pi * (x + x_sx)) * cos(2. * pi * (y + y_sx))) / 5000.;
  }

  T duy_dx(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy) {
    return -(7. * pi * cos(2. * pi * (y + y_sy)) * sin(2. * pi * (x + x_sy))) / 5000.;
  }

  T duy_dy(T t, T x, T y, T x_sx, T x_sy, T y_sx, T y_sy) {
    return -(7. * pi * cos(2. * pi * (x + x_sy)) * sin(2. * pi * (y + y_sy))) / 5000.;
  }
};

void prepareGeometry(SuperGeometry<T, 2>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0, 1 );

  superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareOmegas(LinElaUnitConverter<T, LDESCRIPTOR> converter,
                   T theta,
                   std::vector<T>& allOmegas)
{
  T omega_11 = 1. / (converter.getLatticeShearModulus() / theta + 0.5);
  T omega_d = 1. / (2. * converter.getLatticeShearModulus()/(1. - theta) + 0.5);
  T omega_s = 1. / (2. * (converter.getLatticeShearModulus() + converter.getLatticeLambda()) / (1. + theta) + 0.5);

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

void prepareLattice(SuperLattice<T, LDESCRIPTOR>& lLattice,
                    SuperGeometry<T, 2>& superGeometry,
                    LinElaUnitConverter<T, LDESCRIPTOR> converter,
                    std::vector<T>& allOmegas,
                    ForceField2D<T, T>& forceSol,
                    T theta)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  lLattice.defineDynamics<BoolakeeLinearElasticity>(superGeometry.getMaterialIndicator({1}));

  auto bulkIndicator = superGeometry.getMaterialIndicator({ 1 });

  std::vector<T>          iniPop = {0., 0., 0., 0., 0., 0., 0., 0.};
  std::vector<T>          iniMom = {0., 0.};
  std::vector<T>          iniStress = {0., 0., 0.};

  AnalyticalConst2D<T, T> initialPopulationF(iniPop);
  AnalyticalConst2D<T, T> initialMomentsF(iniMom);
  AnalyticalConst2D<T, T> initialStressF(iniStress);

  // dx, dt, theta, mü, lambda, kappa, epsilon
  T magic[8] = {converter.getConversionFactorLength(),
                converter.getConversionFactorTime(),
                theta,
                converter.getLatticeShearModulus(),
                converter.getLatticeLambda(),
                converter.getDampingFactor(),
                converter.getCharPhysVelocity(),
                converter.getEpsilon()};

  lLattice.setParameter<descriptors::MAGIC_SOLID>(magic);

  lLattice.setParameter<descriptors::OMEGA_SOLID>(allOmegas);
  lLattice.defineField<FORCE>(        superGeometry, 1, forceSol);
  lLattice.defineField<POPULATION>(   superGeometry, 1, initialPopulationF);
  lLattice.defineField<MOMENT_VECTOR>(superGeometry, 1, initialMomentsF);
  lLattice.defineField<SIGMA_SOLID>(  superGeometry, 1, initialStressF);

  /// Make the lattice ready for simulation
  lLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

std::vector<T> getResults( SuperLattice<T, LDESCRIPTOR>& lLattice,
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
    SuperLatticeCuboid2D<T, LDESCRIPTOR>   cuboid(lLattice);
    SuperLatticeRank2D<T, LDESCRIPTOR>     rank(lLattice);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();
  }

  SuperLatticePhysVelocity2D<T, LDESCRIPTOR>          velocityPhys(lLattice, converter);
  SuperLatticeVelocity2D<T, LDESCRIPTOR>              velocity(lLattice);
  SuperLatticeDensity2D<T, LDESCRIPTOR>               density(lLattice);
  SuperLatticeFpop2D<T, LDESCRIPTOR>                  population(lLattice);
  SuperLatticeField2D<T, LDESCRIPTOR, FORCE>          forceField(lLattice);

  ManufacturedSolutionU2D<T, LDESCRIPTOR>             dispSol(converter);
  SuperLatticeFfromAnalyticalF2D<T, LDESCRIPTOR>      dispSolLattice(dispSol, lLattice);

  ManufacturedSolutionStress2D<T, LDESCRIPTOR>        stressSol(converter);
  SuperLatticeFfromAnalyticalF2D<T, LDESCRIPTOR>      stressSolLattice(stressSol, lLattice);

  // Fields for error calc
  SuperLatticeField2D<T, LDESCRIPTOR, DISP_SOLID>     moments(lLattice);
  SuperLatticeField2D<T, LDESCRIPTOR, SIGMA_SOLID>    stress(lLattice);

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

  SuperRelativeErrorL2Norm2D<T>   relUErrorL2Norm(lLattice, moments, dispSol, indicatorF);
  SuperRelativeErrorLinfNorm2D<T> relUErrorLinfNorm(lLattice, moments, dispSol, indicatorF);
  SuperRelativeErrorL2Norm2D<T>   relStressErrorL2Norm( lLattice, stress, stressSol, indicatorF );
  SuperRelativeErrorLinfNorm2D<T> relStressErrorLinfNorm( lLattice, stress, stressSol, indicatorF );

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

  T theta = T(1) / descriptors::invCs2<T,LDESCRIPTOR>();

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

  SuperLattice<T, LDESCRIPTOR> lLattice(converter, superGeometry);

  ForceField2D<T, T> forceSol(theta, converter);

  std::vector<T> allOmegas = {0., 0., 0., 0., 0., 0.};
  prepareOmegas(converter, theta, allOmegas);

  prepareLattice(lLattice, superGeometry, converter, allOmegas, forceSol, theta);

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
      errors = getResults(lLattice, converter, superGeometry, iT, maxIt);
      timer.printStep();
      if (singleton::mpi().isMainProcessor()) {
        fout << N << ";" << converter.getPhysTime(iT) << ";" << errors[0] << ";" << errors[1] << ";" <<errors[2] << ";" << errors[3] << std::endl;
      }
    }

    lLattice.collideAndStream();
  }
  fout.close();
  timer.stop();
}