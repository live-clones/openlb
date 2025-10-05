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


using LDESCRIPTOR = D2Q8<tag::MRT,
                         FORCE,
                         DISP_SOLID,
                         SIGMA_SOLID,
                         MAGIC_SOLID,
                         OMEGA_SOLID,
                         BOUNDARY_COORDS_X,
                         BOUNDARY_COORDS_Y,
                         SOLID_DISTANCE_FIELD,
                         NEUMANN_SOLID_C,
                         NEUMANN_SOLID_NORMAL_X,
                         NEUMANN_SOLID_NORMAL_Y>;


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

    T omega_11 = 1. / (mu / theta + 0.5);
    T omega_d  = 1. / (2 * mu / (1 - theta) + 0.5);
    T omega_s  = 1. / (2 * (mu + lambda) / (1 + theta) + 0.5);

    T tau_11 = 1. / omega_11 - 0.5;
    T tau_s = 1. / omega_d - 0.5;
    T tau_p = 1. / omega_s - 0.5;

    #ifdef usetaucalc
      T tau_12 = (1. / tau_11 + 4. * tau_11) / 8.;
      T tau_22 = -(12. * pow(tau_11, 3)
              + 28. * pow(tau_11, 2) * tau_p + tau_11 * (-5. + 32. * pow(tau_p, 2))
              + tau_p * (-9. + 16. * tau_12 * tau_p + 64. * pow(tau_p, 2))) / (3. * (1. + 4. * pow(tau_11, 2) - 8. * tau_12 * tau_p));
    #else
      T tau_12 = 0.5;
      T tau_22 = 0.5;
    #endif

    T tau_21 = tau_12;

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


    output[0] = bx(0, x, y, dx, mu_phys, lambda_phys, charU) * dt / kappa
              + (d1 * d2bx_dx2(0, x, y, dx, mu_phys, lambda_phys, charU)
              + d2 * d2by_dxdy(0, x, y, dx, mu_phys, lambda_phys, charU)
              + d3 * d2bx_dy2(0, x, y, dx, mu_phys, lambda_phys, charU)) * dx * dx * dt / kappa;

    output[1] = by(0, x, y, dx, mu_phys, lambda_phys, charU) * dt / kappa
              + (d1 * d2by_dy2(0, x, y, dx, mu_phys, lambda_phys, charU)
              + d2 * d2bx_dxdy(0, x, y, dx, mu_phys, lambda_phys, charU)
              + d3 * d2by_dx2(0, x, y, dx, mu_phys, lambda_phys, charU)) * dx * dx * dt / kappa;

    return true;
  };

  T bx(T t, T x, T y, T dx, T mu, T lambda, T charU)
  {
    return   0.0036 * mu * pi * pi * x * x * util::sin(2 * pi * x * y)
           + 0.0036 * pi * pi * y * y * (lambda + 2 * mu) * util::sin(2 * pi * x * y)
           + 0.0028 * pi * x * (lambda + mu) * util::sin(2 * pi * y);
  }

  T d2bx_dx2(T t, T x, T y, T dx, T mu, T lambda, T charU)
  {
    return   pi * pi * (-0.0144 * mu * pi * pi * x * x * y * y * util::sin(2 * pi * x * y)
           + 0.0288 * mu * pi * x * y * util::cos(2 * pi * x * y)
           + 0.0072 * mu * util::sin(2 * pi * x * y)
           - 0.0144 * pi * pi * y * y * y * y * (lambda + 2 * mu) * util::sin(2 * pi * x * y));
  }

  T d2by_dxdy(T t, T x, T y, T dx, T mu, T lambda, T charU)
  {
    return    pi * pi * (-0.0112 * pi * x * (lambda + 2 * mu) * util::sin(2 * pi * y)
           + (lambda + mu) * (-0.0144 * pi * pi * x * x * y * y * util::sin(2 * pi * x * y)
           + 0.0288 * pi * x * y * util::cos(2 * pi * x * y)
           + 0.0072 * util::sin(2 * pi * x * y)));
  }

  T d2bx_dy2(T t, T x, T y, T dx, T mu, T lambda, T charU)
  {
    return   pi * pi * (-0.0144 * mu * pi * pi * x * x * x * x * util::sin(2 * pi * x * y)
           - 0.0144 * pi * pi * x * x * y * y * (lambda + 2 * mu) * util::sin(2 * pi * x * y)
           + 0.0288 * pi * x * y * (lambda + 2 * mu) * util::cos(2 * pi * x * y)
           - 0.0112 * pi * x * (lambda + mu) * util::sin(2 * pi * y)
           + 0.0072 * (lambda + 2 * mu) * util::sin(2 * pi * x * y));
  }

  T by(T t, T x, T y, T dx, T mu, T lambda, T charU)
  {
    return - 0.0014 * mu * cos(2 * pi * y)
           + 0.0028 * pi * pi * (lambda + 2 * mu) * (x * x + 1) * cos(2 * pi * y)
           - pi * (lambda + mu) * (-0.0036 * pi * x * y * sin(2 * pi * x * y)
           + 0.0018 * cos(2 * pi * x * y));
  }

  T d2by_dy2(T t, T x, T y, T dx, T mu, T lambda, T charU)
  {
    return   pi*pi*(0.0056*mu*cos(2*pi*y)
           - 0.0112*pi*pi*(lambda + 2*mu)*(x*x + 1)*cos(2*pi*y)
           - pi*x*x*(lambda + mu)*(0.0144*pi*x*y*sin(2*pi*x*y)
           - 0.0216*cos(2*pi*x*y)));
  }

  T d2bx_dxdy(T t, T x, T y, T dx, T mu, T lambda, T charU)
  {
    return   pi*pi*(-0.0144*mu*pi*pi*x*x*x*y*sin(2*pi*x*y)
           + 0.0216*mu*pi*x*x*cos(2*pi*x*y)
           - 0.0144*pi*pi*x*y*y*y*(lambda + 2*mu)*sin(2*pi*x*y)
           + 0.0216*pi*y*y*(lambda + 2*mu)*cos(2*pi*x*y)
           + 0.0056*(lambda + mu)*cos(2*pi*y));
  }

  T d2by_dx2(T t, T x, T y, T dx, T mu, T lambda, T charU)
  {
    return   pi*pi*(-pi*y*y*(lambda + mu)*(0.0144*pi*x*y*sin(2*pi*x*y)
           - 0.0216*cos(2*pi*x*y))
           + 0.0056*(lambda + 2*mu)*cos(2*pi*y));
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

    // lattice units
    output[0] = 9. / 10000. * util::sin(2. * pi * x * y);
    output[1] = 7. / 10000. * util::cos(2. * pi * y) * (x * x + 1);

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

    T latticeFactor = dt / (kappa * dx);

    // xx
    output[0] = latticeFactor * ((2. * mu + lambda) * dux_dx(0, x, y)
                                          + lambda  * duy_dy(0, x, y));
    // xy
    output[1] = latticeFactor * (mu * (  dux_dy(0, x, y)
                                       + duy_dx(0, x, y)  ));

    // yx
    output[2] = latticeFactor * (mu * (  dux_dy(0, x, y)
              + duy_dx(0, x, y)  ));

    // yy
    output[3] = latticeFactor * ((2. * mu + lambda) * duy_dy(0, x, y)
                                          + lambda  * dux_dx(0, x, y));

    return true;
  };

  T dux_dx(T t, T x, T y) {
    return 0.0018 * pi * y * util::cos(2. * pi * x * y);
  }

  T dux_dy(T t, T x, T y) {
    return 0.0018 * pi * x * util::cos(2 * pi * x * y);
  }

  T duy_dx(T t, T x, T y) {
    return 0.0014 * x * util::cos(2 * pi * y);
  }

  T duy_dy(T t, T x, T y) {
    return -0.0014 * pi * (x * x + 1.) * util::sin(2 * pi * y);
  }
};

void prepareGeometry(SuperGeometry<T, 2>& superGeometry,
                    LinElaUnitConverter<T,LDESCRIPTOR> converter,
                    /*
                    std::shared_ptr<IndicatorF2D<T>> ellipse1,
                    std::shared_ptr<IndicatorF2D<T>> ellipse2,
                    std::shared_ptr<IndicatorF2D<T>> ellipse3)
                    */
                   IndicatorEllipse2D<T>& ellipse1,
                    IndicatorEllipse2D<T>& ellipse2,
                    IndicatorEllipse2D<T>& ellipse3)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0, 5 );

  superGeometry.rename( 5, 1, ellipse1 );
  superGeometry.rename( 1, 6, ellipse2 );
  superGeometry.rename( 1, 7, ellipse3 );

  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareOmegas( LinElaUnitConverter<T, LDESCRIPTOR>  converter,
                    T                                    theta,
                    std::vector<T>&                      allOmegas ) {
    T omega_11  = 1. / (      converter.getLatticeShearModulus()                                 /       theta  + 0.5);
    T omega_d   = 1. / (2. *  converter.getLatticeShearModulus()                                 / (1. - theta) + 0.5);
    T omega_s   = 1. / (2. * (converter.getLatticeShearModulus() + converter.getLatticeLambda()) / (1. + theta) + 0.5);

    T tau_11    =  1. / omega_11 - 0.5;
    T tau_s     =  1. / omega_d - 0.5;
    T tau_p     =  1. / omega_s - 0.5;

    #ifdef usetaucalc
      T tau_12 = (1. / tau_11 + 4. * tau_11) / 8.;
      T tau_22 = -(12. * pow(tau_11, 3)
                 + 28. * pow(tau_11, 2) * tau_p + tau_11 * (-5. + 32. * pow(tau_p, 2))
                 + tau_p * (-9. + 16. * tau_12 * tau_p + 64. * pow(tau_p, 2))) / (3. * (1. + 4. * pow(tau_11, 2) - 8. * tau_12 * tau_p));
    #else
      T tau_12 = 0.5;
      T tau_22 = 0.5;
    #endif

    T tau_21 = tau_12;

    T omega_12  = 1. / (tau_12 + 0.5);
    T omega_21  = 1. / (tau_21 + 0.5);
    T omega_22  = 1. / (tau_22 + 0.5);

    allOmegas   = {omega_11, omega_s, omega_d, omega_12, omega_21, omega_22};

  return;
}

void prepareLattice(SuperLattice<T, LDESCRIPTOR>&       lLattice,
                    SuperGeometry<T, 2>&                superGeometry,
                    LinElaUnitConverter<T, LDESCRIPTOR> converter,
                    std::vector<T>&                     allOmegas,
                    T                                   theta,
                    IndicatorEllipse2D<T>& ellipse1,
                    IndicatorEllipse2D<T>& ellipse2,
                    IndicatorEllipse2D<T>& ellipse3,
                    int                                 bulkNum = 1)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;
  auto bulkIndicator = superGeometry.getMaterialIndicator({ bulkNum });

  lLattice.defineDynamics<BoolakeeLinearElasticityBoundary>( bulkIndicator );

  setBoolakeeNeumannBoundary<T,LDESCRIPTOR,BoolakeeNeumannPostProcessor<T,LDESCRIPTOR>>(lLattice, superGeometry.getMaterialIndicator( 5 ), bulkIndicator, ellipse1, true);
  setBoolakeeNeumannBoundary<T,LDESCRIPTOR,BoolakeeNeumannPostProcessor<T,LDESCRIPTOR>>(lLattice, superGeometry.getMaterialIndicator( 6 ), bulkIndicator, ellipse2, false);
  setBoolakeeNeumannBoundary<T,LDESCRIPTOR,BoolakeeNeumannPostProcessor<T,LDESCRIPTOR>>(lLattice, superGeometry.getMaterialIndicator( 7 ), bulkIndicator, ellipse3, false);

  {
    auto& communicator = lLattice.getCommunicator(stage::PostCollide());
    communicator.template requestField<descriptors::NEUMANN_SOLID_C>();
    communicator.template requestField<descriptors::NEUMANN_SOLID_NORMAL_X>();
    communicator.template requestField<descriptors::NEUMANN_SOLID_NORMAL_Y>();
    communicator.template requestField<descriptors::SOLID_DISTANCE_FIELD>();
    communicator.template requestField<descriptors::BOUNDARY_COORDS_X>();
    communicator.template requestField<descriptors::BOUNDARY_COORDS_Y>();
    communicator.requestOverlap(1);
    communicator.exchangeRequests();
  }

  std::vector<T>          iniPop = {0., 0., 0., 0., 0., 0., 0., 0.};
  std::vector<T>          iniMom = {0., 0.};
  std::vector<T>          iniStress = {0., 0., 0.};

  AnalyticalConst2D<T, T> initialPopulationF(iniPop);
  AnalyticalConst2D<T, T> initialMomentsF(iniMom);
  AnalyticalConst2D<T, T> initialStressF(iniStress);

  // dx, dt, theta, mü, lambda, kappa, uChar, epsilon
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

  T neumann_constants[3] = {(2. * (1. - theta) * (converter.getLatticeBulkModulus() - converter.getLatticeShearModulus())) / (theta * (1. - theta - 4. * converter.getLatticeShearModulus())),
                            2. * converter.getLatticeShearModulus() / (theta - 2. * converter.getLatticeShearModulus()),
                            4. * converter.getLatticeShearModulus() / (1. - theta - 4. * converter.getLatticeShearModulus())};

  lLattice.setParameter<descriptors::NEUMANN_SOLID_C>(neumann_constants);

  ForceField2D<T, T> forceSol(theta, converter);
  lLattice.defineField<FORCE>(superGeometry, bulkNum, forceSol);
  lLattice.defineField<POPULATION>(superGeometry, bulkNum, initialPopulationF);

  lLattice.defineField<MOMENT_VECTOR>(superGeometry, bulkNum, initialMomentsF);
  lLattice.defineField<SIGMA_SOLID>(superGeometry, bulkNum, initialStressF);

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

  SuperVTMwriter2D<T> vtmWriter("ellipseDirichlet");

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

int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  OstreamManager clout(std::cout, "main");
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

  CLIreader args(argc, argv);
  const std::size_t N  = args.getValueOrFallback<std::size_t>("--res", 40);

  // length of domain
  T plateLength = 1.5;

  T charL = 1.;
  T dx = plateLength / N;

  T simTime = 60.0;
  T dt = dx * dx;
  T charU = 1.;

  T nu = 0.7;

  T kappa = 1.;

  T ELattice = 0.1;
  T E = ELattice * (dx * dx * kappa) / dt; // phys

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

  T origin = converter.getConversionFactorLength() / 3.;
  IndicatorCuboid2D<T> cuboid({plateLength * charL, plateLength * charL}, {origin, origin});
  CuboidDecomposition<T, 2> cuboidDecomposition(cuboid, converter.getConversionFactorLength(), noOfCuboids);

  Vector<T,2> center1( .75 * charL, .75 * charL );
  Vector<T,2> center2( .95 * charL, .6 * charL);
  Vector<T,2> center3( .5 * charL, .9 * charL );
  IndicatorEllipse2D<T> ellipse1(center1, .693 * charL, .548 * charL, -20);
  IndicatorEllipse2D<T> ellipse2(center2, .1 * charL, .134 * charL);
  IndicatorEllipse2D<T> ellipse3(center3, .187 * charL, .1 * charL);

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  /// Instantiation of a superGeometry
  SuperGeometry<T, 2> superGeometry(cuboidDecomposition, loadBalancer);
  prepareGeometry(superGeometry, converter, ellipse1, ellipse2, ellipse3);

  std::vector<T> allOmegas = {0., 0., 0., 0., 0., 0.};
  prepareOmegas(converter, theta, allOmegas);

  SuperLattice<T, LDESCRIPTOR> lLattice(converter, superGeometry);
  prepareLattice(lLattice, superGeometry, converter, allOmegas, theta, ellipse1, ellipse2, ellipse3);

  int maxIt = simTime / converter.getConversionFactorTime();

  clout << "Awaiting " << maxIt << " Time steps. Starting simulation..." << std::endl;
  util::Timer<T> timer(converter.getLatticeTime(simTime),
                       superGeometry.getStatistics().getNvoxel());
  timer.start();

  int numDataPoints = 400.;
  std::vector<T> errors = {};
  std::ofstream fout;
  std::string dataFile = singleton::directories().getLogOutDir() + std::to_string(N) +  "_Err.dat";
  fout.open(dataFile.c_str(), std::ios::trunc);
  fout << "N;it;l2UErr;lInfUErr;l2StressErr;lInfStressErr" << std::endl;

  for (int iT = 0; iT <= maxIt; ++iT) {

    timer.update(iT);

    if (iT % (maxIt / numDataPoints) == 0) {
      errors = getResults(lLattice, converter, superGeometry, iT, maxIt);
      timer.printStep();
      fout << N << ";" << converter.getPhysTime(iT) << ";" << errors[0] << ";" << errors[1] << ";" <<errors[2] << ";" << errors[3] << std::endl;
    }
    lLattice.collideAndStream();
  }
  fout.close();
  timer.stop();
}