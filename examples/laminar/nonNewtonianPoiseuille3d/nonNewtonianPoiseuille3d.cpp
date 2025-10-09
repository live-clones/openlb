/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Shota Ito
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

/* 3D force-driven Poiseuille flow with different non-Newtonian viscosity models
 * Implemented options are Newtonian, PowerLaw, Casson, and Carreau-Yasuda models.
 * The first three models are validated for the same pressure drop and dynamic
 * viscosity used for the unit-conversion. The pressure drop is computed according
 * the analytical solution for the Newtonian case.
 * Carreau-Yasuda model is provided as a popular alternative model but still
 * remains untested yet.
*/

#include <olb.h>

using namespace olb;
using namespace olb::names;
using MyCase = Case<NavierStokes, Lattice<double, descriptors::D3Q19<descriptors::FORCE,descriptors::OMEGA>>>;
enum class ViscosityModel: int {
  NEWTONIAN       = 0,
  POWER_LAW       = 1,
  CASSON          = 2,
  CARREAU_YASUDA  = 3
};

namespace olb::parameters {

  struct DYNAMIC_VISCOSITY  : public descriptors::FIELD_BASE<1> {};
  struct POWER_LAW_EXPONENT : public descriptors::FIELD_BASE<1> {};
  struct K0                 : public descriptors::FIELD_BASE<1> {};
  struct K1                 : public descriptors::FIELD_BASE<1> {};
  struct N_CY               : public descriptors::FIELD_BASE<1> {};
  struct MODEL_CONSTANT_A   : public descriptors::FIELD_BASE<1> {};
  struct CHAR_TIME_CONSTANT : public descriptors::FIELD_BASE<1> {};
  struct MU_ZERO            : public descriptors::FIELD_BASE<1> {};
  struct MU_INF             : public descriptors::FIELD_BASE<1> {};
  struct N_PL               : public descriptors::FIELD_BASE<1> {};
  struct CONSISTENCY_INDEX  : public descriptors::FIELD_BASE<1> {};
  struct PRESSURE_DROP      : public descriptors::FIELD_BASE<1> {};
  struct LENGTH             : public descriptors::FIELD_BASE<1> {};
  struct DIAMETER           : public descriptors::FIELD_BASE<1> {};
  struct RADIUS             : public descriptors::FIELD_BASE<1> {};
  struct EOC                : public descriptors::TYPED_FIELD_BASE<bool,1> {};
  struct VISCOSITY_MODEL    : public descriptors::TYPED_FIELD_BASE<ViscosityModel,1> {};
  struct AXIS               : public descriptors::FIELD_BASE<0,1> {};

} // namespace olb::parameters

Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters)
{
  using T                   = MyCase::value_t;
  const T       physDeltaX  = parameters.get<parameters::PHYS_DELTA_X>();
  const T       length      = parameters.get<parameters::LENGTH>();
  const T       radius      = parameters.get<parameters::RADIUS>();
  
  Vector<T, 3> center0(T(0), radius, radius);
  Vector<T, 3> center1(length + 0.5 * physDeltaX, radius, radius);
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  IndicatorLayer3D<T> extendedDomain(pipe, physDeltaX);

  Mesh<T, MyCase::d> mesh(extendedDomain, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

struct Newtonian
{
  using T           = MyCase::value_t;
  using DESCRIPTOR  = MyCase::descriptor_t_of<NavierStokes>;
  using dynamics    = ForcedBGKdynamics<T,DESCRIPTOR>;

private:
  T _Re, _viscosity, _diameter;
public:

  Newtonian( MyCase& myCase ) {
    auto& parameters = myCase.getParameters();
    _Re = parameters.get<parameters::REYNOLDS>();
  }

  static T charVelocity() {
    return  _Re * _viscosity / _diameter;
  }

  static std::shared_ptr<AnalyticalF3D<T,T>> getLatticeVelocityProfile(UnitConverter<T,DESCRIPTOR> const& converter) {
    return std::make_shared<CirclePoiseuille3D<T>>(origin, axis, converter.getCharLatticeVelocity(), radius);
  }

  static std::shared_ptr<AnalyticalF3D<T,T>> getPhysVelocityProfile(UnitConverter<T,DESCRIPTOR> const& converter) {
    return std::make_shared<CirclePoiseuille3D<T>>(origin, axis, converter.getCharPhysVelocity(), radius);
  }

  static void setModelParameter(UnitConverter<T,DESCRIPTOR> const& converter, SuperLattice<T,DESCRIPTOR>& lattice){};
};

struct PowerLaw
{
  using T           = MyCase::value_t;
  using DESCRIPTOR  = MyCase::descriptor_t_of<NavierStokes>;
  using dynamics    = PowerLawForcedBGKdynamics<T,DESCRIPTOR>;

  static T charVelocity() {
    return ((n)/(n+1))*(util::pow(pressureDrop/(2.*consistencyIndex),1./n))*util::pow(radius,(n+1)/n);
  }

  static std::shared_ptr<AnalyticalF3D<T,T>> getLatticeVelocityProfile(UnitConverter<T,DESCRIPTOR> const& converter) {
    return std::make_shared<CirclePowerLaw3D<T>>(origin, axis, converter.getCharLatticeVelocity(), radius, n);
  }

  static std::shared_ptr<AnalyticalF3D<T,T>> getPhysVelocityProfile(UnitConverter<T,DESCRIPTOR> const& converter) {
    return std::make_shared<CirclePowerLaw3D<T>>(origin, axis, converter.getCharPhysVelocity(), radius, n);
  }

  static void setModelParameter(UnitConverter<T,DESCRIPTOR> const& converter, SuperLattice<T,DESCRIPTOR>& lattice){
    T conversionConsistency = util::pow(converter.getPhysDeltaT(),n-2)*util::pow(converter.getPhysDeltaX(),2.);
    lattice.setParameter<powerlaw::M>(consistencyIndex / (conversionConsistency * physRho));
    lattice.setParameter<powerlaw::N>(n);
    const T nuMin = 2.9686e-3;
    const T nuMax = 3.1667;
    lattice.setParameter<powerlaw::OMEGA_MIN>(1./(nuMax*descriptors::invCs2<T,DESCRIPTOR>() + 0.5));
    lattice.setParameter<powerlaw::OMEGA_MAX>(1./(nuMin*descriptors::invCs2<T,DESCRIPTOR>() + 0.5));
  };
};

struct Casson
{
  using T           = MyCase::value_t;
  using DESCRIPTOR  = MyCase::descriptor_t_of<NavierStokes>;
  using dynamics    = CassonForcedBGKdynamics<T,DESCRIPTOR>;

  static T charVelocity() {
    return  Re * viscosity / diameter;
  }

  static std::shared_ptr<AnalyticalF3D<T,T>> getLatticeVelocityProfile(UnitConverter<T,DESCRIPTOR> const& converter) {
    return std::make_shared<CircleCasson3D<T>>(origin, axis, radius, k1*k1, pressureDrop, k0*k0, 1.0/converter.getConversionFactorVelocity());
  }

  static std::shared_ptr<AnalyticalF3D<T,T>> getPhysVelocityProfile(UnitConverter<T,DESCRIPTOR> const& converter) {
    return std::make_shared<CircleCasson3D<T>>(origin, axis, radius, k1*k1, pressureDrop, k0*k0, 1.0);
  }

  static void setModelParameter(UnitConverter<T,DESCRIPTOR> const& converter, SuperLattice<T,DESCRIPTOR>& lattice) {
    const T conversionK1 = util::sqrt( converter.getConversionFactorViscosity());
    const T conversionK0 = conversionK1 / util::sqrt(converter.getConversionFactorTime());
    lattice.setParameter<visco::K_ZERO>(k0 / (util::sqrt(physRho) * conversionK0));
    lattice.setParameter<visco::K_ONE>(k1 / (util::sqrt(physRho) * conversionK1));
  };
};

struct CarreauYasuda
{
  using T           = MyCase::value_t;
  using DESCRIPTOR  = MyCase::descriptor_t_of<NavierStokes>;
  using dynamics    = CarreauYasudaForcedBGKdynamics<T,DESCRIPTOR>;

  static T charVelocity() {
    return  Re * mu_zero / (diameter * physRho);
  }

  static std::shared_ptr<AnalyticalF3D<T,T>> getLatticeVelocityProfile(UnitConverter<T,DESCRIPTOR> const& converter) {
    return std::make_shared<CirclePowerLaw3D<T>>(origin, axis, converter.getCharLatticeVelocity(), radius, n_PL);
  }

  static std::shared_ptr<AnalyticalF3D<T,T>> getPhysVelocityProfile(UnitConverter<T,DESCRIPTOR> const& converter) {
    return std::make_shared<CirclePowerLaw3D<T>>(origin, axis, converter.getCharPhysVelocity(), radius, n_PL);
  }

  static void setModelParameter(UnitConverter<T,DESCRIPTOR> const& converter, SuperLattice<T,DESCRIPTOR>& lattice) {
    lattice.setParameter<visco::N>(n_CY);
    lattice.setParameter<visco::A>(a);
    lattice.setParameter<visco::LAMBDA>(lambda / converter.getConversionFactorTime());
    const T conversionMU = converter.getConversionFactorViscosity();
    lattice.setParameter<visco::MU_ZERO>(mu_zero / (physRho * conversionMU));
    lattice.setParameter<visco::MU_INF>(mu_inf / (physRho * conversionMU));
  }
};                          // PowerLaw index used for comparison


// template<typename MODEL>
// using DYNAMICS = MODEL::dynamics;

void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;
  using T             = MyCase::value_t;
  auto&   geometry    = myCase.getGeometry();
  auto&   parameters  = myCase.getParameters();
  const T physDeltaX  = parameters.get<parameters::PHYS_DELTA_X>();
  const T length      = parameters.get<parameters::LENGTH>();
  const T radius      = parameters.get<parameters::RADIUS>();

  Vector<T, 3> center0(-physDeltaX * 0.2, radius, radius);
  Vector<T, 3> center1(length, radius, radius);
  center0[0] -= 3.*physDeltaX;
  center1[0] += 3.*physDeltaX;
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  geometry.rename(0, 2);
  geometry.rename(2, 1, pipe);
  geometry.clean();
  geometry.innerClean();
  geometry.checkForErrors();
  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

template<typename MODEL>
void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;
  using T           = MyCase::value_t;
  auto& geometry    = myCase.getGeometry();
  auto& lattice     = myCase.getLattice(NavierStokes {});
  auto& parameters  = myCase.getParameters();
  auto& converter   = lattice.getUnitConverter();
  const T length    = parameters.get<parameters::LENGTH>();
  const T radius    = parameters.get<parameters::RADIUS>();

  lattice.defineDynamics<DYNAMICS<MODEL>>(geometry.getMaterialIndicator({1}));

  Vector<T, 3> center0(0, radius, radius);
  Vector<T, 3> center1(length, radius, radius);
  center0[0] -= 0.5*converter.getPhysDeltaX();
  center1[0] += 0.5*converter.getPhysDeltaX();
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  setBouzidiBoundary<T, DESCRIPTOR, BouzidiPostProcessor>(lattice, geometry, 2, pipe);

  Vector<T,3> forceTerm = {0,0,0};
  T conversionF_force = converter.getConversionFactorLength() / util::pow(converter.getConversionFactorTime(), 2.);
  forceTerm[0] = pressureDrop / physRho / conversionF_force;

  AnalyticalConst3D<T,T> force (forceTerm);
  lattice.defineField<FORCE>(geometry, 1, force);
  lattice.defineField<FORCE>(geometry, 2, force);

  clout << "Prepare Lattice ... OK" << std::endl;
}

template<typename MODEL>
void setInitialValues(MyCase& myCase) {
  OstreamManager clout(std::cout, "setInitialValues");
  clout << "Prepare setInitialValues ..." << std::endl;
  using T           = MyCase::value_t;
  auto& geometry    = myCase.getGeometry();
  auto& lattice     = myCase.getLattice(NavierStokes {});
  auto& parameters  = myCase.getParameters();
  auto& converter   = lattice.getUnitConverter();  

  AnalyticalConst3D<T,T> rho(1.);
  std::shared_ptr<AnalyticalF3D<T,T>> profileU = MODEL::getLatticeVelocityProfile(converter);

  lattice.defineRhoU(geometry, 1, rho, *profileU);
  lattice.iniEquilibrium(geometry, 1, rho, *profileU);
  lattice.defineRhoU(geometry, 2, rho, *profileU);
  lattice.iniEquilibrium(geometry, 2, rho,*profileU);

  lattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  MODEL::setModelParameter(converter, lattice);

  lattice.initialize();

  clout << "Prepare setInitialValues ... OK" << std::endl;
}

// Compute error norms
template<typename MODEL, typename T>
Vector<T,3> error( MyCase& myCase)
{
  OstreamManager clout( std::cout,"error" );
  auto& sGeometry   = myCase.getGeometry();
  auto& sLattice    = myCase.getLattice(NavierStokes {});
  auto& parameters  = myCase.getParameters();
  auto& converter   = lattice.getUnitConverter();
  std::vector<T> errors;
  int tmp[]= { };
  T result[2]= { };

  lattice.setProcessingContext(ProcessingContext::Evaluation);
  std::shared_ptr<AnalyticalF3D<T,T>> uSol = MODEL::getPhysVelocityProfile(converter);
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> u( lattice,converter );
  auto indicatorF = geometry.getMaterialIndicator(1);

  SuperAbsoluteErrorL1Norm3D<T> absVelocityErrorNormL1(u, *uSol, indicatorF);
  absVelocityErrorNormL1(result, tmp);
  clout << "velocity-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm3D<T> relVelocityErrorNormL1(u, *uSol, indicatorF);
  relVelocityErrorNormL1(result, tmp);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;
  errors.push_back(result[0]);

  SuperAbsoluteErrorL2Norm3D<T> absVelocityErrorNormL2(u, *uSol, indicatorF);
  absVelocityErrorNormL2(result, tmp);
  clout << "velocity-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm3D<T> relVelocityErrorNormL2(u, *uSol, indicatorF);
  relVelocityErrorNormL2(result, tmp);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;
  errors.push_back(result[0]);

  SuperAbsoluteErrorLinfNorm3D<T> absVelocityErrorNormLinf(u, *uSol, indicatorF);
  absVelocityErrorNormLinf(result, tmp);
  clout << "velocity-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm3D<T> relVelocityErrorNormLinf(u, *uSol, indicatorF);
  relVelocityErrorNormLinf(result, tmp);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;
  errors.push_back(result[0]);

  return errors;
}

// Output to console and files
template<typename MODEL>
void getResults(MyCase& myCase, std::size_t iT, util::Timer<MyCase::value_t>& timer, bool converged)
{
  OstreamManager clout( std::cout,"getResults" );
  using T               = MyCase::value_t;
  auto&       lattice   = myCase.getLattice(NavierStokes {});
  const auto& converter = sLattice.getUnitConverter();
  auto&       geometry  = myCase.getGeometry();
  const bool lastTimeStep = (iT + 1 == converter.getLatticeTime( maxPhysT ));
  const int statIter = converter.getLatticeTime( maxPhysT / 10. );

  SuperVTMwriter3D<T> vtmWriter( "nonNewtonianPoiseuille3d" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( lattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( lattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  std::shared_ptr<AnalyticalF3D<T,T>> uSol = MODEL::getPhysVelocityProfile(converter);
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> analyticalVelocityLattice(*uSol, lattice);
  analyticalVelocityLattice.getName() = "analytical solution";
  vtmWriter.addFunctor(analyticalVelocityLattice);

  SuperLatticeField3D<T,DESCRIPTOR,OMEGA> omega(lattice);
  omega.getName() = "omega";
  vtmWriter.addFunctor(omega);

  const int vtmIter  = converter.getLatticeTime( maxPhysT/100. );

  if ( iT==0 ) {
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( lattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( lattice );

    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }

  // Writes the vtm files and profile text file
  if ( iT%vtmIter==0 || lastTimeStep ) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write( iT );
  }

  // Writes output on the console
  if ( iT%statIter==0 || lastTimeStep ) {
    timer.update( iT );
    timer.printStep();
    lattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
    error<MODEL>( geometry, lattice);
  }
}

template<typename MODEL, typename T>
std::vector<T> simulatePoiseuilleWith()
{
  OstreamManager clout( std::cout,"simulatePoiseuille" );
  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    (T)   N,
    (T)   tau,
    (T)   diameter,
    (T)   MODEL::charVelocity(),
    (T)   viscosity,
    (T)   physRho
  );
  converter.print();
  converter.write("nonNewtonianPoiseuille3d");

  // === 2nd Step: Prepare Geometry ===

  Vector<T, 3> center0(0, radius, radius);
  Vector<T, 3> center1(length + 0.5 * converter.getPhysDeltaX(), radius, radius);
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  IndicatorLayer3D<T> extendedDomain(pipe, converter.getPhysDeltaX());

  // Instantiation of a cuboidDecomposition with weights
  CuboidDecomposition3D<T> cuboidDecomposition( extendedDomain, converter.getPhysDeltaX(), noOfCuboids);

  cuboidDecomposition.setPeriodicity({true, false, false});

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  // Instantiation of a geometry
  const int overlap = 3;
  SuperGeometry<T,3> geometry(cuboidDecomposition, loadBalancer, overlap);

  prepareGeometry(geometry, converter.getPhysDeltaX());

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T,DESCRIPTOR> lattice(converter, geometry);

  //prepareLattice and setBoundaryConditions
  prepareLattice<MODEL>(lattice, geometry);

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), geometry.getStatistics().getNvoxel() );
  timer.start();

  std::vector<T> errors;
  for ( std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {
    lattice.collideAndStream();
    getResults<MODEL>( lattice, iT, geometry, timer );
  }
  timer.stop();
  timer.printSummary();

  return error<MODEL>(geometry, lattice);
}

template<typename MODEL>
void eocPoiseuilleWith() {
  OstreamManager clout( std::cout,"eocPoiseuille" );
  std::vector<int> res {41, 61, 81, 101};
  std::vector<std::vector<T>> errors;
  for (int iter = 0; iter < (int) res.size(); ++iter) {
    N = res[iter];
    errors.push_back(simulatePoiseuilleWith<MODEL>());
  }

  std::vector<T> eoc_L1, eoc_L2, eoc_Linf;
  for (int i = 0; i < (int) res.size() - 1; ++i) {
    eoc_L1.push_back( (util::log(errors.at(i).at(0)) - util::log(errors.at(i+1).at(0)))/(util::log(res[i]) - util::log(res[i+1])) );
    eoc_L2.push_back( (util::log(errors.at(i).at(1)) - util::log(errors.at(i+1).at(1)))/(util::log(res[i]) - util::log(res[i+1])) );
    eoc_Linf.push_back( (util::log(errors.at(i).at(2)) - util::log(errors.at(i+1).at(2)))/(util::log(res[i]) - util::log(res[i+1])) );
  }

  for (int i = 0; i < (int) res.size() -1; ++i) {
    clout << "Error L1: " << errors.at(i).at(0) << std::endl;
    clout << "Error L2: " << errors.at(i).at(1) << std::endl;
    clout << "Error Linf: " << errors.at(i).at(2) << std::endl;
  }
  clout << "EOC with N: " << res << std::endl;
  clout << "EOC L1: " << eoc_L1 << std::endl;
  clout << "EOC L2: " << eoc_L2 << std::endl;
  clout << "EOC Linf: " << eoc_Linf << std::endl;
}

int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  initialize( &argc, &argv );
  
  /// === Step 2: Set Parameters ===
  CLIreader args(argc, argv);
  const std::string viscosityModel = args.getValueOrFallback<std::string>("--model", "PowerLaw");
  MyCase::ParametersD myCaseParameters;
  using T = MyCase::value_t;
  {
    using namespace olb::parameters;
    myCaseParameters.set<EOC>(false);
    myCaseParameters.set<parameters::VISCOSITY_MODEL>(ViscosityModel::NEWTONIAN);
    myCaseParameters.set<RESOLUTION>(31);
    myCaseParameters.set<PHYS_DELTA_X>([&] {
      return myCaseParameters.get<PHYS_CHAR_LENGTH>() / myCaseParameters.get<RESOLUTION>();
    });
    myCaseParameters.set<DOMAIN_EXTENT>({3e-3,1e-3,1e-3});
    myCaseParameters.set<LENGTH>([&] {
      return myCaseParameters.get<DOMAIN_EXTENT>()[0];
    });
    myCaseParameters.set<DIAMETER>([&] {
      return myCaseParameters.get<DOMAIN_EXTENT>()[1];
    });
    myCaseParameters.set<RADIUS>([&] {
      return myCaseParameters.get<DIAMETER>() / 2.;
    });
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.64265);
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1060);
    myCaseParameters.set<MAX_PHYS_T>(.3);
    myCaseParameters.set<AXIS>({1.,0.,0.});
    myCaseParameters.set<REYNOLDS>(100);
    myCaseParameters.set<DYNAMIC_VISCOSITY>(0.0040755);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>([&] {
      return myCaseParameters.get<DYNAMIC_VISCOSITY>() / myCaseParameters.get<PHYS_CHAR_DENSITY>();
    });
    myCaseParameters.set<PRESSURE_DROP>([&] {
      return (16.*util::pow(myCaseParameters.get<DYNAMIC_VISCOSITY>(),2.)*myCaseParameters.get<REYNOLDS>())/
        (myCaseParameters.get<PHYS_CHAR_DENSITY>()*util::pow(myCaseParameters.get<DOMAIN_EXTENT>()[1],3.));
    });

    // PowerLaw model parameter
    myCaseParameters.set<POWER_LAW_EXPONENT>(0.65);
    myCaseParameters.set<CONSISTENCY_INDEX>([&] {
      T dynamicViscosity  = myCaseParameters.get<DYNAMIC_VISCOSITY>();
      T n                 = myCaseParameters.get<POWER_LAW_EXPONENT>();
      T pressureDrop      = myCaseParameters.get<PRESSURE_DROP>();
      T radius            = myCaseParameters.get<RADIUS>();
      return  util::pow(myCaseParameters.get<DYNAMIC_VISCOSITY>(),myCaseParameters.get<POWER_LAW_EXPONENT>())
              * util::pow((myCaseParameters.get<RADIUS>()*myCaseParameters.get<PRESSURE_DROP>()/8.)
                          * util::pow(myCaseParameters.get<POWER_LAW_EXPONENT>()/(myCaseParameters.get<POWER_LAW_EXPONENT>()+1),myCaseParameters.get<POWER_LAW_EXPONENT>())
                          , 1-myCaseParameters.get<POWER_LAW_EXPONENT>());  
    });  // consistency index in Pa s^n, SI unit

    // Casson model parameter
    myCaseParameters.set<K0>(0.07);  // k0 constant (Pa)^0.5, k0^2 corresponds to yield stress
    myCaseParameters.set<K1>([&] {
      return  util::sqrt(myCaseParameters.get<DYNAMIC_VISCOSITY>());  
    });  // k1 constant (Pa*s)^0.5, k1^2 corresponds to Casson viscosity

    // CarreauYasuda model parameter
    myCaseParameters.set<N_CY>(0.9);  // CY fluid index
    myCaseParameters.set<MODEL_CONSTANT_A>(1.5);  // model constant
    myCaseParameters.set<CHAR_TIME_CONSTANT>(3.313);  // characteristic time constant lambda (s)
    myCaseParameters.set<MU_ZERO>([&] {
      return myCaseParameters.get<DYNAMIC_VISCOSITY>()*1.45;  
    });  // zero viscosity (Pa*s)
    myCaseParameters.set<MU_INF>([&] {
      return myCaseParameters.get<DYNAMIC_VISCOSITY>()*10.;  
    });  // infinity viscosity (Pa*s)
    myCaseParameters.set<N_PL>(0.708);  // PL fluid index
  }

  myCaseParameters.fromCLI(argc, argv);

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);
  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  bool eoc = myCaseParameters.get<parameters::EOC>();
  using MODEL = Newtonian;
  switch ( myCaseParameters.get<parameters::VISCOSITY_MODEL>() ) {
    case ViscosityModel::NEWTONIAN:
    default:
      break;
    case ViscosityModel::POWER_LAW:
      using MODEL = PowerLaw;
      break;
    case ViscosityModel::CASSON:
      using MODEL = Casson;
      break;
    case ViscosityModel::CARREAU_YASUDA:
      using MODEL = CarreauYasuda;
      break;
  }

  /// === Step 6: Prepare Lattice ===
  prepareLattice<MODEL>(myCase);

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues(myCase);

  if (!eoc) {
    if (viscosityModel == "Newtonian") {
      simulatePoiseuilleWith<Newtonian>();
    } else if (viscosityModel == "PowerLaw") {
      simulatePoiseuilleWith<PowerLaw>();
    } else if (viscosityModel == "Casson") {
      simulatePoiseuilleWith<Casson>();
    } else if (viscosityModel == "CarreauYasuda") {
      simulatePoiseuilleWith<CarreauYasuda>();
    } else {
      throw std::runtime_error(viscosityModel + " is not a valid viscosity model");
    }
  } else {
    if (viscosityModel == "Newtonian") {
      eocPoiseuilleWith<Newtonian>();
    } else if (viscosityModel == "PowerLaw") {
      eocPoiseuilleWith<PowerLaw>();
    } else if (viscosityModel == "Casson") {
      eocPoiseuilleWith<Casson>();
    } else if (viscosityModel == "CarreauYasuda") {
      eocPoiseuilleWith<CarreauYasuda>();
    } else {
      throw std::runtime_error(viscosityModel + " is not a valid viscosity model");
    }
  }
}
