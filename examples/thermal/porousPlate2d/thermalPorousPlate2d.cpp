/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2008 Orestis Malaspinas
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

/* rayleighBenard2d.cpp:
 * Rayleigh-Benard convection rolls in 2D, simulated with
 * the thermal LB model by Z. Guo e.a., between a hot plate at
 * the bottom and a cold plate at the top.
 */


#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::descriptors;

namespace olb::parameters {
  struct BC_ID  : public descriptors::FIELD_BASE<1> { };
}

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D2Q9<descriptors::FORCE>>,
  AdvectionDiffusion,  Lattice<double, descriptors::D2Q5<descriptors::VELOCITY>>
>;

//#define TemperatureBoundary
//#define RegularizedTemperatureBoundary
// #define RegularizedHeatFluxBoundary

enum class TemperatureBoundary {
  TemperatureDirichlet = 1,
  RegularizedTemperature = 2,
  RegularizedHeatFlux = 3
};

// analytical solution from point light source in infinte domain
// appliacation from R3 to R1.
// effective for x in R3, only the distance to (0,0) is needed.
// documentation e.g. Biomedical Optics, Lihong V. Wang Hsin-I Wu
template <typename T, typename S>
class AnalyticalVelocityPorousPlate2D : public AnalyticalF2D<T, S> {
private:
  T _Re;
  T _u0;
  T _v0;
  T _ly;
public:
  AnalyticalVelocityPorousPlate2D(T Re, T u0, T v0, T ly) : AnalyticalF2D<T, S>(2),
    _Re(Re), _u0(u0), _v0(v0), _ly(ly)
  {
    this->getName() = "AnalyticalVelocityPorousPlate2D";
  };

  bool operator()(T output[2], const S x[2]) override
  {
    output[0] = _u0*((util::exp(_Re* x[1] / _ly) - 1) / (util::exp(_Re) - 1));
    output[1] = _v0;
    return true;
  };
};

template <typename T, typename S>
class AnalyticalTemperaturePorousPlate2D : public AnalyticalF2D<T, S> {
private:
  T _Re;
  T _Pr;
  T _ly;
  T _T0;
  T _deltaT;
public:
  AnalyticalTemperaturePorousPlate2D(T Re, T Pr, T ly, T T0, T deltaT) : AnalyticalF2D<T, S>(1),
    _Re(Re), _Pr(Pr), _ly(ly), _T0(T0), _deltaT(deltaT)
  {
    this->getName() = "AnalyticalTemperaturePorousPlate2D";
  };

  bool operator()(T output[1], const S x[2]) override
  {
    output[0] = _T0 + _deltaT*((util::exp(_Pr*_Re*x[1] / _ly) - 1) / (util::exp(_Pr*_Re) - 1));
    return true;
  };
};

template <typename T, typename S>
class AnalyticalHeatFluxPorousPlate2D : public AnalyticalF2D<T, S> {
private:
  T _Re;
  T _Pr;
  T _deltaT;
  T _ly;
  T _lambda;
public:
  AnalyticalHeatFluxPorousPlate2D(T Re, T Pr, T deltaT, T ly,T lambda) : AnalyticalF2D<T, S>(2),
    _Re(Re), _Pr(Pr), _deltaT(deltaT), _ly(ly), _lambda(lambda)
  {
    this->getName() = "AnalyticalHeatFluxPorousPlate2D";
  };

  bool operator()(T output[2], const S x[2]) override
  {
    output[0] = 0;
    output[1] = - _lambda * _Re * _Pr * _deltaT / _ly * (util::exp(_Pr * _Re * x[1] / _ly))/(util::exp(_Pr * _Re) - 1);
    return true;
  };
};

void error(MyCase& myCase)
{
  OstreamManager clout( std::cout, "error" );
  using T = MyCase::value_t_of<NavierStokes>;
  using NSEDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ADEDESCRIPTOR = MyCase::descriptor_t_of<AdvectionDiffusion>;
  auto& geometry = myCase.getGeometry();
  auto& NSElattice = myCase.getLattice(NavierStokes{});
  auto& ADElattice = myCase.getLattice(AdvectionDiffusion{});
  const auto& converter = NSElattice.getUnitConverter();
  const T Re = converter.getReynoldsNumber();

  int input[1] = { };
  T result[1]  = { };

  auto indicatorF = geometry.getMaterialIndicator({1, 2, 3});

  T u_Re = Re * converter.getPhysViscosity() / converter.getCharPhysLength();
  AnalyticalVelocityPorousPlate2D<T,T> uSol(Re, converter.getCharPhysVelocity(), u_Re, converter.getCharPhysLength());
  SuperLatticePhysVelocity2D<T,NSEDESCRIPTOR> u(NSElattice,converter);

  SuperAbsoluteErrorL2Norm2D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, input);
  clout << "velocity-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm2D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, input);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;

  AnalyticalTemperaturePorousPlate2D<T,T> tSol(Re, Pr, converter.getCharPhysLength(), converter.getCharPhysLowTemperature(), converter.getCharPhysTemperatureDifference());
  SuperLatticePhysTemperature2D<T,NSEDESCRIPTOR,ADEDESCRIPTOR> t(ADElattice,converter);

  SuperAbsoluteErrorL2Norm2D<T> absTemperatureErrorNormL2(t, tSol, indicatorF);
  absTemperatureErrorNormL2(result, input);
  clout << "temperature-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm2D<T> relTemperatureErrorNormL2(t, tSol, indicatorF);
  relTemperatureErrorNormL2(result, input);
  clout << "; temperature-L2-error(rel)=" << result[0] << std::endl;

  AnalyticalHeatFluxPorousPlate2D<T,T> HeatFluxSol(Re, Pr, converter.getCharPhysTemperatureDifference(), converter.getCharPhysLength(), converter.getThermalConductivity());
  SuperLatticePhysHeatFlux2D<T,NSEDESCRIPTOR,ADEDESCRIPTOR> HeatFlux(ADElattice,converter);

  SuperAbsoluteErrorL2Norm2D<T> absHeatFluxErrorNormL2(HeatFlux, HeatFluxSol, indicatorF);
  absHeatFluxErrorNormL2(result, input);
  clout << "heatFlux-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm2D<T> relHeatFluxErrorNormL2(HeatFlux, HeatFluxSol, indicatorF);
  relHeatFluxErrorNormL2(result, input);
  clout << "; heatFlux-L2-error(rel)=" << result[0] << std::endl;
}

/// @brief Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t_of<NavierStokes>;
  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const T physDeltaX  = parameters.get<parameters::PHYS_DELTA_X>();
  const Vector origin{0., 0.};
  IndicatorCuboid2D<T> cuboid(extent, origin);

  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({true,false});
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t_of<NavierStokes>;
  auto& geometry  = myCase.getGeometry();
  auto& parameters    = myCase.getParameters();

  geometry.rename(0,2);
  geometry.rename(2,1,{0,1});

  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const T physDeltaX  = parameters.get<parameters::PHYS_DELTA_X>();
  const Vector origin{-physDeltaX, 0.};
  IndicatorCuboid2D<T> bottom(extent, origin);

  geometry.rename(2,3,1,bottom);

  geometry.clean();
  geometry.innerClean();
  geometry.checkForErrors();

  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase) {
  OstreamManager clout(std::cout,"prepareLattice");

  using T = MyCase::value_t_of<NavierStokes>;
  auto& geometry      = myCase.getGeometry();
  auto& parameters        = myCase.getParameters();
  auto& NSElattice    = myCase.getLattice(NavierStokes{});
  auto& ADElattice    = myCase.getLattice(AdvectionDiffusion{});
  using NSEDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ADEDESCRIPTOR = MyCase::descriptor_t_of<AdvectionDiffusion>;

  const T physCharLength          = parameters.get<parameters::PHYS_CHAR_LENGTH>();
  const T latticeCharVelocity     = parameters.get<parameters::LATTICE_CHAR_VELOCITY>();
  const int N                     = parameters.get<parameters::RESOLUTION>();
  const T physDeltaX              = parameters.get<parameters::PHYS_DELTA_X>();
  const T physViscosity           = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T physDensity             = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T physThermalConductivity = parameters.get<parameters::PHYS_THERMAL_CONDUCTIVITY>();
  const T tau                     = parameters.get<parameters::RELAXATION_TIME>();
  const T Ra                      = parameters.get<parameters::RAYLEIGH>();
  const T Pr                      = parameters.get<parameters::PRANDTL>();
  const T Tcold                   = parameters.get<parameters::T_COLD>();
  const T Thot                    = parameters.get<parameters::T_HOT>();

  NSElattice.setUnitConverter<ThermalUnitConverter<T,NSEDESCRIPTOR,ADEDESCRIPTOR>>(
    (T) physDeltaX, // physDeltaX
    (T) physDeltaX * physDensity / physViscosity * (tau - 0.5) / 3 / N, // physDeltaT
    (T) physCharLength, // charPhysLength
    (T) util::sqrt( 9.81 * Ra * physViscosity * physViscosity / Pr / 9.81 / (Thot - Tcold) / util::pow(physDensity, 3) * (Thot - Tcold) * physDensity ), // charPhysVelocity
    (T) physViscosity, // physViscosity
    (T) physDensity, // physDensity
    (T) physThermalConductivity, // physThermalConductivity
    (T) Pr * physThermalConductivity / physViscosity / physDensity, // physSpecificHeatCapacity
    (T) Ra * physViscosity * physViscosity / Pr / 9.81 / (Thot - Tcold) / util::pow(physDensity, 3), // physThermalExpansionCoefficient
    (T) Tcold, // charPhysLowTemperature
    (T) Thot // charPhysHighTemperature
  );
  const auto& NSEconverter = NSElattice.getUnitConverter();
  NSEconverter.print();

  ADElattice.setUnitConverter(NSEconverter);
  const auto& ADEconverter = ADElattice.getUnitConverter();

  T ADEomega  = NSEconverter.getLatticeThermalRelaxationFrequency();
  T NSEomega  = ADEconverter.getLatticeRelaxationFrequency();

  ADElattice.defineDynamics<AdvectionDiffusionBGKdynamics>(geometry, 1);
  ADElattice.defineDynamics<AdvectionDiffusionBGKdynamics>(geometry, 2);
  ADElattice.defineDynamics<AdvectionDiffusionBGKdynamics>(geometry, 3);
  NSElattice.defineDynamics<ForcedBGKdynamics>(geometry, 1);
  NSElattice.defineDynamics<ForcedBGKdynamics>(geometry, 2);
  NSElattice.defineDynamics<ForcedBGKdynamics>(geometry, 3);

  /// sets boundary
  boundary::set<boundary::LocalVelocity>(NSElattice, geometry, 2);
  boundary::set<boundary::LocalVelocity>(NSElattice, geometry, 3);
  const int bcId = parameters.get<parameters::BC_ID>();
  TemperatureBoundary bc = static_cast<TemperatureBoundary>(bcId);

  switch (bc):
    case: TemperatureBoundary::TemperatureDirichlet {
      boundary::set<boundary::AdvectionDiffusionDirichlet>(ADElattice, geometry, 2);
      boundary::set<boundary::AdvectionDiffusionDirichlet>(ADElattice, geometry, 3);
      break;
    }
    case: TemperatureBoundary::RegularizedTemperature {
      boundary::set<boundary::RegularizedTemperature>(ADElattice, geometry.getMaterialIndicator(2));
      boundary::set<boundary::RegularizedTemperature>(ADElattice, geometry.getMaterialIndicator(3));
      break;
    }
    case: TemperatureBoundary::RegularizedHeatFlux {
      T heatFlux[2];
      T input[2] = {0.,1.};
      AnalyticalHeatFluxPorousPlate2D<T,T> HeatFluxSol(Re, Pr, ADEconverter.getCharPhysTemperatureDifference(), ADEconverter.getCharPhysLength(), ADEconverter.getThermalConductivity());
      HeatFluxSol(heatFlux, input);
      T temp = ADEconverter.getLatticeSpecificHeatCapacity(ADEconverter.getPhysSpecificHeatCapacity())*(ADEconverter.getLatticeThermalRelaxationTime() - 0.5) / ADEconverter.getLatticeThermalRelaxationTime();
      heatFlux[0] = ADEconverter.getLatticeHeatFlux(heatFlux[0]) / temp;
      heatFlux[1] = ADEconverter.getLatticeHeatFlux(heatFlux[1]) / temp;
      AnalyticalConst2D<T,T> heatFluxC(heatFlux[0], heatFlux[1]);
      boundary::set<boundary::RegularizedHeatFlux>(ADElattice, geometry, 2);
      ADElattice.defineU(geometry, 2, heatFluxC);
      boundary::set<boundary::RegularizedTemperature>(ADElattice, geometry, 3);
      ADElattice.setParameter<descriptors::OMEGA>(ADEomega);
      break;
    }

  ADElattice.setParameter<descriptors::OMEGA>(ADEomega);
  NSElattice.setParameter<descriptors::OMEGA>(NSEomega);

  T boussinesqForcePrefactor = 9.81 / ADEconverter.getConversionFactorVelocity() * ADEconverter.getConversionFactorTime() *
                               ADEconverter.getCharPhysTemperatureDifference() * ADEconverter.getPhysThermalExpansionCoefficient();

  auto& coupling = myCase.setCouplingOperator(
    "NavierStokesAdvectionDiffusionCoupling",
    NavierStokesAdvectionDiffusionCoupling{},
    names::NavierStokes{}, NSElattice,
    names::AdvectionDiffusion{}, ADElattice
  );

  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::T0>(
    ADEconverter.getLatticeTemperature(Tcold)
  );
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::FORCE_PREFACTOR>(
    boussinesqForcePrefactor * Vector<T,2>{0.0,1.0}
  );

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t_of<NavierStokes>;
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  const T Tcold = parameters.get<parameters::T_COLD>();
  const T Thot = parameters.get<parameters::T_HOT>();

  auto& NSElattice = myCase.getLattice(NavierStokes{});
  auto& NSEconverter = NSElattice.getUnitConverter();
  auto& ADElattice = myCase.getLattice(AdvectionDiffusion{});
  auto& ADEconverter = ADElattice.getUnitConverter();

  /// for each material set the defineRhoU and the Equilibrium
  std::vector<T> zero(2,T());
  AnalyticalConst2D<T,T> u(zero);
  AnalyticalConst2D<T,T> rho(1.);
  AnalyticalConst2D<T,T> force(zero);

  T u_Re = NSEconverter.getLatticeVelocity( NSEconverter.getReynoldsNumber() * NSEconverter.getPhysViscosity() / NSEconverter.getCharPhysLength() );

  AnalyticalConst2D<T,T> u_top(NSEconverter.getCharLatticeVelocity(), u_Re);
  AnalyticalConst2D<T,T> u_bot(0.0, u_Re);

  NSElattice.defineRhoU(geometry, 1, rho, u);
  NSElattice.iniEquilibrium(geometry, 1, rho, u);
  NSElattice.defineField<FORCE>(geometry, 1, force);
  NSElattice.defineRhoU(geometry, 2, rho, u_top);
  NSElattice.iniEquilibrium(geometry, 2, rho, u_top);
  NSElattice.defineField<FORCE>(geometry, 2, force);
  NSElattice.defineRhoU(geometry, 3, rho, u_bot);
  NSElattice.iniEquilibrium(geometry, 3, rho, u_bot);
  NSElattice.defineField<FORCE>(geometry, 3, force);

  AnalyticalConst2D<T,T> Cold(ADEconverter.getLatticeTemperature(Tcold));
  AnalyticalConst2D<T,T> Hot(ADEconverter.getLatticeTemperature(Thot));

  ADElattice.defineRho(geometry, 1, Cold);
  ADElattice.iniEquilibrium(geometry, 1, Cold, u);
  ADElattice.defineField<descriptors::VELOCITY>(geometry, 1, u);

  ADElattice.defineField<descriptors::VELOCITY>(geometry, 2, u_top);
  ADElattice.defineRho(geometry, 2, Hot);
  ADElattice.iniEquilibrium(geometry, 2, Hot, u_top);
  ADElattice.defineField<descriptors::VELOCITY>(geometry, 2, u);

  ADElattice.defineField<descriptors::VELOCITY>(geometry, 3, u_bot);
  ADElattice.defineRho(geometry, 3, Cold);
  ADElattice.iniEquilibrium(geometry, 3, Cold, u_bot);
  ADElattice.defineField<descriptors::VELOCITY>(geometry, 3, u);

  NSElattice.initialize();
  ADElattice.initialize();
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{
  // Nothing to do here
}

void getResults(MyCase& myCase)
{
  OstreamManager clout(std::cout,"getResults");
  using T = MyCase::value_t_of<NavierStokes>;
  using NSEDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ADEDESCRIPTOR = MyCase::descriptor_t_of<AdvectionDiffusion>;
  auto& NSElattice = myCase.getLattice(NavierStokes{});
  auto& ADElattice = myCase.getLattice(AdvectionDiffusion{});
  const auto& converter = NSElattice.getUnitConverter();

  // TODO such times should also be parameters
  const int statIter = converter.getLatticeTime(10.0);
  const int vtkIter = converter.getLatticeTime(100.);

  SuperVTMwriter2D<T> vtkWriter("thermalPorousPlate2d");
  SuperLatticePhysVelocity2D<T, NSEDESCRIPTOR> velocity(NSElattice, converter);
  SuperLatticePhysPressure2D<T, NSEDESCRIPTOR> pressure(NSElattice, converter);
  SuperLatticePhysTemperature2D<T, NSEDESCRIPTOR, ADEDESCRIPTOR> temperature(ADElattice, converter);
  SuperLatticePhysHeatFlux2D<T, NSEDESCRIPTOR, ADEDESCRIPTOR> heatflux(ADElattice, converter);

  T u_Re = Re * converter.getPhysViscosity() / converter.getCharPhysLength();
  AnalyticalVelocityPorousPlate2D<T,T> uSol(Re, converter.getCharPhysVelocity(), u_Re, converter.getCharPhysLength());
  SuperLatticeFfromAnalyticalF2D<T,NSEDESCRIPTOR> uSolLattice(uSol,NSElattice);
  AnalyticalTemperaturePorousPlate2D<T,T> TSol(Re, Pr, converter.getCharPhysLength(), converter.getCharPhysLowTemperature(), converter.getCharPhysTemperatureDifference());
  SuperLatticeFfromAnalyticalF2D<T,ADEDESCRIPTOR> TSolLattice(TSol,ADElattice);
  AnalyticalHeatFluxPorousPlate2D<T,T> HeatFluxSol(Re, Pr, converter.getCharPhysTemperatureDifference(), converter.getCharPhysLength(), converter.getThermalConductivity());
  SuperLatticeFfromAnalyticalF2D<T,ADEDESCRIPTOR> HeatFluxSolLattice(HeatFluxSol,ADElattice);

  vtkWriter.addFunctor( pressure );
  vtkWriter.addFunctor( velocity );
  vtkWriter.addFunctor( temperature );
  vtkWriter.addFunctor( heatflux );
  vtkWriter.addFunctor( uSolLattice );
  vtkWriter.addFunctor( TSolLattice );
  vtkWriter.addFunctor( HeatFluxSolLattice );


  if (iT == 0) {
    /// Writes the converter log file
    // writeLogFile(converter,"thermalPorousPlate2d");
    T tmpIn[2] = {0.,1.};
    T tmpOut[2];
    HeatFluxSol(tmpOut,tmpIn);
    clout << converter.getLatticeHeatFlux(tmpOut[0]) << " " << converter.getLatticeHeatFlux(tmpOut[1]) << std::endl;
    clout << tmpOut[0] << " " << tmpOut[1] << std::endl;

    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, NSEDESCRIPTOR> cuboid(NSElattice);
    SuperLatticeRank2D<T, NSEDESCRIPTOR> rank(NSElattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
  }

  /// Writes the VTK files
  if (iT % vtkIter == 0 || converged) {
    NSElattice.getStatistics().print(iT,converter.getPhysTime(iT));
    ADElattice.setProcessingContext(ProcessingContext::Evaluation);
    NSElattice.setProcessingContext(ProcessingContext::Evaluation);
    timer.print(iT);
    error(geometry, NSElattice, ADElattice, converter, Re);
    vtkWriter.write(iT);

    ///writes Jpeg
    //SuperEuklidNorm2D<T, DESCRIPTOR> normVel(velocity);
    BlockReduction2D2D<T> planeReduction(temperature, N, BlockDataSyncMode::ReduceOnly);
    // write output of velocity as JPEG
    heatmap::plotParam<T> jpeg_Param;
    jpeg_Param.maxValue = Thot;
    jpeg_Param.minValue = Tcold;
    heatmap::write(planeReduction, iT, jpeg_Param);
  }

}

/*
int main(int argc, char *argv[])
{

  /// === 1st Step: Initialization ===
  OstreamManager clout(std::cout,"main");
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

  /// === 2nd Step: Prepare Geometry ===
  std::vector<T> extend(2,T());
  extend[0] = lx;
  extend[1] = ly;
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extend, origin);

  /// Instantiation of a cuboidDecomposition with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidDecomposition2D<T> cuboidDecomposition(cuboid, converter.getPhysDeltaX(), noOfCuboids);
  cuboidDecomposition.setPeriodicity({true, false});

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  /// Instantiation of a geometry
  geometry<T,2> geometry(cuboidDecomposition, loadBalancer);

  prepareGeometry(geometry, converter);

  /// === 3rd Step: Prepare Lattice ===

  SuperLattice<T, ADEDESCRIPTOR> ADElattice(dummy_converter, geometry);
  SuperLattice<T, NSEDESCRIPTOR> NSElattice(converter, geometry);

  prepareLattice(converter, NSElattice, ADElattice, geometry );

  T boussinesqForcePrefactor = 9.81 / converter.getConversionFactorVelocity() * converter.getConversionFactorTime() *
                               converter.getCharPhysTemperatureDifference() * converter.getPhysThermalExpansionCoefficient();
  SuperLatticeCoupling coupling(
    NavierStokesAdvectionDiffusionCoupling{},
    names::NavierStokes{}, NSElattice,
    names::Temperature{},  ADElattice);
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::T0>(
    converter.getLatticeTemperature(Tcold));
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::FORCE_PREFACTOR>(
    boussinesqForcePrefactor * Vector<T,2>{0.0,1.0});

  /// === 4th Step: Main Loop with Timer ===
  util::Timer<T> timer(converter.getLatticeTime(maxPhysT), geometry.getStatistics().getNvoxel() );
  timer.start();

  util::ValueTracer<T> converge(converter.getLatticeTime(1.0),epsilon);
  for (iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {
    if (converge.hasConverged()) {
      clout << "Simulation converged." << std::endl;
      getResults(converter, NSElattice, ADElattice, iT, geometry, timer, converge.hasConverged());
      break;
    }

    /// === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(converter, NSElattice, ADElattice, iT, geometry);

    /// === 6th Step: Collide and Stream Execution ===
    NSElattice.collideAndStream();
    coupling.apply();
    ADElattice.collideAndStream();

    /// === 7th Step: Computation and Output of the Results ===
    getResults(converter, NSElattice, ADElattice, iT, geometry, timer, converge.hasConverged());
    converge.takeValue(NSElattice.getStatistics().getAverageEnergy());
  }

  timer.stop();
  timer.printSummary();
}
*/

int main(int argc, char* argv[]) {
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<DOMAIN_EXTENT        >({1, 1});
    myCaseParameters.set<PHYS_CHAR_LENGTH     >( 1.0 );
    myCaseParameters.set<LATTICE_CHAR_VELOCITY>( 0.1 );
    myCaseParameters.set<PHYS_CHAR_VISCOSITY  >( 1e-3 );
    myCaseParameters.set<RESOLUTION           >(20);
    myCaseParameters.set<RELAXATION_TIME      >(1.);
    myCaseParameters.set<REYNOLDS             >(20);
    myCaseParameters.set<RAYLEIGH             >(100);
    myCaseParameters.set<PRANDTL              >(0.71);
    myCaseParameters.set<MAX_PHYS_T           >(1e4);
    myCaseParameters.set<T_HOT                >(274.15);
    myCaseParameters.set<T_COLD               >(273.15);
    myCaseParameters.set<BC_ID                >(1);
  }
  myCaseParameters.set<PHYS_DELTA_X>([&]() -> MyCase::value_t {
    return {myCaseParameters.get<PHYS_CHAR_LENGTH>() / myCaseParameters.get<RESOLUTION>()};
  });

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