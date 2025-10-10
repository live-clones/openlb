/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Florian Kaiser
 *                2008 Orestis Malaspinas
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
  NavierStokes, Lattice<float, descriptors::D2Q9<descriptors::FORCE>>,
  Temperature,  Lattice<float, descriptors::D2Q5<descriptors::VELOCITY>>
>;

// #define AdvectionDiffusionDirichletBC
// #define RegularizedTemperatureBC
#define RegularizedHeatFluxBC

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

void computeError(MyCase& myCase)
{
  OstreamManager clout( std::cout, "computeError" );
  using T = MyCase::value_t_of<NavierStokes>;
  using NSEDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ADEDESCRIPTOR = MyCase::descriptor_t_of<Temperature>;
  auto& geometry = myCase.getGeometry();
  auto& NSElattice = myCase.getLattice(NavierStokes{});
  auto& ADElattice = myCase.getLattice(Temperature{});
  const auto& converter = NSElattice.getUnitConverter();
  const T Re = converter.getReynoldsNumber();
  const T Pr = converter.getPrandtlNumber();

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
  const Vector bottomExtent{extent[0], 0.};
  const Vector origin{0., 0.};
  IndicatorCuboid2D<T> bottom(bottomExtent, origin);

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
  auto& parameters    = myCase.getParameters();
  auto& NSElattice    = myCase.getLattice(NavierStokes{});
  auto& ADElattice    = myCase.getLattice(Temperature{});
  using NSEDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ADEDESCRIPTOR = MyCase::descriptor_t_of<Temperature>;

  const T physCharLength          = parameters.get<parameters::PHYS_CHAR_LENGTH>();
  const int N                     = parameters.get<parameters::RESOLUTION>();
  const T physDeltaX              = parameters.get<parameters::PHYS_DELTA_X>();
  const T physCharVelocity        = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T physDensity             = parameters.get<parameters::PHYS_CHAR_DENSITY>();
  const T gravitationalConstant   = parameters.get<parameters::GRAVITATIONAL_ACC>();
  const T physThermalConductivity = parameters.get<parameters::PHYS_THERMAL_CONDUCTIVITY>();
  const T tau                     = parameters.get<parameters::LATTICE_RELAXATION_TIME>();
  const T Ra                      = parameters.get<parameters::RAYLEIGH>();
  const T Re                      = parameters.get<parameters::REYNOLDS>();
  const T Pr                      = parameters.get<parameters::PRANDTL>();
  const T Tcold                   = parameters.get<parameters::T_COLD>();
  const T Thot                    = parameters.get<parameters::T_HOT>();

  const T physViscosity = physCharVelocity * physCharLength / Re;
  const T physDeltaT = physDeltaX * physDeltaX * (1. / descriptors::invCs2<T, NSEDESCRIPTOR>()) * (tau - 0.5) / physViscosity; // physDeltaT
  const T physThermalExpansionCoefficient = Ra * physViscosity * physViscosity / (Pr * gravitationalConstant * (Thot - Tcold) * physCharLength * physCharLength * physCharLength);
  const T physSpecificHeatCapacity = Pr * physThermalConductivity / (physViscosity * physDensity);

  NSElattice.setUnitConverter<ThermalUnitConverter<T,NSEDESCRIPTOR,ADEDESCRIPTOR>>(
    (T) physDeltaX,
    (T) physDeltaT,
    (T) physCharLength,
    (T) physCharVelocity,
    (T) physViscosity,
    (T) physDensity,
    (T) physThermalConductivity,
    (T) physSpecificHeatCapacity,
    (T) physThermalExpansionCoefficient,
    (T) Tcold,
    (T) Thot
  );
  const auto& converter = NSElattice.getUnitConverter();
  converter.print();

  ADElattice.setUnitConverter(converter);

  T ADEomega  = converter.getLatticeThermalRelaxationFrequency();
  T NSEomega  = converter.getLatticeRelaxationFrequency();

  ADElattice.defineDynamics<AdvectionDiffusionBGKdynamics>(geometry, 1);
  ADElattice.defineDynamics<AdvectionDiffusionBGKdynamics>(geometry, 2);
  ADElattice.defineDynamics<AdvectionDiffusionBGKdynamics>(geometry, 3);
  NSElattice.defineDynamics<ForcedBGKdynamics>(geometry, 1);
  NSElattice.defineDynamics<ForcedBGKdynamics>(geometry, 2);
  NSElattice.defineDynamics<ForcedBGKdynamics>(geometry, 3);

  /// sets boundary
  boundary::set<boundary::LocalVelocity>(NSElattice, geometry, 2);
  boundary::set<boundary::LocalVelocity>(NSElattice, geometry, 3);

  #ifdef AdvectionDiffusionDirichletBC
    boundary::set<boundary::AdvectionDiffusionDirichlet>(ADElattice, geometry, 2);
    boundary::set<boundary::AdvectionDiffusionDirichlet>(ADElattice, geometry, 3);
  #endif

  #ifdef RegularizedTemperatureBC
    boundary::set<boundary::RegularizedTemperature>(ADElattice, geometry, 2);
    boundary::set<boundary::RegularizedTemperature>(ADElattice, geometry, 3);
  #endif
  #ifdef RegularizedHeatFluxBC
    T heatFlux[2];
    T input[2] = {0.,1.};
    AnalyticalHeatFluxPorousPlate2D<T,T> HeatFluxSol(Re, Pr, converter.getCharPhysTemperatureDifference(), converter.getCharPhysLength(), converter.getThermalConductivity());
    HeatFluxSol(heatFlux, input);
    T temp = converter.getLatticeSpecificHeatCapacity(converter.getPhysSpecificHeatCapacity())*(converter.getLatticeThermalRelaxationTime() - 0.5) / converter.getLatticeThermalRelaxationTime();
    heatFlux[0] = converter.getLatticeHeatFlux(heatFlux[0]) / temp;
    heatFlux[1] = converter.getLatticeHeatFlux(heatFlux[1]) / temp;
    AnalyticalConst2D<T,T> heatFluxC(heatFlux[0], heatFlux[1]);
    boundary::set<boundary::RegularizedHeatFlux>(ADElattice, geometry, 2);
    ADElattice.defineU(geometry, 2, heatFluxC);
    boundary::set<boundary::RegularizedTemperature>(ADElattice, geometry, 3);
    ADElattice.setParameter<descriptors::OMEGA>(ADEomega);
  #endif


  ADElattice.setParameter<descriptors::OMEGA>(ADEomega);
  NSElattice.setParameter<descriptors::OMEGA>(NSEomega);

  const T boussinesqForcePrefactor = 9.81 / converter.getConversionFactorVelocity()
                                   * converter.getConversionFactorTime()
                                   * converter.getCharPhysTemperatureDifference()
                                   * converter.getPhysThermalExpansionCoefficient();

  auto& coupling = myCase.setCouplingOperator(
    "Boussinesq",
    NavierStokesAdvectionDiffusionCoupling{},
    names::NavierStokes{}, NSElattice,
    names::Temperature{},  ADElattice);
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::T0>(
    converter.getLatticeTemperature(Tcold));
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::FORCE_PREFACTOR>(
    boussinesqForcePrefactor * Vector<T,2>{0.0,1.0});

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t_of<NavierStokes>;
  auto& geometry = myCase.getGeometry();
  auto& NSElattice = myCase.getLattice(NavierStokes{});
  auto& converter = NSElattice.getUnitConverter();
  auto& ADElattice = myCase.getLattice(Temperature{});

  const T Tcold = converter.getCharPhysLowTemperature();
  const T Thot = converter.getCharPhysHighTemperature();
  /// for each material set the defineRhoU and the Equilibrium
  std::vector<T> zero(2,T());
  AnalyticalConst2D<T,T> u(zero);
  AnalyticalConst2D<T,T> rho(1.);
  AnalyticalConst2D<T,T> force(zero);

  T u_Re = converter.getLatticeVelocity( converter.getReynoldsNumber() * converter.getPhysViscosity() / converter.getCharPhysLength() );

  AnalyticalConst2D<T,T> u_top(converter.getCharLatticeVelocity(), u_Re);
  AnalyticalConst2D<T,T> u_bot(0.0, u_Re);

  NSElattice.defineRhoU(geometry, 1, rho, u);
  NSElattice.iniEquilibrium(geometry, 1, rho, u);
  NSElattice.defineField<descriptors::FORCE>(geometry, 1, force);
  NSElattice.defineRhoU(geometry, 2, rho, u_top);
  NSElattice.iniEquilibrium(geometry, 2, rho, u_top);
  NSElattice.defineField<descriptors::FORCE>(geometry, 2, force);
  NSElattice.defineRhoU(geometry, 3, rho, u_bot);
  NSElattice.iniEquilibrium(geometry, 3, rho, u_bot);
  NSElattice.defineField<descriptors::FORCE>(geometry, 3, force);

  AnalyticalConst2D<T,T> Cold(converter.getLatticeTemperature(Tcold));
  AnalyticalConst2D<T,T> Hot(converter.getLatticeTemperature(Thot));

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

void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT,
                bool converged)
{
  OstreamManager clout(std::cout,"getResults");

  using T               = MyCase::value_t_of<NavierStokes>;
  using NSEDESCRIPTOR   = MyCase::descriptor_t_of<NavierStokes>;
  using ADEDESCRIPTOR   = MyCase::descriptor_t_of<Temperature>;
  auto& NSElattice      = myCase.getLattice(NavierStokes{});
  auto& ADElattice      = myCase.getLattice(Temperature{});
  auto& parameters      = myCase.getParameters();
  const auto& converter = NSElattice.getUnitConverter();
  const T Re            = converter.getReynoldsNumber();
  const T Pr            = converter.getPrandtlNumber();
  const std::size_t N   = converter.getResolution();
  const T Tcold = converter.getCharPhysLowTemperature();
  const T Thot = converter.getCharPhysHighTemperature();
  const int statIter = converter.getLatticeTime(parameters.get<parameters::PHYS_STAT_ITER_T>());
  const int vtkIter = converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T>());

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
    computeError(myCase);
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

  if (iT % statIter == 0 || converged) {
    timer.print(iT);
  }

}

void simulate(MyCase& myCase) {
  OstreamManager clout(std::cout,"simulate");

  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& NSElattice = myCase.getLattice(NavierStokes{});
  auto& converter  = NSElattice.getUnitConverter();

  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(
    parameters.get<parameters::MAX_PHYS_T>()
  );

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  const T epsilon       = parameters.get<parameters::EPSILON>();
  const T convCheckTime = parameters.get<parameters::CONV_ITER>();
  util::ValueTracer<T> converge(converter.getLatticeTime(convCheckTime), epsilon);
  for (std::size_t iT=0; iT < iTmax; ++iT) {
    if (converge.hasConverged()) {
      clout << "Simulation converged." << std::endl;
      getResults(myCase, timer, iT, converge.hasConverged());
      break;
    }

    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    myCase.getLattice(NavierStokes{}).collideAndStream();
    myCase.getOperator("Boussinesq").apply();
    myCase.getLattice(Temperature{}).collideAndStream();


    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, timer, iT, converge.hasConverged());
    converge.takeValue(NSElattice.getStatistics().getAverageEnergy());
  }

  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[]) {
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<DOMAIN_EXTENT            >({1, 1});
    myCaseParameters.set<PHYS_CHAR_LENGTH         >( 1.0 );
    myCaseParameters.set<LATTICE_CHAR_VELOCITY    >( 0.1 );
    myCaseParameters.set<PHYS_CHAR_VISCOSITY      >( 1e-3 );
    myCaseParameters.set<PHYS_CHAR_DENSITY        >( 1. );
    myCaseParameters.set<PHYS_THERMAL_CONDUCTIVITY>( 0.03 );
    myCaseParameters.set<RESOLUTION               >(20);
    myCaseParameters.set<LATTICE_RELAXATION_TIME  >(1.);
    myCaseParameters.set<REYNOLDS                 >(20);
    myCaseParameters.set<RAYLEIGH                 >(100);
    myCaseParameters.set<PRANDTL                  >(0.71);
    myCaseParameters.set<MAX_PHYS_T               >(1e4);
    myCaseParameters.set<T_HOT                    >(274.15);
    myCaseParameters.set<T_COLD                   >(273.15);
    myCaseParameters.set<EPSILON                  >(1e-7);
    myCaseParameters.set<CONV_ITER                >(1.);
    myCaseParameters.set<PHYS_STAT_ITER_T         >(10.0);
    myCaseParameters.set<PHYS_VTK_ITER_T          >(100.);
  }
  myCaseParameters.set<parameters::PHYS_DELTA_X>([&]() -> MyCase::value_t {
    return {myCaseParameters.get<parameters::PHYS_CHAR_LENGTH>() / myCaseParameters.get<parameters::RESOLUTION>()};
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
