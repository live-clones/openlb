/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2020 Louis Kronberg, Stephan Simonis
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

#pragma once

#include <olb.h>

using namespace olb;
//using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::names;

using MyCase = Case<
  AdvectionDiffusion, Lattice<double, descriptors::D2Q5<descriptors::VELOCITY>>
>;

namespace olb::parameters {

  struct RUNS : public descriptors::FIELD_BASE<1> {};
  struct START_RESOLUTION : public descriptors::FIELD_BASE<1> {};

}

template <typename T, typename TDESCRIPTOR>
class AdePhysTemp1D : public AnalyticalF2D<T,T> {

protected:
  MyCase myCase;
  T t;
public:
  AdePhysTemp1D(T time, MyCase& _myCase) : AnalyticalF2D<T,T>(1), myCase(_myCase),
    t(time) {};

  bool operator()(T output[], const T input[]) override
  {
    using T = MyCase::value_t;
    auto& lattice = myCase.getLattice(AdvectionDiffusion{});
    const auto& converter = lattice.getUnitConverter();
    const T mue = converter.getPhysDiffusivity();
    const T uMag = converter.getCharPhysVelocity();
    const T res = converter.getResolution();

    T x = input[0];
    T gf = res/(res+1.);

    // initial condition (pseudo 1D)
    output[0] = util::sin(M_PI*(x-uMag*t)*gf) * util::exp(-mue*util::pow(M_PI,2)*t*util::pow(gf,2));

    return true;
  };
};

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  const T physCharLength = params.get<parameters::PHYS_CHAR_LENGTH>();
  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  Vector origin(-physCharLength/2, -physCharLength/2);

  IndicatorCuboid2D<T> cuboid(extent, origin);

#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif

  const T physDeltaX = physCharLength / params.get<parameters::RESOLUTION>();
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, noOfCuboids);
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({true,true});
  return mesh;
}

void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();

  geometry.rename(0,1);
  geometry.communicate();
  geometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  geometry.innerClean();
  geometry.checkForErrors();
  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase)
{
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<AdvectionDiffusion>;
  using namespace olb::parameters;

  auto& ADlattice = myCase.getLattice(AdvectionDiffusion{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const T physDeltaX = parameters.get<PHYS_CHAR_LENGTH>() / parameters.get<RESOLUTION>();
  const T physDeltaT = physDeltaX * physDeltaX;
  const T physLength = parameters.get<PHYS_CHAR_LENGTH>();
  const T mue = parameters.get<PHYS_DIFFUSIVITY>();
  const T physVelocity = parameters.get<PECLET>() * mue / physLength;
  const T physDensity = parameters.get<PHYS_CHAR_DENSITY>();

  ADlattice.setUnitConverter<AdeUnitConverter<T,DESCRIPTOR>>(
    physDeltaX,    // physDeltaX
    physDeltaT,    // physDeltaT (diffusive scaling)
    physLength,    // charPhysLength
    physVelocity,  // charPhysVelocity from Peclet
    mue,           // physDiffusivity
    physDensity    // physDensity,
  );

  const auto& converter = ADlattice.getUnitConverter();
  converter.print();

  ADlattice.defineDynamics<AdvectionDiffusionBGKdynamics>(geometry.getMaterialIndicator({1}));
}

void setInitialValues(MyCase& myCase) {
  using T = myCase::value_t;
  auto& ADlattice = myCase.getLattice();
  auto& geometry = myCase.getGeometry();
  const auto& converter = ADlattice.getUnitConverter();

  //TODO make 1D?
  AnalyticalConst1D<T,T> u0( converter.getCharLatticeVelocity());
  AdePhysTemp1D<T> Tinit( 0.0, converter );

  auto bulkIndicator = geometry.getMaterialIndicator({0,1});

  ADlattice.defineField<descriptors::VELOCITY>( bulkIndicator, u0 );
  ADlattice.defineRho( bulkIndicator, Tinit );
  ADlattice.iniEquilibrium( bulkIndicator, Tinit, u0 );

  ADlattice.setParameter<descriptors::OMEGA>( converter.getLatticeAdeRelaxationFrequency() );

  /// Make the lattice ready for simulation
  ADlattice.initialize();
}

T errorOverLine( SuperLattice<T, TDESCRIPTOR>& ADlattice,
                 SuperGeometry<T,2>& superGeometry,
                 int iT,
                 AdeUnitConverter<T, TDESCRIPTOR> converter )
{
  OstreamManager clout(std::cout,"error");

  T lx = converter.getCharPhysLength();
  T dist = converter.getPhysDeltaX();

  SuperLatticeDensity2D<T, TDESCRIPTOR> temperature(ADlattice);
  AnalyticalFfromSuperF2D<T> aTemp( temperature, true, 1);
  AdePhysTemp1D<T> temperatureSol( converter.getPhysTime(iT), converter);

  T result[1] = {0.};    // will contain result of analytical solution
  T resultSim[1] = {0.}; // will contain result of lattice evaluation
  T input[2] = {T(), 0}; // lattice node at which to evaluate the temparature

  int ndatapoints = converter.getResolution(); // number of data points on line
  T tempMSE = .0;                              // mean squared error
  T tempSquaredSum = .0;                       // sum of analytical temperature

  CSV<T> csvWriterSim("simulation" + std::to_string(converter.getPhysTime(iT)));
  CSV<T> csvWriterAna("analytical" + std::to_string(converter.getPhysTime(iT)));
  for (int i = 0; i <= ndatapoints; i++) {
    input[0] = -T(.5)*lx + i*dist;
    aTemp(resultSim, input);
    temperatureSol(result, input);

    tempMSE += (result[0]  - resultSim[0]) * (result[0]  - resultSim[0]);
    // remove division by zero
    tempSquaredSum += (result[0]+T(1)) * (result[0]+T(1));
    csvWriterSim.writeDataFile( input[0], resultSim[0], 16);
    csvWriterAna.writeDataFile(input[0], result[0], 16);
  }

  CSV<T> csvWriterErr;
  csvWriterErr.writeDataFile(converter.getPhysTime(iT), util::sqrt(tempMSE/tempSquaredSum), "averageL2RelError", 16);
  return util::sqrt(tempMSE/tempSquaredSum);
}


T getResults( SuperLattice<T, TDESCRIPTOR>& ADlattice,
              int iT,
              int statIter,
              AdeUnitConverter<T, TDESCRIPTOR> converter,
              SuperGeometry<T,2>& superGeometry,
              util::Timer<T>& timer)
{
  OstreamManager clout(std::cout,"getResults");
  SuperVTMwriter2D<T> vtkWriter("advectionDiffusion1d");

  AdePhysTemp1D<T> temperatureSol( converter.getPhysTime(iT), converter );

  T avg = 0;

  if (iT == 0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, TDESCRIPTOR> cuboid(ADlattice);
    SuperLatticeRank2D<T, TDESCRIPTOR> rank(ADlattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();

    SuperLatticeDensity2D<T, TDESCRIPTOR> temperature(ADlattice);
    SuperLatticeFfromAnalyticalF2D<T, TDESCRIPTOR> solution(temperatureSol, ADlattice);

    vtkWriter.addFunctor( temperature );
    vtkWriter.addFunctor( solution );
    vtkWriter.write( iT );

    timer.update(iT);
    timer.printStep();
    ADlattice.getStatistics().print( iT, converter.getPhysTime(iT) );

    /// Output data for plot of inital temperatre value of lattice
    CSV<T> csvWriterInit("initalValue");
    T result[1] = {0.};
    T input[2] = {T(), 0};
    for (int i = 0; i <= converter.getResolution(); i++) {
      input[0] = -T(.5)*converter.getCharPhysLength() + i*converter.getPhysDeltaX();
      temperatureSol(result, input);
      csvWriterInit.writeDataFile(input[0],result[0], 16);
    }
  }
  else if (iT % statIter == 0 ) {
    /// Writes the VTK files
    SuperLatticeDensity2D<T, TDESCRIPTOR> temperature(ADlattice);
    SuperLatticeFfromAnalyticalF2D<T, TDESCRIPTOR> solution(temperatureSol, ADlattice);
    ADlattice.setProcessingContext(ProcessingContext::Evaluation);

    vtkWriter.addFunctor( temperature );
    vtkWriter.addFunctor( solution );
    vtkWriter.write( iT );

    /// ADLattice statistics console output
    timer.update(iT);
    timer.printStep();
    ADlattice.getStatistics().print(iT, converter.getPhysTime( iT ));

    // compute relative error
    avg = errorOverLine(ADlattice, superGeometry, iT, converter);
    clout << "Relative L2-error norm over centerline: "  << avg << std::endl;
  }

  return avg;
}

void simulate(int N, int statIter, T mue, T peclet, T physLength)
{

  OstreamManager clout(std::cout,"simulate");
  clout << "Executing the simulation with N=" << std::to_string(N) << std::endl;

  AdeUnitConverter<T,TDESCRIPTOR> converter(
    physLength/N,            // physDeltaX
    util::pow(physLength/N, 2),    // physDeltaT (diffusive scaling)
    physLength,              // charPhysLength
    peclet*mue/physLength,   // charPhysVelocity from Peclet
    mue,                     // physDiffusivity
    1                        // physDensity,
  );

  converter.print();

  // switch outdirectory only if there are multiple simulation runs
  if (runs > 1) {
    singleton::directories().setOutputDir("./tmp/N" + std::to_string(N) + "/");
    /// file clean-up
    CSV<T> csvWriterErr("averageL2RelError");
    csvWriterErr.clearFile();
  }

  /// === 2nd Step: Prepare Geometry ===
  std::vector<T> extend(2,T());
  std::vector<T> origin(2,T());

  extend[0] = converter.getCharPhysLength();
  extend[1] = converter.getCharPhysLength();

  origin[0] = - converter.getCharPhysLength()/2;
  origin[1] = - converter.getCharPhysLength()/2;

  IndicatorCuboid2D<T> cuboid(extend, origin);

  /// Instantiation of an empty cuboidDecomposition
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidDecomposition2D<T> cuboidDecomposition(cuboid, converter.getPhysDeltaX(), noOfCuboids);
  cuboidDecomposition.setPeriodicity({true, true});

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  /// Instantiation of a superGeometry
  SuperGeometry<T,2> superGeometry(cuboidDecomposition, loadBalancer, 2);
  prepareGeometry(superGeometry, cuboid);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice<T, TDESCRIPTOR> ADlattice(converter, superGeometry);

  prepareLattice(ADlattice, superGeometry, converter);

  /// === 4th Step: Main Loop with Timer ===
  util::Timer<T> timer( superGeometry.getStatistics().getNvoxel() );
  timer.start();

  T pulseDiffBound = 1e-3;
  int timeCount = util::ceil( -1./converter.getPhysDeltaT() * util::log(pulseDiffBound)/(converter.getPhysDiffusivity()*util::pow(M_PI,2)) /statIter );
  int iTmax = timeCount * statIter;

  int iT;
  T simulationAverage = .0;
  for (iT = 0; iT < iTmax; ++iT) {
    simulationAverage += getResults(ADlattice, iT, statIter, converter, superGeometry, timer);

    /// === 6th Step: Collide and Stream Execution ===
    ADlattice.collideAndStream();
  }

  simulationAverage /= timeCount;

  // this outputs into ./tmp/gnuplotData/data/averageSimL2RelErr
  singleton::directories().setOutputDir("./tmp/");
  CSV<T> csvWriterErr;
  csvWriterErr.writeDataFile(N, simulationAverage, "averageSimL2RelErr", 16);

  clout << "Simulation Average Relative L2 Error: " << simulationAverage << std::endl;

  ADlattice.setProcessingContext(ProcessingContext::Evaluation);
  timer.stop();
  timer.printSummary();
}
