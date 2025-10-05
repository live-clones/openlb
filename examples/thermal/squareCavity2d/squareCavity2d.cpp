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

/* squareCavity2d.cpp
 * The reference is the paper in "Gaedtke, M., Wachter, S., Raedle, M., Nirschl, H., & Krause, M. J. (2018).
 * Application of a lattice Boltzmann method combined with a Smagorinsky turbulence model to spatially resolved heat flux inside a refrigerated vehicle.
 * Computers & Mathematics with Applications, 76(10), 2315-2329."
 */

// natural convection of air in a square cavity in 2D


#include <olb.h>
#include "squareCavity2d.h"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = double;

#define TURBULENT

#ifdef TURBULENT
typedef D2Q9<FORCE,TAU_EFF> NSDESCRIPTOR;
typedef D2Q5<VELOCITY,TAU_EFF> TDESCRIPTOR;
#else
typedef D2Q9<FORCE> NSDESCRIPTOR;
typedef D2Q5<VELOCITY> TDESCRIPTOR;
#endif


#ifdef TURBULENT
const int statisticsIntervall = 10; // take the turbulent statistics every 10 time steps after convergence
const int statisticsEnsembles = 20; // take 20 ensembles for the turbulent statistics
#endif

/// Compute the nusselt number at the left wall
T computeNusselt(SuperGeometry<T,2>& superGeometry,
                 SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                 SuperLattice<T, TDESCRIPTOR>& ADlattice)
{
  int voxel = 0, material = 0;
  T T_x = 0, T_xplus1 = 0, T_xplus2 = 0;
  T q = 0;

  for (int iC = 0; iC < NSlattice.getLoadBalancer().size(); iC++) {
    int ny = NSlattice.getBlock(iC).getNy();
    int iX = 0;
    for (int iY = 0; iY < ny; ++iY) {
      material = superGeometry.getBlockGeometry(iC).getMaterial(iX,iY);

      T_x = ADlattice.getBlock(iC).get(iX,iY).computeRho();
      T_xplus1 = ADlattice.getBlock(iC).get(iX+1,iY).computeRho();
      T_xplus2 = ADlattice.getBlock(iC).get(iX+2,iY).computeRho();

      if ( material == 2 ) {
        q += (3.0*T_x - 4.0*T_xplus1 + 1.0*T_xplus2)/2.0*N;
        voxel++;
      }
    }
  }

#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(q, MPI_SUM);
  singleton::mpi().reduceAndBcast(voxel, MPI_SUM);
#endif

  return q / (T)voxel;
}

/// Stores geometry information in form of material numbers
void prepareGeometry(SuperGeometry<T,2>& superGeometry,
                     ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter)
{

  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0,4);

  std::vector<T> extend(2,T());
  extend[0] = lx;
  extend[1] = lx;
  std::vector<T> origin(2,T());
  origin[0] = converter.getPhysLength(1);
  origin[1] = 0.5*converter.getPhysLength(1);
  IndicatorCuboid2D<T> cuboid2(extend, origin);

  superGeometry.rename(4,1,cuboid2);

  std::vector<T> extendwallleft(2,T(0));
  extendwallleft[0] = converter.getPhysLength(1);
  extendwallleft[1] = lx;
  std::vector<T> originwallleft(2,T(0));
  originwallleft[0] = 0.0;
  originwallleft[1] = 0.0;
  IndicatorCuboid2D<T> wallleft(extendwallleft, originwallleft);

  std::vector<T> extendwallright(2,T(0));
  extendwallright[0] = converter.getPhysLength(1);
  extendwallright[1] = lx;
  std::vector<T> originwallright(2,T(0));
  originwallright[0] = lx+1.5*converter.getPhysLength(1);
  originwallright[1] = 0.0;
  IndicatorCuboid2D<T> wallright(extendwallright, originwallright);

  superGeometry.rename(4,2,1,wallleft);
  superGeometry.rename(4,3,1,wallright);


  /// Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;

}

void prepareLattice( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                     SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                     SuperLattice<T, TDESCRIPTOR>& ADlattice,
                     SuperGeometry<T,2>& superGeometry )
{

  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  T omega  =  converter.getLatticeRelaxationFrequency();
  T Tomega  =  converter.getLatticeThermalRelaxationFrequency();

  #ifdef TURBULENT
    NSlattice.defineDynamics<ExternalTauEffLESForcedBGKdynamics<T,NSDESCRIPTOR,momenta::AdvectionDiffusionBulkTuple>>(superGeometry.getMaterialIndicator({1, 2, 3}));
    ADlattice.defineDynamics<ExternalTauEffLESBGKadvectionDiffusionDynamics>(superGeometry.getMaterialIndicator({1, 2, 3}));
  #else
    NSlattice.defineDynamics<ForcedBGKdynamics>(superGeometry.getMaterialIndicator({1, 2, 3}));
    ADlattice.defineDynamics<AdvectionDiffusionBGKdynamics>(superGeometry.getMaterialIndicator({1, 2, 3}));
  #endif

  boundary::set<boundary::BounceBack>(ADlattice, superGeometry, 4);
  boundary::set<boundary::BounceBack>(NSlattice, superGeometry, 4);

  /// sets boundary
  boundary::set<boundary::AdvectionDiffusionDirichlet>(ADlattice, superGeometry.getMaterialIndicator({2, 3}));
  boundary::set<boundary::LocalVelocity>(NSlattice, superGeometry.getMaterialIndicator({2, 3}));

  /// define initial conditions
  AnalyticalConst2D<T,T> rho(1.);
  AnalyticalConst2D<T,T> u0(0.0, 0.0);
  AnalyticalConst2D<T,T> T_cold(converter.getLatticeTemperature(Tcold));
  AnalyticalConst2D<T,T> T_hot(converter.getLatticeTemperature(Thot));
  AnalyticalConst2D<T,T> T_mean(converter.getLatticeTemperature(Tmean));

  /// for each material set Rho, U and the Equilibrium
  NSlattice.defineRhoU(superGeometry.getMaterialIndicator({1, 2, 3}), rho, u0);
  NSlattice.iniEquilibrium(superGeometry.getMaterialIndicator({1, 2, 3}), rho, u0);

  ADlattice.defineRho(superGeometry, 1, T_mean);
  ADlattice.iniEquilibrium(superGeometry, 1, T_mean, u0);
  ADlattice.defineRho(superGeometry, 2, T_hot);
  ADlattice.iniEquilibrium(superGeometry, 2, T_hot, u0);
  ADlattice.defineRho(superGeometry, 3, T_cold);
  ADlattice.iniEquilibrium(superGeometry, 3, T_cold, u0);

#ifdef TURBULENT
  AnalyticalConst2D<T,T> tauNS(1./omega);
  AnalyticalConst2D<T,T> tauAD(1./Tomega);

  NSlattice.defineField<descriptors::TAU_EFF>( superGeometry.getMaterialIndicator({1, 2, 3}), tauNS );
  ADlattice.defineField<descriptors::TAU_EFF>( superGeometry.getMaterialIndicator({1, 2, 3}), tauAD );
#endif

  NSlattice.setParameter<descriptors::OMEGA>(omega);
  ADlattice.setParameter<descriptors::OMEGA>(Tomega);

  /// Make the lattice ready for simulation
  NSlattice.initialize();
  ADlattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                        SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                        SuperLattice<T, TDESCRIPTOR>& ADlattice,
                        int iT, SuperGeometry<T,2>& superGeometry)
{

  // nothing to do here

}

void getResults( Vector<T, 18>& output,
                 ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                 SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                 SuperLattice<T, TDESCRIPTOR>& ADlattice, int iT,
                 SuperGeometry<T,2>& superGeometry,
                 util::Timer<T>& timer,
                 bool converged)
{

  OstreamManager clout(std::cout,"getResults");

  SuperVTMwriter2D<T> vtkWriter("squareCavity2d");
  SuperLatticePhysVelocity2D<T, NSDESCRIPTOR> velocity(NSlattice, converter);
  SuperLatticePhysPressure2D<T, NSDESCRIPTOR> pressure(NSlattice, converter);
  SuperLatticePhysTemperature2D<T, NSDESCRIPTOR, TDESCRIPTOR> temperature(ADlattice, converter);
  vtkWriter.addFunctor( pressure );
  vtkWriter.addFunctor( velocity );
  vtkWriter.addFunctor( temperature );

  AnalyticalFfromSuperF2D<T> interpolation(velocity, true);

  const int statIter = 2000.;

  if (iT == 0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, NSDESCRIPTOR> cuboid(NSlattice);
    SuperLatticeRank2D<T, NSDESCRIPTOR> rank(NSlattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
  }

  /// Writes the VTK files
  if (iT % statIter == 0 || converged) {
    NSlattice.setProcessingContext(ProcessingContext::Evaluation);
    ADlattice.setProcessingContext(ProcessingContext::Evaluation);

    timer.update(iT);
    timer.printStep();

    /// NSLattice statistics console output
    NSlattice.getStatistics().print(iT,converter.getPhysTime(iT));
    /// ADLattice statistics console output
    ADlattice.getStatistics().print(iT,converter.getPhysTime(iT));

    vtkWriter.write(iT);

    BlockReduction2D2D<T> planeReduction(temperature, 600, BlockDataSyncMode::ReduceOnly);
    BlockGifWriter<T> gifWriter;
    gifWriter.write(planeReduction, Tcold-1, Thot+1, iT, "temperature");

    SuperEuklidNorm2D<T, NSDESCRIPTOR> normVel( velocity );
    BlockReduction2D2D<T> planeReduction2(normVel, 600, BlockDataSyncMode::ReduceOnly);
    BlockGifWriter<T> gifWriter2;
    gifWriter2.write( planeReduction2, iT, "velocity" );

  }

  if ( converged ) {
    T nusselt = computeNusselt(superGeometry, NSlattice, ADlattice);

    /// Initialize vectors for data output
    T xVelocity[2] = { T() };
    T outputVelX[2] = { T() };
    T yVelocity[2] = { T() };
    T outputVelY[2] = { T() };
    const int outputSize = 512;
    Vector<T, outputSize> velX;
    Vector<T, outputSize> posX;
    Vector<T, outputSize> velY;
    Vector<T, outputSize> posY;

    /// loop for the resolution of the cavity at x = lx/2 in yDirection and vice versa
    for (int n = 0; n < outputSize; ++n) {
      T yPosition[2] = { lx / 2, lx * n / (T) outputSize };
      T xPosition[2] = { lx * n / (T) outputSize, lx / 2 };

      /// Interpolate xVelocity at x = lx/2 for each yPosition
      interpolation(xVelocity, yPosition);
      interpolation(yVelocity, xPosition);
      /// Store the interpolated values to compare them among each other in order to detect the maximum
      velX[n] = xVelocity[0];
      posY[n] = yPosition[1];
      velY[n] = yVelocity[1];
      posX[n] = xPosition[0];

      /// Initialize output with the corresponding velocities and positions at the origin
      if (n == 0) {
        outputVelX[0] = velX[0];
        outputVelX[1] = posY[0];
        outputVelY[0] = velY[0];
        outputVelY[1] = posX[0];
      }
      /// look for the maximum velocity in xDirection and the corresponding position in yDirection
      if (n > 0 && velX[n] > outputVelX[0]) {
        outputVelX[0] = velX[n];
        outputVelX[1] = posY[n];
      }
      /// look for the maximum velocity in yDirection and the corresponding position in xDirection
      if (n > 0 && velY[n] > outputVelY[0]) {
        outputVelY[0] = velY[n];
        outputVelY[1] = posX[n];
      }
    }

    Vector<T, 5> parameterCollection = {outputVelX[0], outputVelY[0], outputVelX[1], outputVelY[1], nusselt};
    compareToStudy(output, parameterCollection);
  }
}

void simulate(Vector<T, 18>& output)
{
  OstreamManager clout(std::cout,"simulate");

  ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> converter(
    (T) lx / N,
    #ifdef TURBULENT
    // (T) 2.*0.056/charU*lx/N,
    (T) 2. * 0.01 / charU * lx / N,
    #else
    (T) (tau - 0.5) / descriptors::invCs2<T,NSDESCRIPTOR>() * util::pow((lx/N),2) / kinematicViscosity,
    #endif
    (T) lx,
    (T) charU,
    (T) kinematicViscosity,
    (T) density,
    (T) thermalDiffusivity,
    (T) Pr * thermalDiffusivity / kinematicViscosity / density,
    (T) thermalExpansion,
    (T) Tcold,
    (T) Thot
  );
  converter.print();


  UnitConverter<T,TDESCRIPTOR> dummy_converter(
    (T)   converter.getPhysDeltaX(),      // deltaX
    (T)   converter.getPhysDeltaT(),      // deltaT
    (T)   converter.getCharPhysLength(),  // from converter
    (T)   converter.getCharPhysVelocity(),// from converter
    (T)   converter.getPhysViscosity(),   // from converter
    (T)   converter.getPhysDensity()      // from converter
  );
  physLambda = converter.getPhysThermalDiffusivity();
  charL = converter.getCharPhysLength();

  /// === 2nd Step: Prepare Geometry ===
  std::vector<T> extend(2,T());
  extend[0] = lx + 2*converter.getPhysLength(1);
  extend[1] = lx + converter.getPhysLength(1);
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extend, origin);

  /// Instantiation of an empty cuboidDecomposition
  CuboidDecomposition2D<T> cuboidDecomposition(cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize());

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  /// Instantiation of a superGeometry
  SuperGeometry<T,2> superGeometry(cuboidDecomposition, loadBalancer);

  prepareGeometry(superGeometry, converter);

  /// === 3rd Step: Prepare Lattice ===

  SuperLattice<T, TDESCRIPTOR> ADlattice(dummy_converter, superGeometry);
  SuperLattice<T, NSDESCRIPTOR> NSlattice(converter, superGeometry);

  //prepareLattice and setBoundaryCondition
  prepareLattice(converter, NSlattice, ADlattice, superGeometry);

  T boussinesqForcePrefactor = 9.81 / converter.getConversionFactorVelocity() * converter.getConversionFactorTime() *
                               converter.getCharPhysTemperatureDifference() * converter.getPhysThermalExpansionCoefficient();

#ifdef TURBULENT
  const T preFactor = smagoConst*smagoConst
                    * descriptors::invCs2<T,NSDESCRIPTOR>()*descriptors::invCs2<T,NSDESCRIPTOR>()
                    * 2*util::sqrt(2);

  SuperLatticeCoupling coupling(
    SmagorinskyBoussinesqCoupling{},
    names::NavierStokes{}, NSlattice,
    names::Temperature{},  ADlattice);
  coupling.setParameter<SmagorinskyBoussinesqCoupling::T0>(
    converter.getLatticeTemperature(Tmean));
  coupling.setParameter<SmagorinskyBoussinesqCoupling::FORCE_PREFACTOR>(
    boussinesqForcePrefactor * Vector<T,2>{0.0,1.0});
  coupling.setParameter<SmagorinskyBoussinesqCoupling::SMAGORINSKY_PREFACTOR>(preFactor);
  coupling.setParameter<SmagorinskyBoussinesqCoupling::PR_TURB>(PrTurb);
  coupling.setParameter<SmagorinskyBoussinesqCoupling::OMEGA_NSE>(
    converter.getLatticeRelaxationFrequency());
  coupling.setParameter<SmagorinskyBoussinesqCoupling::OMEGA_ADE>(
    converter.getLatticeThermalRelaxationFrequency());
#else
  SuperLatticeCoupling coupling(
    NavierStokesAdvectionDiffusionCoupling{},
    names::NavierStokes{}, NSlattice,
    names::Temperature{},  ADlattice);
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::T0>(
    converter.getLatticeTemperature(Tcold));
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::FORCE_PREFACTOR>(
    boussinesqForcePrefactor * Vector<T,2>{0.0,1.0});
#endif

#ifdef TURBULENT
  SuperVTMwriter2D<T> vtkWriter("squareCavity2d");

  SuperLatticePhysTemperature2D<T,NSDESCRIPTOR, TDESCRIPTOR> sTemp(ADlattice,converter);
  SuperLatticePhysVelocity2D<T,NSDESCRIPTOR> sVel(NSlattice,converter);

  SuperLatticeTimeAveragedF2D<T> sAveragedTemp(sTemp);
  SuperLatticeTimeAveragedF2D<T> sAveragedVel(sVel);
  SuperLatticeTimeAveragedCrossCorrelationF2D<T> sAveragedTempVelCross(sTemp,sVel);
  SuperLatticeTimeAveragedCrossCorrelationF2D<T> sAveragedVelVelCross(sVel,sVel);
#endif

  /// === 4th Step: Main Loop with Timer ===
  util::Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  util::ValueTracer<T> converge(6,epsilon);
  bool converged = false;
  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {

    if (converge.hasConverged() && !converged) {
      converged = true;
      clout << "Simulation converged." << std::endl;
      clout << "Time " << iT << "." << std::endl;

      getResults(output, converter, NSlattice, ADlattice, iT, superGeometry, timer, converge.hasConverged());
    }

    /// === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(converter, NSlattice, ADlattice, iT, superGeometry);

    /// === 6th Step: Collide and Stream Execution ===
    coupling.apply();

    NSlattice.collideAndStream();
    ADlattice.collideAndStream();

    /// === 7th Step: Computation and Output of the Results ===
    if ( !converged ) {
      getResults(output, converter, NSlattice, ADlattice, iT, superGeometry, timer, converge.hasConverged());
    }
    if (!converged && iT % 1000 == 0) {
      ADlattice.setProcessingContext(ProcessingContext::Evaluation);
      converge.takeValue(computeNusselt(superGeometry, NSlattice, ADlattice),true);
    }
#ifdef TURBULENT
    if (converged && iT % statisticsIntervall == 0) {
      NSlattice.communicate();
      ADlattice.communicate();
      NSlattice.setProcessingContext(ProcessingContext::Evaluation);
      ADlattice.setProcessingContext(ProcessingContext::Evaluation);
      sAveragedTemp.addEnsemble();
      sAveragedVel.addEnsemble();
      sAveragedTempVelCross.addEnsemble();
      sAveragedVelVelCross.addEnsemble();
      if ( sAveragedTemp.getEnsembles() >= statisticsEnsembles ) {
        break;
      }
    }
    #endif
    if (converged) break;
  }
#ifdef TURBULENT
  vtkWriter.write(sAveragedTemp);
  vtkWriter.write(sAveragedVel);
  vtkWriter.write(sTemp);
  vtkWriter.write(sVel);
#endif

  timer.stop();
  timer.printSummary();
  return;
}

void setRayleighDependentParameters(T& Ra)
{
  lx = util::pow(Ra * kinematicViscosity * kinematicViscosity / Pr / gravitationalConstant / (Thot - Tcold) / thermalExpansion, (T) 1 / 3);  // length of the square
  charU = 1.0 / lx /( Pr * thermalDiffusivity / kinematicViscosity / 1.0 * 1.0 / thermalDiffusivity) * charUMultiplier();
}

void defineResolutions(std::vector<int>& convergenceRes, int startRes, int nRes) {
  for (int i = 0; i < nRes; i++) {
      convergenceRes[i] = startRes * (1 << i);
  }
}

int main(int argc, char *argv[])
{

  /// === 1st Step: Initialization ===
  OstreamManager clout(std::cout,"main");
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  CLIreader args(argc, argv);

  std::string fName( "squareCavity2d.xml" );
  XMLreader config( fName );

#ifndef TURBULENT
  tau                   = config["Application"]["Discretization"]["LatticeRelaxationTime"].get<T>();
#endif
  smagoConst            = config["Application"]["PhysParameters"]["SmagorinskyConstant"].get<T>();

  maxPhysT              = config["Application"]["PhysParameters"]["PhysMaxTime"].get<T>();
  Pr                    = config["Application"]["PhysParameters"]["Prandtl"].get<T>();
  PrTurb                = config["Application"]["PhysParameters"]["PrandtlTurb"].get<T>();
  thermalExpansion      = config["Application"]["PhysParameters"]["ThermalExpansion"].get<T>();
  thermalDiffusivity    = config["Application"]["PhysParameters"]["ThermalDiffusivity"].get<T>();
  density               = config["Application"]["PhysParameters"]["PhysDensity"].get<T>();
  heatCapacity          = config["Application"]["PhysParameters"]["HeatCapacity"].get<T>();
  kinematicViscosity    = config["Application"]["PhysParameters"]["KinematicViscosity"].get<T>();
  gravitationalConstant = config["Application"]["PhysParameters"]["GravitationalConstant"].get<T>();
  charUConfig           = config["Application"]["PhysParameters"]["CharLatticeVelocity"].get<T>();
  Thot                  = config["Application"]["PhysParameters"]["HotTemperature"].get<T>();
  Tcold                 = config["Application"]["PhysParameters"]["ColdTemperature"].get<T>();
  Tmean                 = (Thot + Tcold) / 2.0;

  epsilon               = config["Application"]["ConvergenceCheck"]["Residuum"].get<T>();

  bool convergenceTest  = config["Application"]["ConvergenceTest"]["Mode"].get<bool>();
  clout << "convergence test " << convergenceTest << std::endl;
  Vector<T, 18> output;
  if (convergenceTest)
  {
    Ra = config["Application"]["PhysParameters"]["Rayleigh"].get<T>();
    const int startRes = config["Application"]["ConvergenceTest"]["StartResolution"].get<int>();
    const int numRes = config["Application"]["ConvergenceTest"]["NumberResolutions"].get<int>();
    std::vector<int> convergenceRes(numRes);
    defineResolutions(convergenceRes, startRes, numRes);
    chooseComparisonParameters(Ra);
    if (comparable) {
      std::fstream fs;
      setRayleighDependentParameters(Ra);
      if (singleton::mpi().isMainProcessor()) {
        fs.open("convergence.csv",
          std::fstream::in | std::fstream::out | std::fstream::app);
        fs << "N;Ra;xVel_sim;xVel_ref;xVel_err;yVel_sim;yVel_ref;yVel_err;yVel_sim/xVel_sim;yVel_ref/xVel_ref;yVel/xVel_err;xmax_sim;xmax_ref;xmax_err;ymax_sim;ymax_ref;ymax_err;Nusselt_sim;Nusselt_ref;Nusselt_err" << std::endl;
      }
      for (int testId = 0; testId < numRes; testId++)
      {
        N = convergenceRes[testId];
        simulate(output);
        // if (singleton::mpi().isMainProcessor()) {
        //  fs.open("convergence.csv",
        //     std::fstream::in | std::fstream::out | std::fstream::app);
          fs  <<     N  << ";"
          <<        Ra  << ";"
          << output[0]  << ";"
          << output[1]  << ";"
          << output[2]  << ";"
          << output[3]  << ";"
          << output[4]  << ";"
          << output[5]  << ";"
          << output[6]  << ";"
          << output[7]  << ";"
          << output[8]  << ";"
          << output[9]  << ";"
          << output[10] << ";"
          << output[11] << ";"
          << output[12] << ";"
          << output[13] << ";"
          << output[14] << ";"
          << output[15] << ";"
          << output[16] << ";"
          << output[17] << ";"
          << std::endl;
        }
      // }
      fs.close();
      } else {
        clout << "No comparable values for Ra = " << Ra << " found in De Vahl Davis (1983)." << std::endl;
      }
  }
  else
  {
    clout << "no convergence study." << std::endl;
    Ra  = args.getValueOrFallback<std::size_t>("--Ra", 1000);
    N  = args.getValueOrFallback<std::size_t>("--N", 64);
    chooseComparisonParameters(Ra);
    setRayleighDependentParameters(Ra);
    simulate(output);
  }
}
