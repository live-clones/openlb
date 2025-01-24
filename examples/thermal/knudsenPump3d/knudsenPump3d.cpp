/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Anas Selmi
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

/* knudsenPump3d.cpp:
 * example of Knudsen pump effect to test thermal creep boundary
 * Left wall: Cold
 * Right wall: Hot
 * Thermal creep: from Cold to hot
 *
 * G. Karniadakis, A. Beskok and N. Aluru, Microflows and Nanoflows: Fundamentals and Simulation, 2005.
 */

#include "analyticalSolutionThermalCreep.h"
#include "olb3D.h"
#include "olb3D.hh"


using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;

using TDESCRIPTOR = D3Q7<VELOCITY>;
using NSDESCRIPTOR =
    D3Q19<NORMAL, VELOCITY2, AVERAGE_VELOCITY, BOUZIDI_SLIP_CREEP>;

// simple: simple cylinder with temperature difference
// buffer: added Buffers with constant temperature on the left and right
typedef enum { simple, buffer } Geometry;
Geometry geometryType = simple;
// Parameters for the simulation setup
const T   lx        = 2e-6; // length of the channel
const T   diameter  = 1e-6;
const T   ly        = 1e-6;
const T   lz        = ly;
const T   radius    = diameter / 2;
const T   radiusBig = 1e-6;
T         length    = lx + 2e-6;
const int N         = 21;         // resolution of the model
const T   Ra        = 1e6;        // Rayleigh number
const T   Pr        = 0.71;       // Prandtl number
const T   maxPhysT  = 5e-7 * 0.4; // max. simulation time in s, SI unit
const T   epsilon   = 1.e-5;      // precision of the convergence (residuum)
T         Kn        = 0.053;
T         dp        = 0.;
//const T gamma = 1.4;        // specific heat ratio
const T PI                = (T)3.141592653589793238463;
const T R_specific        = 287;
const T dynamic_viscosity = 18.13e-6; // [Pa.s]
const T DT                = 2 * 0.5 / (lx * 1e6);
const T Tcold             = 293;        // temperature of the left wall in K
const T Thot              = Tcold + DT; // temperature of the right wall in K

void prepareGeometry(
    SuperGeometry<T, 3>&                                superGeometry,
    ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR>& converter)
{

  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;
  const T dx = converter.getPhysDeltaX();
  if (geometryType == simple) {
    length = lx;
  }
  Vector<T, 3>                     center0(-length / 2 - dx * 0.5, 0, 0);
  Vector<T, 3>                     center1(-lx / 2 - dx * 0.5, 0, 0);
  Vector<T, 3>                     center2(lx / 2 + dx * 0.5, 0, 0);
  Vector<T, 3>                     center3(length / 2 + dx * 0.5, 0, 0);
  std::shared_ptr<IndicatorF3D<T>> pipe =
      std::make_shared<IndicatorCylinder3D<T>>(center1, center2, radius);
  std::shared_ptr<IndicatorF3D<T>> leftContainer =
      std::make_shared<IndicatorCylinder3D<T>>(center0, center1,
                                               radiusBig + dx / 2);
  std::shared_ptr<IndicatorF3D<T>> rightContainer =
      std::make_shared<IndicatorCylinder3D<T>>(center2, center3,
                                               radiusBig + dx / 2);

  IndicatorIdentity3D<T> pump(leftContainer + pipe + rightContainer);
  // Sets material number for fluid and boundary
  superGeometry.rename(0, 2);
  if (geometryType == buffer) {
    superGeometry.rename(2, 1, pump);
  }
  else {
    superGeometry.rename(2, 1, pipe);
  }

  // set material number for right and left walls

  // Set material number for left wall
  if (geometryType == buffer) {

    IndicatorLayer3D<T> leftC(leftContainer, dx);
    superGeometry.rename(2, 3, leftC);
    // Set material number for right wall

    IndicatorLayer3D<T> rightC(rightContainer, dx);
    superGeometry.rename(2, 4, rightC);
    center1[0] -= -4 * dx;
    center2[0] += 4 * dx;
    //IndicatorCylinder3D<T> pipeextend(center1, center2, radius+dx);
    //superGeometry.rename(3,2,pipeextend);
    //superGeometry.rename(4,2,pipeextend);
  }
  else if (geometryType == simple) {
    Vector<T, 3> origin(-lx / 2, 0, 0);
    Vector<T, 3> extend = origin;
    origin[0]           = -lx / 2 - 3 * dx;
    extend[0]           = -lx / 2;
    IndicatorCylinder3D<T> leftWall(origin, extend, radius);
    superGeometry.rename(2, 3, leftWall);

    origin[0] = lx / 2 - 3 * converter.getPhysDeltaX();
    extend[0] = lx / 2 + 3 * converter.getPhysDeltaX();
    IndicatorCylinder3D<T> rightWall(origin, extend, radius);
    superGeometry.rename(2, 4, rightWall);
  }

  /*
  // setting corners to Matr nbr 2
  center1[0] -= 4*converter.getPhysDeltaX();
  center2[0] += 4*converter.getPhysDeltaX();
  IndicatorCylinder3D<T> pipe2(center1, center2, radius);
  superGeometry.rename(3,2,pipe2);
  superGeometry.rename(4,2,pipe2);
  */

  superGeometry.clean();
  //superGeometry.rename(1, 13, *leftContainer);
  //superGeometry.rename(1, 14, *rightContainer);

  /// Removes all not needed boundary voxels outside the surface
  //superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(
    ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR>& converter,
    SuperLattice<T, NSDESCRIPTOR>&                      NSlattice,
    SuperLattice<T, TDESCRIPTOR>& ADlattice, SuperGeometry<T, 3>& superGeometry)
{

  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  T Tomega  = converter.getLatticeThermalRelaxationFrequency();
  T NSomega = converter.getLatticeRelaxationFrequency();

  /// define lattice Dynamics
  clout << "defining dynamics" << std::endl;

  auto walls = superGeometry.getMaterialIndicator({2, 3, 4});
  auto bulk  = superGeometry.getMaterialIndicator(1);
  ADlattice.defineDynamics<AdvectionDiffusionBGKdynamics>(
      superGeometry.getMaterialIndicator({1}));
  NSlattice.defineDynamics<BGKdynamics>(
      superGeometry.getMaterialIndicator({1}));
  //ADlattice.defineDynamics<ConstRhoBGKdynamics>(superGeometry.getMaterialIndicator({13,14}));
  //NSlattice(names::NavierStokes());
  //ADlattice(Temperature());
  Vector<T, 3>                     center0(-length / 2, 0, 0);
  Vector<T, 3>                     center1(-lx / 2, 0, 0);
  Vector<T, 3>                     center2(lx / 2, 0, 0);
  Vector<T, 3>                     center3(length / 2, 0, 0);
  std::shared_ptr<IndicatorF3D<T>> pipe =
      std::make_shared<IndicatorCylinder3D<T>>(center1, center2, radius);
  std::shared_ptr<IndicatorF3D<T>> leftContainer =
      std::make_shared<IndicatorCylinder3D<T>>(center0, center1, radiusBig);
  std::shared_ptr<IndicatorF3D<T>> rightContainer =
      std::make_shared<IndicatorCylinder3D<T>>(center2, center3, radiusBig);
  IndicatorIdentity3D<T> pump(leftContainer + pipe + rightContainer);

  clout << "defining ADE boundaries" << std::endl;
  /// set Temperature boundaries
  setBouzidiBoundary<T, TDESCRIPTOR, BouzidiAdeDirichletPostProcessor>(
      ADlattice, walls, bulk, pump);
  //setAdvectionDiffusionTemperatureBoundary(ADlattice, superGeometry, 13);
  //setAdvectionDiffusionTemperatureBoundary(ADlattice, superGeometry, 14);

  /// define initial conditions
  clout << "defining ADE initial conditions" << std::endl;

  AnalyticalConst3D<T, T>  u0(0.0, 0.0, 0.0);
  AnalyticalConst3D<T, T>  T_cold(converter.getLatticeTemperature(Tcold));
  AnalyticalConst3D<T, T>  T_hot(converter.getLatticeTemperature(Thot));
  T                        TCold = converter.getLatticeTemperature(Tcold);
  T                        THot  = converter.getLatticeTemperature(Thot);
  AnalyticalLinear3D<T, T> T_scaled((THot - TCold) / lx, 0, 0,
                                    (TCold + THot) / 2);

  ADlattice.defineRho(superGeometry, 1, T_scaled);
  ADlattice.iniEquilibrium(superGeometry, 1, T_scaled, u0);

  setBouzidiAdeDirichlet(ADlattice, superGeometry, 2, T_scaled);
  ADlattice.defineRho(superGeometry, 2, T_scaled);
  ADlattice.iniEquilibrium(superGeometry, 2, T_scaled, u0);

  setBouzidiAdeDirichlet(ADlattice, superGeometry, 3, T_cold);
  ADlattice.defineRho(superGeometry, 3, T_cold);
  ADlattice.iniEquilibrium(superGeometry, 3, T_cold, u0);

  setBouzidiAdeDirichlet(ADlattice, superGeometry, 4, T_hot);
  ADlattice.defineRho(superGeometry, 4, T_hot);
  ADlattice.iniEquilibrium(superGeometry, 4, T_hot, u0);
  /*
  //setBouzidiAdeDirichlet(ADlattice, superGeometry, 13, T_cold);
  ADlattice.defineRho(superGeometry, 13, T_cold);
  ADlattice.iniEquilibrium(superGeometry, 13, T_cold, u0);

  //setBouzidiAdeDirichlet(ADlattice, superGeometry, 14, T_hot);
  ADlattice.defineRho(superGeometry, 14, T_hot);
  ADlattice.iniEquilibrium(superGeometry, 14, T_hot, u0);*/
  clout << "defining NS boundaries" << std::endl;
  /// Set NS boundaries

  const T dx = converter.getPhysDeltaX();
  center1[0] -= 8 * dx;
  center2[0] += 8 * dx;
  IndicatorCylinder3D<T> pipe2(center1, center2, radius);

  clout << "Applied Bouzidi general slip boundary with Kn = " << Kn
        << std::endl;
  setBouzidiBoundary<T, NSDESCRIPTOR, BouzidiSlipVelocityPostProcessor>(
      NSlattice, superGeometry, 2, pipe2);
  //setBouzidiBoundary<T, NSDESCRIPTOR, BouzidiPostProcessor>(
  //    NSlattice, superGeometry, 2, pipe);
  setBouzidiSlipVelocity(NSlattice, superGeometry, 2, pipe2);
  NSlattice.setParameter<descriptors::BOUZIDI_TUNER>(T(0));
  NSlattice.setParameter<descriptors::LAMBDA>(Kn *
                                              converter.getCharPhysLength());
  NSlattice.setParameter<descriptors::CONVERSION_FACTOR_VELOCITY>(
      converter.getConversionFactorVelocity());
  NSlattice.setParameter<descriptors::CONVERSION_FACTOR_LENGTH>(
      converter.getConversionFactorLength());
  NSlattice.setParameter<descriptors::CHAR_LENGTH>(
      converter.getCharPhysLength());

  setBouzidiBoundary(NSlattice, superGeometry, 3, *leftContainer);
  setBouzidiBoundary(NSlattice, superGeometry, 4, *rightContainer);
  //setBouzidiBoundary<T, NSDESCRIPTOR, BouzidiSlipVelocityPostProcessor> (NSlattice, superGeometry, 3, *leftContainer);
  //setBouzidiSlipVelocity(NSlattice, superGeometry, 3, *leftContainer);

  //setBouzidiBoundary<T, NSDESCRIPTOR, BouzidiSlipVelocityPostProcessor> (NSlattice, superGeometry, 4, *rightContainer);
  //setBouzidiSlipVelocity(NSlattice, superGeometry, 4, *rightContainer);

  T rho_i = 1.;
  T rho_o = 1.;

  //AnalyticalLinear3D<T, T> rho( (rho_o - rho_i)/lx, 0, 0, rho_i);
  T p_m = converter.getLatticePressure(converter.getCharPhysPressure());
  //AnalyticalConst3D<T,T> rho(p_m * descriptors::invCs2<T,NSDESCRIPTOR>() + 1);
  AnalyticalConst3D<T, T> rho(1.);
  NSlattice.defineRhoU(superGeometry.getMaterialIndicator({1, 2, 3, 4}), rho,
                       u0);
  NSlattice.iniEquilibrium(superGeometry.getMaterialIndicator({1, 2, 3, 4}),
                           rho, u0);

  ADlattice.setParameter<descriptors::OMEGA>(Tomega);
  NSlattice.setParameter<descriptors::OMEGA>(NSomega);

  for (int iC = 0; iC < NSlattice.getLoadBalancer().size(); ++iC) {
    auto& blockLattice = NSlattice.getBlock(iC);
    blockLattice.forSpatialLocations([&](LatticeR<NSDESCRIPTOR::d> latticeR) {
      Vector<T,3> physR{};
      superGeometry.getBlockGeometry(iC).getPhysR(physR, latticeR);
      Vector<T,3>    origin = superGeometry.getBlockGeometry(iC).getOrigin();
      Vector<int,3>  originLat(converter.getLatticeLength(origin[0]),
                               converter.getLatticeLength(origin[1]),
                               converter.getLatticeLength(origin[2]));
      blockLattice.get(latticeR).template setField<descriptors::COORDINATE>(
          latticeR + originLat);
    });
  }
  /// Make the lattice ready for simulation
  NSlattice.initialize();
  ADlattice.initialize();

  clout << "Prepare Lattice: Communicate ..." << std::endl;

  {
    auto& communicator = NSlattice.getCommunicator(stage::PostPostProcess());
    communicator.requestField<descriptors::POPULATION>();
    communicator.requestField<descriptors::NORMAL>();
    communicator.requestField<descriptors::VELOCITY>();
    communicator.requestField<descriptors::BOUZIDI_DISTANCE>();
    communicator.requestField<descriptors::BOUZIDI_VELOCITY>();
    communicator.requestField<descriptors::BOUZIDI_SLIP_CREEP>();
    communicator.requestOverlap(NSlattice.getOverlap());
    communicator.exchangeRequests();
  }
  {
    auto& communicator = ADlattice.getCommunicator(stage::PostPostProcess());
    communicator.requestField<descriptors::POPULATION>();
    communicator.requestField<descriptors::VELOCITY>();
    communicator.requestOverlap(ADlattice.getOverlap());
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

void getResults(ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR>& converter,
                SuperLattice<T, NSDESCRIPTOR>&                      NSlattice,
                SuperLattice<T, TDESCRIPTOR>& ADlattice, int iT,
                SuperGeometry<T, 3>& superGeometry, util::Timer<T>& timer,
                bool converged)
{
  // update temperature field for boundaries
  //SuperLatticePhysTemperature3D<T, NSDESCRIPTOR, TDESCRIPTOR> temperature(ADlattice, converter);
  OstreamManager      clout(std::cout, "getResults");
  SuperVTMwriter3D<T> vtkWriter("knudsenPump3d");

  if (iT == 0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization

    SuperGeometryF<T, 3>       geometry(superGeometry);
    SuperLatticeCuboid3D<T, NSDESCRIPTOR> cuboid(NSlattice);
    SuperLatticeRank3D<T, NSDESCRIPTOR>   rank(NSlattice);
    vtkWriter.write(geometry);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
  }

  const int statIter = converter.getLatticeTime(maxPhysT / 10);
  const int saveIter = converter.getLatticeTime(maxPhysT / 20);

  if (iT % statIter == 0 || converged) {
    /// Timer console output
    timer.update(iT);
    timer.printStep();

    /// Lattice statistics console output
    NSlattice.getStatistics().print(iT, converter.getPhysTime(iT));
    ADlattice.getStatistics().print(iT, converter.getPhysTime(iT));
  }

  /// Writes the VTK files and prints statistics
  if (iT % saveIter == 0 || converged) {
    ADlattice.setProcessingContext(ProcessingContext::Evaluation);
    NSlattice.setProcessingContext(ProcessingContext::Evaluation);
    ThermalCreep3d<T, NSDESCRIPTOR, TDESCRIPTOR> uSol(6.4, R_specific, DT, Kn,
                                                      lx, converter);

    SuperVTMwriter3D<T>        vtkWriter("knudsenPump3d");
    SuperLatticePhysVelocity3D velocity(NSlattice, converter);
    SuperLatticePhysPressure3D pressure(NSlattice, converter);
    SuperLatticePhysTemperature3D<T, NSDESCRIPTOR, TDESCRIPTOR> temperature(
        ADlattice, converter);
    SuperLatticeField3D<T, NSDESCRIPTOR, descriptors::NORMAL> normal(NSlattice);
    SuperGeometryF<T, 3> materials(superGeometry);
    //SuperLatticeField3D<T, NSDESCRIPTOR, descriptors::COORDINATE> coords(NSlattice);
    //SuperLatticeField3D<T, NSDESCRIPTOR, descriptors::VELOCITY2> slip(NSlattice);
    //SuperLatticeField3D<T, NSDESCRIPTOR, descriptors::AVERAGE_VELOCITY> creep(NSlattice);
    //SuperLatticeField3D<T, NSDESCRIPTOR, descriptors::LOCATION> tang2(NSlattice);
    SuperLatticeFfromAnalyticalF3D<T, NSDESCRIPTOR> analyticalVelocityLattice(
        uSol, NSlattice);
    //SuperLatticeField3D<T, NSDESCRIPTOR, descriptors::BOUZIDI_DISTANCE> qBouzidi(NSlattice);
    //SuperLatticeField3D<T, NSDESCRIPTOR, descriptors::BOUZIDI_VELOCITY> vBouzidi(NSlattice);
    //SuperLatticeField3D<T, NSDESCRIPTOR, descriptors::BOUZIDI_SLIP_CREEP> cBouzidi(NSlattice);

    vtkWriter.addFunctor(pressure);
    vtkWriter.addFunctor(velocity);
    vtkWriter.addFunctor(temperature);
    vtkWriter.addFunctor(normal, "normal vector");
    vtkWriter.addFunctor(materials);
    //vtkWriter.addFunctor(coords, "Point coordinates");
    //vtkWriter.addFunctor(slip, "slip velocity");
    //vtkWriter.addFunctor(creep, "creep velocity");
    //vtkWriter.addFunctor(tang2, "tang2");
    vtkWriter.addFunctor(analyticalVelocityLattice, "Analyt sol");
    //vtkWriter.addFunctor(qBouzidi, "bouzidiDist");
    //vtkWriter.addFunctor(vBouzidi, "bouzidiVelocity");
    //vtkWriter.addFunctor(cBouzidi, "bouzidiCreep");
    //task(vtkWriter, iT);
    vtkWriter.write(iT);

    if (iT > 0) {
      NSlattice.communicate();
      SuperLatticePhysVelocity3D velocity(NSlattice, converter);
      AnalyticalFfromSuperF3D<T> intpolateVelocity(velocity, true, 1);

      std::string filename = "centerVelocity" + std::to_string(iT);
      T           numerical[3] {};
      T           analytical[3] {};
      T           point[3] {};
      T           D  = converter.getLatticeLength(diameter);
      T           dx = converter.getPhysDeltaX();

      Gnuplot<T> gplot(filename);
      point[0] = 0.5e-6;
      point[2] = 0.;
      for (int iY = -D / 2; iY < D / 2; ++iY) {
        point[1] = (T)converter.getPhysLength(iY);
        uSol(analytical, point);
        intpolateVelocity(numerical, point);
        gplot.setData(iY * dx * 2e6, {analytical[0], numerical[0]},
                      {"analytical", "numerical"});
      }
      // Create PNG file
      gplot.writePNG();
    }
  }
  const bool lastTimeStep =
      (converged || (iT + 1 == converter.getLatticeTime(maxPhysT)));
  if (lastTimeStep) {
    T L  = converter.getLatticeLength(length / 2);
    T dx = 1. / T(converter.getResolution());
    T point[3] {};
    point[1] = 0;
    point[2] = 0;

    SuperLatticePhysPressure3D<T, NSDESCRIPTOR> pressure(NSlattice, converter);
    SuperLatticePhysTemperature3D<T, NSDESCRIPTOR, TDESCRIPTOR> temperature(
        ADlattice, converter);

    AnalyticalFfromSuperF3D<T> pressure_intp(pressure, true, 1);
    AnalyticalFfromSuperF3D<T> temperature_intp(temperature, true, 1);
    point[0] = 0;
    T pmean[3] {};
    pressure_intp(pmean, point);
    T cnt = 6. * dynamic_viscosity * dynamic_viscosity * R_specific /
            ((radius * radius * pmean[0]) * (1 + 8 * Kn));
    T Dp = cnt * (Thot - Tcold);

    AnalyticalLinear3D<T, T> p_sol(Dp / lx, 0, 0, pmean[0]);
    T                        analytical[3] {};
    T                        numerical[3] {};
    T                        p_analyt = 0;
    T                        p_sim    = 0;
    T                        temp[3];
    static Gnuplot<T>        gplot("pressure_along x_axis");
    CSV<T>                   csvWriter("pressure_Temperature_alongX",
                                       {"x", "pressure", "temperature", "cnst=dp/dT"});
    for (int iX = -L; iX <= L; ++iX) {
      point[0] = (T)converter.getPhysLength(iX);
      p_sol(analytical, point);
      pressure_intp(numerical, point);
      temperature_intp(temp, point);
      p_analyt = 2 * (analytical[0] - pmean[0]) / Dp;
      p_sim    = 2 * (numerical[0] - pmean[0]) / Dp;
      gplot.setData(iX * dx, {p_analyt, p_sim}, {"analytical", "numerical"});
      csvWriter.writeDataFile(iX * dx, {numerical[0], temp[0], cnt});
    }
    gplot.writePNG();
    clout << "The mean pressure is " << pmean[0] << std::endl;
    clout << "The viscosity is " << dynamic_viscosity << " and Kn = " << Kn
          << std::endl;
    clout << "The constant to change temperature to pressure is " << cnt
          << std::endl;
  }
}

int main(int argc, char* argv[])
{
  /// === 1st Step: Initialization ===
  OstreamManager clout(std::cout, "main");
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

  // Computing physical parameters for unitConverter
  T Taverage            = (Tcold + Thot) / 2;
  T lambda              = Kn * diameter;
  T kinematic_viscosity = 0.5 * lambda * sqrt(8 * R_specific * Taverage / PI);
  T rho                 = dynamic_viscosity / kinematic_viscosity;
  T p_m                 = rho * R_specific * Taverage;
  T maxVelocity =
      0.75 * dynamic_viscosity * R_specific / p_m * (Thot - Tcold) / lx;
  T dx = diameter / N;
  T dt = (0.25 / 3.) * (dx * dx / kinematic_viscosity);
  //dt = 1e-11;
  clout << "The maximum thermal creep velocity is " << maxVelocity << std::endl;
  clout << "The mean free path is " << lambda << std::endl;
  clout << "The kinematic viscosity is" << kinematic_viscosity << std::endl;
  clout << "The mean pressure is " << p_m << std::endl;
  clout << "The mean temperature is " << Taverage << std::endl;
  ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> converter(
      (T)diameter / N, // physDeltaX
      (T)dt, // physDeltaT = charLatticeVelocity / charPhysVelocity * physDeltaX
      (T)diameter,            // charPhysLength
      (T)maxVelocity,         // charPhysVelocity
      (T)kinematic_viscosity, // physViscosity
      (T)rho,                 // physDensity
      (T)0.015,               // physThermalConductivity
      (T)Pr * 0.015 /
          dynamic_viscosity, // physSpecificHeatCapacity = Pr* physThermalConductivity / dynamic_viscosity [J / kg K]
      (T)1 / Taverage,       // physThermalExpansionCoefficient (ideal gas)
      (T)Tcold,              // charPhysLowTemperature
      (T)Thot,               // charPhysHighTemperature
      (T)p_m                 // charPhysPressure
  );
  converter.print();

  /// === 2nd Step: Prepare Geometry ===
  //Vector<T,3> center0(0, radius, radius);
  //Vector<T,3> center1(lx+0.5*converter.getPhysDeltaX(), radius, radius);
  //IndicatorCylinder3D<T> pipe(center0, center1, radius);
  if (geometryType == simple) {
    length = lx;
  }
  Vector<T, 3> center0(-length / 2, 0, 0);
  Vector<T, 3> center1(length / 2 + 0.5 * converter.getPhysDeltaX(), 0, 0);
  T            radius_sim = radius;
  if (geometryType == buffer) {
    radius_sim = radiusBig + converter.getPhysDeltaX() / 2;
  }
  IndicatorCylinder3D<T> pipe(center0, center1, radius_sim);
  IndicatorLayer3D<T>    extendedDomain(pipe, converter.getPhysDeltaX());

#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  /// Instantiation of a cuboidGeometry with weights
  CuboidGeometry3D<T> cuboidGeometry(extendedDomain, converter.getPhysDeltaX(),
                                     noOfCuboids);

  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  const int           overlap = 4;
  SuperGeometry<T, 3> superGeometry(cuboidGeometry, loadBalancer, overlap);

  prepareGeometry(superGeometry, converter);
  superGeometry.communicate();

  /// === 3rd Step: Prepare Lattice ===

  SuperLattice<T, TDESCRIPTOR>  ADlattice(superGeometry);
  SuperLattice<T, NSDESCRIPTOR> NSlattice(superGeometry);

  util::Timer<T> timer(converter.getLatticeTime(maxPhysT),
                       superGeometry.getStatistics().getNvoxel());
  getResults(converter, NSlattice, ADlattice, 0, superGeometry, timer, false);

  prepareLattice(converter, NSlattice, ADlattice, superGeometry);

  // Only coupling for velocity needed here, coupling for temperature done through boundary
  SuperLatticeCoupling coupling(ThermalCreepBouzidiCoupling {},
                                names::NavierStokes {}, NSlattice,
                                names::Temperature {}, ADlattice);
  coupling.restrictTo(superGeometry.getMaterialIndicator({1, 2}));
  coupling.setParameter<ThermalCreepBouzidiCoupling::R>(R_specific);
  coupling.setParameter<ThermalCreepBouzidiCoupling::MU>(dynamic_viscosity);
  coupling
      .setParameter<ThermalCreepBouzidiCoupling::CONVERSION_FACTOR_PRESSURE>(
          converter.getConversionFactorPressure());
  coupling.setParameter<ThermalCreepBouzidiCoupling::CHAR_PHYS_PRESSURE>(
      converter.getCharPhysPressure());
  coupling
      .setParameter<ThermalCreepBouzidiCoupling::CONVERSION_FACTOR_TEMPERATURE>(
          converter.getConversionFactorTemperature());
  coupling.setParameter<descriptors::CONVERSION_FACTOR_VELOCITY>(
      converter.getConversionFactorVelocity());
  coupling.setParameter<descriptors::CONVERSION_FACTOR_LENGTH>(
      converter.getConversionFactorLength());
  Vector<int, 3> N_xyz(converter.getLatticeLength(lx),
                       converter.getLatticeLength(ly),
                       converter.getLatticeLength(lz));
  coupling.setParameter<ThermalCreepBouzidiCoupling::N_XYZ>(N_xyz);
  clout << "The vector N_XYZ is " << N_xyz << std::endl;
  /// === 4th Step: Main Loop with Timer ===
  timer.start();

  util::ValueTracer<T> converge(converter.getLatticeTime(50.), epsilon);
  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {
    if (converge.hasConverged()) {
      clout << "Simulation converged." << std::endl;

      getResults(converter, NSlattice, ADlattice, iT, superGeometry, timer,
                 converge.hasConverged());

      clout << "Time " << iT << "." << std::endl;

      break;
    }
    //setBoundaryValues(NSlattice, ADlattice, superGeometry, converter);
    getResults(converter, NSlattice, ADlattice, iT, superGeometry, timer,
               converge.hasConverged());

    /// === 6th Step: Collide and Stream Execution ===
    ADlattice.collideAndStream();
    NSlattice.collideAndStream();

    NSlattice.getCommunicator(stage::PreCoupling()).communicate();
    ADlattice.getCommunicator(stage::PreCoupling()).communicate();
    coupling.execute();

    NSlattice.getCommunicator(stage::PostCoupling()).communicate();
    ADlattice.getCommunicator(stage::PostCoupling()).communicate();

    /// === 7th Step: Computation and Output of the Results ===
    //getResults
    converge.takeValue(ADlattice.getStatistics().getAverageEnergy(), true);
  }

  NSlattice.setProcessingContext(ProcessingContext::Evaluation);
  ADlattice.setProcessingContext(ProcessingContext::Evaluation);

  timer.stop();
  timer.printSummary();
}
