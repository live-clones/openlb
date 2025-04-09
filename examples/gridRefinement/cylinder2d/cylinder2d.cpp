/*  Lattice Boltzmann sample, written in C++, using the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
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

#include "olb2D.h"
#include "olb2D.hh"

using namespace olb;
using namespace olb::descriptors;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = D2Q9<>;

struct GeometryParameters {
  T L;
  T lengthX;
  T lengthY;
  T centerCylinderX;
  T centerCylinderY;
  T radiusCylinder;

  GeometryParameters(int N) {
    L = 0.1 / N;
    lengthX = 2.2;
    lengthY = 0.41 + L;
    centerCylinderX = 0.2;
    centerCylinderY = 0.2 + L / 2.0;
    radiusCylinder = 0.05;
  }
};


void prepareGeometry(const UnitConverter<T,DESCRIPTOR>& converter,
                     SuperGeometry<T,2>& sGeometry,
                     std::shared_ptr<IndicatorF2D<T>> circle,
                     GeometryParameters geomParams)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  Vector<T,2> extend(geomParams.lengthX, geomParams.lengthY);
  Vector<T,2> origin;

  sGeometry.rename(0,2);
  sGeometry.rename(2,1,{1,1});

  // Set material number for inflow
  extend[0] = 2.*geomParams.L;
  origin[0] = -geomParams.L;
  IndicatorCuboid2D<T> inflow( extend, origin );
  sGeometry.rename( 2,3,1,inflow );

  // Set material number for outflow
  origin[0] = geomParams.lengthX - geomParams.L;
  IndicatorCuboid2D<T> outflow( extend, origin );
  sGeometry.rename( 2,4,1,outflow );
  // Set material number for cylinder
  // sGeometry.rename( 1,5, circle );

  // Removes all not needed boundary voxels outside the surface
  sGeometry.clean();
  sGeometry.checkForErrors();

  sGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareGeometryMiddle(const UnitConverter<T,DESCRIPTOR>& converter,
                         SuperGeometry<T,2>& sGeometry)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  sGeometry.rename(0,1);

  sGeometry.clean();
  sGeometry.checkForErrors();

  sGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareGeometryFine(const UnitConverter<T,DESCRIPTOR>& converter,
                         SuperGeometry<T,2>& sGeometry,
                         std::shared_ptr<IndicatorF2D<T>> circle,
                         GeometryParameters geomParams)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  Vector<T,2> extend(geomParams.lengthX, geomParams.lengthY);
  Vector<T,2> origin;

  sGeometry.rename(0,1);
  sGeometry.rename(1,5, circle);

  sGeometry.clean();
  sGeometry.checkForErrors();

  sGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(SuperLattice<T,DESCRIPTOR>& sLattice,
                    const UnitConverter<T,DESCRIPTOR>& converter,
                    SuperGeometry<T,2>& sGeometry,
                    std::shared_ptr<IndicatorF2D<T>> circle)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  auto bulkIndicator = sGeometry.getMaterialIndicator({1});
  sLattice.defineDynamics<BGKdynamics>(bulkIndicator);

  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(sLattice, sGeometry, 2);

  // Setting of the boundary conditions
  // if boundary conditions are chosen to be local
  //boundary::set<boundary::LocalVelocity>(sLattice, sGeometry, 3);
  //boundary::set<boundary::LocalPressure>(sLattice, sGeometry, 4);

  //if boundary conditions are chosen to be interpolatedy, 3);
  boundary::set<boundary::InterpolatedVelocity>(sLattice, sGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, sGeometry, 4);

  // if boundary conditions are chosen to be Zou He type
  //boundary::set<boundary::ZouHeVelocity>(sLattice, sGeometry, 3);
  //boundary::set<boundary::ZouHePressure>(sLattice, sGeometry, 4);

  // Material=5 -->bouzidi / bounce back
  #ifdef BOUZIDI
  setBouzidiBoundary(sLattice, sGeometry, 5, *circle);
  #else
  boundary::set<boundary::BounceBack>(sLattice, sGeometry, 5);
  #endif

  // Initial conditions
  AnalyticalConst2D<T,T> rhoF( 1 );
  std::vector<T> velocity( 2,T( 0 ) );
  AnalyticalConst2D<T,T> uF( velocity );

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, rhoF, uF );
  sLattice.iniEquilibrium( bulkIndicator, rhoF, uF );

  sLattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues(SuperLattice<T,DESCRIPTOR>& sLattice,
                       const UnitConverter<T,DESCRIPTOR>& converter,
                       std::size_t iT,
                       SuperGeometry<T,2>& sGeometry,
                       size_t iTmaxStart,
                       T maxLatticeU,
                       T maxPhysU,
                       T iniPhysU)
{

  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  std::size_t iTupdate = 20;

  if (iT % iTupdate == 0 && iT<= iTmaxStart) {
    PolynomialStartScale<T,T> StartScale( iTmaxStart, T( maxLatticeU - iniPhysU ) );

    // Creates and sets the Poiseuille inflow profile using functors
    T iTvec[1] = {T( iT )};
    T frac[1] = {};
    StartScale( frac,iTvec );
    T maxVelocity = converter.getCharLatticeVelocity()*3./2.*frac[0];
    T distance2Wall = converter.getPhysDeltaX()/2.;
    Poiseuille2D<T> poiseuilleU( sGeometry, 3, maxVelocity, distance2Wall );

    sLattice.defineU( sGeometry, 3, poiseuilleU );

    clout << "startup step=" << iT << "/" << iTmaxStart << "; t=" << converter.getPhysTime(iT)
          << "; ux_LU=" << maxVelocity << "; maxLatticeU=" << maxLatticeU
          << "; ux_PU=" << converter.getPhysVelocity( maxVelocity ) << "; ux_max_PU=" << maxPhysU << std::endl;

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

void writeResults(SuperLattice<T, DESCRIPTOR>& sLattice,
                  const UnitConverter<T,DESCRIPTOR>& converter,
                  std::size_t iT,
                  SuperGeometry<T,2>& sGeometry,
                  std::string name,
                  std::string vtmName="",
                  bool vtk=false,
                  bool images=false)
{
  OstreamManager clout(std::cout, name);

  if ( iT == 0 && vtk ) {
    SuperVTMwriter2D<T> vtmWriter(vtmName, 0);
    if (name == "level0") {
      SuperLatticeRank2D rank( sLattice );
      vtmWriter.write(rank);
    }
    SuperGeometryF<T,DESCRIPTOR::d> geometryF(sGeometry);
    geometryF.getName() = name + "_geometry";
    vtmWriter.write(geometryF);
    SuperLatticeCuboid2D<T,DESCRIPTOR> cuboidF(sLattice);
    cuboidF.getName() = name + "_cuboid";
    vtmWriter.write(cuboidF);
    vtmWriter.createMasterFile();
  }

  sLattice.setProcessingContext(ProcessingContext::Evaluation);
  // Background not implemented for 2D
  SuperVTMwriter2D<T> vtmWriter(vtmName);
  SuperLatticePhysVelocity2D velocityF(sLattice, converter);
  SuperLatticePhysPressure2D pressureF(sLattice, converter);
  SuperLatticeRefinementMetricKnudsen2D qualityF(sLattice, converter);
  SuperLatticePhysField2D<T,DESCRIPTOR,fields::refinement::PREV_RHO> prev_rho( sLattice, 1. );
  prev_rho.getName() = "PREV_RHO";
  vtmWriter.addFunctor(prev_rho);
  vtmWriter.addFunctor(qualityF);
  vtmWriter.addFunctor(velocityF);
  vtmWriter.addFunctor(pressureF);
  vtmWriter.write(iT);

  if ( images ) {
    SuperEuklidNorm2D<T,DESCRIPTOR> normVel( velocityF );
    BlockReduction2D2D<T> planeReductionU( normVel );
    // write output as JPEG
    heatmap::plotParam<T> uParams;
    uParams.name = "velocity_"+name;
    heatmap::write( planeReductionU, iT, uParams) ;
    
    BlockReduction2D2D<T> planeReductionP( pressureF );
    heatmap::plotParam<T> pParams;
    pParams.name = "pressure_"+name;
    heatmap::write( planeReductionP, iT, pParams );
  }
}

int main(int argc, char* argv[])
{
  olbInit(&argc, &argv);

  CLIreader args(argc, argv);
  std::string outdir        = args.getValueOrFallback<std::string>( "--outdir", "./tmp" );
  const int N               = args.getValueOrFallback<int>( "--res",        11);
  const int Re              = args.getValueOrFallback<int>( "--Re",         100);
  const int maxPhysT        = args.getValueOrFallback<int>( "--maxPhysT",   4);
  const T tMaxInit          = args.getValueOrFallback<T>  ( "--iniPhysT",   2);
  size_t iTmax              = args.getValueOrFallback<int>( "--iTmax",      0);
  size_t nout               = args.getValueOrFallback( "--nout",            5);     // minimum number of vtk outputs
  size_t iTout              = args.getValueOrFallback( "--iTout",           0);     // iterations for vtk outputs
  T physTout                = args.getValueOrFallback( "--physTout",        0);     // timestep for vtk outputs
  T maxPhysU                = args.getValueOrFallback<T>  ( "--maxPhysU",   1.);    // in m/s
  T iniPhysU                = args.getValueOrFallback<T>  ( "--iniPhysU",   .0);
  bool allVtk               = args.contains( "--allVtk" );  // vtk for all substeps
  bool subImages          = args.contains( "--subImages" );  // vtk for all substeps

  std::stringstream outdir_mod;
  outdir_mod << outdir;
  outdir_mod << "_u" << maxPhysU << "_Re" << Re << "_res" << N << "_iTmax" << iTmax << "_iTout" << iTout;
  singleton::directories().setOutputDir( outdir_mod.str()+"/" );  // set output directory
  std::ofstream fileStream( outdir_mod.str() + "/output.txt" );
  DoubleBuffer doubleBuffer( std::cout.rdbuf(), fileStream.rdbuf() );
  std::streambuf* originalCoutBuffer = std::cout.rdbuf( &doubleBuffer );
  OstreamManager clout( std::cout, "main" );
  clout << "outdir specified to " << outdir_mod.str() << std::endl;

  GeometryParameters geomParams(N);

  const UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> converterLevel0(
    int {N},                                   // resolution: number of voxels per charPhysL
    (T)   0.51,                                // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   2.0*geomParams.radiusCylinder,       // charPhysLength: reference length of simulation geometry
    (T)   maxPhysU,                            // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   0.2*2.*geomParams.radiusCylinder/Re, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0                                  // physDensity: physical density in __kg / m^3__
  );
  converterLevel0.print();
  T maxLatticeU                 = converterLevel0.getLatticeVelocity( maxPhysU );

  Vector<T,2> extend( geomParams.lengthX, geomParams.lengthY );
  Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid(extend, origin);

  CuboidGeometry2D<T> cGeometryLevel0(cuboid, geomParams.L);
  cGeometryLevel0.splitFractional(0, 0, {0.025,0.5,0.475});
  cGeometryLevel0.splitFractional(1, 1, {0.15,0.7,0.15});
  cGeometryLevel0.splitFractional(3, 0, {0.05,0.5,0.45});
  cGeometryLevel0.splitFractional(5, 1, {0.15,0.7,0.15});

  Vector<T,2> center(geomParams.centerCylinderX, geomParams.centerCylinderY);
  std::shared_ptr<IndicatorF2D<T>> circle = std::make_shared<IndicatorCircle2D<T>>(center, geomParams.radiusCylinder);

  auto cGeometryLevel1 = cGeometryLevel0;
  IndicatorCuboid2D<T> refineCuboidI({1.0, 0.25}, {0.08,0.08});
  cGeometryLevel1.remove(refineCuboidI);
  cGeometryLevel1.refine(2);
  cGeometryLevel1.print();

  auto cGeometryLevel2 = cGeometryLevel1;
  cGeometryLevel2.remove(*circle);
  cGeometryLevel2.refine(2);
  cGeometryLevel2.print();

  // Adjust weights for balancing
  for (int iC=0; iC < cGeometryLevel0.size(); ++iC) {
    auto& cuboid = cGeometryLevel0.get(iC);
    cuboid.setWeight(cuboid.getLatticeVolume());

    auto origin = cGeometryLevel0.get(iC).getOrigin();

    for (int jC=0; jC < cGeometryLevel1.size(); ++jC) {
      if (cGeometryLevel1.get(jC).getOrigin() == origin) {
        cuboid.setWeight(cuboid.getWeight() + 2*cGeometryLevel1.get(jC).getLatticeVolume());
      }
    }
    for (int jC=0; jC < cGeometryLevel2.size(); ++jC) {
      if (cGeometryLevel2.get(jC).getOrigin() == origin) {
        cuboid.setWeight(cuboid.getWeight() + 4*cGeometryLevel2.get(jC).getLatticeVolume());
      }
    }
  }

  HeuristicLoadBalancer<T> loadBalancerLevel0(cGeometryLevel0);
  SuperGeometry<T,2> sGeometryLevel0(cGeometryLevel0, loadBalancerLevel0);
  prepareGeometry(converterLevel0, sGeometryLevel0, circle, geomParams);
  SuperLattice<T,DESCRIPTOR> sLatticeLevel0(cGeometryLevel0,
                                            loadBalancerLevel0,
                                            3,
                                            converterLevel0);
  prepareLattice(sLatticeLevel0, converterLevel0, sGeometryLevel0, circle);

  RefinedLoadBalancer<T,2> loadBalancerLevel1(cGeometryLevel0,
                                              loadBalancerLevel0,
                                              cGeometryLevel1);
  auto converterLevel1 = convectivelyRefineUnitConverter(converterLevel0, 2);
  converterLevel1.print();
  SuperGeometry<T,2> sGeometryLevel1(cGeometryLevel1, loadBalancerLevel1);
  prepareGeometryMiddle(converterLevel1, sGeometryLevel1);
  SuperLattice<T,DESCRIPTOR> sLatticeLevel1(cGeometryLevel1,
                                            loadBalancerLevel1,
                                            3,
                                            converterLevel1);
  prepareLattice(sLatticeLevel1, converterLevel1, sGeometryLevel1, circle);

  auto coarseToFineLevel1 = refinement::lagrava::makeCoarseToFineCoupler(
    sLatticeLevel0, sGeometryLevel0,
    sLatticeLevel1, sGeometryLevel1);
  auto fineToCoarseLevel1 = refinement::lagrava::makeFineToCoarseCoupler(
    sLatticeLevel0, sGeometryLevel0,
    sLatticeLevel1, sGeometryLevel1);

  RefinedLoadBalancer<T,2> loadBalancerLevel2(cGeometryLevel1,
                                              loadBalancerLevel1,
                                              cGeometryLevel2);
  auto converterLevel2 = convectivelyRefineUnitConverter(converterLevel1, 2);
  converterLevel2.print();
  SuperGeometry<T,2> sGeometryLevel2(cGeometryLevel2, loadBalancerLevel2);
  prepareGeometryFine(converterLevel2, sGeometryLevel2, circle, geomParams);
  SuperLattice<T,DESCRIPTOR> sLatticeLevel2(cGeometryLevel2,
                                            loadBalancerLevel2,
                                            3,
                                            converterLevel2);
  prepareLattice(sLatticeLevel2, converterLevel2, sGeometryLevel2, circle);

  auto coarseToFineLevel2 = refinement::lagrava::makeCoarseToFineCoupler(
    sLatticeLevel1, sGeometryLevel1,
    sLatticeLevel2, sGeometryLevel2);
  auto fineToCoarseLevel2 = refinement::lagrava::makeFineToCoarseCoupler(
    sLatticeLevel1, sGeometryLevel1,
    sLatticeLevel2, sGeometryLevel2);

  clout << "Initializing PREV_RHO ..." << std::endl;
  coarseToFineLevel1->initialize_prev(meta::id<refinement::lagrava::FullTimeCoarseToFineO>{});
  coarseToFineLevel2->initialize_prev(meta::id<refinement::lagrava::FullTimeCoarseToFineO>{});
  clout << "Initializing PREV_RHO ... OK" << std::endl;

  // iTmax depends on maximum physical time. If iTmax is provided in command line, it is an upper bound
  if ( iTmax == 0 ) iTmax = converterLevel0.getLatticeTime( maxPhysT );
  else iTmax = std::min( iTmax, converterLevel0.getLatticeTime( maxPhysT ) );
  size_t iTmaxStart = converterLevel0.getLatticeTime( tMaxInit );
  // nout is the minimum number of vtk outputs --> take max between nout and nout derived from iTout or physTout
  size_t nout_from_iTout = 0, nout_from_physTout = 0;
  if ( iTout != 0 ) { nout_from_iTout = size_t( iTmax / iTout ); nout = std::max( nout, nout_from_iTout ); }
  if ( physTout != 0 ) { nout_from_physTout = size_t( iTmax / physTout ); nout = std::max( nout, nout_from_physTout ); }
  size_t iTvtk = std::size_t( std::max( T(iTmax / nout), T(1) ) );
  clout << "Set nout to " << nout << ", so iTvtk=" << iTvtk << std::endl;
  
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer(converterLevel0.getLatticeTime(maxPhysT),
                           sGeometryLevel0.getStatistics().getNvoxel()
                       + 2*sGeometryLevel1.getStatistics().getNvoxel()
                       + 4*sGeometryLevel2.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT < iTmax; ++iT) {
    if (iT % iTvtk == 0) {
      writeResults(sLatticeLevel0, converterLevel0, iT, sGeometryLevel0, std::to_string(iT)+"_00_level0_preBC", "00_level0_preBC", true);
      writeResults(sLatticeLevel1, converterLevel1, iT, sGeometryLevel1, std::to_string(iT)+"_00_level1_preBC", "00_level1_preBC", true);
      writeResults(sLatticeLevel2, converterLevel2, iT, sGeometryLevel2, std::to_string(iT)+"_00_level2_preBC", "00_level2_preBC", true);

      setBoundaryValues( sLatticeLevel0, converterLevel0, iT, sGeometryLevel0, iTmaxStart, maxLatticeU, maxPhysU, iniPhysU );
      // writeResults(sLatticeLevel0, converterLevel0, iT, sGeometryLevel0, std::to_string(iT)+"_01_level0_postBC", "01_level0_postBC", allVtk, subImages);
      sLatticeLevel0.collideAndStream();
      // writeResults(sLatticeLevel0, converterLevel0, iT, sGeometryLevel0, std::to_string(iT)+"_02_level0_postCS", "02_level0_postCS", allVtk, subImages);
  
        sLatticeLevel1.collideAndStream();
        // writeResults(sLatticeLevel1, converterLevel1, iT, sGeometryLevel1, std::to_string(iT)+"_03_level1_postCS1", "03_level1_postCS1", allVtk, subImages);
        coarseToFineLevel1->apply(meta::id<refinement::lagrava::HalfTimeCoarseToFineO>{});
        writeResults(sLatticeLevel1, converterLevel1, iT, sGeometryLevel1, std::to_string(iT)+"_04_level1_postC2F1", "04_level1_postC2F1", allVtk, subImages);
  
          sLatticeLevel2.collideAndStream();
          // writeResults(sLatticeLevel2, converterLevel2, iT, sGeometryLevel2, std::to_string(iT)+"_05_level2_postCS1", "05_level2_postCS1", allVtk, subImages);
          coarseToFineLevel2->apply(meta::id<refinement::lagrava::HalfTimeCoarseToFineO>{});
          // writeResults(sLatticeLevel2, converterLevel2, iT, sGeometryLevel2, std::to_string(iT)+"_06_level2_postC2F1", "06_level2_postC2F1", allVtk, subImages);
  
          sLatticeLevel2.collideAndStream();
          // writeResults(sLatticeLevel2, converterLevel2, iT, sGeometryLevel2, std::to_string(iT)+"_07_level2_postCS2", "07_level2_postCS2", allVtk, subImages);
          coarseToFineLevel2->apply(meta::id<refinement::lagrava::FullTimeCoarseToFineO>{});
          // writeResults(sLatticeLevel2, converterLevel2, iT, sGeometryLevel2, std::to_string(iT)+"_08_level2_postC2F2", "08_level2_postC2F2", allVtk, subImages);
  
        fineToCoarseLevel2->apply(meta::id<refinement::lagrava::FineToCoarseO>{});
        writeResults(sLatticeLevel1, converterLevel1, iT, sGeometryLevel1, std::to_string(iT)+"_09_level1_postF2C1", "09_level1_postF2C1", allVtk, subImages);
        sLatticeLevel1.collideAndStream();
        // writeResults(sLatticeLevel1, converterLevel1, iT, sGeometryLevel1, std::to_string(iT)+"_10_level1_postCS2", "10_level1_postCS2", allVtk, subImages);
        coarseToFineLevel1->apply(meta::id<refinement::lagrava::FullTimeCoarseToFineO>{});
        writeResults(sLatticeLevel1, converterLevel1, iT, sGeometryLevel1, std::to_string(iT)+"_11_level1_postC2F2", "11_level1_postC2F2", allVtk, subImages);
  
          sLatticeLevel2.collideAndStream();
          // writeResults(sLatticeLevel2, converterLevel2, iT, sGeometryLevel2, std::to_string(iT)+"_12_level2_postCS3", "12_level2_postCS3", allVtk, subImages);
          coarseToFineLevel2->apply(meta::id<refinement::lagrava::HalfTimeCoarseToFineO>{});
          // writeResults(sLatticeLevel2, converterLevel2, iT, sGeometryLevel2, std::to_string(iT)+"_13_level2_postC2F3", "13_level2_postC2F3", allVtk, subImages);
  
          sLatticeLevel2.collideAndStream();
          // writeResults(sLatticeLevel2, converterLevel2, iT, sGeometryLevel2, std::to_string(iT)+"_14_level2_postCS4", "14_level2_postCS4", allVtk, subImages);
          coarseToFineLevel2->apply(meta::id<refinement::lagrava::FullTimeCoarseToFineO>{});
          // writeResults(sLatticeLevel2, converterLevel2, iT, sGeometryLevel2, std::to_string(iT)+"_15_level2_postC2F4", "15_level2_postC2F4", allVtk, subImages);
  
        fineToCoarseLevel2->apply(meta::id<refinement::lagrava::FineToCoarseO>{});
        // writeResults(sLatticeLevel1, converterLevel1, iT, sGeometryLevel1, std::to_string(iT)+"_16_level1_postF2C2", "16_level1_postF2C2", allVtk, subImages);
  
      fineToCoarseLevel1->apply(meta::id<refinement::lagrava::FineToCoarseO>{});
      // writeResults(sLatticeLevel0, converterLevel0, iT, sGeometryLevel0, std::to_string(iT)+"_17_level0_postF2C", "17_level0_postF2C", allVtk, subImages);
  
    } else {
      setBoundaryValues( sLatticeLevel0, converterLevel0, iT, sGeometryLevel0, iTmaxStart, maxLatticeU, maxPhysU, iniPhysU );
      sLatticeLevel0.collideAndStream();

      sLatticeLevel1.collideAndStream();
      coarseToFineLevel1->apply(meta::id<refinement::lagrava::HalfTimeCoarseToFineO>{});

      {
        sLatticeLevel2.collideAndStream();
        coarseToFineLevel2->apply(meta::id<refinement::lagrava::HalfTimeCoarseToFineO>{});

        sLatticeLevel2.collideAndStream();
        coarseToFineLevel2->apply(meta::id<refinement::lagrava::FullTimeCoarseToFineO>{});

        fineToCoarseLevel2->apply(meta::id<refinement::lagrava::FineToCoarseO>{});
      }

      sLatticeLevel1.collideAndStream();
      coarseToFineLevel1->apply(meta::id<refinement::lagrava::FullTimeCoarseToFineO>{});

      {
        sLatticeLevel2.collideAndStream();
        coarseToFineLevel2->apply(meta::id<refinement::lagrava::HalfTimeCoarseToFineO>{});

        sLatticeLevel2.collideAndStream();
        coarseToFineLevel2->apply(meta::id<refinement::lagrava::FullTimeCoarseToFineO>{});

        fineToCoarseLevel2->apply(meta::id<refinement::lagrava::FineToCoarseO>{});
      }

      fineToCoarseLevel1->apply(meta::id<refinement::lagrava::FineToCoarseO>{});
    }

    if (iT % std::max(int(iTvtk/2), 1) == 0) {
      timer.update(iT);
      timer.printStep();
      sLatticeLevel0.getStatistics().print(iT, converterLevel0.getPhysTime(iT));
    }
  }

  singleton::pool().wait();

  timer.stop();
  timer.printSummary();
  std::cout.rdbuf(originalCoutBuffer);
}
