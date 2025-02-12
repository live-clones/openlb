/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
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

#include "olb3D.h"
#include "olb3D.hh"

using namespace olb;
using namespace olb::descriptors;

using T = float;
using DESCRIPTOR = D3Q19<>;

void prepareGeometry(const UnitConverter<T,DESCRIPTOR>& converter,
                     SuperGeometry<T,3>& sGeometry,
                     Vector<T,3> extendDomain,
                     Vector<T,3> originDomain)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;
  
  IndicatorCuboid3D<T> cuboidDomain(extendDomain, originDomain);
  IndicatorLayer3D<T> domainLayer(cuboidDomain, converter.getPhysDeltaX());
  Vector<T,3> origin = domainLayer.getMin();
  Vector<T,3> extend = domainLayer.getMax() - domainLayer.getMin();

  sGeometry.rename( 0, 2, domainLayer );
  sGeometry.rename( 2, 1, cuboidDomain );
  sGeometry.clean();

  // Vector<T,3> origin = sGeometry.getStatistics().getMinPhysR(2);
  origin[1] += converter.getPhysDeltaX()/2.;
  origin[2] += converter.getPhysDeltaX()/2.;

  // Vector<T,3> extend = sGeometry.getStatistics().getMaxPhysR(2);
  extend[1] = extend[1]-origin[1]-converter.getPhysDeltaX()/2.;
  extend[2] = extend[2]-origin[2]-converter.getPhysDeltaX()/2.;

  // Set material number for inflow
  origin[0] = sGeometry.getStatistics().getMinPhysR(2)[0]-converter.getPhysDeltaX();
  extend[0] = 2*converter.getPhysDeltaX();
  IndicatorCuboid3D<T> inflow(extend, origin);
  sGeometry.rename(2,3, inflow);

  // Set material number for outflow
  origin[0] = sGeometry.getStatistics().getMaxPhysR(2)[0]-converter.getPhysDeltaX();
  extend[0] = 2*converter.getPhysDeltaX();
  IndicatorCuboid3D<T> outflow(extend, origin);
  sGeometry.rename(2,4,outflow);

  sGeometry.clean();
  sGeometry.checkForErrors();

  sGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareGeometryFine(const UnitConverter<T,DESCRIPTOR>& converter,
                         SuperGeometry<T,3>& sGeometry,
                         STLreader<T>& foilBody,
                         STLreader<T>& foilTail)
{
  OstreamManager clout( std::cout,"prepareGeometryFine");
  clout << "Prepare Geometry ..." << std::endl;

  sGeometry.rename(0,1);

  // IndicatorSphere3D<T> sphereI({0.55, 0.25, 0.25}, 0.1);
  sGeometry.rename(1,2, foilBody);
  sGeometry.rename(1,2, foilTail);

  sGeometry.clean();
  sGeometry.checkForErrors();

  sGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(SuperLattice<T,DESCRIPTOR>& sLattice,
                    const UnitConverter<T,DESCRIPTOR>& converter,
                    SuperGeometry<T,3>& sGeometry)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  auto bulkIndicator = sGeometry.getMaterialIndicator({1});
  sLattice.defineDynamics<BGKdynamics>(bulkIndicator);

  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(sLattice, sGeometry, 2);

  //if interpolated boundary conditions are chosen
  boundary::set<boundary::InterpolatedVelocity>(sLattice, sGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, sGeometry, 4);

  // Initial conditions
  AnalyticalConst3D<T,T> rhoF( 1 );
  Vector<T,3> velocityV;
  AnalyticalConst3D<T,T> uF(velocityV);

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, rhoF, uF );
  sLattice.iniEquilibrium( bulkIndicator, rhoF, uF );

  sLattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(SuperLattice<T, DESCRIPTOR>& sLattice,
                       const UnitConverter<T,DESCRIPTOR>& converter,
                       std::size_t iT,
                       SuperGeometry<T,3>& sGeometry,
                       size_t iTmaxStart,
                       T maxLatticeU,
                       T maxPhysU,
                       T iniPhysU)
{
  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  std::size_t iTupdate = 30;

  if (iT % iTupdate == 0 && iT <= iTmaxStart) {
    PolynomialStartScale<T,std::size_t> StartScale( iTmaxStart, T( maxLatticeU - iniPhysU ) );

    // Creates and sets the Poiseuille inflow profile using functors
    size_t iTvec[1] = { iT };
    T ux_i[1] = {};
    StartScale( ux_i, iTvec );
    AnalyticalConst3D<T,T> ux( ux_i[0] + iniPhysU );
    AnalyticalConst3D<T,T> uy( 0. );
    AnalyticalConst3D<T,T> uz( 0. );
    AnalyticalComposed3D<T,T> u( ux, uy, uz );
    sLattice.defineU( sGeometry, 3, u );

    clout << "startup step=" << iT << "/" << iTmaxStart << "; t=" << converter.getPhysTime(iT)
          << "; ux_LU=" << ux_i[0] << "; maxLatticeU=" << maxLatticeU
          << "; ux_PU=" << converter.getPhysVelocity( ux_i[0] ) << "; ux_max_PU=" << maxPhysU << std::endl;

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

void writeResults(SuperLattice<T, DESCRIPTOR>& sLattice,
                  const UnitConverter<T,DESCRIPTOR>& converter,
                  std::size_t iT,
                  SuperGeometry<T,3>& sGeometry,
                  std::string name)
{
  OstreamManager clout(std::cout, name);

  if (iT == 0) {
    SuperVTMwriter3D<T> vtmWriter(name, 0);
    if (name == "level0") {
      SuperLatticeRank3D rank( sLattice );
      vtmWriter.write(rank);
    }
    SuperGeometryF<T,DESCRIPTOR::d> geometryF(sGeometry);
    geometryF.getName() = name + "_geometry";
    vtmWriter.write(geometryF);
    vtmWriter.createMasterFile();
  }

  sLattice.setProcessingContext(ProcessingContext::Evaluation);
  sLattice.scheduleBackgroundOutputVTK([&,name,iT](auto task) {
    SuperVTMwriter3D<T> vtmWriter(name);
    SuperLatticePhysVelocity3D velocityF(sLattice, converter);
    SuperLatticePhysPressure3D pressureF(sLattice, converter);
    SuperLatticeRefinementMetricKnudsen3D qualityF(sLattice, converter);
    vtmWriter.addFunctor(qualityF);
    vtmWriter.addFunctor(velocityF);
    vtmWriter.addFunctor(pressureF);
    task(vtmWriter, iT);
  });

  SuperLatticePhysVelocity3D velocityF(sLattice, converter);
  SuperEuklidNorm3D<T> normVel( velocityF );
  BlockReduction3D2D<T> planeReductionU( normVel, Vector<T,3>({0, 0, 1}) );
  // write output as JPEG
  heatmap::plotParam<T> uScale;
  uScale.minValue = 0;
  uScale.name = "velocity_"+name;
  heatmap::write( planeReductionU, iT, uScale) ;
  
  // SuperLatticePhysPressure3D pressureF(sLattice, converter);
  // BlockReduction3D2D<T> planeReductionP( pressureF, Vector<T,3>({0, 0, 1}) );
  // heatmap::plotParam<T> pScale;
  // pScale.name = "pressure_"+name;
  // heatmap::write( planeReductionP, iT, pScale );
  
  // SuperLatticeYplus3D<T, DESCRIPTOR> yPlusF( sLattice, converter, sGeometry, foilBody, 5 );
  // SuperLatticeRefinementMetricKnudsen3D qualityF(sLattice, converter);
  // SuperRoundingF3D<T, T> roundedQuality ( qualityF, RoundingMode::NearestInteger );
  // SuperDiscretizationF3D<T> discretizationF ( roundedQuality, 0., 2. );
  // BlockReduction3D2D<T> planeReductionQ( discretizationF, Vector<T,3>({0, 0, 1}) );
  // heatmap::plotParam<T> jpeg_scale;
  // jpeg_scale.colour = "blackbody";
  // jpeg_scale.name = "quality_"+name;
  // heatmap::write( planeReductionQ, iT, jpeg_scale );
}

int main(int argc, char* argv[])
{
  olbInit(&argc, &argv);
  OstreamManager clout(std::cout, "main");

  CLIreader args(argc, argv);
  const int N               = args.getValueOrFallback<int>( "--res",        31);
  const int Re              = args.getValueOrFallback<int>( "--Re",         100);
  const int maxPhysT        = args.getValueOrFallback<int>( "--maxPhysT",   4);
  const T tMaxInit          = args.getValueOrFallback<T>  ( "--iniPhysT",   2);
  size_t iTmax              = args.getValueOrFallback<int>( "--iTmax",      0);
  T maxPhysU                = args.getValueOrFallback<T>  ( "--maxPhysU",   .2);    // in m/s
  T iniPhysU                = args.getValueOrFallback<T>  ( "--iniPhysU",   .00);
  T lengthDomain            = args.getValueOrFallback<T>  ( "--lx",         2.5);     // in m
  const T heightDomain      = args.getValueOrFallback<T>  ( "--ly",         0.5);     // in m
  const T depthDomain       = args.getValueOrFallback<T>  ( "--lz",         0.2);   // in m

  const UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> converterL0(
    int {N},             // resolution: number of voxels per charPhysL
    (T)   0.505,         // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   0.2,           // charPhysLength: reference length of simulation geometry
    (T)   maxPhysU,      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   0.2*0.2/Re,    // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0            // physDensity: physical density in __kg / m^3__
  );
  T maxLatticeU                 = converterL0.getLatticeVelocity( maxPhysU );
  converterL0.print();
  converterL0.write("du93");

  Vector<T,3> extendDomain = {lengthDomain, heightDomain, depthDomain};
  Vector<T,3> originDomain = {-lengthDomain/3, -heightDomain/2, -depthDomain/2};
  IndicatorCuboid3D<T> domainI(extendDomain,originDomain);
  IndicatorLayer3D<T> domainLayerI(domainI, converterL0.getPhysDeltaX());
  CuboidGeometry3D<T> cGeometryL0(domainLayerI, converterL0.getPhysDeltaX());
  // split 0th cuboid along dimension 0, creating cuboids 0,1,2 with ratios .15, .5, .35
  std::vector<T> xSplit = {0.3/2.5,1.2/2.5,1/2.5};
  cGeometryL0.splitFractional(0, 0, xSplit);
  // split middle cuboid (1) along dimension 1, inserting 2 back at 1, splitting 1 into 2,3,4
  std::vector<T> ySplit = {0.3,0.4,0.3};
  cGeometryL0.splitFractional(1, 1, ySplit);
  // // split centre cuboid (3) along dimension 0, inserting 4 back at 3, splitting 3 into 4,5
  // cGeometryL0.splitFractional(3, 0, {0.5,0.5});
  // // split centre cuboids (3,5) along dimension 1
  // cGeometryL0.splitFractional(4, 1, {0.5,0.5});  // 5 back at 4, 4 into 5,6
  // cGeometryL0.splitFractional(4, 1, {0.5,0.5});  // 5,6 back at 4,5, 4 into 6,7
  cGeometryL0.print();

  auto cGeometryL1 = cGeometryL0;
  Vector<T,3> originRefinement = {xSplit[0]*lengthDomain, ySplit[0]*heightDomain, 0};
  originRefinement += originDomain;
  Vector<T,3> extendRefinement = {xSplit[1]*lengthDomain*T(.99), ySplit[1]*heightDomain*T(.99), depthDomain};
  clout << "originRefinement=" << originRefinement << "; extendRefinement=" << extendRefinement << std::endl;
  IndicatorCuboid3D<T> toBeRefinedI( extendRefinement,
                                     originRefinement );
  cGeometryL1.remove(toBeRefinedI);
  cGeometryL1.refine(2);
  cGeometryL1.print();

  // split 0th cuboid along dimension 0, creating cuboids 0,1,2 with ratios .15, .5, .35
  xSplit = {0.15,0.5,0.35};
  cGeometryL1.splitFractional(0, 0, xSplit);
  // split middle cuboid (1) along dimension 1, inserting 2 back at 1, splitting 1 into 2,3,4
  std::vector<T> ySplit = {0.2,0.6,0.2};
  cGeometryL1.splitFractional(1, 1, ySplit);
  cGeometryL1.print();

  auto cGeometryL2 = cGeometryL1;
  Vector<T,3> originRefinement2 = {xSplit[0]*extendRefinement[0], ySplit[0]*extendRefinement[1], extendRefinement[2]};
  originRefinement += originDomain;
  Vector<T,3> extendRefinement2 = {xSplit[1]*extendRefinement[0]*T(.99), ySplit[1]*extendRefinement[1]*T(.99), extendRefinement[2]};
  clout << "originRefinement=" << originRefinement << "; extendRefinement=" << extendRefinement << std::endl;
  IndicatorCuboid3D<T> toBeRefinedI( extendRefinement,
                                     originRefinement );
  cGeometryL1.remove(toBeRefinedI);
  cGeometryL1.refine(2);
  cGeometryL1.print();

  // Adjust weights for balancing
  for (int iC=0; iC < cGeometryL0.size(); ++iC) {
    auto& cuboid = cGeometryL0.get(iC);
    cuboid.setWeight(cuboid.getLatticeVolume());
    auto origin = cGeometryL0.get(iC).getOrigin();
    for (int jC=0; jC < cGeometryL1.size(); ++jC) {
      for (int kC=0; kC < cGeometryL2.size(); ++kC) {

      }
      if (cGeometryL1.get(jC).getOrigin() == origin) {
        cuboid.setWeight(cuboid.getWeight() + 2*cGeometryL1.get(jC).getLatticeVolume());
      }
    }
  }

  HeuristicLoadBalancer<T> loadBalancerLevel0(cGeometryL0);
  SuperGeometry<T,3> sGeometryLevel0(cGeometryL0, loadBalancerLevel0);
  prepareGeometry(converterL0, sGeometryLevel0, extendDomain, originDomain);
  SuperLattice<T,DESCRIPTOR> sLatticeLevel0(cGeometryL0,
                                            loadBalancerLevel0,
                                            3,
                                            converterL0);
  prepareLattice(sLatticeLevel0, converterL0, sGeometryLevel0);

  RefinedLoadBalancer<T,3> loadBalancerLevel1(cGeometryL0,
                                              loadBalancerLevel0,
                                              cGeometryL1);
  auto converterL1 = convectivelyRefineUnitConverter(converterL0, 2);
  converterL1.print();
  SuperGeometry<T,3> sGeometryLevel1(cGeometryL1, loadBalancerLevel1);
  STLreader<T> foilBody( "DU93W210TET05_2deg_body_small.stl", converterL0.getPhysDeltaX() );
  STLreader<T> foilTail( "DU93W210TET05_2deg_tail_small.stl", converterL0.getPhysDeltaX() );
  prepareGeometryFine(converterL1, sGeometryLevel1, foilBody, foilTail);
  SuperLattice<T,DESCRIPTOR> sLatticeLevel1(cGeometryL1,
                                            loadBalancerLevel1,
                                            3,
                                            converterL1);
  prepareLattice(sLatticeLevel1, converterL1, sGeometryLevel1);

  clout << "Coupling coarse to fine grid(s) ..." << std::endl;
  auto coarseToFine = refinement::lagrava::makeCoarseToFineCoupler(
    sLatticeLevel0, sGeometryLevel0,
    sLatticeLevel1, sGeometryLevel1);
  auto fineToCoarse = refinement::lagrava::makeFineToCoarseCoupler(
    sLatticeLevel0, sGeometryLevel0,
    sLatticeLevel1, sGeometryLevel1);
  clout << "Coupling coarse to fine grid(s) ... OK" << std::endl;

  // iTmax depends on maximum physical time. If iTmax is provided in command line, it is an upper bound
  if ( iTmax == 0 ) iTmax = converterL0.getLatticeTime( maxPhysT );
  else iTmax = std::min( iTmax, converterL0.getLatticeTime( maxPhysT ) );
  size_t iTmaxStart = converterL0.getLatticeTime( tMaxInit );

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer(converterL0.getLatticeTime(maxPhysT),
                           sGeometryLevel0.getStatistics().getNvoxel()
                       + 2*sGeometryLevel1.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT < converterL0.getLatticeTime(maxPhysT); ++iT) {
    if (iT % converterL0.getLatticeTime(1.0) == 0) {
      writeResults(sLatticeLevel0, converterL0, iT, sGeometryLevel0, "level0");
      writeResults(sLatticeLevel1, converterL1, iT, sGeometryLevel1, "level1");
    }

    setBoundaryValues( sLatticeLevel0, converterL0, iT, sGeometryLevel0, iTmaxStart, maxLatticeU, maxPhysU, iniPhysU );
    sLatticeLevel0.collideAndStream();

    sLatticeLevel1.collideAndStream();
    coarseToFine->apply(meta::id<refinement::lagrava::HalfTimeCoarseToFineO>{});

    sLatticeLevel1.collideAndStream();
    coarseToFine->apply(meta::id<refinement::lagrava::FullTimeCoarseToFineO>{});

    fineToCoarse->apply(meta::id<refinement::lagrava::FineToCoarseO>{});

    if (iT % converterL0.getLatticeTime(0.1) == 0) {
      timer.update(iT);
      timer.printStep();
      sLatticeLevel0.getStatistics().print(iT, converterL0.getPhysTime(iT));
    }
  }

  singleton::pool().wait();

  timer.stop();
  timer.printSummary();
}
