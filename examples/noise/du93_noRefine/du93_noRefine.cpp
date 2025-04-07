/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Philipp Spelten
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

/* du93_3d.cpp:
 * TODO: change description
 */


#include "olb3D.h"
#include "olb3D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = float;
using DESCRIPTOR = D3Q19<>;

T* uAverage = NULL;

typedef enum {localpressure, localvelocity, interpressure, intervelocity, zerogradient, interconvection} OutletType;

#define BOUZIDI

// Stores data from stl file in geometry in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter, //std::shared_ptr<IndicatorF3D<T>> domain,
                      STLreader<T>& foilBody,
                      STLreader<T>& foilTail,
                      SuperGeometry<T,3>& sGeometry,
                      Vector<T,3> extendDomain,
                      Vector<T,3> originDomain,
                      bool withDampingLayer,
                      T boundaryDepth )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl
        << "Materials:" << std::endl
        << "1 = Fluid if no damping layer, else damping layer" << std::endl
        << "2 = Check (should be empty!)" << std::endl
        << "3 = Inflow" << std::endl
        << "4 = Outflow" << std::endl
        << "5 = Airfoil Bounce Back or Bouzidi" << std::endl
        << "6 = Airfoil porous trailing edge" << std::endl
        << "7 = Fluid if 1 is damping layer" << std::endl;
  
  IndicatorCuboid3D<T> cuboidDomain(extendDomain, originDomain);
  IndicatorLayer3D<T> domainLayer(cuboidDomain, converter.getPhysDeltaX());
  Vector<T,3> origin = domainLayer.getMin();
  Vector<T,3> extend = domainLayer.getMax() - domainLayer.getMin();

  sGeometry.rename( 0, 2 );

  sGeometry.rename( 2, 1, {1, 0, 0} );
  // Set material number for inflow
  origin[0] = domainLayer.getMin()[0]-converter.getPhysDeltaX()/2;
  extend[0] = converter.getPhysDeltaX();
  IndicatorCuboid3D<T> inflow( extend, origin );
  sGeometry.rename( 2, 3, 1, inflow );
  // Set material number for outflow
  origin[0] = domainLayer.getMax()[0]-converter.getPhysDeltaX()/2;
  extend[0] = converter.getPhysDeltaX();
  IndicatorCuboid3D<T> outflow( extend,origin );
  sGeometry.rename( 2, 4, outflow );
  sGeometry.rename( 1, 5, foilBody );
  sGeometry.rename( 1, 6, foilTail );

  // clout << "cleaning geometry" << std::endl;
  // sGeometry.clean();
  // sGeometry.getStatistics().print();

  if ( withDampingLayer ) {
    // fluid domain part to fluid
    T bd_pu = converter.getPhysLength( boundaryDepth );
    origin = domainLayer.getMin();
    extend = domainLayer.getMax() - domainLayer.getMin();
    origin[0] += bd_pu;
    extend[0] -= 2*bd_pu;
    IndicatorCuboid3D<T> fluid_domain( extend, origin );
    // clout << "1 to 7 in fluid_domain" << std::endl;
    sGeometry.rename( 1, 7, fluid_domain );  // all nodes except boundary length
    // sGeometry.getStatistics().print();
  }

  // Removes all not needed boundary voxels outside the surface
  // sGeometry.clean();
  sGeometry.checkForErrors();
  sGeometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void setUAverage( SuperLattice<T,DESCRIPTOR>& sLattice,
                  SuperGeometry<T,3>& sGeometry,
                  int fluidNumber
                  )
{
  sLattice.setProcessingContext(ProcessingContext::Evaluation);  // important for synchronization (?) on GPU

  // update outflow boundary value (adaptive convection boundary for smaller domains)
  SuperLatticeVelocity3D<T, DESCRIPTOR> velocity( sLattice );
  std::unique_ptr<SuperIndicatorF<T,3>> fluidIndicator = sGeometry.getMaterialIndicator({fluidNumber});
  SuperSum3D<T> sum( velocity, fluidIndicator );
  int input[1];
  T output[3];
  sum(output, input);
  *uAverage = output[0] / sGeometry.getStatistics().getNvoxel( 1 );

  sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
}

// Set up the geometry of the simulation
void prepareLattice( SuperLattice<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     STLreader<T>& foilBody,
                     STLreader<T>& foilTail,
                     SuperGeometry<T,3>& sGeometry,
                     bool porousTE,
                     T Kin,
                     T initialUx,
                     bool withDampingLayer,
                     T boundaryDepth,
                     T dampingStrength,
                     OutletType outlet,
                     Vector<T,3> extendDomain,
                     Vector<T,3> originDomain,
                     std::string name )
{
  OstreamManager clout( std::cout,"prepareLattice_"+name );
  clout << "Prepare Lattice ..." << std::endl;

  IndicatorCuboid3D<T> cuboidDomain(extendDomain, originDomain);
  IndicatorLayer3D<T> domainLayer(cuboidDomain, converter.getPhysDeltaX());

  T tmp = T();
  uAverage = &tmp;  //*uAverage = 0.0;
  int fluidNumber = 1;
  if ( withDampingLayer ) fluidNumber = 7;
  setUAverage( sLattice, sGeometry, fluidNumber );
  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  clout << "Setting bulk dynamics in fluid and PML in damping layer" << std::endl;
  auto bulkIndicator = sGeometry.getMaterialIndicator({1});
  if ( withDampingLayer ) {
    boundary::set<boundary::SpongeLayer>( sLattice, sGeometry, 1 );
    Vector<T,3> extend = domainLayer.getMax() - domainLayer.getMin();
    DampingTerm<3,T,DESCRIPTOR> sigma( converter, boundaryDepth, extend, dampingStrength );
    sLattice.defineField<descriptors::DAMPING>( bulkIndicator, sigma );

    bulkIndicator = sGeometry.getMaterialIndicator({7});
    sLattice.defineDynamics<BGKdynamics>(bulkIndicator);

    bulkIndicator = sGeometry.getMaterialIndicator({1,7});
  } else {
    sLattice.defineDynamics<BGKdynamics>(bulkIndicator);
  }

  // Material=3 -->inlet
  clout << "Setting LocalVelocity inlet on 3" << std::endl;
  if ( sGeometry.getStatistics().getNvoxel(3) > 0 ) boundary::set<boundary::LocalVelocity>(sLattice, sGeometry, 3);

  // Material=4 -->outlet
  if ( sGeometry.getStatistics().getNvoxel(4) > 0 ) {
    clout << "Setting outlet on 4" << std::endl;
    switch ( outlet ) {
      case localpressure:   boundary::set<boundary::LocalPressure>(sLattice, sGeometry, 4); break;
      case localvelocity:   boundary::set<boundary::LocalVelocity>(sLattice, sGeometry, 4); break;
      case interpressure:   boundary::set<boundary::InterpolatedPressure>(sLattice, sGeometry, 4); break;
      case intervelocity:   boundary::set<boundary::InterpolatedVelocity>(sLattice, sGeometry, 4); break;
      case zerogradient:    setZeroGradientBoundary(sLattice, sGeometry, 4); break;
      case interconvection: boundary::set<boundary::InterpolatedConvection>(sLattice, sGeometry, 4); break;
    }
  }

  // Material=5 -->bouzidi / bounce back
  clout << "Setting solid boundary on 5" << std::endl;
  #ifdef BOUZIDI
  if ( sGeometry.getStatistics().getNvoxel(5) > 0 ) setBouzidiBoundary<T,DESCRIPTOR>(sLattice, sGeometry, 5, foilBody);
  #else
  if ( sGeometry.getStatistics().getNvoxel(5) > 0 ) setBounceBackBoundary(sLattice, sGeometry, 5);
  #endif

  if ( porousTE ) {
    clout << "Setting porous dynamics on 6" << std::endl;
    // Material=6 --> porous media
    sLattice.defineDynamics<PorousBGKdynamics>(sGeometry, 6);  // velocity will always be multiplied by d
    T tau = converter.getLatticeRelaxationTime();
    T nu = (tau-0.5)/3.;
    T h = converter.getPhysDeltaX();
    T d = 1. - (h*h*nu*tau/Kin);
    clout << "Lattice Porosity: " << d << "(1 = permeable , 0 = not permeable)" << std::endl;
    clout << "Kmin: " << h*h*nu*tau << std::endl;
    if (Kin < h*h*nu*tau) {
      clout << "WARNING: Chosen K is too small!" << std::endl;
      Kin = h*h*nu*tau;  // exit(1);
    }
    AnalyticalConst3D<T,T> porosity(d);
    sLattice.defineField<POROSITY>(sGeometry.getMaterialIndicator({6}), porosity);
  } else {
    clout << "Setting solid boundary on 6" << std::endl;
    #ifdef BOUZIDI
    setBouzidiBoundary<T,DESCRIPTOR>(sLattice, sGeometry, 6, foilTail);
    #else
    setBounceBackBoundary(sLattice, sGeometry, 5);
    #endif
  }

  // Initial conditions
  AnalyticalConst3D<T,T> rhoF( 1 );
  Vector<T,3> velocityV;
  velocityV[0] = initialUx;
  AnalyticalConst3D<T,T> uF(velocityV);
  AnalyticalConst3D<T,T> ux(velocityV[0]);
  AnalyticalConst3D<T,T> uy(velocityV[1]);
  AnalyticalConst3D<T,T> uz(velocityV[2]);

  if ( withDampingLayer ) {
    sLattice.defineField<descriptors::UX>( bulkIndicator, ux );
    sLattice.defineField<descriptors::UY>( bulkIndicator, uy );
    sLattice.defineField<descriptors::UZ>( bulkIndicator, uz );
    sLattice.defineField<descriptors::DENSITY>( bulkIndicator, rhoF );
  }

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, rhoF, uF );
  sLattice.iniEquilibrium( bulkIndicator, rhoF, uF );

  sLattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( SuperLattice<T, DESCRIPTOR>& sLattice,
                        const UnitConverter<T,DESCRIPTOR>& converter,
                        size_t iT,
                        SuperGeometry<T,3>& sGeometry,
                        size_t iTmaxStart,
                        T maxLatticeU,
                        T maxPhysU,
                        T iniPhysU,
                        bool withDampingLayer,
                        OutletType outlet,
                        std::string name )
{
  OstreamManager clout( std::cout,"setBoundaryValues_"+name );

  // No of time steps for smooth start-up#
  int iTupdate = int(iTmaxStart/50/10)*10;  // --> about 50 updates, rounded to 10 iterations

  if ( iT%iTupdate == 0 && iT <= iTmaxStart ) {
    // Smooth start curve, sinus
    // SinusStartScale<T,int> StartScale(iTmaxStart, T( maxLatticeU ));

    // Smooth start curve, polynomial
    PolynomialStartScale<T, size_t> StartScale( iTmaxStart, T( maxLatticeU - iniPhysU ) );

    // Creates and sets the Poiseuille inflow profile using functors
    size_t iTvec[1] = { iT };
    T ux_i[1] = {};
    StartScale( ux_i, iTvec );
    AnalyticalConst3D<T,T> ux( ux_i[0] + iniPhysU );
    AnalyticalConst3D<T,T> uy( 0. );
    AnalyticalConst3D<T,T> uz( 0. );
    AnalyticalComposed3D<T,T> u( ux, uy, uz );
    sLattice.defineU( sGeometry, 3, u );

    if ( withDampingLayer ) {
      auto bulkIndicator = sGeometry.getMaterialIndicator({1,7});
      sLattice.defineField<descriptors::UX>( bulkIndicator, ux );
    }

    clout << "startup step=" << iT << "/" << iTmaxStart << "; t=" << converter.getPhysTime(iT)
          << "; ux_LU=" << ux_i[0] << "; maxLatticeU=" << maxLatticeU
          << "; ux_PU=" << converter.getPhysVelocity( ux_i[0] ) << "; ux_max_PU=" << maxPhysU << std::endl;
  }

  // AnalyticalConst3D<T,T> rho_out( 1. );
  // AnalyticalConst3D<T,T> ux_out( maxLatticeU );
  // AnalyticalConst3D<T,T> uy_out( 0. );
  // AnalyticalConst3D<T,T> uz_out( 0. );
  // AnalyticalComposed3D<T,T> u_out( ux_out, uy_out, uz_out );
  // switch ( outlet )
  // {
  // case localpressure:
  // case interpressure:
  //   sLattice.defineRho( sGeometry, 4, rho_out );
  //   break;
  // case localvelocity:
  // case intervelocity:
  //   sLattice.defineU( sGeometry, 4, u_out );
  //   break;
  // case zerogradient:
  // case interconvection:
  //   break;
  // }

  sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
}

// Computes the pressure drop between the voxels before and after the airfoil
void getResults( SuperLattice<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry<T,3>& sGeometry, util::Timer<T>& timer,
                 STLreader<T>& foilBody,
                 size_t iTmax, size_t iTmaxStart, size_t iTvtk,
                 bool withDampingLayer,
                 std::string name )
{
  OstreamManager clout( std::cout, "getResults_"+name+", iT=" + std::to_string(iT) );

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

  const size_t statIter = int( iTvtk/10 );
  const size_t iTcheck  = 500;
  bool lastIteration = false;

  if ( iT%iTcheck == 0 ) {
    if ( sLattice.getStatistics().getAverageRho() > 2. ) {
      lastIteration = true;
    }
  }

  // Writes output on the console
  if ( iT%statIter == 0 ) {
    SuperLatticeYplus3D<T, DESCRIPTOR> yPlus( sLattice, converter, sGeometry, foilBody, 5 );
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
    
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT, converter.getPhysTime( iT ) );

    // Drag, lift, pressure drop
    AnalyticalFfromSuperF3D<T> intpolatePressure( pressure, true );
    SuperLatticePhysDrag3D<T,DESCRIPTOR> drag( sLattice, sGeometry.getMaterialIndicator({5,6}), converter );

    olb::Vector<T, 3> point1V = sGeometry.getStatistics().getCenterPhysR( 5 );
    olb::Vector<T, 3> point2V = sGeometry.getStatistics().getCenterPhysR( 5 );
    T point1[3] = {};
    T point2[3] = {};
    for ( int i = 0; i<3; i++ ) {
      point1[i] = point1V[i];
      point2[i] = point2V[i];
    }
    point1[0] = std::min(sGeometry.getStatistics().getMinPhysR( 5 )[0], sGeometry.getStatistics().getMinPhysR( 6 )[0]) - converter.getPhysDeltaX();
    point2[0] = std::max(sGeometry.getStatistics().getMaxPhysR( 5 )[0], sGeometry.getStatistics().getMaxPhysR( 6 )[0]) + converter.getPhysDeltaX();

    T p1, p2;
    intpolatePressure( &p1,point1 );
    intpolatePressure( &p2,point2 );

    clout << name << ": pressure1=" << p1;
    clout << "; pressure2=" << p2;

    T pressureDrop = p1-p2;
    clout << "; pressureDrop=" << pressureDrop;

    T dragA[3];
    int input1[0];
    drag( dragA, input1 );
    clout << "; drag=" << dragA[0] << "; lift=" << dragA[1] << std::endl;

    int input[4] = {};
    SuperMax3D<T> yPlusMaxF( yPlus, sGeometry, 1 );
    T yPlusMax[1];
    yPlusMaxF( yPlusMax,input );
    clout << "yPlusMax=" << yPlusMax[0] << std::endl;

    if ( p1 != p1 ) lastIteration = true;

  }

  // Writes the vtk files
  if ( iT%iTvtk == 0 || lastIteration ) {
    sLattice.scheduleBackgroundOutputVTK([&,name,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriter(name);
      SuperLatticePhysVelocity3D velocityF(sLattice, converter);
      SuperLatticePhysPressure3D pressureF(sLattice, converter);
      vtmWriter.addFunctor(velocityF);
      vtmWriter.addFunctor(pressureF);
      task(vtmWriter, iT);
    });
    
    // SuperLatticePhysVelocity3D velocityF(sLattice, converter);
    // SuperEuklidNorm3D<T> normVel( velocityF );
    // BlockReduction3D2D<T> planeReductionU( normVel, Vector<T,3>({0, 0, 1}) );
    // // write output as JPEG
    // heatmap::plotParam<T> uScale;
    // uScale.minValue = 0;
    // uScale.name = "velocity_"+name;
    // heatmap::write( planeReductionU, iT, uScale) ;
    
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

  if ( iT%50 == 0 ) {
    int fluidNumber = 1;
    if ( withDampingLayer ) fluidNumber = 7;
    setUAverage( sLattice, sGeometry, fluidNumber );
    clout << "uAverage(" << iT << ")_" << name << "= " << *uAverage << std::endl;
  }

  if ( lastIteration ) {
    clout << "Stopping earyl after " << iT << " iterations due to too high average pressure. Sorry..." << std::endl;
    timer.stop();
    timer.printSummary();
    exit(1);
  }

  sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
}

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  // display messages from every single mpi process
  //clout.setMultiOutput(true);
  CLIreader args(argc, argv);
  std::string outdir        = args.getValueOrFallback<std::string>( "--outdir", "./tmp" );
  std::string foilName      = args.getValueOrFallback<std::string>( "--foilname", "DU93W210TET05");
  const int angle           = args.getValueOrFallback( "--angle",         2);   // in deg
  T lengthDomain            = args.getValueOrFallback( "--lx",            6);   // in m
  const T heightDomain      = args.getValueOrFallback( "--ly",            2);   // in m
  const T depthDomain       = args.getValueOrFallback( "--lz",            0.2); // in m
  const size_t res          = args.getValueOrFallback( "--res",           50);  // voxel/m (dx_LU/m)
  size_t boundaryDepth      = args.getValueOrFallback( "--bd",            20);  // depth of damping layer in LU
  T maxPhysT                = args.getValueOrFallback( "--tmax",          10);  // in s
  size_t iTmax              = args.getValueOrFallback( "--imax",          0);
  size_t nout               = args.getValueOrFallback( "--nout",          5);   // minimum number of vtk outputs
  size_t iout               = args.getValueOrFallback( "--iout",          0);   // iterations for vtk outputs
  T tout                    = args.getValueOrFallback( "--tout",          0);   // timestep for vtk outputs
  T maxPhysU                = args.getValueOrFallback( "--umax",          .01);   // in m/s
  T Re                      = args.getValueOrFallback( "--Re",            0);
  T tau                     = args.getValueOrFallback( "--tau",           0);   // previously tau=0.53 fixed
  const T Kin               = args.getValueOrFallback( "--permeability",  1e-8);
  const T dampingStrength   = args.getValueOrFallback( "--dampingStrength", .5);
  const T tMaxInit          = args.getValueOrFallback( "--tmaxinit",      2);
  T charL                   = args.getValueOrFallback( "--charL",         1);
  T iniPhysU                = args.getValueOrFallback( "--iniPhysU",      0);
  const int outlet          = args.getValueOrFallback( "--outlet",        3);  // default: interpolated pressure outlet
  const bool debug          = args.contains("--debug");
  const bool porousTE       = !args.contains("--no-porous");                   // --no-porous = no porous material in trailing edge
  const bool withDampingLayer   = !args.contains("--no-damping");              // --no-damping = no damping layer around domain

  const T cs_LU             = 1 / std::sqrt(3.0);
  T charV                   = 4 * maxPhysU;
  T viscosity;
  if ( tau == 0 ) {
    if ( Re == 0 ) {
      viscosity             = 1.48e-5; // 1e-2;
      Re                    = charV * charL / viscosity;
    } else viscosity        = charV * charL / Re;
    tau                     = viscosity / (cs_LU*cs_LU) + 0.5;
  }
  else {
    if ( Re == 0 ) Re       = 200;
    viscosity               = charV * charL / Re;
  }
  if ( debug ) {
    iTmax                   = 100;
  }

  std::stringstream outdir_mod;
  outdir_mod << outdir;
  if ( !porousTE ) outdir_mod << "_noPorous";
  else outdir_mod << "_" << Kin << "porous";
  outdir_mod << "_" << angle << "deg_u" << maxPhysU << "_Re" << Re << "_" << lengthDomain << "x" << heightDomain << "x" << depthDomain << "_res" << res;
  if ( withDampingLayer ) outdir_mod << "_bd" << boundaryDepth << "x" << dampingStrength;
  switch ( outlet ) {
    case 1: outdir_mod << "_outletLocalPressure";         break;
    case 2: outdir_mod << "_outletLocalVelocity";         break;
    case 3: outdir_mod << "_outletInterpolatedPressure";  break;
    case 4: outdir_mod << "_outletInterpolatedVelocity";  break;
    case 5: outdir_mod << "_outletZeroGradient";          break;
    case 6: outdir_mod << "_outletInterpolatedConvection";break;
  }

  singleton::directories().setOutputDir( outdir_mod.str()+"/" );  // set output directory
  std::ofstream fileStream( outdir_mod.str() + "/output.txt" );
  DoubleBuffer doubleBuffer( std::cout.rdbuf(), fileStream.rdbuf() );
  std::streambuf* originalCoutBuffer = std::cout.rdbuf( &doubleBuffer );
  OstreamManager clout( std::cout, "main" );
  clout << "outdir specified to " << outdir_mod.str() << std::endl;

  OutletType outlettype = intervelocity;
  switch ( outlet ) {
    case 1:   outlettype = localpressure;     clout << "Outlet boundary condition type specified to local pressure."            << std::endl; break;
    case 2:   outlettype = localvelocity;     clout << "Outlet boundary condition type specified to local velocity."            << std::endl; break;
    case 3:   outlettype = interpressure;     clout << "Outlet boundary condition type specified to interpolated pressure."     << std::endl; break;
    case 4:   outlettype = intervelocity;     clout << "Outlet boundary condition type specified to interpolated velocity."     << std::endl; break;
    case 5:   outlettype = zerogradient;      clout << "Outlet boundary condition type specified to zero gradient."             << std::endl; break;
    case 6:   outlettype = interconvection;   clout << "Outlet boundary condition type specified to interpolated convection."   << std::endl; break;
  }

  if ( porousTE ) clout << "Calculating with porous trailing edge" << std::endl;
  else clout << "Calculating without trailing edge" << std::endl;

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converterL0(
    (size_t)  res,            // resolution: number of voxels per charPhysL
    (T)       tau,            // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)       charL,          // charPhysLength: reference length of simulation geometry
    (T)       charV,          // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)       viscosity,      // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)       1.204           // physDensity: physical density in __kg / m^3__; air at 101.325 kPa and 20 °C
  );
  T maxLatticeU                 = converterL0.getLatticeVelocity( maxPhysU );
  // Prints the converter log as console output
  converterL0.print();
  // Writes the converter log in a file
  converterL0.write("du93_3d");

  // === 2nd Step: Prepare Geometry ===

  // Instantiation of the STLreader class
  // file name, voxel size in meter, stl unit in meter, outer voxel no., inner voxel no.
  std::string foilBodyFilename = foilName + "_" + std::to_string(angle) + "deg_body.stl";
  std::string foilTailFilename = foilName + "_" + std::to_string(angle) + "deg_tail.stl";
  clout << "Loading airfoil body from " << foilBodyFilename << " and airfoil tail from " << foilTailFilename << std::endl;
  STLreader<T> foilBody( foilBodyFilename, converterL0.getPhysDeltaX() );
  STLreader<T> foilTail( foilTailFilename, converterL0.getPhysDeltaX() );

  // Instantiation of a cuboidGeometry with weights
// #ifdef PARALLEL_MODE_MPI
//   const int noOfCuboids = singleton::mpi().getSize();
// #else
//   const int noOfCuboids = 7;
// #endif
  // setup domain
  if ( withDampingLayer ) lengthDomain += converterL0.getPhysLength( boundaryDepth );

  Vector<T,3> extendDomain = {lengthDomain, heightDomain, depthDomain};
  Vector<T,3> originDomain = {-lengthDomain/3, -heightDomain/2, -depthDomain/2};
  IndicatorCuboid3D<T> cuboidDomainL0(extendDomain, originDomain);
  IndicatorLayer3D<T> domainL0(cuboidDomainL0, converterL0.getPhysDeltaX());
  CuboidGeometry3D<T> cGeometryL0(domainL0, converterL0.getPhysDeltaX());
  cGeometryL0.setPeriodicity( {false, true, true} );
  cGeometryL0.print();

  // // === 3rd Step: Prepare Lattice ===
  HeuristicLoadBalancer<T> loadBalancerL0(cGeometryL0);
  SuperGeometry<T,3> sGeometryL0(cGeometryL0, loadBalancerL0);
  prepareGeometry( converterL0, foilBody, foilTail, sGeometryL0, extendDomain, originDomain, withDampingLayer, boundaryDepth );
  SuperLattice<T,DESCRIPTOR> sLatticeL0( cGeometryL0, loadBalancerL0, 3, converterL0);
  prepareLattice( sLatticeL0, converterL0, foilBody, foilTail, sGeometryL0, porousTE, Kin, iniPhysU, withDampingLayer, boundaryDepth, dampingStrength, outlettype, extendDomain, originDomain, "level0" );

  // === 3a-rd Step: calculate iterations from input ===
  // iTmax depends on maximum physical time. If iTmax is provided in command line, it is an upper bound
  if ( iTmax == 0 ) iTmax = converterL0.getLatticeTime( maxPhysT );
  else iTmax = std::min( iTmax, converterL0.getLatticeTime( maxPhysT ) );
  size_t iTmaxStart = converterL0.getLatticeTime( tMaxInit );
  // nout is the minimum number of vtk outputs --> take max between nout and nout derived from iout or tout
  size_t nout_from_iout = 0, nout_from_tout = 0;
  if ( iout != 0 ) { nout_from_iout = size_t( iTmax / iout ); nout = std::max( nout, nout_from_iout ); }
  if ( tout != 0 ) { nout_from_tout = size_t( iTmax / tout ); nout = std::max( nout, nout_from_tout ); }
  size_t iTvtk = size_t( iTmax / nout );
  clout << "Set nout to " << nout << ", so iTvtk=" << iTvtk << std::endl;

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timerL0( iTmax, sGeometryL0.getStatistics().getNvoxel() );
  timerL0.start();

  for (size_t iT = 0; iT < iTmax; ++iT) {
    // === Computation and Output of the Results ===
    getResults(sLatticeL0, converterL0, iT, sGeometryL0, timerL0, foilBody, iTmax, iTmaxStart, iTvtk, withDampingLayer, "level0");
    
    // === Definition of Initial and Boundary Conditions ===
    setBoundaryValues( sLatticeL0, converterL0, iT, sGeometryL0, iTmaxStart, maxLatticeU, maxPhysU, iniPhysU, withDampingLayer, outlettype, "level0" );
    sLatticeL0.collideAndStream();
    
    if (iT % converterL0.getLatticeTime(0.1) == 0) {
      timerL0.update(iT);
      timerL0.printStep();
      sLatticeL0.getStatistics().print(iT, converterL0.getPhysTime(iT));
    }

  }

  timerL0.stop();
  timerL0.printSummary();
  std::cout.rdbuf(originalCoutBuffer);
}
