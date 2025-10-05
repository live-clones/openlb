/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2017 Albert Mink, Christopher McHardy
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

/* cube3d.cpp:
 * A 3D implementation of radiation being emitted from a square or one side of a cube and
 * spreading in a bigger cuboid through media. The model offers 35 different cases of
 * optical parameters to describe the participating medium and two methods to
 * choose from.
 * The theoretical background and validation of said methods are detailed in
 * [A. Mink, C. McHardy, L. Bressel, C. Rauh and M. J. Krause. “Radiative transfer lattice Boltzmann
 * methods: 3D models and their performance in different regimes of radiative transfer”. In: Journal of
 * Quantitative Spectroscopy & Radiative Transfer, Volume 243 (2020). DOI: 10.1016/j.jqsrt.2019.106810.]
*/

#include "olb3D.h"
#include "olb3D.hh"   // include full template code

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;

typedef double T;

// select between different approaches: the mesoscopic method (default) or the macroscopic one based on the P1 approximation
//#define MINK  //macroscopic approach
// select between a directed or diffuse (default) boundary condition
#define DIRECTED  //directed boundary condition

using DESCRIPTOR = descriptors::D3Q27DescriptorLebedev;

struct INTENSITY : public descriptors::FIELD_BASE<1> { };

// PostProcessor to implement a directed boundary at the emittor
template<typename T, typename DESCRIPTOR, int discreteNormalX, int discreteNormalY, int discreteNormalZ>
struct RtlbmDirectedBoundaryPostProcessor3D
{
static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::list<INTENSITY>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform{
    T intensity = parameters.template get<INTENSITY>();
    for(int iPop = 0; iPop<DESCRIPTOR::q; iPop++){
      const auto c = descriptors::c<DESCRIPTOR>(iPop);
      T k = c[0]*discreteNormalX + c[1]*discreteNormalY + c[2]*discreteNormalZ;
      T norm_c = std::sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
      T norm_n = std::sqrt(discreteNormalX*discreteNormalX + discreteNormalY*discreteNormalY + discreteNormalZ*discreteNormalZ);
      T cos_theta = k / (norm_c * norm_n);
      if (util::nearZero(cos_theta + 1.)) {
        cell[iPop] = (1 - descriptors::t<T,DESCRIPTOR>(0))*intensity - descriptors::t<T,DESCRIPTOR>(iPop);
      } else {
        cell[iPop] = - descriptors::t<T,DESCRIPTOR>(iPop);
      }
    }
  }
};

// Set RtlbmDirectedBoundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR>
void setRtlbmDirectedBoundary(BlockLattice<T,DESCRIPTOR>& _block, BlockIndicatorF3D<T>& indicator)
{
  using namespace boundaryhelper;
  const auto& blockGeometryStructure = indicator.getBlockGeometry();
  std::vector<int> discreteNormal(4,0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (indicator(iX, iY, iZ)) {
      discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);
      if (discreteNormal[1]!=0 || discreteNormal[2]!=0 || discreteNormal[3]!=0) {
        _block.addPostProcessor(
          typeid(stage::PostCollide), {iX,iY,iZ},
          promisePostProcessorForNormal<T, DESCRIPTOR, RtlbmDirectedBoundaryPostProcessor3D>(
            Vector<int,3>(discreteNormal.data()+1)
          )
        );
        _block.template defineDynamics<NoCollideDynamics>({iX, iY, iZ});
      } else {
        _block.template defineDynamics<BounceBack>({iX, iY, iZ});
      }
    }
  });
}

// Initialising the setRtlbmDirectedBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setRtlbmDirectedBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setRtlbmDirectedBoundary");

  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); iCloc++) {
    setRtlbmDirectedBoundary<T,DESCRIPTOR>(sLattice.getBlock(iCloc), indicator->getBlockIndicatorF(iCloc));
  }
}

// Store all geometric information in material numbers
SuperGeometry<T,3> prepareGeometryCube( T conversionFactorLength, T boundaryShift )
{
  OstreamManager clout( std::cout, "prepareGeometryCube" );
  clout << "Prepare cubeGeometry ..." << std::endl;

  Vector<T,3> origin = {0+boundaryShift,-0.5,-0.5};
  Vector<T,3> extent = {1-boundaryShift,1,1};
  IndicatorCuboid3D<T> span(extent, origin);
  origin = {0+boundaryShift,0,0};

  // select between a finite or infinite source of light
  // finite light beam source
  IndicatorCuboid3D<T> emittor(2*conversionFactorLength, 0.2, 0.2, origin);
  // infinite light beam source
  //IndicatorCuboid3D<T> emittor(2*conversionFactorLength, 1.0-2*conversionFactorLength, 1.0-2*conversionFactorLength, origin);

  //CuboidDecomposition<T,3>* cuboid = new CuboidDecomposition<T,3>( span, conversionFactorLength, 1 );
  CuboidDecomposition<T,3>* cuboid = new CuboidDecomposition<T,3>( span, conversionFactorLength, 2*singleton::mpi().getSize() );
  HeuristicLoadBalancer<T>* loadBalancer = new HeuristicLoadBalancer<T>( *cuboid );

  SuperGeometry<T,3> superGeometry( *cuboid, *loadBalancer, 2 );
  // set material number, 0 outside, 1 domain, 2 outflow, 3 inflow
  superGeometry.rename( 0, 2, span );
  superGeometry.rename( 2, 1, {1, 1, 1} ); // last 3 digits form overlap
  superGeometry.rename( 2, 3, emittor );
  // clean voxels with material number 0, inside and outside geometry
  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();
  clout << "Prepare cubeGeometry ... OK" << std::endl;
  return superGeometry;
}

// Set up the geometry of the simulation in the macroscopic case
void prepareLatticeMink( SuperLattice<T,DESCRIPTOR>& sLattice,
                     RadiativeUnitConverter<T, DESCRIPTOR>const& converter,
                     SuperGeometry<T,3>& superGeometry,
                     T inletDirichlet)
{
  OstreamManager clout( std::cout, "prepareLatticeMink" );
  clout << "Prepare Lattice ..." << std::endl;

  T latticeRelaxationFrequency = converter.getLatticeRelaxationFrequency();
  T latticeSink = 3.*converter.getLatticeAbsorption()*(converter.getLatticeAbsorption()+converter.getLatticeScattering()) / 8.;
  T latticeAbsorption = converter.getLatticeAbsorption();
  T latticeScattering = converter.getLatticeScattering();

  sLattice.defineDynamics<NoDynamics<T, DESCRIPTOR>>( superGeometry, 0 );
  //sLattice.defineDynamics<PoissonDynamics<T, DESCRIPTOR>>( superGeometry, 1 );
  sLattice.defineDynamics<P1Dynamics<T, DESCRIPTOR>>( superGeometry, 1 );
  sLattice.defineDynamics<EquilibriumBoundaryFirstOrder>( superGeometry, 2 );
  sLattice.defineDynamics<EquilibriumBoundaryFirstOrder>( superGeometry, 3 );

  #ifdef DIRECTED
    setRtlbmDirectedBoundary<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator({3}));
  #endif

  sLattice.setParameter<descriptors::OMEGA>( latticeRelaxationFrequency );
  sLattice.setParameter<collision::P1::ABSORPTION>( latticeAbsorption );
  sLattice.setParameter<collision::P1::SCATTERING>( latticeScattering );
  sLattice.setParameter<collision::Poisson::SINK>( latticeSink );
  sLattice.setParameter<INTENSITY>( inletDirichlet );


  clout << "Prepare Lattice ... OK" << std::endl;
  return;
}

// Set up the geometry of the simulation in the mesoscopic case
void prepareLatticeMcHardy( SuperLattice<T,DESCRIPTOR>& sLattice,
                     RadiativeUnitConverter<T, DESCRIPTOR>const& converter,
                     SuperGeometry<T,3>& superGeometry,
                     T inletDirichlet,
                     T anisotropyFactor )
{
  OstreamManager clout( std::cout, "prepareLatticeMcHardy" );
  clout << "Prepare Lattice ..." << std::endl;

  T latticeRelaxationFrequency = converter.getLatticeRelaxationFrequency();

  int const q = descriptors::q<DESCRIPTOR>() - 1;
  std::array<std::array<T,q>, q> phi;
  T solution [q*(q +1)/2];
  computeAnisotropyMatrix<DESCRIPTOR>(1e-4, anisotropyFactor, solution, phi);

  T anisoMatrix[(q+1)*(q+1)] {};
  for ( int m = 0; m < q; m++ ) {
    for ( int n = 0; n < q; n++ ) {
      anisoMatrix[m+1+n+1] = phi[m][n];
    }
  }

  T latticeAbsorption = converter.getLatticeAbsorption();
  T latticeScattering = converter.getLatticeScattering();

  sLattice.defineDynamics<NoDynamics<T, DESCRIPTOR>>( superGeometry, 0 );
  // select between the mesoscopic method with or without a Runge Kutta scheme
  sLattice.defineDynamics<RTLBMdynamicsMcHardyRK<T, DESCRIPTOR>>( superGeometry, 1 );
  // sLattice.defineDynamics<RTLBMdynamicsMcHardy<T, DESCRIPTOR>>( superGeometry, 1 );
  sLattice.defineDynamics<EquilibriumBoundaryFirstOrder>( superGeometry, 2 );
  sLattice.defineDynamics<EquilibriumBoundaryFirstOrder>( superGeometry, 3 );

  #ifdef DIRECTED
    setRtlbmDirectedBoundary<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator({3}));
  #endif

  sLattice.setParameter<descriptors::OMEGA>( latticeRelaxationFrequency );
  sLattice.setParameter<Light::ANISOMATRIX>( anisoMatrix );
  sLattice.setParameter<Light::ABSORPTION>( latticeAbsorption );
  sLattice.setParameter<Light::SCATTERING>( latticeScattering );
  sLattice.setParameter<INTENSITY>( inletDirichlet );



  clout << "Prepare Lattice ... OK" << std::endl;
  return;
}


void setBoundaryValues( SuperLattice<T,DESCRIPTOR>& sLattice,
                        UnitConverter<T,DESCRIPTOR> const& converter,
                        SuperGeometry<T,3>& superGeometry,
                        T inletDirichlet )
{
  // since adv-diffusion model is used, the velocity it set to 0
  AnalyticalConst3D<T,T> u0( 0.0, 0.0, 0.0 );        // 3D -> 3D
  AnalyticalConst3D<T,T> rho0( 0.0 );                // 3D -> 1D
  AnalyticalConst3D<T,T> rho1( inletDirichlet );                // 3D -> 1D

  // initialize media with density from analytical solution
  // at iT=0 the error is given by the maschinen genauigkeit
  sLattice.iniEquilibrium( superGeometry, 1, rho0, u0 );
  sLattice.iniEquilibrium( superGeometry, 2, rho0, u0 );
  sLattice.iniEquilibrium( superGeometry, 3, rho1, u0 );
  sLattice.defineRho( superGeometry, 2, rho0 );
  sLattice.defineRho( superGeometry, 3, rho1 );
  sLattice.initialize();
}

int main( int argc, char *argv[] )
{
  // ===== 1st Step: Initialization =====
  initialize( &argc, &argv );
  std::string fName("cube3d.xml");
  XMLreader config(fName);

  int RESOLUTION;
  T LATTICERELAXATIONTIME, ANISOTROPYFACTOR;
  std::string NAME;
  config["Application"]["Name"].read(NAME);
  config["Application"]["Discretization"]["Resolution"].read<int>(RESOLUTION);
  config["Application"]["Discretization"]["LatticeRelaxationTime"].read(LATTICERELAXATIONTIME);
  T maxPhysT;
  config["Application"]["PhysParam"]["maxTime"].read(maxPhysT);
  config["Application"]["AnisotropyFactor"].read(ANISOTROPYFACTOR);

  std::string caseName{"case"};
  if ( argc < 2 ) {
    std::cerr << "Error: Please provide the case number (1 to 35) as a first argument and optionally the resolution as a second argument." << std::endl;
    return 1;
  }
  else if ( argc < 3 && argc > 1 ) {
    caseName.append(argv[1]);
  }
  else if ( argc < 4 && argc > 2 ) {
    caseName.append(argv[1]);
    RESOLUTION = std::atoi(argv[2]);
  }
  else if ( argc > 4 ) {
    std::cerr << "Error: Too many arguments provided." << std::endl;
    return 1;
  }
#ifdef MINK
  singleton::directories().setOutputDir( "./"+std::to_string(RESOLUTION)+caseName+"_mink/" );
#else
  singleton::directories().setOutputDir( "./"+std::to_string(RESOLUTION)+caseName+"_mcHardy/" );
#endif
  OstreamManager clout( std::cout, "main" );
  clout << caseName << std::endl;


  T ABSORPTION, SCATTERING, MCVALUE, TOTENERGY;
  std::stringstream xmlAbsorption(config["Application"][caseName].getAttribute( "absorption"));
  xmlAbsorption >> ABSORPTION;
  std::stringstream xmlScattering(config["Application"][caseName].getAttribute( "scattering"));
  xmlScattering >> SCATTERING;
  std::stringstream xmlMCBoundary(config["Application"][caseName].getAttribute( "mcvalue"));
  xmlMCBoundary >> MCVALUE;
  std::stringstream xmlTotalEnergy(config["Application"][caseName].getAttribute( "totalEnergy"));
  xmlMCBoundary >> TOTENERGY;

  //LATTICERELAXATIONTIME = 1./(RESOLUTION*DESCRIPTOR<T>::invCs2*(1./(3*(ABSORPTION+SCATTERING))));
  LATTICERELAXATIONTIME = 1;// 1./(RESOLUTION*DESCRIPTOR<T>::invCs2*(1./(3*(ABSORPTION+SCATTERING)))+0.5);
  clout << "omega = .... " << LATTICERELAXATIONTIME << std::endl;
  RadiativeUnitConverter<T,DESCRIPTOR> const converter( RESOLUTION, LATTICERELAXATIONTIME, ABSORPTION, SCATTERING, ANISOTROPYFACTOR );
  clout << "omega = " << converter.getLatticeRelaxationTime() << std::endl;
#ifdef MINK
  clout << "working with diffuse approximation" << std::endl;
  T latticeSink = converter.getLatticeAbsorption()/converter.getLatticeDiffusion()/8.;
  clout << "latticeSink=" << latticeSink << std::endl;
#else
  clout << "working with direct discretization" << std::endl;
#endif
  converter.print();
  converter.write();


  // ===== 2nd Step: Prepare Geometry =====
  SuperGeometry<T,3> superGeometry( prepareGeometryCube(converter.getConversionFactorLength(), 0. ) );
  //SuperGeometry<T,3> superGeometry( prepareGeometryCube(converter.getConversionFactorLength(), 1./converter.getExtinction() ) );


  // ===== 3rd Step: Prepare Lattice =====
  SuperLattice<T,DESCRIPTOR> sLattice(superGeometry);
  #ifdef MINK
  prepareLatticeMink( sLattice,
                  converter,
                  superGeometry,
                  1.0 );
  #else
    prepareLatticeMcHardy( sLattice,
                  converter,
                  superGeometry,
                  1.0,
                  ANISOTROPYFACTOR );
  #endif


  // ===== 4th Step: Main Loop with Timer =====
  util::Timer<double> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  // ===== 5th Step: Definition of Initial and Boundary Conditions =====
  setBoundaryValues(sLattice, converter, superGeometry, 1.0);
  // setBoundaryValues(sLattice, converter, superGeometry, MCVALUE);

  SuperVTMwriter3D<T> vtmWriter("cube3d");
  SuperGeometryF3D<T> geometry(superGeometry);
  vtmWriter.write( geometry );
  vtmWriter.createMasterFile();
  SuperLatticeDensity3D<T,DESCRIPTOR> density( sLattice );
  SuperLatticeFlux3D<T,DESCRIPTOR> flux( sLattice );
  vtmWriter.addFunctor(density);
  vtmWriter.addFunctor(flux);
  vtmWriter.addFunctor(geometry);

  util::ValueTracer<T> converge( 160, 1e-8);
  clout << "iT convergence criteria: " << 160 << std::endl;

  for ( int iT = 0; iT >= -1  /*converter.getLatticeTime( maxPhysT )*/; ++iT ) {

#ifdef MINK
    if ( iT % (2*RESOLUTION) == 0 ) {
#else
    if ( iT % (RESOLUTION/5) == 0 ) {
#endif
      sLattice.getStatistics().print( iT, converter.getPhysTime( iT ) );
      timer.print( iT );
      timer.printStep();
      vtmWriter.write(iT);
    }

    converge.takeValue( sLattice.getStatistics().getAverageRho(), true );
    // ===== 5.5th Step: Check for Convergence =====
    if ( converge.hasConverged() && iT > 4*RESOLUTION ) {
      clout << "Simulation converged. -- " << iT << std::endl;
      vtmWriter.write(iT);
      clout << "------" << iT << std::endl;

      // write and save results in a text file
      AnalyticalFfromSuperF3D<T> analytDen( density, true, 1 );
      std::string str = singleton::directories().getLogOutDir() + std::to_string(RESOLUTION) + caseName + ".csv";
      AnalyticalFfromSuperF3D<T> analytFlu( flux, true, 1 );
      clout << str << std::endl;
      FILE * pFile;
      const char * fileName = str.data();
      pFile  = fopen(fileName, "w");
      //fprintf(pFile, "%i\n", iT);
      //fprintf(pFile, "%s, %s, %s, %s\n", "position x", "0.0", "0.25", "0.375");
      for ( int nZ = 0; nZ <= 100; ++nZ ) {
        double position1[3] = {1.0*double(nZ)/100, 0, 0};
        double position2[3] = {1.0*double(nZ)/100, 0.25, 0};
        double position3[3] = {1.0*double(nZ)/100, 0.375, 0};
        double light1[1] = {0};
        double fluxx1[3] = {0.,0.,0.};
        double light2[1] = {0};
        double light3[1] = {0};
        analytDen( light1, position1 );
        analytFlu( fluxx1, position1 );
        analytDen( light2, position2 );
        analytDen( light3, position3 );
        if (singleton::mpi().getRank()==0) {
          printf("%4.3f, %8.6e, %8.6e, %8.6e, %8.6e, %8.6e, %8.6e\n", position1[0], light1[0], light2[0], light3[0],
                 fluxx1[0], fluxx1[1], fluxx1[2]);
          fprintf(pFile, "%4.3f, %8.6e, %8.6e, %8.6e, %8.6e, %8.6e, %8.6e\n", position1[0], light1[0], light2[0], light3[0],
                  fluxx1[0], fluxx1[1], fluxx1[2]);
        }
      }
      fclose(pFile);
      break;
    }
    // ===== 6th Step: Collide and Stream Execution =====
    sLattice.collideAndStream();
  }

  timer.stop();
  timer.printSummary();

  return 0;
}