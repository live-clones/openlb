/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023 Adrian Kummerlaender, Dennis Teutscher
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net
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

#include <olb3D.h>
#include <olb3D.hh>

//#define DYNAMIC_XML
using namespace olb;

using T = float;
using DESCRIPTOR = descriptors::D3Q19<descriptors::POROSITY>;
using BulkDynamics_ = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Porous<momenta::BulkTuple>,
  equilibria::SecondOrder,
  collision::SmagorinskyEffectiveOmega<collision::BGK>
>;

void prepareGeometry(const UnitConverter<T,DESCRIPTOR>& converter,
                     SuperGeometry<T,3>& sGeometry,
                     IndicatorF3D<T>& bounding)
{
  sGeometry.rename(0,2);
  sGeometry.rename(2,1,{1,1,1});
  // Set material number for floor
  {
    auto origin = bounding.getMin()-converter.getPhysDeltaX();
    auto extent = bounding.getMax() - bounding.getMin();
    extent[2] = 2*converter.getPhysDeltaX();
    extent[0] = bounding.getMax()[0]+0.1*bounding.getMax()[0];
    extent[1] = bounding.getMax()[1]+0.1*bounding.getMax()[0];
    IndicatorCuboid3D<T> floor(extent, origin);
    sGeometry.rename(2,4,floor);
    sGeometry.rename(1,4,floor);
  }
  //// Set material number for roof
  {
    auto origin = bounding.getMin()-converter.getPhysDeltaX();
    auto extent = bounding.getMax() - bounding.getMin();
    origin[2] = bounding.getMax()[2]-converter.getPhysDeltaX();
    extent[2] = 1*converter.getPhysDeltaX();
    extent[0] = bounding.getMax()[0]+0.1*bounding.getMax()[0];
    extent[1] = bounding.getMax()[1]+0.1*bounding.getMax()[0];
    IndicatorCuboid3D<T> roof(extent, origin);
    sGeometry.rename(2,5,roof);
  }
  // Set material number for right wall
  {
    auto origin = bounding.getMin()-converter.getPhysDeltaX();
    auto extent = bounding.getMax() - bounding.getMin();
    extent[0] = bounding.getMax()[0]+0.1*bounding.getMax()[0];
    extent[1] = converter.getPhysDeltaX();
    extent[2] = bounding.getMax()[2]+0.1*bounding.getMax()[0];
    IndicatorCuboid3D<T> right(extent, origin);
    sGeometry.rename(2,5,right);
  }
  // Set material number for left wall
  {
    auto origin = bounding.getMin()-converter.getPhysDeltaX();
    auto extent = bounding.getMax() - bounding.getMin();
    origin[1] = bounding.getMax()[1]-converter.getPhysDeltaX();
    extent[0] = bounding.getMax()[0]+0.1*bounding.getMax()[0];
    extent[1] = converter.getPhysDeltaX();
    extent[2] = bounding.getMax()[2]+0.1*bounding.getMax()[0];
    IndicatorCuboid3D<T> left(extent, origin);
    sGeometry.rename(2,5,left);
  }

  // Removes all not needed boundary voxels outside the surface
  sGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  sGeometry.innerClean();
  sGeometry.checkForErrors();
  sGeometry.print();
}

void prepareLattice(SuperLattice<T,DESCRIPTOR>& sLattice,
                    SuperGeometry<T,3>& sGeometry,
                    const UnitConverter<T,DESCRIPTOR>& converter, IndicatorF3D<T>& volume)
{
  const T omega = converter.getLatticeRelaxationFrequency();
  // Material=1,2 -->bulk dynamics
  #ifdef DYNAMIC_XML
  defineDynamicsAndSetParameterFromXml<T,DESCRIPTOR>("porousCity3d.xml",sGeometry,sLattice);
  #else
  sLattice.defineDynamics<NoDynamics>(sGeometry,0);
  sLattice.defineDynamics<BulkDynamics_>(sGeometry,1);
  sLattice.defineDynamics<BulkDynamics_>(sGeometry,2);
  sLattice.setParameter<collision::LES::SMAGORINSKY>(0.3);
  #endif
  // Material 4 -->bounceback(ground)
  boundary::set<boundary::BounceBack>(sLattice, sGeometry, 4);
  // Material 5 --> slip(sides and roof)
  boundary::set<boundary::FullSlip>(sLattice, sGeometry, 5);

  {
    boundary::set<boundary::InterpolatedVelocity>(sLattice, sGeometry, 2);
    AnalyticalConst3D<T,T> rhoF(1);
    AnalyticalConst3D<T,T> uF( T( 0 ), T( 0 ), T( 0 ) );
    sLattice.defineRhoU(sGeometry.getMaterialIndicator(2), rhoF, uF);
    sLattice.iniEquilibrium(sGeometry.getMaterialIndicator(2), rhoF, uF);
  }
  sLattice.setParameter<descriptors::OMEGA>(omega);
  {
    auto bulkIndicator = sGeometry.getMaterialIndicator({1,3,4,5});
    AnalyticalConst3D<T,T> rhoF( T( 1 ) );
    AnalyticalConst3D<T,T> uF( T( 0 ), T( 0 ), T( 0 ) );
    sLattice.defineRhoU(bulkIndicator, rhoF, uF);
    sLattice.iniEquilibrium(bulkIndicator, rhoF, uF);
  }

  //Defining porosity of fluid (1 = permeable )
  AnalyticalConst3D<T,T> initialPorosityF(1);
  sLattice.defineField<descriptors::POROSITY>(sGeometry.getMaterialIndicator({0,1,2,3,4,5}), initialPorosityF);

  //Defining porosity for geometry (0 not permeable)
  AnalyticalConst3D<T,T> solidPorosityF(0);
  SuperIndicatorFfromIndicatorF3D<T> discreteVolume(volume, sGeometry);
  sLattice.defineField<descriptors::POROSITY>(discreteVolume, solidPorosityF);
  // Make the lattice ready for simulation
  sLattice.initialize();

}

void getResults(SuperLattice<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT, util::Timer<T>* timer,
                 const T logT, const T maxPhysT, const T vtkSave,
                 std::string filenameVtk,
                 const int timerPrintMode,
                 SuperGeometry<T,3>& sGeometry)
{
  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter3D<T> vtmWriter( "porousCity3d" );

  //Writes geometry and porosity field
  if (iT == 0) {
    SuperLatticeCuboid3D<T,DESCRIPTOR> cuboidF(sLattice);
    SuperLatticeRank3D<T,DESCRIPTOR> rankF(sLattice);
    SuperLatticeField3D<T,DESCRIPTOR,descriptors::POROSITY> porosityF(sLattice);

    vtmWriter.write(cuboidF);
    vtmWriter.write(rankF);
    vtmWriter.write(porosityF);
    vtmWriter.createMasterFile();
  }
  //Get statistics
  if ( iT%converter.getLatticeTime( logT )==0) {
    sLattice.getStatistics().print( iT, converter.getPhysTime( iT ) );
    timer->print( iT,timerPrintMode );
  }

  // Writes the VTK
  if ( iT%converter.getLatticeTime( vtkSave )==0 && iT>0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    sLattice.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtkWriter("porousCity3d");
      SuperLatticePhysVelocity3D<T,DESCRIPTOR> velocity(sLattice, converter);
      vtkWriter.addFunctor(velocity);
      task(vtkWriter, iT);
    });
  }

  if ( iT ==converter.getLatticeTime(30)) {
    SuperLatticePhysVelocity3D<T,DESCRIPTOR> velocity(sLattice, converter);
    SuperEuklidNorm3D<T> normVel( velocity );
    auto line3d = Line3D<T>().originAt({-70,400,0}).parallelTo({0,0,1});
    LineLattice3D<T> lineLattice(sGeometry.getCuboidGeometry(),line3d);
    BlockReduction3D1D<T> line(normVel,lineLattice,BlockDataSyncMode::ReduceAndBcast,BlockDataReductionMode::Discrete);
    int n = line.getN();
    static Gnuplot<T> gplot( "inflow" );
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    for(int i =0; i < n; i++)
    {
      T res[3]{};
      line(res,i);
      auto physR= line.getPhysR(i);
      clout << res[0]<<" "<<res[1]<<" "<<res[2]<<std::endl;
      gplot.setData(res[0],T(physR[2]),"olb");

    }
    gplot.writePNG();
  }
}

void setBoundaryValues(const UnitConverter<T,DESCRIPTOR>& converter,
                       SuperGeometry<T,3>& sGeometry,
                       SuperLattice<T,DESCRIPTOR>& sLattice,
                       int iT,T z_ref, T z_0, T d, T kappa)
{
  OstreamManager clout(std::cout, "boundary");

  const auto maxStartT =  converter.getLatticeTime(10);
  const auto startIterT = converter.getLatticeTime(0.1);

  if (iT < maxStartT) {
    SinusStartScale<T,int> startScale(maxStartT, T(1));
    int iTvec[1]= {iT};
    T frac[1]= {};
    startScale( frac,iTvec );
    //AnalyticalWindProfileF3D<T> uF(converter.getCharLatticeVelocity()*frac[0],z_0,z_ref,kappa,d,0,2);

    AnalyticalWindProfilePowerLawF3D<T> uF(converter.getCharLatticeVelocity()*frac[0],z_ref,0.3f,0.05,0.3f);


    sLattice.defineU(sGeometry, 2, uF);
    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}


int main(int argc, char* argv[]) {
  // === 1st Step: Initialization ===
  olbInit(&argc, &argv);
  OstreamManager clout(std::cout, "main");
  std::string fName( "porousCity3d.xml" );
  XMLreader config( fName );

  std::string olbdir, outputdir;
  config["Application"]["OlbDir"].read( olbdir );
  config["Output"]["OutputDir"].read( outputdir );
  singleton::directories().setOlbDir( olbdir );
  singleton::directories().setOutputDir( outputdir );

  //Get converter information from XML
  UnitConverter<T,DESCRIPTOR>* converter = createUnitConverter<T,DESCRIPTOR>( config );

  CLIreader args(argc, argv);
  clout << "To utilize additional STL files, specify the desired file using the --volume flag." << std::endl;
  clout << "Example: ./porousCity3d --volume example.stl"<< std::endl;
  const std::string volumeFile = args.getValueOrFallback<std::string>("--volume", "kit_campus.stl");
  // Prints the converter log as console output
  converter->print();

  //Load specific parameters from XML
  T logT = config["Output"]["Log"]["SaveTime"].get<T>();
  T vtkSave = config["Output"]["VisualizationVTK"]["SaveTime"].get<T>();
  T maxPhysT = config["Application"]["PhysParameters"]["PhysMaxTime"].get<T>();
  int timerPrintMode = config["Output"]["Timer"]["PrintMode"].get<int>();
  std::string filenameVtk = config["Output"]["VisualizationVTK"]["Filename"].get<std::string>();

  //atmLayer paramters
  T z_ref = config["Application"]["PhysParameters"]["atmReferenceHeight"].get<T>();
  T z_0 = config["Application"]["PhysParameters"]["atmRoughnessLength"].get<T>();
  T d = config["Application"]["PhysParameters"]["atmGroundNormalDisplacement"].get<T>();
  T kappa = config["Application"]["PhysParameters"]["karmanConstant"].get<T>();

  // === 2nd Step: Prepare Geometry ===
  STLreader<T> volume(volumeFile, converter->getPhysDeltaX(), 1,olb::RayMode::Robust);
  auto origin =volume.getMin();
  origin[0] = -0.1*volume.getMax()[0];
  origin[1] = -0.1*volume.getMax()[0];
  origin[2] = -converter->getCharPhysLength();
  auto extent = volume.getMax() - origin;
  extent[0] += 0.1*volume.getMax()[0];
  extent[1] += 0.1*volume.getMax()[0];
  extent[2] *=2.;

  IndicatorCuboid3D<T> bounding(extent, origin);

  CuboidGeometry3D<T> cGeometry(bounding,converter->getPhysDeltaX(), singleton::mpi().getSize());

  BlockLoadBalancer<T> loadBalancer(cGeometry);
  SuperGeometry<T,3> sGeometry(cGeometry, loadBalancer);

  prepareGeometry(*converter, sGeometry, bounding);

  SuperLattice<T,DESCRIPTOR> sLattice(sGeometry);

  // === 3rd Step: Prepare Lattice ===
  prepareLattice(sLattice, sGeometry, *converter, volume);

  // === 4th Step: Main Loop with Timer ===
  util::Timer<T>* timer = util::createTimer<T>(config, *converter, sGeometry.getStatistics().getNvoxel());
  timer->start();

  for (std::size_t iT=0; iT <= converter->getLatticeTime(maxPhysT); ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(*converter,sGeometry,sLattice,iT, z_ref, z_0, d, kappa);
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, *converter, iT, timer, logT, maxPhysT, vtkSave, filenameVtk, timerPrintMode, sGeometry);
  }

  sLattice.setProcessingContext(ProcessingContext::Evaluation);

  timer->stop();
  timer->printSummary();
  delete converter;
  delete timer;
  return 0;
}
