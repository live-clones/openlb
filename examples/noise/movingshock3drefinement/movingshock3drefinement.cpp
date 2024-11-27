/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006, 2007, 2012 Jonas Latt, Mathias J. Krause,
 *  Louis Kronberg, Christian Vorwerk, Bastian Schäffauer
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

/* movingshock3d.cpp:
 * The implementation of a backward facing step. It is furthermore
 * shown how to use checkpointing to save the state of the
 * simulation regularly.
 * The geometry of the step is based on the experiment described in
 * [Armaly, B.F., Durst, F., Pereira, J. C. F. and Schönung, B. Experimental
 * and theoretical investigation of backward-facing step flow. 1983.
 * J. Fluid Mech., vol. 127, pp. 473-496, DOI: 10.1017/S0022112083002839]
 */

#include "olb3D.h"
#include "olb3D.hh"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <getopt.h>

using namespace olb;
using namespace olb::descriptors;
using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = D3Q27<>;

using BulkDynamics = BGKdynamics<T,DESCRIPTOR>;
typedef enum {periodic, local, damping} BoundaryType;

// Stores geometry information in form of material numbers
void prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperGeometry<T,3>& superGeometry,
                     std::shared_ptr<IndicatorF3D<T>> domain,
                     BoundaryType boundarytype,
                     int res,
                     bool debug,
                     int boundary_depth,
                     Vector<T,3> domain_lengths
                     )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl
        << "Materials:" << std::endl
        << "1 = Fluid" << std::endl
        << "2 = Check (should be empty)" << std::endl
        << "3 = Outer Boundaries (far field)" << std::endl;

  superGeometry.rename( 0, 2 );   // all nodes to temporary type

  T bd_pu = converter.getPhysLength(boundary_depth);
  if ( boundarytype == damping ) {
    bd_pu = converter.getPhysLength(boundary_depth);
  } else {
    bd_pu = converter.getPhysDeltaX();
  }
  Vector<T,3> extend( domain_lengths[0]-2*bd_pu, domain_lengths[1]-2*bd_pu, domain_lengths[2]-2*bd_pu );
  Vector<T,3> origin( -domain_lengths[0]/2+bd_pu, -domain_lengths[1]/2+bd_pu, -domain_lengths[2]/2+bd_pu );
  std::shared_ptr<IndicatorF3D<T>> fluid_domain = std::make_shared<IndicatorCuboid3D<T>>( extend, origin );
  // fluid_domain = cuboid2;
  if ( boundarytype == periodic ) {
    superGeometry.rename(2, 1);  // all nodes are fluid
  } else {
    superGeometry.rename(2, 1, fluid_domain);  // all nodes except boundary length
  }

  // defining boundary conditions
  if ( boundarytype != periodic ) {
    superGeometry.rename( 2, 3, domain-fluid_domain );   // remaining to outer/reflecting/far-field boundaries
  }

  superGeometry.communicate();

  // Removes all not needed boundary voxels outside the surface
  // superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  // superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}


template <typename T>
class AcousticPulse3D : public AnalyticalF3D<T,T> {
protected:
  T rho0;
  T amplitude;
  T alpha;
  UnitConverter<T, DESCRIPTOR> _converter;
public:
  AcousticPulse3D(UnitConverter<T, DESCRIPTOR> converter, T rho0, T amplitude, T alpha )
      : AnalyticalF3D<T,T>(1), rho0(rho0), amplitude(amplitude), alpha(alpha),
      _converter(converter) {};

  bool operator()(T output[], const T input[]) override
  {
    T x = input[0]*50+0.1;
    T y = input[1]*50+0.1;
    T z = input[2]*50+0.1;
    output[0] = rho0+amplitude*util::exp(-alpha*(x*x+y*y+z*z));
    return true;
  };
};

template <typename T>
class DampingTerm3D : public AnalyticalF3D<T,T> {
protected:
  Vector<T,3> x_0, Lx;
  int A=4;
  int n=4;
  T _boundary_depth_pu;
  UnitConverter<T, DESCRIPTOR> _converter;
  T _ds;
public:
  DampingTerm3D(UnitConverter<T, DESCRIPTOR> converter, int boundary_depth_lu, Vector<T,3> domain_lengths, T damping_strength )
      : AnalyticalF3D<T,T>(1), _converter(converter), _ds(damping_strength)
      {
        _boundary_depth_pu = converter.getPhysLength(boundary_depth_lu);
        for (size_t d=0; d<3; d++) {
          Lx[d] = domain_lengths[d]/2;  // Lx is usually half the domain
          x_0[d] = Lx[d] - _boundary_depth_pu;  // subtract boundary depth from domain length
          // x_0[d] /= Lx[d];  // normalize x_0 (Lx is still needed in operator() to normalize x; Lx will be replaced by 1)
        }
      };

  bool operator()(T output[], const T input[]) override
  {
    // normalize x to [0,1]
    Vector<T,3> x, distance_from_border;
    bool is_boundary = false;
    for (size_t d=0; d<3; d++) {
      x[d] = ( std::abs(input[d]) - (Lx[d] - _boundary_depth_pu) ) / _boundary_depth_pu;  // x_0 becomes 0
      x[d] = std::max(x[d], 0.);  // if x[d]<0, ignore for X; if x[d]>1, set to one to avoid negative values of sigma
      if ( x[d] > 0 ) {
        is_boundary = true;
      }
      distance_from_border[d] = Lx[d] - input[d];
    }
    if ( !is_boundary ) { output[0] = 0; return true; }

    T X=0.;
    for ( size_t d=0; d<3; d++ ) {
      if ( distance_from_border[d] == *std::min_element(std::begin(distance_from_border), std::end(distance_from_border)) ) {
        X = x[d];
      }
    }
    T sigma;
    sigma = A * ( ( (std::pow(X, n))*(1-X)*(n+1)*(n+2) ) / std::pow(_boundary_depth_pu,n+2) );  // x_0 left out
    sigma /= 9.8304 / std::pow(_boundary_depth_pu, 6);  // calculated as max of function (at X=0.8) for A=n=4
    sigma = std::max(sigma, 0.);
    sigma = std::min(sigma, 1.);  // hard set 0<sigma<1
    output[0] = _ds*sigma;
    return true;
  };
};

// Set up the geometry of the simulation
void prepareLattice(UnitConverter<T,DESCRIPTOR> const& converter,
                    SuperLattice<T,DESCRIPTOR>& sLattice,
                    SuperGeometry<T,3>& superGeometry,
                    T rho0,
                    T Ma,
                    T amplitude,
                    T alpha,
                    BoundaryType boundarytype,
                    int boundary_depth,
                    T damping_strength,
                    Vector<T,3> domain_lengths
                    )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 --> bulk dynamics
  auto bulkIndicator = superGeometry.getMaterialIndicator({1});
  sLattice.defineDynamics<BulkDynamics>(bulkIndicator);

  // Material=3 --> dynamics depend..
  clout << "boundarytype=" << boundarytype << std::endl;
  if ( boundarytype != periodic ) {
    sLattice.defineDynamics<BulkDynamics>(superGeometry.getMaterialIndicator({3}));
  }
  if ( boundarytype == local ) {
    setLocalPressureBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 3);
    setLocalVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 3);
  } else if ( boundarytype == damping ) {
    bulkIndicator = superGeometry.getMaterialIndicator({1, 3});
    setDampingBoundary<T,DESCRIPTOR>(sLattice, superGeometry, 3);
  }

  // Initial conditions
  AnalyticalConst3D<T,T> *ux, *uy, *uz;
  ux = new AnalyticalConst3D<T,T>( Ma );
  uy = new AnalyticalConst3D<T,T>( 0. );
  uz = new AnalyticalConst3D<T,T>( 0. );
  AnalyticalComposed3D<T,T> u( *ux, *uy, *uz );

  T ndatapoints = converter.getResolution();
  T dist = converter.getPhysDeltaX();
  
  AcousticPulse3D<T> pressureProfile( converter, rho0, amplitude, alpha );
  
  Gnuplot<T> gplot_hline_p( "shock_hline");
  Gnuplot<T> gplot_diag_p( "shock_diag");
  gplot_hline_p.setLabel("distance [m]", "density [LU]");
  gplot_diag_p.setLabel("distance [m]", "density [LU]");
  for (int n = 0; n <= int(ndatapoints/2); n++) {
    T input_hline[3] =  {n*dist, 0, 0};
    T output_hline[3];
    pressureProfile(output_hline, input_hline);
    gplot_hline_p.setData(input_hline[0], output_hline[0]);
    T input_diag[3] =  {n*dist, n*dist, 0};
    T output_diag[3];
    pressureProfile(output_diag, input_diag);
    gplot_diag_p.setData(input_diag[0], output_diag[0]);
  }
  gplot_hline_p.writePNG(-1, -1, "shock_hline");
  gplot_diag_p.writePNG(-1, -1, "shock_diag");

  //Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, pressureProfile, u );
  sLattice.iniEquilibrium( bulkIndicator, pressureProfile, u );

  // Define fields
  sLattice.setParameter<descriptors::OMEGA>( omega );
  AnalyticalConst3D<T,T> rhoField( rho0 );
  sLattice.defineField<descriptors::DENSITY>( superGeometry, 3, rhoField );
  AnalyticalConst3D<T,T> uxField( Ma );
  sLattice.defineField<descriptors::UX>( superGeometry, 3, uxField );
  AnalyticalConst3D<T,T> uyField( 0. );
  sLattice.defineField<descriptors::UY>( superGeometry, 3, uyField );
  AnalyticalConst3D<T,T> uzField( 0. );
  sLattice.defineField<descriptors::UZ>( superGeometry, 3, uzField );  

  if ( boundarytype == damping ) {
    // output damping layer parameter
    DampingTerm3D<T> sigma_plot( converter, boundary_depth, domain_lengths, damping_strength );
    Gnuplot<T> gplot_hline( "sigma_hline");
    Gnuplot<T> gplot_diag( "sigma_diag");
    gplot_hline.setLabel("distance [m]", "density [LU]");
    gplot_diag.setLabel("distance [m]", "density [LU]");
    for (int n = 0; n <= int(ndatapoints/2); n++) {
      T input_hline[3] =  {n*dist, 0, 0};
      T output_hline[3];
      sigma_plot(output_hline, input_hline);
      gplot_hline.setData(input_hline[0], output_hline[0]);
      T input_diag[3] =  {n*dist, n*dist, 0};
      T output_diag[3];
      sigma_plot(output_diag, input_diag);
      gplot_diag.setData(input_diag[0], output_diag[0]);
    }
    gplot_hline.writePNG(-1, -1, "sigma_hline");
    gplot_diag.writePNG(-1, -1, "sigma_diag");

    DampingTerm3D<T> sigma( converter, boundary_depth, domain_lengths, damping_strength );
    sLattice.defineField<descriptors::DAMPING>( superGeometry, 3, sigma );
  }
  
  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a sinusoidal pressure
void setBoundaryValues( UnitConverter<T,DESCRIPTOR> const& converter,
                        SuperLattice<T,DESCRIPTOR>& sLattice, int iT,
                        SuperGeometry<T,3>& superGeometry,
                        T rho0, T Ma,
                        BoundaryType boundarytype )
{
  if ( boundarytype != periodic and boundarytype != damping ) {
    OstreamManager clout( std::cout,"setBoundaryValues" );
    AnalyticalConst3D<T,T> ux( Ma );
    AnalyticalConst3D<T,T> uy( 0. );
    AnalyticalConst3D<T,T> uz( 0. );
    AnalyticalComposed3D<T,T> u( ux, uy, uz );
    AnalyticalConst3D<T,T> rho( rho0 );

    sLattice.defineRhoU( superGeometry, 3, rho, u );
    sLattice.iniEquilibrium( superGeometry.getMaterialIndicator({3}), rho, u );

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }
}

T L2Norm(SuperLattice<T,DESCRIPTOR>& sLattice,
        SuperGeometry<T,3>& superGeometry,
        int iT,
        UnitConverter<T,DESCRIPTOR> const& converter) {
  OstreamManager clout(std::cout, "L2Norm");

  T result[3] = {T(), T(), T()};
  int tmp[] = {int()};

  SuperLatticePhysPressure3D<T,DESCRIPTOR> pressure3D(sLattice, converter);
  AnalyticalConst<3,T,T> rho03D( 0. );

  auto indicatorF = superGeometry.getMaterialIndicator({1});
  SuperAbsoluteErrorL2Norm3D <T> absPressureErrorL2Norm(pressure3D, rho03D, indicatorF);

  absPressureErrorL2Norm(result, tmp);
  T l2_abs = result[0];

  return l2_abs;
}

// write data to termimal and file system
void getResults(SuperLattice<T,DESCRIPTOR>& sLattice,
                UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer,
                SuperPlaneIntegralFluxVelocity3D<T>& velocityFlux,
                SuperPlaneIntegralFluxPressure3D<T>& pressureFlux,
                T rho0, T amplitude,
                Gnuplot<T>& gplot_l2_abs, T Lp0,
                size_t iout, T tmax, size_t imax
                )
{
#ifdef PLATFORM_GPU_CUDA
  gpu::cuda::device::synchronize();
#endif
  OstreamManager clout( std::cout,"getResults" );
  // clout << "getResults startet at iT=" << iT << ", iout=" << iout << std::endl;
  SuperVTMwriter3D<T> vtmWriter( "movingshock3d" );

  if ( iT==0 ) {
    // Writes geometry, cuboid no. and rank no. to file system
//    SuperLatticeGeometry3D<T,DESCRIPTOR> geometry( sLattice, superGeometry );
//    SuperLatticeCuboid3D<T,DESCRIPTOR> cuboid( sLattice );
//    SuperLatticeRank3D<T,DESCRIPTOR> rank( sLattice );
//    vtmWriter.write( geometry );
//    vtmWriter.write( cuboid );
//    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Writes to console 10 times
  if ( iT%iout==0 ) {
    // sLattice.setProcessingContext(ProcessingContext::Evaluation);
    // velocityFlux.print();
    // pressureFlux.print();
    // write to terminal
    timer.update( iT );
    timer.printStep();
    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }

  if ( iT%std::max(int(iout/10), 1)==0 ) {
    gplot_l2_abs.setData(T(iT), L2Norm(sLattice, superGeometry, iT,
                                       converter) / Lp0 );
  }

  // often get LU density plot along x
  if ( iT%iout == 0 ) {
    SuperLatticePhysPressure3D<T,DESCRIPTOR> pressure( sLattice, converter );
    SuperLatticePhysVelocity3D<T,DESCRIPTOR> velocity( sLattice, converter );

    std::stringstream ss;
    ss << std::setw(4) << std::setfill('0') << iT;
    Gnuplot<T> gplot_hline( "pressure_hline_" + ss.str());
    Gnuplot<T> gplot_vline( "pressure_vline_" + ss.str());
    Gnuplot<T> gplot_diagonal( "pressure_diagonal_" + ss.str());
    Gnuplot<T> gplot_udiagonal( "velocity_diagonal_" + ss.str());
    gplot_hline.setLabel("distance [m]", "density [LU]");
    gplot_vline.setLabel("distance [m]", "density [LU]");
    gplot_diagonal.setLabel("distance [m]", "density [LU]");
    gplot_udiagonal.setLabel("distance [m]", "velocity [LU]");
    // Vectors for simulated solution
    T densities_hline[1] = {T()};
    T densities_vline[1] = {T()};
    T densities_diagonal[1] = {T()};
    T velocities_diagonal[1] = {T()};
    T dist = converter.getPhysDeltaX();
    T ndatapoints = converter.getResolution(); // number of data points on line
    // CSV<T> csvWriterConcentration;
    // save concentration along the middle of the PRF
    for (int n = -int(ndatapoints/2); n <= int(ndatapoints/2); n++) {
      T input_hline[3] =  {n*dist, 0, 0};
      AnalyticalFfromSuperF3D<T> interpolation_hline( pressure, true, true );
      interpolation_hline(densities_hline, input_hline);
    // csvWriterConcentration.writeDataFile(input_hline[0], densities_hline[0], "simulation" , 16);
      gplot_hline.setData(input_hline[0], densities_hline[0]);

      T input_vline[3] =  {0, n*dist, 0};
      AnalyticalFfromSuperF3D<T> interpolation_vline( pressure, true, true );
      interpolation_vline(densities_vline, input_vline);
      gplot_vline.setData(input_vline[1], densities_vline[0]);

      T input_diagonal[3] =  {n*dist, n*dist, n*dist};
      AnalyticalFfromSuperF3D<T> interpolation_diagonal( pressure, true, true );
      interpolation_diagonal(densities_diagonal, input_diagonal);
      // csvWriterConcentration.writeDataFile(input_diagonal[0], densities_diagonal[0], "simulation" , 16);
      gplot_diagonal.setData(input_diagonal[0], densities_diagonal[0]);
      AnalyticalFfromSuperF3D<T> interpolation_udiagonal( velocity, true, true );
      interpolation_udiagonal(velocities_diagonal, input_diagonal);
      gplot_udiagonal.setData(input_diagonal[0], velocities_diagonal[0]);
    }
    T ymin(converter.getPhysPressure(-amplitude/200));
    T ymax(converter.getPhysPressure(+amplitude/200));
    gplot_hline.setYrange(ymin, ymax);
    gplot_vline.setYrange(ymin, ymax);
    gplot_diagonal.setYrange(ymin, ymax);
    // plot is generated
    gplot_hline.writePNG(-1, -1, "pressure_hline");
    gplot_vline.writePNG(-1, -1, "pressure_vline");
    // gplot_hline.writePDF("densities_hline");
    gplot_diagonal.writePNG(-1, -1, "pressure_diagonal");
    gplot_udiagonal.writePNG(-1, -1, "velocity_diagonal");
    // gplot_diagonal.writePDF("densities_diagonal");

    // get VTK and images
    if ( iT%(iout*2)==0) {
      clout << "vtmWriter startet" << std::endl;
      SuperLatticePhysVelocity3D<T,DESCRIPTOR> velocity( sLattice, converter );
      SuperLatticePhysPressure3D<T,DESCRIPTOR> pressure( sLattice, converter );
      // VTK
      vtmWriter.addFunctor( velocity );
      vtmWriter.addFunctor( pressure );
      vtmWriter.write( iT );

      // pressure image
      BlockReduction3D2D<T> pressureReduction( pressure, Vector<T,3>({0, 0, 1}) );
      heatmap::plotParam<T> jpeg_ParamP;
      jpeg_ParamP.maxValue = converter.getPhysPressure(+amplitude/200);
      jpeg_ParamP.minValue = converter.getPhysPressure(-amplitude/200);
      jpeg_ParamP.colour = "rainbow";
      jpeg_ParamP.fullScreenPlot = true;
      heatmap::write(pressureReduction, iT, jpeg_ParamP);

      // velocity image
      SuperEuklidNorm3D<T> normVel( velocity );
      // BlockReduction3D2D<T> planeReduction( normVel, origin, u, v, 600, BlockDataSyncMode::ReduceOnly );
      BlockReduction3D2D<T> planeReduction( normVel, Vector<T,3>({0, 0, 1}) );
      heatmap::plotParam<T> plotParam;
      jpeg_ParamP.maxValue = converter.getCharPhysVelocity()*1.1;
      jpeg_ParamP.minValue = converter.getCharPhysVelocity()*0.9;
      jpeg_ParamP.colour = "rainbow";
      jpeg_ParamP.fullScreenPlot = true;
      heatmap::write(planeReduction, iT, plotParam);
    }
  }

  // Saves lattice data
  if ( iT%imax==0 && iT>0 ) {
    //clout << "Checkpointing the system at t=" << iT << std::endl;
    //sLattice.save( "movingshock3d.checkpoint" );
    // The data can be reloaded using
    //     sLattice.load("movingshock3d.checkpoint");
  }
}


int main( int argc, char* argv[], char *envp[] )
{
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  if ( singleton::mpi().isMainProcessor() ) {
    std::copy(argv, argv + argc, std::ostream_iterator<char *>(std::cout, " "));
    printf("\n");
  }
  std::copy(argv, argv + argc, std::ostream_iterator<char *>(std::cout, " "));
  printf("\n");
  singleton::directories().setOutputDir( "./tmp/" );  // set output directory
  OstreamManager clout( std::cout, "main" );

  // === get Command Line Arguments
  static const struct option long_options[] =
      {
          { "rho0",         required_argument, 0,  1  },
          { "res",          required_argument, 0,  2  },
          { "Ma",           required_argument, 0,  3  },
          { "targetMa",     required_argument, 0,  4  },
          { "amplitude",    required_argument, 0,  5  },
          { "length",       required_argument, 0, 'l' },
          { "tau",          required_argument, 0, 't' },
          { "imax",         required_argument, 0,  6  },
          { "tmax",         required_argument, 0,  7  },
          { "iout",         required_argument, 0,  8  },
          { "nout",         required_argument, 0,  9  },
          { "Re",           required_argument, 0,  10 },
          { "boundarycondition", required_argument, 0, 'b' },
          { "debug",        no_argument,       0,  11 },
          { "boundarydepth",   required_argument, 0,  12 },
          { "damping_strength", required_argument, 0, 13 },
          0
  };
  // initialize arguments
  // Parameters for the simulation setup
  const T b_shock           = 3;
  const T alpha_shock       = log(2)/(b_shock*b_shock);

  int res       = 0;
  size_t imax   = 0;
  size_t iout   = 0;
  size_t nout   = 0;
  T tmax        = 0.;
  T rho0        = 0.;
  T Ma          = -100.;
  T targetMa    = 0.;
  T amplitude   = 0.;
  T tau         = 0.;
  T lengthDomain = 0.;
  T Re          = 0.;
  int BC_type   = 0;
  int SC_type   = 0;
  int boundary_depth = 0;
  T damping_strength = 0.;
  bool debug    = false;

  while (1)
  {
    int index = -1;
    struct option * opt = 0;
    int result = getopt_long(argc, argv, "abc:d", long_options, &index);
    if (result == -1) break; /* end of list */
    switch (result)
    {
    case 1: /* same as index==0 */
      rho0 = T(atof(optarg));
      clout << "'rho0' was specified to " << optarg << std::endl;
      break;
    case 2: /* same as index==1 */
      res = atoi(optarg);
      clout << "'res' was specified to " << optarg << std::endl;
      break;
    case 3: /* same as index==1 */
      Ma = T(atof(optarg));
      clout << "'Ma' was specified to " << optarg << std::endl;
      break;
    case 4: /* same as index==1 */
      targetMa = T(atof(optarg));
      clout << "'targetMa' was specified to " << optarg << std::endl;
      break;
    case 5: /* same as index==1 */
      amplitude = T(atof(optarg));
      clout << "'amplitude' was specified to " << optarg << std::endl;
      break;
    case 6: /* same as index==1 */
      imax = atoi(optarg);
      clout << "'imax' was specified to " << optarg << std::endl;
      break;
    case 7: /* same as index==1 */
      tmax = T(atof(optarg));
      clout << "'tmax' was specified to " << optarg << std::endl;
      break;
    case 8: /* same as index==1 */
      iout = atoi(optarg);
      clout << "'iout' was specified to " << optarg << std::endl;
      break;
    case 9: /* same as index==1 */
      nout = atoi(optarg);
      clout << "'nout' was specified to " << optarg << std::endl;
      break;
    case 10: /* same as index==1 */
      Re = T(atof(optarg));
      clout << "'Re' was specified to " << optarg << std::endl;
      break;
    case 'l': /* same as index==2 */
      lengthDomain = T(atof(optarg));
      clout << "'l'/'length' was specified to " << optarg << std::endl;
      break;
    case 't': /* same as index==3 */
      tau = T(atof(optarg));
      clout << "'t'/'tau' was specified to " << optarg << std::endl;
      break;
    case 'b':
      BC_type = atoi(optarg);
      clout << "'b'/'boundarycondition' was specified to " << optarg << std::endl;
      break;
    case 11:
      debug = true;
      clout << "'debug' was activated" << std::endl;
      break;
    case 12:
      boundary_depth = atoi(optarg);
      clout << "'boundarydepth' was set to " << optarg << std::endl;
      break;
    case 13:
      damping_strength = T(atof(optarg));
      clout << "'damping_strength' was set to " << optarg << std::endl;
      break;
    case 14:
      SC_type = atoi(optarg);
      clout << "'sc_type' was set to " << optarg << std::endl;
      break;
    case 0: /* all parameter that do not appear in the optstring */
      opt = (struct option *)&(long_options[index]);
      clout << "'" << opt->name << "' was specified to ";
      if (opt->has_arg == required_argument) {
        clout << " to " << optarg;
      }
      clout << ".\n";
      break;
    default: /* unknown */
      break;
    }
  }
  /* print all other parameters */
  while (optind < argc)
  {
    clout << "other parameter: <" << argv[optind++] << ">\n";
  }
  /* default values */
  BoundaryType boundarytype;
  switch ( BC_type ) {
    case 1: boundarytype = periodic; clout << "Boundary condition type specified to periodic." << std::endl; break;
    case 2: boundarytype = local; clout << "Boundary condition type specified to local." << std::endl; break;
    case 3: boundarytype = damping; clout << "Boundary condition type specified to damping." << std::endl; break;
    default: boundarytype = local; clout << "Boundary condition type not specified. Default to local." << std::endl; break;
  }
  if ( res == 0 ) {
    if ( debug )            res             = 20;
    else                    res             = 81;
  }
  if ( rho0 == 0 )          rho0            = 1.;
  if ( imax == 0 )          imax            = 10000;
  if ( nout == 0 )          nout            = 20;
  if ( boundary_depth == 0 ) boundary_depth = 20;
  if ( tmax == 0 ) {
    if ( debug )            tmax            = 1;
    else                    tmax            = T(10)/300.;  // should be n*100[LU]
  }          
  if ( Ma == -100 )         Ma              = 0.25;
  if ( targetMa == 0 )      targetMa        = 0.5;
  if ( amplitude == 0 )     amplitude       = 0.1;
  if ( damping_strength == 0 ) damping_strength = 1.;
  
  if ( lengthDomain == 0 )  lengthDomain    = 1.;
  T heightDomain            = 2*lengthDomain;
  T depthDomain             = 2*lengthDomain;

  const T charL             = lengthDomain;
  T charV                   = std::max(Ma*1.1, 0.1);  // /targetMa;
  const T cs_LU             = 1 / std::sqrt(3.0);
  T viscosity;
  if ( tau == 0 ) {
    if ( Re == 0 ) {
      viscosity             = 1.48e-5; // 1e-2;
      Re                    = charV * charL / viscosity;
    } else viscosity        = charV * charL / Re;
    tau                     = viscosity / (cs_LU*cs_LU) + 0.5;
  }
  else {
    if ( Re == 0 ) Re       = 20000;
    viscosity               = charV * charL / Re;
  }
  if ( debug ) {
    Re                      = 200;
    charV                   = 0.25;
    heightDomain            = lengthDomain;
    depthDomain             = lengthDomain;
    viscosity               = charV * charL / Re;
    tau                     = viscosity / (cs_LU*cs_LU) + 0.5;
  }

  clout << "Input to unitconverter: res=" << res << ", tau=" << tau << ", lengthDomain=" << lengthDomain << ", charV=" << charV << ", viscosity=" << viscosity << ", rho0=" << rho0 << std::endl;
  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> converter(
    (size_t) res,               // resolution: number of voxels per charPhysL
    (T)      tau,               // relaxation time
    (T)      lengthDomain,      // charPhysLength: reference length of simulation geometry
    (T)      charV,             // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)      viscosity,         // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)      rho0               // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("movingshock3d");

  // === 2nd Step: Prepare Geometry ===

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 6;
#endif

  // setup domain
  Vector<T,3> domain_lengths = {lengthDomain, heightDomain, depthDomain};
  Vector<T,3> originDomain( - lengthDomain/2, - heightDomain/2, - depthDomain/2 );  //
  Vector<T,3> extendDomain( lengthDomain, heightDomain, depthDomain );  // size of the domain
  std::shared_ptr<IndicatorF3D<T>> domain =
      std::make_shared<IndicatorCuboid3D<T>>( extendDomain, originDomain );

  CuboidGeometry3D<T> cuboidGeometry(*(domain),
                                     converter.getConversionFactorLength(),
                                     noOfCuboids
                                     );
  if ( boundarytype == periodic || boundarytype == damping ) {
    cuboidGeometry.setPeriodicity(true, true, true);
  } else {
    cuboidGeometry.setPeriodicity(false, false, false);
  }

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,3> superGeometry( cuboidGeometry, loadBalancer );

  T boundary_depth_pu = converter.getPhysLength(boundary_depth);
  clout << "Setup: debug=" << debug << "; boundary_depth=" << boundary_depth << "; bd_depth_pu=" << boundary_depth_pu << "; damping_strength=" << damping_strength << "; overlap=" << superGeometry.getOverlap() << std::endl;
  prepareGeometry( converter, superGeometry, domain, boundarytype, res, debug, boundary_depth, domain_lengths );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T,DESCRIPTOR> sLattice( superGeometry, debug );

  //prepare Lattice and set boundaryConditions
  prepareLattice( converter, sLattice, superGeometry, rho0, Ma, amplitude, alpha_shock, boundarytype, boundary_depth, damping_strength, domain_lengths );

  // instantiate reusable functors
  SuperPlaneIntegralFluxVelocity3D<T> velocityFlux(
      sLattice,
      converter,
      superGeometry,
      { lengthDomain/T(2), heightDomain/T(2), depthDomain/T(2) },
      Vector<T, 3>( 0., 1., 1. )
  );

  SuperPlaneIntegralFluxPressure3D<T> pressureFlux(
      sLattice,
      converter,
      superGeometry,
      { lengthDomain/T(2), heightDomain/T(2), depthDomain/T(2) },
      Vector<T, 3>( 0., 1., 1. )
  );

  Gnuplot<T> gplot_l2_abs("l2_absolute");//, Gnuplot<T>::LOGLOG, Gnuplot<T>::OFF);
  gplot_l2_abs.setLabel("time []", "absolute L2 norm []");
  T Lp0 = L2Norm( sLattice, superGeometry, 0, converter );

  clout << "tmax=" << tmax << "[PU]=" << converter.getLatticeTime( tmax ) << "[PU], while imax=" << imax << " as input." << std::endl;
  imax = std::min( converter.getLatticeTime( tmax ), imax );
  tmax = converter.getPhysTime( imax );
  if ( debug ) { imax=100; iout=10;  }
  if ( iout == 0 )          iout            = imax / nout;
  clout << "tmax=" << tmax << ", while imax=" << imax << " recalculated. iout=" << iout << std::endl;

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer(converter.getLatticeTime( tmax ),
                       superGeometry.getStatistics().getNvoxel()
                       );
  timer.start();

  size_t iT = 0;
  while ( iT < imax ) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( converter, sLattice, iT, superGeometry, rho0, Ma, boundarytype );
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults(sLattice, converter, iT, superGeometry, timer,
               velocityFlux, pressureFlux,
               rho0, amplitude,
               gplot_l2_abs, Lp0,
               iout, tmax, imax
               );
    // test density failure
    if (sLattice.getStatistics().getAverageRho() > 2.) {
      clout << "breaking simulation loop because density is too high: "
            << sLattice.getStatistics().getAverageRho() << std::endl;
      break;
    }
    // if (sLattice.getStatistics().getAverageEnergy() > 2.) {
    //   clout << "breaking simulation loop because energy is too high: "
    //         << sLattice.getStatistics().getAverageEnergy() << std::endl;
    //   break;
    // }
    iT++;
  }
  clout << "Simulation stopped after " << iT << " iterations." << std::endl;

  gplot_l2_abs.setYrange(1e-3, 1);
  gplot_l2_abs.setLogScale(2);
  gplot_l2_abs.writePNG(-1, -1, "gplot_l2_abs");
  timer.stop();
  timer.printSummary();
}
