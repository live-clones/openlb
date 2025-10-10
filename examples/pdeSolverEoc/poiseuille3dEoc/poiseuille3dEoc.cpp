/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2018 Marc Haußmann, Mathias J. Krause, Jonathan Jeppener-Haltenhoff
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

/* poiseuille3d.cpp:
 * This example examines a 3D Poseuille flow
 * It illustrates the usage of different flow setups, boundary conditions
 * and computation of error norms.
 * The following simulation options can be combined freely:
 * - if the compiler flag ENABLE_MRT is set, mrt collision operators are used
 * - forced/ nonForced flow
 * - different boundary conditions
 * - simulation only or eoc convergence analysis
 */

// the main code of the simulation is in poiseuille3d.h as it is also used by the
// example ../../laminar/poiseuille3d

// set flag in order to use mrt collision operators instead of bgk
//#define ENABLE_MRT

#include "../../laminar/poiseuille3d/case.h"

void simulatePoiseuilleForEOC(MyCase::ParametersD& parameters, Gnuplot<MyCase::value_t>& gplot) {

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(parameters);
  
  /// === Step 4: Create Case ===
  MyCase myCase(parameters, mesh);
  
  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);
  
  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase);
  
  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues(myCase);
  
  /// === Step 8: Simulate ===
  simulate(myCase);

}

int main( int argc, char* argv[] )
{
  OstreamManager clout( std::cout,"main" );

  // === 1st Step: Initialization ===
  initialize( &argc, &argv );
  MyCase::ParametersD myCaseParameters;
  setGetParameters(myCaseParameters, argc, argv);
  myCaseParameters.set<parameters::EOC>(true);
  
  BoundaryType boundaryType = myCaseParameters.get<parameters::BOUNDARY_TYPE>();
  bool forbiddenEOCCombination = (boundaryType == FREE_SLIP) || (boundaryType == PARTIAL_SLIP);
  if ( forbiddenEOCCombination ) std::runtime_error("eoc computation is currently not supported for slip boundary conditions");

  std::string runName = "bc" + std::to_string(int(myCaseParameters.get<parameters::BOUNDARY_TYPE>())) + "_force" + std::to_string(int(myCaseParameters.get<parameters::FLOW_TYPE>()));
  singleton::directories().setOutputDir( "./tmp/" + runName "/" );

  // Initialize gnuplot
  Gnuplot<MyCase::value_t> gplot(
    "Velocity_and_StrainRate_eoc",
    false,
    "set terminal png size 720, 720 font 'Arial,10'",
    Gnuplot<MyCase::value_t>::LOGLOG,
    Gnuplot<MyCase::value_t>::LINREG);
  
  // set the labels for the plot
  gplot.setLabel("Resolution test", "average Error");

  size_t startN = myCaseParameters.get<parameters::EOC_START_RESOLUTION>();
  size_t maxN   = myCaseParameters.get<parameters::EOC_MAX_RESOLUTION>();
  size_t stepN  = myCaseParameters.get<parameters::EOC_RESOLUTION_STEP>();

  // loop over the different simulations
  for(size_t simuN = startN; simuN < maxN; simuN += stepN){
    /// Run the simulations
    clout << "Starting next simulation with N = " << simuN << std::endl;
    myCaseParameters.set<parameters::RESOLUTION>(simuN);
    simulatePoiseuilleForEOC(myCaseParameters, gplot);
  }

  gplot.writePNG(-1,-1,runName);

  return 0;
}
