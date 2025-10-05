/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Felix Schuhmann, Shota Ito
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

/* ellipsoidShape3d.cpp:
 * This example examines a steady flow past an ellipsoid placed inside a rectangular channel.
 * The half-axes are alligned along the x, y and z coordinate directions and periodic boundary treatment
 * is applied in y and z directions.
 *
 * Contains:
 * - Standard channel flow simulation
 * - Flow simulation while computing derivatives with respect to ellipsoid radii in y and z direction using ADf
 * - Shape optimization using LBFGS:
 *   radii in y and z direction are optimized to minimize L2-Norm of dissipation over the entire channel
 *   radius in x direction is adjusted to enforce fixed ellipsoid volume
 */

#include "ellipsoidShape3d.h"

template<typename T>
T objective(Vector<T,2> radiusEllipsoidYZ){
  return simulateEllipsoid3D<T>( radiusEllipsoidYZ );
}

int main( int argc, char* argv[] )
{
  initialize( &argc, &argv );
  OstreamManager clout(std::cout, "main");

  if constexpr (false){
    Vector<S,2> radiusEllipsoidYZ(.05, .05);
    clout << simulateEllipsoid3D<S>( radiusEllipsoidYZ ) << std::endl;
  }

  if constexpr (false) {
    Vector<U,2> radiusEllipsoidYZ(.05, .05);
    radiusEllipsoidYZ[0].setDiffVariable(0);
    radiusEllipsoidYZ[1].setDiffVariable(1);
    clout << simulateEllipsoid3D<U>( radiusEllipsoidYZ ) << std::endl;
  }

  if constexpr (true){
    OptiCaseAD<S,2,VectorHelp> optiCase(
      objective<S>,
      objective<U>);
    OptimizerLBFGS<S,Vector<S,2>> optimizer(
    2, 1.e-16, 10, .01, 10, "StrongWolfe", 20, 1.e-4, true, "", "log",
    true, 0.19, true, 0.01, false, 0., true,
    {OptimizerLogType::value, OptimizerLogType::control, OptimizerLogType::derivative});

    Vector<S,2> startValue(0.08, 0.08);
    optimizer.setControl(startValue);
    optimizer.optimize(optiCase);
  }
}
