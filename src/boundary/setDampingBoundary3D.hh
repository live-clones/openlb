/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023 Dennis Teutscher, Alexander Schulz
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

//This file contains the Damping Boundary
//This is an onLattice boundary
#ifndef SET_DAMPING_BOUNDARY_HH
#define SET_DAMPING_BOUNDARY_HH

#include "setDampingBoundary3D.h"



namespace olb {

///Initialising the setDampingBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setDampingBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,3>& superGeometry, int material)
{
  setDampingBoundary<T, DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material));
}

///Initialising the setDampingBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setDampingBoundary( SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setDampingBoundary");
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setDampingBoundary( sLattice.getBlock(iCloc), indicator->getBlockIndicatorF(iCloc) );
  }
}

/// Set Damping boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR>
void setDampingBoundary( BlockLattice<T,DESCRIPTOR>& block, BlockIndicatorF3D<T>& indicator )
{
  OstreamManager clout(std::cout, "setDampingBoundary");
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    Dynamics<T,DESCRIPTOR>* dynamics = block.template getDynamics<DampingBoundaryDynamics<T,DESCRIPTOR,momenta::BulkTuple,equilibria::SecondOrder>>();
    setBoundary(block, iX,iY,iZ, dynamics);
  });
}
}//namespace olb

#endif
