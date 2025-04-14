/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023
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

//This file contains the Sponge Layer Damping Boundary
//This is an onLattice boundary

#ifndef OLB_SPONGE_LAYER_BC_H
#define OLB_SPONGE_LAYER_BC_H

#include "boundary/dynamics/dampingBoundaryDynamics.h"

namespace olb {

namespace boundary {

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
struct SpongeLayer {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 0;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  return meta::id<SpongeLayerDynamics<T,DESCRIPTOR,momenta::BulkTuple,equilibria::SecondOrder>>();
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  return std::nullopt;
}

};

}


}

#endif
