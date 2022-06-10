/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Adrian Kummerlaender
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

#ifndef NAVIER_STOKES_ADVECTION_DIFFUSION_COUPLING_H
#define NAVIER_STOKES_ADVECTION_DIFFUSION_COUPLING_H

#include "core/operator.h"

namespace olb {

/**
* Class for the coupling between a Navier-Stokes (NS) lattice and an
* Advection-Diffusion (AD) lattice.
*/
struct NavierStokesAdvectionDiffusionCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct FORCE_PREFACTOR : public descriptors::FIELD_BASE<0,1> { };
  struct T0 : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<FORCE_PREFACTOR,T0>;

  int getPriority() const {
    return 0;
  }

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;

    auto& cellNSE = cells.template get<names::NavierStokes>();
    auto& cellADE = cells.template get<names::Temperature>();

    // computation of the bousinessq force
    auto force = cellNSE.template getFieldPointer<descriptors::FORCE>();
    auto forcePrefactor = parameters.template get<FORCE_PREFACTOR>();
    V temperatureDifference = cellADE.computeRho() - parameters.template get<T0>();
    for (unsigned iD = 0; iD < DESCRIPTOR::d; ++iD) {
      force[iD] = forcePrefactor[iD] * temperatureDifference;
    }
    // Velocity coupling
    V u[DESCRIPTOR::d] { };
    cellNSE.computeU(u);
    cellADE.template setField<descriptors::VELOCITY>(u);
  }

};

}

#endif
