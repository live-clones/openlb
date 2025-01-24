/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
 *                2022 Nando Suntoyo, Adrian Kummerlaender, Shota Ito
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

/*
 * Implements the absorbing layer described in H. Xu and P. Sagaut, “Analysis of the absorbing layers for the weakly-compressible lattice Boltzmann methods,” Journal of Computational Physics, vol. 245, pp. 14–42, Jul. 2013, doi: 10.1016/j.jcp.2013.02.051.

*/

#ifndef DAMPING_BOUNDARY_DYNAMICS_H
#define DAMPING_BOUNDARY_DYNAMICS_H

// #include "dynamics/latticeDescriptors.h"
#include "dynamics/dynamics.h"

namespace olb {

namespace boundaryhelper {

//===================================================================================
//================= DampingDynamics =========
//===================================================================================
template<typename T, typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
struct DampingBoundaryDynamics final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
  using parameters = typename meta::list<descriptors::OMEGA>;
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

  std::type_index id() override {
    return typeid(DampingBoundaryDynamics);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<DampingBoundaryDynamics>>();
  }

  using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
    V fEq[DESCRIPTOR::q] { };
    const auto statistic = EquilibriumF().compute(cell, parameters, fEq);
    const V rhoRef = cell.template getField<descriptors::DENSITY>();
    const V uxRef = cell.template getField<descriptors::UX>();
    const V uyRef = cell.template getField<descriptors::UY>();
    const V uzRef = cell.template getField<descriptors::UZ>();
    const V uRef[3] = {uxRef, uyRef, uzRef};
    
    const V omega = parameters.template get<descriptors::OMEGA>();
    const V sigma = cell.template getField<descriptors::DAMPING>();
    
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      const V fEqRef = computeEquilibrium(iPop, rhoRef, uRef);
      cell[iPop] = cell[iPop] + omega * ( fEq[iPop] - cell[iPop]) + sigma * ( fEqRef - fEq[iPop] );
    }
    return statistic;
  };

  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d]) const override any_platform {
    return equilibrium<DESCRIPTOR>::template secondOrder(iPop, rho, u);
  };

  std::string getName() const override {
    return "DampingBoundaryDynamics";
  };
};
}
}
#endif
