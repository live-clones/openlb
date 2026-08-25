/*  This file is part of the OpenLB library
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

/* spongeLayerDynamics.h:
 * Implements the absorbing layer described in equation (24) of
 * H. Xu and P. Sagaut, “Analysis of the absorbing layers for the
 * weakly-compressible lattice Boltzmann methods”, Journal of
 * Computational Physics, vol. 245, pp. 14–42, Jul. 2013,
 * doi: 10.1016/j.jcp.2013.02.051.
*/

#ifndef SPONGE_LAYER_DYNAMICS_H
#define SPONGE_LAYER_DYNAMICS_H

/*Whatever collision operator whithout turbulence model*/
namespace olb {

//===================================================================================
//================= DampingDynamics =========
//===================================================================================
template<typename T, typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
struct SpongeLayerDynamics : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
  using parameters = typename meta::list<descriptors::OMEGA>;
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;
  using CollisionO = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

  template <typename NEW_T>
  using exchange_value_type = SpongeLayerDynamics<NEW_T,DESCRIPTOR,MOMENTA,EQUILIBRIUM,COLLISION>;

  template <typename M>
  using exchange_momenta = SpongeLayerDynamics<T,DESCRIPTOR,M,EQUILIBRIUM,COLLISION>;

  std::type_index id() override {
    return typeid(SpongeLayerDynamics);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<SpongeLayerDynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    V pi[util::TensorVal<DESCRIPTOR>::n] { };
    MomentaF().computeStress(cell, pi);
    const V sigma = cell.template getField<descriptors::DAMPING>();
    const V omega = parameters.template get<descriptors::OMEGA>();
    V fEq[DESCRIPTOR::q] { };
    const auto statistic = EquilibriumF().compute(cell, parameters, fEq);
    const V rhoRef = cell.template getField<descriptors::DENSITY>();
    const V uxRef = cell.template getField<descriptors::UX>();
    const V uyRef = cell.template getField<descriptors::UY>();
    const V uzRef = cell.template getField<descriptors::UZ>();
    const V uRef[3] = {uxRef, uyRef, uzRef};
    V fEqRef[DESCRIPTOR::q] { };
    EquilibriumF().compute(cell, rhoRef, uRef, fEqRef);
    
    V fEqSum=0.0;
    V rho, u[DESCRIPTOR::d];
    MomentaF().computeRhoU(cell, rho, u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      fEqSum=fEqSum+ sigma*(fEqRef[iPop]-fEq[iPop]);
    }
    V rhoRef2 = rho+0.5*fEqSum;
    V uRef2[DESCRIPTOR::d];
    for (int iDim=0; iDim < DESCRIPTOR::d; ++iDim) {
      uRef2[iDim]=u[iDim];
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        uRef2[iDim]+=V {0.5}*descriptors::c<DESCRIPTOR>(iPop,iDim)*sigma*(fEqRef[iPop]-fEq[iPop])
        /rhoRef2;
      }
    }
    V fEqRef2[DESCRIPTOR::q] { };
    EquilibriumF().compute(cell, rhoRef2, uRef2, fEqRef2);
    CollisionO().apply(cell, parameters);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] += omega * ( fEqRef2[iPop] - fEq[iPop] ) + sigma * ( fEqRef[iPop] - fEqRef2[iPop] );
    }
    return statistic;
  };
  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    EquilibriumF().compute(cell, rho, u, fEq);
  };


  std::string getName() const override {
    return "SpongeLayerDynamics";
  };

};

}

#endif
