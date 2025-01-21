/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Tim Bingert
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

#ifndef PHASE_FIELD_COUPLING_H
#define PHASE_FIELD_COUPLING_H

#include "core/operator.h"

namespace olb {

// =========================================================================//
// =============Hybrid Allen-Cahn + signed distance function================//
// =========================================================================//

struct initialPsi {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::INTERFACE_WIDTH>;

  int getPriority() const {
    return 1;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::value_t;
    auto w = parameters.template get<descriptors::INTERFACE_WIDTH>();

    V phi = cell.template getFieldComponent<descriptors::STATISTIC>(0);
    V psi = 0;
    if (phi > (0.99)) {
      psi += -4.595120;
    }
    else if (phi < (0.01)) {
      psi += 4.595120;
    }
    else {
      psi += util::log(1/phi-1);
    }
    psi *= -w/4.;
    cell.template setField<descriptors::SCALAR>(psi/util::sqrt(psi*psi+1.));
    cell.template setField<descriptors::PSI>(psi);
    cell.template setField<descriptors::PSI0>(psi);
  }

};

struct normGradPsi {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::INTERFACE_WIDTH>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;

    auto psi = cell.template getField<descriptors::PSI>();
    auto psi0 = cell.template getField<descriptors::PSI0>();
    auto w = parameters.template get<descriptors::INTERFACE_WIDTH>();
    Vector<V,DESCRIPTOR::d> gradPsi{};
    V normGradPsi = 1.;
    if (fabs(psi) <= 3.*w) {
      V a = psi - cell.neighbor({-1,0}).template getField<descriptors::PSI>();
      V b = cell.neighbor({1,0}).template getField<descriptors::PSI>() - psi;
      V c = psi - cell.neighbor({0,-1}).template getField<descriptors::PSI>();
      V d = cell.neighbor({0,1}).template getField<descriptors::PSI>() - psi;
      V gradX, gradY;
      if (psi0 > 0) {
        if (a<0) a = 0;
        if (b>0) b = 0;
        if (c<0) c = 0;
        if (d>0) d = 0;
        gradX = a;
        if (b*b>a*a) gradX = b;
        gradY = c;
        if (d*d>c*c) gradY = d;
        normGradPsi = util::sqrt(gradX*gradX+gradY*gradY);
      }
      else if (psi0 < 0) {
        if (a>0) a = 0;
        if (b<0) b = 0;
        if (c>0) c = 0;
        if (d<0) d = 0;
        gradX = a;
        if (b*b>a*a) gradX = b;
        gradY = c;
        if (d*d>c*c) gradY = d;
        normGradPsi = util::sqrt(gradX*gradX+gradY*gradY);
      }
      else {}
      cell.template setField<descriptors::NORMGRADPSI>(normGradPsi);
    } else {
      cell.template setField<descriptors::NORMGRADPSI>(normGradPsi);
    }
  }
};

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal>
struct normGradPsiBoundary2D {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::INTERFACE_WIDTH>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    auto psi = cell.template getField<descriptors::PSI>();
    auto psi0 = cell.template getField<descriptors::PSI0>();
    auto w = parameters.template get<descriptors::INTERFACE_WIDTH>();
    Vector<T,DESCRIPTOR::d> gradPsi{};
    T normGradPsi = 1.;
    if (fabs(psi) <= 3.*w) {
      T a = psi - cell.neighbor({-1,0}).template getField<descriptors::PSI>();
      T b = cell.neighbor({1,0}).template getField<descriptors::PSI>() - psi;
      T c = psi - cell.neighbor({0,-1}).template getField<descriptors::PSI>();
      T d = cell.neighbor({0,1}).template getField<descriptors::PSI>() - psi;
      if (xNormal == 1) b = a;
      else if (xNormal == -1) a = b;
      if (yNormal == 1) d = c;
      else if (yNormal == -1) c = d;
      T gradX, gradY;
      if (psi0 > 0) {
        if (a<0) a = 0;
        if (b>0) b = 0;
        if (c<0) c = 0;
        if (d>0) d = 0;
        gradX = a;
        if (b*b>a*a) gradX = b;
        gradY = c;
        if (d*d>c*c) gradY = d;
        normGradPsi = util::sqrt(gradX*gradX+gradY*gradY);
      }
      else if (psi0 < 0) {
        if (a>0) a = 0;
        if (b<0) b = 0;
        if (c>0) c = 0;
        if (d<0) d = 0;
        gradX = a;
        if (b*b>a*a) gradX = b;
        gradY = c;
        if (d*d>c*c) gradY = d;
        normGradPsi = util::sqrt(gradX*gradX+gradY*gradY);
      }
      else {}
      cell.template setField<descriptors::NORMGRADPSI>(normGradPsi);
    } else {
      cell.template setField<descriptors::NORMGRADPSI>(normGradPsi);
    }
  }
};

struct psiEvolve {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  struct DELTAT    : public descriptors::FIELD_BASE<1> { };
  using parameters = meta::list<DELTAT,descriptors::INTERFACE_WIDTH>;

  int getPriority() const {
    return 1;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::value_t;
    auto dt = parameters.template get<DELTAT>();
    auto w = parameters.template get<descriptors::INTERFACE_WIDTH>();
    auto s = cell.template getField<descriptors::SCALAR>();
    auto psi = cell.template getField<descriptors::PSI>();
    V psi_new = psi;
    if (fabs(psi) >= 3.*w) cell.template setField<descriptors::PSI>(psi_new);
    else {
      auto normGradPsi = cell.template getField<descriptors::NORMGRADPSI>();
      psi_new = psi-dt*s*(normGradPsi-1.);
      cell.template setField<descriptors::PSI>(psi_new);
    }
  }
};

struct AllenCahnNonLocalHelper {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::INTERFACE_WIDTH>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::value_t;
    V phi = cell.template getFieldComponent<descriptors::STATISTIC>(0);
    auto w = parameters.template get<descriptors::INTERFACE_WIDTH>();

    V psi = cell.template getField<descriptors::PSI>();
    V interfaceIndicator = 0;
    if (fabs(psi) < 3.*w) interfaceIndicator = 1;
    V top = (1.-interfaceIndicator)*phi*(phi-1.)*(phi-0.5);
    V bottom = (1.-interfaceIndicator)*(1.-phi)*phi;

    cell.template setField<descriptors::TOP>(top);
    cell.template setField<descriptors::BOTTOM>(bottom);
  }
};

struct AllenCahnPostProcessor {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  struct SIGMA       : public descriptors::FIELD_BASE<1> { };
  struct W           : public descriptors::FIELD_BASE<1> { };
  struct TAUS        : public descriptors::FIELD_BASE<3> { };
  struct RHOS        : public descriptors::FIELD_BASE<2> { };
  struct NONLOCALITY : public descriptors::FIELD_BASE<1> { };
  using parameters = meta::list<SIGMA,W,TAUS,RHOS,NONLOCALITY>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELL::template value_t<names::NavierStokes>::descriptor_t;

    auto& cellNS = cells.template get<names::NavierStokes>();
    auto& cellAC = cells.template get<names::Component1>();

    V phi = cellAC.template getFieldComponent<descriptors::STATISTIC>(0);
    Vector<V,DESCRIPTOR::d> gradPhi{};
    V laplacePhi = 0;
    //TODO: use lattice gradient schemes here from Cahn-Hilliard implementation
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      const V phi_i = cellAC.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getFieldComponent<descriptors::STATISTIC>(0);
      gradPhi += phi_i * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
      laplacePhi += 2*(phi_i - phi) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
    }

    auto sigma = parameters.template get<SIGMA>();
    auto w = parameters.template get<W>();
    auto rho_v = parameters.template get<RHOS>()[0];
    auto rho_l = parameters.template get<RHOS>()[1];
    V rho = rho_v + (rho_l-rho_v)*phi;
    auto tau_v = parameters.template get<TAUS>()[0];
    auto tau_l = parameters.template get<TAUS>()[1];
    V tau = tau_v + (tau_l-tau_v)*phi;
    auto gamma = parameters.template get<NONLOCALITY>();
    auto tau_mobil = parameters.template get<TAUS>()[2];
    V M = (tau_mobil-0.5)/descriptors::invCs2<V,DESCRIPTOR>();

    // Computation and storage of forces
    Vector<V,DESCRIPTOR::d> forceNS{};
    V u[DESCRIPTOR::d] {};
    cellNS.computeU(u);
    V k = 1.5*sigma*w;
    V beta = 12*sigma/w;
    V mu = 4*beta*phi*(phi-1.)*(phi-0.5)-k*laplacePhi;
    forceNS += mu*gradPhi;
    auto externalBlockForce = cellNS.template getField<descriptors::EXTERNAL_FORCE>();
    cellNS.template setField<descriptors::FORCE>(externalBlockForce + forceNS/rho);
    cellNS.template setField<descriptors::NABLARHO>((rho_l-rho_v)*gradPhi);
    cellNS.template setField<descriptors::TAU_EFF>(tau);
    cellNS.template setField<descriptors::RHO>(rho);

    V psi = cellAC.template getField<descriptors::PSI>();
    V interfaceIndicator = 0.;
    if (fabs(psi) < 3.*w) interfaceIndicator = 1;
    Vector<V,DESCRIPTOR::d> forceAC{};
    Vector<V,DESCRIPTOR::d> old_phiU{};
    Vector<V,DESCRIPTOR::d> n{};
    Vector<V,DESCRIPTOR::d> phiU{};
    V gradPhiSqr = 0;
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      gradPhiSqr += gradPhi[iD]*gradPhi[iD];
      phiU[iD] = phi*u[iD];
    }
    if (gradPhiSqr >= 1e-28) {
      n += gradPhi/util::sqrt(gradPhiSqr);
    }
    old_phiU = cellAC.template getField<descriptors::OLD_PHIU>();
    cellAC.template setField<descriptors::OLD_PHIU>(phiU);
    V lambda = 4*phi*(1.-phi)/w;
    V D_N = 4*beta/k*(interfaceIndicator-1.)*(phi*(phi-1.)*(phi-0.5)-gamma*(1.-phi)*phi);
    V source_old = cellAC.template getField<descriptors::SOURCE_OLD>();
    cellAC.template setField<descriptors::SOURCE_OLD>(M*D_N);
    forceAC += (phiU - old_phiU) + interfaceIndicator*lambda*n/descriptors::invCs2<V,DESCRIPTOR>();
    cellAC.template setField<descriptors::FORCE>(forceAC);
    V source = 1.5*M*D_N-0.5*source_old;
    cellAC.template setField<descriptors::SOURCE>(source);
    cellAC.template setField<descriptors::VELOCITY>(u);
  }
};

struct LiangPostProcessor {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  struct SIGMA       : public descriptors::FIELD_BASE<1> { };
  struct W           : public descriptors::FIELD_BASE<1> { };
  struct TAUS        : public descriptors::FIELD_BASE<2> { };
  struct RHOS        : public descriptors::FIELD_BASE<2> { };
  using parameters = meta::list<SIGMA,W,TAUS,RHOS>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELL::template value_t<names::NavierStokes>::descriptor_t;

    auto& cellNS = cells.template get<names::NavierStokes>();
    auto& cellAC = cells.template get<names::Component1>();

    V phi = cellAC.template getFieldComponent<descriptors::STATISTIC>(0);
    Vector<V,DESCRIPTOR::d> gradPhi{};
    V laplacePhi = 0;
    //TODO: use lattice gradient schemes here from Cahn-Hilliard implementation
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      const V phi_i = cellAC.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getFieldComponent<descriptors::STATISTIC>(0);
      gradPhi += phi_i * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
      laplacePhi += 2*(phi_i - phi) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
    }

    auto sigma = parameters.template get<SIGMA>();
    auto w = parameters.template get<W>();
    auto rho_v = parameters.template get<RHOS>()[0];
    auto rho_l = parameters.template get<RHOS>()[1];
    V rho = rho_v + (rho_l-rho_v)*phi;
    auto tau_v = parameters.template get<TAUS>()[0];
    auto tau_l = parameters.template get<TAUS>()[1];
    V tau = tau_v + (tau_l-tau_v)*phi;

    // Computation and storage of forces
    Vector<V,DESCRIPTOR::d> forceNS{};
    V u[DESCRIPTOR::d] {};
    cellNS.computeU(u);
    V k = 1.5*sigma*w;
    V beta = 12*sigma/w;
    V mu = 4*beta*phi*(phi-1.)*(phi-0.5)-k*laplacePhi;
    forceNS += mu*gradPhi;
    auto externalBlockForce = cellNS.template getField<descriptors::EXTERNAL_FORCE>();
    cellNS.template setField<descriptors::FORCE>(externalBlockForce*(rho-rho_l)/rho + forceNS/rho);
    cellNS.template setField<descriptors::NABLARHO>((rho_l-rho_v)*gradPhi);
    cellNS.template setField<descriptors::TAU_EFF>(tau);
    cellNS.template setField<descriptors::RHO>(rho);

    Vector<V,DESCRIPTOR::d> forceAC{};
    Vector<V,DESCRIPTOR::d> old_phiU{};
    Vector<V,DESCRIPTOR::d> n{};
    Vector<V,DESCRIPTOR::d> phiU{};
    V gradPhiSqr = 0;
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      gradPhiSqr += gradPhi[iD]*gradPhi[iD];
      phiU[iD] = phi*u[iD];
    }
    if (gradPhiSqr >= 1e-28) {
      n += gradPhi/util::sqrt(gradPhiSqr);
    }
    old_phiU = cellAC.template getField<descriptors::OLD_PHIU>();
    cellAC.template setField<descriptors::OLD_PHIU>(phiU);
    V lambda = 4*phi*(1.-phi)/w;
    forceAC += (phiU - old_phiU) + lambda*n/descriptors::invCs2<V,DESCRIPTOR>();
    cellAC.template setField<descriptors::FORCE>(forceAC);
    cellAC.template setField<descriptors::VELOCITY>(u);
  }
};

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal>
struct GeometricPhaseFieldWallProcessor2D {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::THETA>;

  int getPriority() const {
    return 1;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    auto theta = parameters.template get<descriptors::THETA>();
    Vector<int,DESCRIPTOR::d> tangent{yNormal*(-1),xNormal*1};
    Vector<int,DESCRIPTOR::d> opp_tang{tangent[0]*(-1),tangent[1]*(-1)};

    auto phi_1 = cell.neighbor({-xNormal,-yNormal}).template getFieldComponent<descriptors::STATISTIC>(0);
    auto phi_1r = cell.neighbor({-xNormal,-yNormal}).neighbor(tangent).template getFieldComponent<descriptors::STATISTIC>(0);
    auto phi_1l = cell.neighbor({-xNormal,-yNormal}).neighbor(opp_tang).template getFieldComponent<descriptors::STATISTIC>(0);
    auto phi_2r = cell.neighbor({-2*xNormal,-2*yNormal}).neighbor(tangent).template getFieldComponent<descriptors::STATISTIC>(0);
    auto phi_2l = cell.neighbor({-2*xNormal,-2*yNormal}).neighbor(opp_tang).template getFieldComponent<descriptors::STATISTIC>(0);

    T dphi_1 = ( phi_1r - phi_1l ) / 2.;
    T dphi_2 = ( phi_2r - phi_2l ) / 2.;
    T tau_x_dphi = 1.5*dphi_1 - 0.5*dphi_2;

    auto phi = cell.template getField<descriptors::STATISTIC>();
    phi[0] = phi_1 + tan( M_PI/2. - theta ) * abs(tau_x_dphi);
    cell.template setField<descriptors::STATISTIC>(phi);
  }
};


};

#endif
