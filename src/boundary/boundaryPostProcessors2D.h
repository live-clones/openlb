/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  Generic version of the collision, which modifies the particle
 *  distribution functions, by Orestis Malaspinas.
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

#ifndef FD_BOUNDARIES_2D_H
#define FD_BOUNDARIES_2D_H

#include "core/postProcessing.h"

#include "core/operator.h"

namespace olb {

/**
* This class computes the skordos BC
* on a flat wall in 2D but with a limited number of terms added to the
* equilibrium distributions (i.e. only the Q_i : Pi term)
*/
template <typename T, typename DESCRIPTOR, int direction, int orientation>
class StraightFdBoundaryProcessor2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const { return 0; }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform;

private:
  template <int deriveDirection, typename CELL>
  void interpolateGradients(CELL& blockLattice, T velDeriv[DESCRIPTOR::d]) const any_platform;
};

/// PostProcessors for the chemical potential boundary condition in the free energy model.
/// The chemical potentials on the boundary are set equal to the chemical potential on the
/// fluid cell normal to the boundary. This is necessary because the coupling between the
/// lattices requires the calculation of the gradient of the chemical potential.
///
/// It would be preferable if these were implemented as a lattice coupling that ran
/// between the chemical potential and force lattice couplings. However there is no
/// access to the discrete normals in lattice couplings.
template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y>
struct FreeEnergyInletMomentum2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform {

    cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y}).template getField<descriptors::CHEM_POTENTIAL>());

    T rho0 = cell.computeRho();
    T rho1 = cell.neighbor({-NORMAL_X,-NORMAL_Y}).computeRho();

    cell.template setField<descriptors::CHEM_POTENTIAL>(
      cell.template getField<descriptors::CHEM_POTENTIAL>() + (rho1 / rho0 - 1) / descriptors::invCs2<T,DESCRIPTOR>());

  }

};

template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y>
struct FreeEnergyInletOrderParameter2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform {

    cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y}).template getField<descriptors::CHEM_POTENTIAL>());

  }

};

template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y>
struct FreeEnergyOutletMomentum2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform {

    T rhoBoundaryNew, rhoBoundaryOld, rhoBulk, u[2];

    rhoBoundaryOld = cell.computeRho();

    cell.neighbor({-NORMAL_X,-NORMAL_Y}).computeRhoU(rhoBulk, u);

    T uPerp = 0;

    Vector<T,2> normalVec({NORMAL_X,NORMAL_Y});

    if (normalVec[0] == 0) {
      uPerp = normalVec[1] * u[1];
    } else if (normalVec[1] == 0) {
            uPerp = normalVec[0] * u[0];
    }

    rhoBoundaryNew = (rhoBoundaryOld + uPerp * rhoBulk) / (1. + uPerp);
    cell.defineRho(rhoBoundaryNew);

    cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y}).template getField<descriptors::CHEM_POTENTIAL>());

    cell.template setField<descriptors::CHEM_POTENTIAL>(
      cell.template getField<descriptors::CHEM_POTENTIAL>() + (rhoBulk / rhoBoundaryNew - 1) / descriptors::invCs2<T,DESCRIPTOR>());

  }

};

template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y>
struct FreeEnergyOutletOrderParameter2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform {

    T rhoBoundaryNew, rhoBoundaryOld, rhoBulk, u[2];

    rhoBoundaryOld = cell.computeRho();

    cell.neighbor({-NORMAL_X,-NORMAL_Y}).computeRhoU(rhoBulk, u);

    T uPerp = 0;

    Vector<T,2> normalVec({NORMAL_X,NORMAL_Y});

    if (normalVec[0] == 0) {
      uPerp = normalVec[1] * u[1];
    } else if (normalVec[1] == 0) {
            uPerp = normalVec[0] * u[0];
    }

    rhoBoundaryNew = (rhoBoundaryOld + uPerp * rhoBulk) / (1. + uPerp);
    cell.defineRho(rhoBoundaryNew);

    cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y}).template getField<descriptors::CHEM_POTENTIAL>());

  }

};


/// PostProcessor for the wetting boundary condition in the free energy model. This is
/// required to set rho on the boundary (using the denisty of the neighbouring cell in
/// direction of inwards facing normal at the boundary), as the coupling between the
/// lattices requires the calculation of a density gradient.
template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y>
struct FreeEnergyWallMomentumProcessor2D {

public:
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<olb::descriptors::ADDEND>;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform{

    auto addend = parameters.template get<descriptors::ADDEND>();

    T rhoBulk = cell.neighbor({-NORMAL_X,-NORMAL_Y}).computeRho();
    T rhoTmp = 0.;

    for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
      rhoTmp += cell[iPop];
    }

    T rhoBoundary = rhoBulk + addend;
    rhoBoundary -= rhoTmp;

    cell[0] = rhoBoundary - 1.;

    cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y}).template getField<descriptors::CHEM_POTENTIAL>());

    cell.template setField<descriptors::CHEM_POTENTIAL>(
      cell.template getField<descriptors::CHEM_POTENTIAL>() + (rhoBulk / rhoBoundary - 1) / descriptors::invCs2<T,DESCRIPTOR>());

  }

};

template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y>
struct FreeEnergyWallOrderParameterProcessor2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<olb::descriptors::ADDEND>;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform{

    auto addend = parameters.template get<descriptors::ADDEND>();

    T rhoBulk = cell.neighbor({-NORMAL_X,-NORMAL_Y}).computeRho();
    T rhoTmp = 0.;

    for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
      rhoTmp += cell[iPop];
    }

    T rhoBoundary = rhoBulk + addend;
    rhoBoundary -= rhoTmp;

    cell[0] = rhoBoundary - 1.;

    cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y}).template getField<descriptors::CHEM_POTENTIAL>());

  }

};

/// PostProcessor for pressure / velocity outflow boundaries in the free energy model.
/// The density / order parameters are prescribed to the outflow nodes such that they
/// obey the local-velocity convective boundary condition given in Lou, Gou, Shi (2013).
template <typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y>
class FreeEnergyConvectiveProcessor2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform;
};

/**
* This class computes a convection BC on a flat wall in 2D
*/
template <typename T, typename DESCRIPTOR, int direction, int orientation>
class StraightConvectionBoundaryProcessor2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  struct PREV_CELL
      : public descriptors::FIELD_BASE<
            util::populationsContributingToVelocity<DESCRIPTOR, direction,
                                                    -orientation>()
                .size()> {};

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void initialize(CELL& cell) any_platform;

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform;
};

/**
* This class computes a slip BC in 2D
*/

template <typename T, typename DESCRIPTOR>
class SlipBoundaryProcessor2D : public LocalPostProcessor2D<T, DESCRIPTOR> {
public:
  SlipBoundaryProcessor2D(int x0_, int x1_, int y0_, int y1_,
                          int discreteNormalX_, int discreteNormalY_);
  int  extent() const override { return 0; }
  int  extent(int whichDirection) const override { return 0; }
  void process(BlockLattice<T, DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice<T, DESCRIPTOR>& blockLattice, int x0_,
                        int x1_, int y0_, int y1_) override;

private:
  int reflectionPop[DESCRIPTOR::q];
  int x0, x1, y0, y1;
};

template <typename T, typename DESCRIPTOR>
class SlipBoundaryProcessorGenerator2D
    : public PostProcessorGenerator2D<T, DESCRIPTOR> {
public:
  SlipBoundaryProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_,
                                   int discreteNormalX_, int discreteNormalY_);
  PostProcessor2D<T, DESCRIPTOR>*          generate() const override;
  PostProcessorGenerator2D<T, DESCRIPTOR>* clone() const override;

private:
  int discreteNormalX;
  int discreteNormalY;
};


/**
* This class computes the skordos BC in 2D on a convex
* corner but with a limited number of terms added to the
* equilibrium distributions (i.e. only the Q_i : Pi term)
*/
template <typename T, typename DESCRIPTOR, int xNormal, int yNormal>
class OuterVelocityCornerProcessor2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 1;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform;
};

} // namespace olb

#endif
