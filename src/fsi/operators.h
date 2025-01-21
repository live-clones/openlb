/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
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

#ifndef FSI_OPERATORS_H
#define FSI_OPERATORS_H

#include <ranges>

#include "fields.h"
#include "elements/concept.h"

#include "elements/couplingFaceF2D.h"

namespace olb {

/// Operator for discretizing FSI elements into HLBM porosities
template <concepts::PorosityElementF PorosityF>
struct InitializePorosityO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::merge<
    typename PorosityF::parameters,
    meta::list<
      fields::array_of<fields::fsi::ELEMENT_TAG>,
      fields::fsi::ELEMENTS_COUNT
    >
  >;

  int getPriority() const {
    return 0;
  }

  template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    const auto physR = cell.template getField<descriptors::LOCATION>();
    const auto nElements = parameters.template get<fields::fsi::ELEMENTS_COUNT>();
    if (cell.template getField<descriptors::POROSITY>() == 1) {
      for (std::size_t iElement=0; iElement < nElements; ++iElement) {
        const auto porosity = PorosityF().compute(parameters, physR, iElement);
        const auto prevPorosity = cell.template getField<descriptors::POROSITY>();
        if (porosity < prevPorosity) {
          const auto newTag = PorosityF().tag(parameters, physR, iElement);
          cell.template setField<fields::fsi::ELEMENT_TAG>(newTag);
          cell.template setField<descriptors::POROSITY>(porosity);
          const auto u = PorosityF().computeU(parameters, physR, iElement);
          cell.template setField<descriptors::VELOCITY>(u);
        }
      }
      if (cell.template getField<descriptors::POROSITY>() == 1) {
        cell.template setField<fields::fsi::ELEMENT_TAG>(0);
        cell.template setField<descriptors::VELOCITY>(0);
      }
    }
  }
};

/// Operator for updating FSI element discretizations into HLBM porosities
template <concepts::PorosityElementF PorosityF>
struct UpdatePorosityO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::merge<
    typename PorosityF::parameters,
    meta::list<fields::array_of<fields::fsi::ELEMENT_TAG>>
  >;

  int getPriority() const {
    return 0;
  }

  template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    const int tag = cell.template getField<fields::fsi::ELEMENT_TAG>();
    if (tag > 0) {
      const int iElement = tag - 1;
      const auto physR = cell.template getField<descriptors::LOCATION>();

      const auto newTag      = PorosityF().tag(parameters, physR, iElement);
      const auto newPorosity = PorosityF().compute(parameters, physR, iElement);
      const auto newU        = PorosityF().computeU(parameters, physR, iElement);

      const bool isPorous = newPorosity < 1;
      cell.template setField<fields::fsi::ELEMENT_TAG>(isPorous * newTag);
      cell.template setField<descriptors::VELOCITY>(isPorous * newU);

      const bool isFrontier = !PorosityF().isInterior(parameters, physR, iElement);
      cell.template setField<descriptors::POROSITY>(isFrontier * newPorosity);
    }
  }
};

/// Operator for calculating per-cell momentum exchange forces in a HLBM context
struct CollectPorousBoundaryForceO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::list<
    fields::array_of<fields::fsi::ELEMENT_PIVOT>,
    fields::converter::PHYS_LENGTH
  >;

  int getPriority() const {
    return 1;
  }

  template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    FieldD<V,DESCRIPTOR,fields::fsi::BOUNDARY_FORCE> force { };
    FieldD<V,DESCRIPTOR,fields::fsi::BOUNDARY_TORQUE> torque { };
    const V porosity = cell.template getField<descriptors::POROSITY>();
    const int tag = cell.template getField<fields::fsi::ELEMENT_TAG>();
    cell.template setField<fields::fsi::REDUCED_ELEMENT_TAG>(0);

    bool isAtSurface = false;
    if (tag > 0 && porosity < 1) {
      Vector<V,DESCRIPTOR::d> u{};
      cell.computeU(u.data());

      for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
        const auto c = descriptors::c<DESCRIPTOR>(iPop);
        auto neighbor = cell.neighbor(c);

        const bool isNeighborFluid      = neighbor.template getField<descriptors::POROSITY>() == 1;
        //const bool isIncreasingPorosity = neighbor.template getField<descriptors::POROSITY>() >  porosity;

        const V pop1 = neighbor[iPop];
        const V pop2 = cell[descriptors::opposite<DESCRIPTOR>(iPop)];
        force -= isNeighborFluid * (pop1*(c-u) + pop2*(c+u));

        isAtSurface |= isNeighborFluid;
      }

      auto physR = cell.template getField<descriptors::LOCATION>();
      auto pivots = parameters.template get<fields::array_of<fields::fsi::ELEMENT_PIVOT>>();
      const Vector<V,DESCRIPTOR::d> pivot{pivots, static_cast<std::size_t>(tag-1)};
      torque = crossProduct((physR - pivot) * parameters.template get<fields::converter::PHYS_LENGTH>(),
                            force);

      if (isAtSurface) {
        cell.template setField<fields::fsi::REDUCED_ELEMENT_TAG>(tag);
        for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
          const auto c = descriptors::c<DESCRIPTOR>(iPop);
          auto neighbor = cell.neighbor(c);
          if (neighbor.template getField<fields::fsi::ELEMENT_TAG>() == 0) {
            neighbor.template setField<fields::fsi::ELEMENT_TAG>(tag);
          }
        }
      } else {
        force = 0;
        torque = 0;
      }
    }

    cell.template setField<fields::fsi::BOUNDARY_FORCE>(force);
    cell.template setField<fields::fsi::BOUNDARY_TORQUE>(torque);
  }

};

/// Operator for calculating per-cell stress in a HLBM context
struct CollectPorousBoundaryStressO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::list<
    fields::array_of<fields::fsi::ELEMENT_PIVOT>,
    fields::array_of<CouplingFaceF2D::NEIGHBORS>,
    fields::converter::PHYS_LENGTH
  >;

  int getPriority() const {
    return 1;
  }

  template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    FieldD<V,DESCRIPTOR,fields::fsi::BOUNDARY_STRESS> stress { };
    const int tag = cell.template getField<fields::fsi::ELEMENT_TAG>();
    cell.template setField<fields::fsi::REDUCED_ELEMENT_TAG>(0);

    bool isAtSurface = false;
    if (tag > 0 && cell.template getField<descriptors::POROSITY>() < 1) {
      for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
        const auto c = descriptors::c<DESCRIPTOR>(iPop);
        auto neighbor = cell.neighbor(c);
        const bool isNeighborFluid = neighbor.template getField<descriptors::POROSITY>() == 1;
        isAtSurface |= isNeighborFluid;
      }

      if (isAtSurface) {
        cell.computeStress(stress.data());
        //stress = strainRateTensorFDM(cell, n);
        cell.template setField<fields::fsi::REDUCED_ELEMENT_TAG>(tag);

        for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
          const auto c = descriptors::c<DESCRIPTOR>(iPop);
          auto neighbor = cell.neighbor(c);
          if (neighbor.template getField<fields::fsi::ELEMENT_TAG>() == 0) {
            neighbor.template setField<fields::fsi::ELEMENT_TAG>(tag);
          }
        }
      } else {
        stress = 0;
      }
    }

    cell.template setField<fields::fsi::BOUNDARY_STRESS>(stress);
  }

};

/// Operator for allowing migration of padding layer
struct GrowPaddingLayerO {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 1;
  }

  template <concepts::Cell CELL>
  void apply(CELL& cell) any_platform
  {
    using DESCRIPTOR = typename CELL::descriptor_t;
    cell.template setField<fields::fsi::REDUCED_ELEMENT_TAG>(0);
    const int tag = cell.template getField<fields::fsi::ELEMENT_TAG>();
    if (tag > 0 && cell.template getField<descriptors::POROSITY>() < 1) {
      bool isAtSurface = false;

      for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
        const auto c = descriptors::c<DESCRIPTOR>(iPop);
        auto neighbor = cell.neighbor(c);
        const bool isNeighborFluid = neighbor.template getField<descriptors::POROSITY>() == 1;
        isAtSurface |= isNeighborFluid;
      }

      if (isAtSurface) {
        for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
          const auto c = descriptors::c<DESCRIPTOR>(iPop);
          auto neighbor = cell.neighbor(c);
          if (neighbor.template getField<fields::fsi::ELEMENT_TAG>() == 0) {
            neighbor.template setField<fields::fsi::ELEMENT_TAG>(tag);
          }
        }
      }
    }
  }
};

/// Operator for integrating per-element momentum exchange forces in a HLBM-FSI context
struct IntegratePorousBoundaryForceO {
  static constexpr OperatorScope scope = OperatorScope::PerBlock;

  using parameters = meta::list<
    fields::array_of<fields::fsi::ELEMENT_FORCE>,
    fields::array_of<fields::fsi::ELEMENT_TORQUE>,
    fields::array_of<fields::fsi::REDUCED_ELEMENT_TAG>,
    fields::fsi::REDUCED_ELEMENTS_COUNT
  >;

  int getPriority() const {
    return 2;
  }

  template <typename BLOCK>
  struct type {
    void setup(BLOCK& blockLattice) { }
    void apply(BLOCK& blockLattice) {
      throw std::runtime_error("IntegratePorousBoundaryForceO not implemented");
    }
  };

  template <typename BLOCK>
  void setup(BLOCK& blockLattice) {
    type<BLOCK>{}.setup(blockLattice);
  }

  template <typename BLOCK>
  void apply(BLOCK& blockLattice) {
    type<BLOCK>{}.apply(blockLattice);
  }

};

template <typename T, typename DESCRIPTOR, Platform PLATFORM>
requires (isPlatformCPU(PLATFORM))
struct IntegratePorousBoundaryForceO::type<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>> {

void setup(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& blockLattice) {
  blockLattice.template getData<OperatorParameters<IntegratePorousBoundaryForceO>>();
}

void apply(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& blockLattice) {
  auto& parameters = blockLattice.template getData<OperatorParameters<IntegratePorousBoundaryForceO>>().parameters;

  const auto& tag = blockLattice.template getField<fields::fsi::REDUCED_ELEMENT_TAG>();
  const auto& boundaryForce = blockLattice.template getField<fields::fsi::BOUNDARY_FORCE>();
  const auto& boundaryTorque = blockLattice.template getField<fields::fsi::BOUNDARY_TORQUE>();

  std::map<int,std::pair< FieldD<T,DESCRIPTOR,fields::fsi::ELEMENT_FORCE>
                        , FieldD<T,DESCRIPTOR,fields::fsi::ELEMENT_TORQUE>>> data;
  // Sum up per-element boundary forces and torques
  blockLattice.forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> latticeR) {
    const std::size_t iCell = blockLattice.getCellId(latticeR);
    if (tag[0][iCell] > 0) {
      auto& [force, torque] = data[tag[0][iCell]];
      force  += boundaryForce.getRow(iCell);
      torque += boundaryTorque.getRow(iCell);
    }
  });

  auto reducedElementTag = parameters.template get<fields::array_of<fields::fsi::REDUCED_ELEMENT_TAG>>();
  auto elementForce  = parameters.template get<fields::array_of<fields::fsi::ELEMENT_FORCE>>();
  FieldD<T,DESCRIPTOR,fields::array_of<fields::fsi::ELEMENT_TORQUE>> elementTorque
    = parameters.template get<fields::array_of<fields::fsi::ELEMENT_TORQUE>>();

  std::vector<int> tags;
  tags.reserve(data.size());
  for (int tag : std::views::keys(data)) {
    tags.emplace_back(tag);
  }
  std::sort(std::begin(tags), std::end(tags));

  for (std::size_t iReduced=0; iReduced < tags.size(); ++iReduced) {
    reducedElementTag[iReduced] = tags[iReduced];
    const auto& [f, t] = data[tags[iReduced]];
    for (unsigned iD=0; iD < elementForce.d; ++iD) {
      elementForce[iD][iReduced] = f[iD];
    }
    for (unsigned iD=0; iD < elementTorque.d; ++iD) {
      elementTorque[iD][iReduced] = t[iD];
    }
  }

  parameters.template set<fields::fsi::REDUCED_ELEMENTS_COUNT>(tags.size());
}

};

/// Operator for integrating per-element stress in a HLBM-FSI context
struct IntegratePorousBoundaryStressO {
  static constexpr OperatorScope scope = OperatorScope::PerBlock;

  using parameters = meta::list<
    fields::array_of<fields::fsi::ELEMENT_STRESS>,
    fields::array_of<fields::fsi::REDUCED_ELEMENT_TAG>,
    fields::fsi::REDUCED_ELEMENTS_COUNT
  >;

  int getPriority() const {
    return 2;
  }

  template <typename BLOCK>
  struct type {
    void setup(BLOCK& blockLattice) { }
    void apply(BLOCK& blockLattice) {
      throw std::runtime_error("IntegratePorousBoundaryStressO not implemented");
    }
  };

  template <typename BLOCK>
  void setup(BLOCK& blockLattice) {
    type<BLOCK>{}.setup(blockLattice);
  }

  template <typename BLOCK>
  void apply(BLOCK& blockLattice) {
    type<BLOCK>{}.apply(blockLattice);
  }

};

template <typename T, typename DESCRIPTOR, Platform PLATFORM>
requires (isPlatformCPU(PLATFORM))
struct IntegratePorousBoundaryStressO::type<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>> {

void setup(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& blockLattice) {
  blockLattice.template getData<OperatorParameters<IntegratePorousBoundaryStressO>>();
}

void apply(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& blockLattice) {
  auto& parameters = blockLattice.template getData<OperatorParameters<IntegratePorousBoundaryStressO>>().parameters;

  const auto& tag = blockLattice.template getField<fields::fsi::REDUCED_ELEMENT_TAG>();
  const auto& boundaryStress = blockLattice.template getField<fields::fsi::BOUNDARY_STRESS>();

  std::map<int,FieldD<T,DESCRIPTOR,fields::fsi::ELEMENT_STRESS>> data;
  // Sum up per-element boundary stress
  blockLattice.forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> latticeR) {
    const std::size_t iCell = blockLattice.getCellId(latticeR);
    if (tag[0][iCell] > 0) {
      auto& stress = data[tag[0][iCell]];
      stress += boundaryStress.getRow(iCell);
    }
  });

  auto reducedElementTag = parameters.template get<fields::array_of<fields::fsi::REDUCED_ELEMENT_TAG>>();
  auto elementStress = parameters.template get<fields::array_of<fields::fsi::ELEMENT_STRESS>>();

  std::vector<int> tags;
  tags.reserve(data.size());
  for (int tag : std::views::keys(data)) {
    tags.emplace_back(tag);
  }
  std::sort(std::begin(tags), std::end(tags));

  for (std::size_t iReduced=0; iReduced < tags.size(); ++iReduced) {
    reducedElementTag[iReduced] = tags[iReduced];
    const auto& f = data[tags[iReduced]];
    for (unsigned iD=0; iD < elementStress.d; ++iD) {
      elementStress[iD][iReduced] = f[iD];
    }
  }

  parameters.template set<fields::fsi::REDUCED_ELEMENTS_COUNT>(tags.size());
}

};


}

#endif
