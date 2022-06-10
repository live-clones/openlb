/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Adrian Kummerlaender
 *
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

#ifndef SUPER_LATTICE_COUPLING_H
#define SUPER_LATTICE_COUPLING_H

#include "operator.h"
#include "superLattice.h"

#include "utilities/typeMap.h"
#include "solver/names.h"

namespace olb {


template <typename COUPLER, typename COUPLEES, Platform PLATFORM>
class ConcreteBlockCouplingO<COUPLEES,PLATFORM,COUPLER,OperatorScope::PerCell>
  final : public AbstractCouplingO<COUPLEES> {
private:
  template <typename VALUED_DESCRIPTOR>
  using ptr_to_lattice = ConcreteBlockLattice<typename VALUED_DESCRIPTOR::value_t,
                                              typename VALUED_DESCRIPTOR::descriptor_t,
                                              PLATFORM>*;

  utilities::TypeIndexedTuple<typename COUPLEES::template map_values<
    ptr_to_lattice
  >> _lattice;

  typename AbstractCouplingO<COUPLEES>::ParametersD _parameters;

  void execute(typename AbstractCouplingO<COUPLEES>::LatticeR latticeR)
  {
    CellID iCell = std::get<0>(_lattice.tuple)->getCellId(latticeR);
    auto cells = _lattice.exchange_values([&](auto name) -> auto {
      return cpu::Cell{*_lattice.get(name), iCell};
    });
    COUPLER().apply(cells);
  }

public:
  template <typename LATTICE>
  ConcreteBlockCouplingO(LATTICE&& lattice):
    _lattice{lattice}
  { }

  std::type_index id() const override {
    return typeid(COUPLER);
  }

  typename AbstractCouplingO<COUPLEES>::AbstractParameters& getParameters() override {
    return _parameters;
  }

  void execute() override
  {
    using loc = typename AbstractCouplingO<COUPLEES>::LatticeR::value_t;
    auto* lattice = std::get<0>(_lattice.tuple);
    #ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for schedule(static) collapse(1)
    #endif
    for (loc iX=0; iX < lattice->getNx(); ++iX) {
      for (loc iY=0; iY < lattice->getNy(); ++iY) {
        if constexpr (AbstractCouplingO<COUPLEES>::descriptor_t::d == 3) {
          for (loc iZ=0; iZ < lattice->getNz(); ++iZ) {
            execute({iX,iY,iZ});
          }
        } else {
          execute({iX,iY});
        }
      }
    }
  }

};

template <typename COUPLER, typename COUPLEES, Platform PLATFORM>
class ConcreteBlockCouplingO<COUPLEES,PLATFORM,COUPLER,OperatorScope::PerCellWithParameters>
  final : public AbstractCouplingO<COUPLEES> {
private:
  template <typename VALUED_DESCRIPTOR>
  using ptr_to_lattice = ConcreteBlockLattice<typename VALUED_DESCRIPTOR::value_t,
                                              typename VALUED_DESCRIPTOR::descriptor_t,
                                              PLATFORM>*;

  utilities::TypeIndexedTuple<typename COUPLEES::template map_values<
    ptr_to_lattice
  >> _lattice;

  typename COUPLER::parameters::template decompose_into<
    AbstractCouplingO<COUPLEES>::ParametersD::template include_fields
  > _parameters;

  void execute(typename AbstractCouplingO<COUPLEES>::LatticeR latticeR)
  {
    CellID iCell = std::get<0>(_lattice.tuple)->getCellId(latticeR);
    auto cells = _lattice.exchange_values([&](auto name) -> auto {
      return cpu::Cell{*_lattice.get(name), iCell};
    });
    COUPLER().apply(cells, _parameters);
  }

public:
  template <typename LATTICE>
  ConcreteBlockCouplingO(LATTICE&& lattice):
    _lattice{lattice}
  { }

  std::type_index id() const override {
    return typeid(COUPLER);
  }

  typename AbstractCouplingO<COUPLEES>::AbstractParameters& getParameters() override {
    return _parameters;
  }

  void execute() override
  {
    using loc = typename AbstractCouplingO<COUPLEES>::LatticeR::value_t;
    auto* lattice = std::get<0>(_lattice.tuple);
    #ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for schedule(static) collapse(1)
    #endif
    for (loc iX=0; iX < lattice->getNx(); ++iX) {
      for (loc iY=0; iY < lattice->getNy(); ++iY) {
        if constexpr (AbstractCouplingO<COUPLEES>::descriptor_t::d == 3) {
          for (loc iZ=0; iZ < lattice->getNz(); ++iZ) {
            execute({iX,iY,iZ});
          }
        } else {
          execute({iX,iY});
        }
      }
    }
  }

};

template <typename COUPLER, typename COUPLEES>
class SuperLatticeCoupling {
private:
  template <typename VALUED_DESCRIPTOR>
  using ptr_to_lattice = SuperLattice<typename VALUED_DESCRIPTOR::value_t,
                                      typename VALUED_DESCRIPTOR::descriptor_t>*;

  utilities::TypeIndexedTuple<typename COUPLEES::template map_values<
    ptr_to_lattice
  >> _lattice;

  std::vector<std::unique_ptr<AbstractCouplingO<COUPLEES>>> _block;

  template <Platform PLATFORM>
  auto constructConcreteBlockCoupling(int iC)
  {
    auto block = _lattice.exchange_values([&](auto name) -> auto {
      using NAME = typename decltype(name)::type;
      using T = typename COUPLEES::template value<NAME>::value_t;
      using DESCRIPTOR = typename COUPLEES::template value<NAME>::descriptor_t;
      return dynamic_cast<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>*>(
        &_lattice.get(name)->getBlock(iC));
    });
    return std::make_unique<ConcreteBlockCouplingO<COUPLEES,PLATFORM,COUPLER,COUPLER::scope>>(block);
  }

public:
  template <typename... MAP>
  SuperLatticeCoupling(COUPLER,
                       MAP&&... args) {
    auto map = std::make_tuple(&args...);
    COUPLEES::keys_t::for_each([&](auto id) {
      using name_t = typename decltype(id)::type;
      constexpr unsigned idx = COUPLEES::keys_t::template index<name_t>();
      _lattice.template set<name_t>(std::get<2*idx+1>(map));
    });

    auto& load = std::get<0>(_lattice.tuple)->getLoadBalancer();
    for (int iC = 0; iC < load.size(); ++iC) {
      Platform reference = std::get<0>(_lattice.tuple)->getBlock(iC).getPlatform();
      _lattice.for_each([&](auto name, auto lattice) {
        if (lattice->getBlock(iC).getPlatform() != reference) {
          throw std::runtime_error("Platforms of coupled block lattices must match");
        }
      });
      callUsingConcretePlatform(reference,
                                [&](auto platform) {
        _block.emplace_back(constructConcreteBlockCoupling<platform.value>(iC));
      });
    }
  }

  void execute()
  {
    auto& load = std::get<0>(_lattice.tuple)->getLoadBalancer();
    for (int iC = 0; iC < load.size(); ++iC) {
      _block[iC]->execute();
    }
  }

  template <typename FIELD>
  void setParameter(typename AbstractCouplingO<COUPLEES>::template FieldD<FIELD>&& field)
  {
    auto& load = std::get<0>(_lattice.tuple)->getLoadBalancer();
    for (int iC = 0; iC < load.size(); ++iC) {
      _block[iC]->getParameters().template set<FIELD>(std::forward<decltype(field)>(field));
    }
  }

};

template <typename COUPLER, typename... MAP>
SuperLatticeCoupling(COUPLER, MAP&&...) -> SuperLatticeCoupling<
  COUPLER,
  typename meta::map<MAP...>::template map_values<descriptors::extract_valued_descriptor_t>
>;


}

#endif
