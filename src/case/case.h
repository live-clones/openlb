/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Adrian Kummerlaender, Shota Ito, Julius Jeßberger, Mathias J. Krause
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

#ifndef CASE_CASE_H
#define CASE_CASE_H

#include "mesh.h"
#include "parameters.h"

#include "utilities/typeMap.h"
#include "utilities/typeIndexedContainers.h"

namespace olb {

template <typename T, typename DESCRIPTOR>
struct Lattice {
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  template <typename VALUE_TYPE>
  using exchange_value_t = Lattice<VALUE_TYPE, DESCRIPTOR>;
};

template <typename MAP>
class PlainCase {
public:
  // Select first descriptor as reference for value type and dimension (for now)
  using reference_discretization_t = MAP::values_t::template get<0>;
  
  // Expose reference value type to use in apps
  using value_t = reference_discretization_t::value_t;

  template <typename VALUE_TYPE>
  using exchange_value_t = PlainCase<typename MAP::template map_values<
    meta::exchange_value_type<VALUE_TYPE>::template type>
  >;

  using ParametersD = olb::ParametersD<typename reference_discretization_t::value_t,
                                       typename reference_discretization_t::descriptor_t>;

private:
  /// Reference to mesh (shared between multiple cases)
  Mesh<typename reference_discretization_t::value_t,
       reference_discretization_t::descriptor_t::d>& _mesh;

  ParametersD& _parameters;

  /// Case-specific geometry
  std::unique_ptr<SuperGeometry<typename reference_discretization_t::value_t,
                                reference_discretization_t::descriptor_t::d>> _geometry;

  template <typename DISCRETIZATION>
  using ptr_to_lattice = std::unique_ptr<SuperLattice<typename DISCRETIZATION::value_t,
                                                      typename DISCRETIZATION::descriptor_t>>;

  /// Case-specific lattices
  utilities::TypeIndexedTuple<typename MAP::template map_values<
    ptr_to_lattice
  >> _lattices;

  std::unordered_map<std::string, std::unique_ptr<ApplicableO>> _operators;

public:
  static constexpr unsigned d = reference_discretization_t::descriptor_t::d;

  template <typename NAME>
  using value_t_of = MAP::template value<NAME>::value_t;
  template <typename NAME>
  using descriptor_t = MAP::template value<NAME>::descriptor_t;

  void constructLatticesFromMesh() {
    MAP::keys_t::for_each([&](auto name) {
      using DISCRETIZATION = MAP::template value<typename decltype(name)::type>;
      using V = DISCRETIZATION::value_t;
      using DESCRIPTOR = DISCRETIZATION::descriptor_t;
      _lattices.template set(name, std::make_unique<SuperLattice<V,DESCRIPTOR>>(
          _mesh.getCuboidDecomposition()
        , _mesh.getLoadBalancer()));
    });
  }

  PlainCase(ParametersD& parameters,
            Mesh<typename reference_discretization_t::value_t,d>& mesh)
    : _parameters{parameters}
    , _mesh{mesh}
    , _geometry{new SuperGeometry<typename reference_discretization_t::value_t,
                                  reference_discretization_t::descriptor_t::d>(
          mesh.getCuboidDecomposition()
        , mesh.getLoadBalancer())}
    , _lattices{}
  {
    constructLatticesFromMesh();
  };

  ParametersD& getParameters() {
    return _parameters;
  }

  auto& getMesh() {
    return _mesh;
  }

  auto& getGeometry() {
    return *_geometry;
  }

  template <typename NAME>
  auto& getLattice(NAME) {
    return *_lattices.template get<NAME>();
  }

  ApplicableO& getOperator(std::string name) {
    return *_operators.at(name);
  }

  template <typename... ARGS>
  auto& setCouplingOperator(std::string name, ARGS&&... args) {
    auto ptr = constructUniqueCoupling(std::forward<ARGS&&>(args)...);
    auto& concrete = *ptr;
    _operators[name].reset(ptr.release());
    return concrete;
  }

};

template <typename... MAP>
using Case = PlainCase<meta::map<MAP...>>;

}

#endif
