/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Adrian Kummerlaender
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

#ifndef CASE_PARAMETERS_H
#define CASE_PARAMETERS_H

#include "io/cliReader.h"

#include <functional>

namespace olb {

namespace parameters {

// Converter-related parameters
struct RESOLUTION : public descriptors::TYPED_FIELD_BASE<int,1> { };
struct PHYS_DELTA_X : public descriptors::FIELD_BASE<1> { };
struct PHYS_DELTA_T : public descriptors::FIELD_BASE<1> { };
struct PHYS_CHAR_VELOCITY : public descriptors::FIELD_BASE<1> { };
struct PHYS_CHAR_VISCOSITY : public descriptors::FIELD_BASE<1> { };
struct PHYS_CHAR_DENSITY : public descriptors::FIELD_BASE<1> { };

struct DOMAIN_EXTENT : public descriptors::FIELD_BASE<0,1> { };
struct RAYLEIGH : public descriptors::FIELD_BASE<1> { };
struct PRANDTL : public descriptors::FIELD_BASE<1> { };
struct MAX_PHYS_T : public descriptors::FIELD_BASE<1> { };

}

/// On-demand allocating parameter field storage
/**
 * To be used in non-critical sections for convenient storage of parameter values
 **/
template <concepts::BaseType T, concepts::Descriptor DESCRIPTOR>
class ParametersD final {
private:
  struct TypeErasedFieldD {
    void* field = nullptr;

    std::string name;
    std::function<void()> deleter;
    std::function<std::string()> value;
    std::function<void(std::string)> set;

    ~TypeErasedFieldD() {
      if (field && deleter) {
        deleter();
      }
    }
  };
  std::unordered_map<std::type_index, TypeErasedFieldD> _map;

  template <concepts::Field FIELD>
  void set(TypeErasedFieldD& typeErasedField,
           const FieldD<T,DESCRIPTOR,FIELD>& value) {
    typeErasedField.field = new FieldD<T,DESCRIPTOR,FIELD>{value};
    typeErasedField.name = fields::name<FIELD>();
    typeErasedField.deleter = [&typeErasedField]() {
      delete static_cast<FieldD<T,DESCRIPTOR,FIELD>*>(typeErasedField.field);
    };
    typeErasedField.value = [&typeErasedField]() -> std::string {
      std::stringstream out;
      out << *static_cast<FieldD<T,DESCRIPTOR,FIELD>*>(typeErasedField.field);
      return out.str();
    };
    typeErasedField.set = [&typeErasedField](std::string str) {
      if (auto v = FieldD<T,DESCRIPTOR,FIELD>::fromString(str)) {
        *static_cast<FieldD<T,DESCRIPTOR,FIELD>*>(typeErasedField.field) = *v;
      }
    };
  };

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  ParametersD() = default;

  template <concepts::Field FIELD>
  void set(FieldD<T,DESCRIPTOR,FIELD>&& value) {
    TypeErasedFieldD& typeErasedField = _map[typeid(FIELD)];
    if (!typeErasedField.field) {
      set<FIELD>(typeErasedField, value);
    } else {
      *static_cast<FieldD<T,DESCRIPTOR,FIELD>*>(typeErasedField.field) = value;
    }
  };

  template <concepts::Field FIELD>
  auto get() {
    TypeErasedFieldD& typeErasedField = _map[typeid(FIELD)];
    if (!typeErasedField.field) {
      set<FIELD>(typeErasedField, FIELD::template getInitialValue<T,DESCRIPTOR>());
    }
    if constexpr (DESCRIPTOR::template size<FIELD>() == 1) {
      return (*static_cast<FieldD<T,DESCRIPTOR,FIELD>*>(typeErasedField.field))[0];
    } else {
      return *static_cast<FieldD<T,DESCRIPTOR,FIELD>*>(typeErasedField.field);
    }
  }

  void print() {
    OstreamManager clout(std::cout, "ParametersD");
    clout << "-- Simulation Parameters --" << std::endl;
    std::size_t nCharName = 0;
    for (auto& [_, typeErasedField] : _map) {
      nCharName = std::max(nCharName, typeErasedField.name.length());
    }
    for (auto& [_, typeErasedField] : _map) {
      clout << typeErasedField.name << std::string(nCharName - typeErasedField.name.length() + 1, ' ')
            << "= " << typeErasedField.value()
            << std::endl;
    }
    clout << "---------------------------" << std::endl;
  }

  void fromCLI(int& argc, char** argv) {
    CLIreader args(argc, argv);
    if (args.contains("--help")) {
      print();
      throw std::runtime_error("Program aborted on user choice");
    }
    bool updated = false;
    for (auto& [_, typeErasedField] : _map) {
      const std::string key = "--" + typeErasedField.name;
      if (args.contains(key)) {
        typeErasedField.set(args.getValueOrFallback<std::string>(key, "[]"));
        updated = true;
      }
    }
    if (updated) {
      print();
    }
  }

};

}

#endif
