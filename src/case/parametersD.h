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

#ifndef CASE_PARAMETERSD_H
#define CASE_PARAMETERSD_H

#include "io/cliReader.h"
#include "parameters.h"

#include <functional>

namespace olb {

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
    typeErasedField.name = parameters::name<FIELD>();
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

  std::optional<CLIreader> _args;

  void tryUpdateFromCLI(TypeErasedFieldD& typeErasedField) {
    if (typeErasedField.field && _args) {
      const std::string key = "--" + typeErasedField.name;
      if (_args->contains(key)) {
        typeErasedField.set(_args->getValueOrFallback<std::string>(key, "[]"));
      }
    }
  }

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
      // Catch potential updates of deferred parameter accesses
      tryUpdateFromCLI(typeErasedField);
    }
    if constexpr (DESCRIPTOR::template size<FIELD>() == 1) {
      return (*static_cast<FieldD<T,DESCRIPTOR,FIELD>*>(typeErasedField.field))[0];
    } else {
      return *static_cast<FieldD<T,DESCRIPTOR,FIELD>*>(typeErasedField.field);
    }
  }

  void print(OstreamManager clout={std::cout, "ParametersD"}) {
    std::size_t nCharName = 0;
    for (auto& [_, typeErasedField] : _map) {
      nCharName = std::max(nCharName, typeErasedField.name.length());
    }
    for (auto& [_, typeErasedField] : _map) {
      clout << typeErasedField.name << std::string(nCharName - typeErasedField.name.length() + 1, ' ')
            << "= " << typeErasedField.value()
            << std::endl;
    }
  }

  void fromCLI(int& argc, char** argv) {
    if (!_args) {
      _args = CLIreader(argc, argv);
    }
    if (_args->contains("--help")) {
      OstreamManager clout(std::cout, "Help");
      clout << std::endl;
      clout << "-- Parameters --" << std::endl;
      print(clout);
      clout << std::endl;
      clout << "You can change any of these parameters via `./app --NAME VALUE`." << std::endl;
      clout << "The defined values are printed at the start of the simulation." << std::endl;
      clout << std::endl;
      throw std::runtime_error("Program aborted on user choice");
    }
    for (auto& [_, typeErasedField] : _map) {
      tryUpdateFromCLI(typeErasedField);
    }
  }

};

}

#endif
