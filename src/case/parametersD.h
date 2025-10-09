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
    std::function<void()> calculate;

    ~TypeErasedFieldD() {
      if (field && deleter) {
        deleter();
      }
    }
  };

  std::unordered_map<std::type_index, TypeErasedFieldD> _map;
  std::optional<CLIreader> _args;

  template <concepts::Field FIELD>
  void setupTypeErasedField(TypeErasedFieldD& typeErasedField) {
    typeErasedField.name = parameters::name<FIELD>();
    typeErasedField.deleter = [&typeErasedField]() {
      delete static_cast<FieldD<T,DESCRIPTOR,FIELD>*>(typeErasedField.field);
    };
    typeErasedField.value = [this]() -> std::string {
      std::stringstream out;
      out << *static_cast<FieldD<T,DESCRIPTOR,FIELD>*>(this->getFieldPtr<FIELD>());
      return out.str();
    };
    typeErasedField.set = [this](std::string str) {
      if (auto v = FieldD<T,DESCRIPTOR,FIELD>::fromString(str)) {
        this->set<FIELD>(*v);
      }
    };
  }

  void tryUpdateFromCLI(TypeErasedFieldD& typeErasedField) {
    if (_args) {
      const std::string key = "--" + typeErasedField.name;
      if (_args->contains(key)) {
        typeErasedField.set(_args->getValueOrFallback<std::string>(key, "[]"));
      }
    }
  }

  // Helper to get the typed field pointer, executing calculate if needed
  template <concepts::Field FIELD>
  void* getFieldPtr() {
    TypeErasedFieldD& typeErasedField = _map[typeid(FIELD)];
    // If not initialized at all, use default value
    if (!typeErasedField.field && !typeErasedField.calculate) {
      set<FIELD>(FIELD::template getInitialValue<T,DESCRIPTOR>());
      tryUpdateFromCLI(typeErasedField); // Allow CLI override of default
    }
    else if (typeErasedField.calculate) {
      typeErasedField.calculate();
      typeErasedField.calculate = nullptr; // Ensure calculate runs only once
    }
    return typeErasedField.field;
  }

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  ParametersD() = default;

  // Set parameter with a concrete value
  template <concepts::Field FIELD>
  void set(FieldD<T,DESCRIPTOR,FIELD> value) {
    TypeErasedFieldD& typeErasedField = _map[typeid(FIELD)];
    if (typeErasedField.field && typeErasedField.deleter) {
      typeErasedField.deleter();
    }
    typeErasedField.calculate = nullptr;
    typeErasedField.field = new FieldD<T,DESCRIPTOR,FIELD>{value};
    setupTypeErasedField<FIELD>(typeErasedField);
  }

  //template <concepts::Field FIELD>
  //requires (   DESCRIPTOR::template size<FIELD>() == 1
  //          && std::is_same_v<typename FIELD::value_t,std::string>)
  //void set(char* text) {
  //  set<FIELD>(std::string(text));
  //}

  // Set dependent parameter with a lambda for lazy evaluation
  template <concepts::Field FIELD>
  void set(concepts::CallableReturning<FieldD<T,DESCRIPTOR,FIELD>> auto&& f) {
    TypeErasedFieldD& typeErasedField = _map[typeid(FIELD)];
    if (typeErasedField.field && typeErasedField.deleter) {
      typeErasedField.deleter();
      typeErasedField.field = nullptr;
    }
    typeErasedField.calculate = [this,f]() {
      this->set<FIELD>(f());
    };
    setupTypeErasedField<FIELD>(typeErasedField);
  }


  template <concepts::Field FIELD>
  auto get() {
    void* ptr = getFieldPtr<FIELD>();
    if constexpr (DESCRIPTOR::template size<FIELD>() == 1) {
      return (*static_cast<FieldD<T,DESCRIPTOR,FIELD>*>(ptr))[0];
    } else {
      return *static_cast<FieldD<T,DESCRIPTOR,FIELD>*>(ptr);
    }
  }

  void print(OstreamManager clout = {std::cout, "ParametersD"}) {
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
