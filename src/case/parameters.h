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
#include "core/meta.h"

namespace olb {

namespace parameters {

struct OVERLAP : public descriptors::TYPED_FIELD_BASE<unsigned,1> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,1>{3};
  }
};

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

/// Returns name of PARAMETER for human consumption
template <typename PARAMETER>
std::string name() {
  auto raw = meta::nice_name<PARAMETER>();
  if (raw.starts_with("olb::")) {
    raw = std::string_view(raw.cbegin() + 5,
                           raw.cend());
  }
  if (raw.starts_with("parameters::")) {
    raw = std::string_view(raw.cbegin() + 12,
                           raw.cend());
  }
  return std::string(raw);
}

}

}

#endif
