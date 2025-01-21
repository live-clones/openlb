/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender
 *                2024 Dennis Teutscher
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

#ifndef DESCRIPTOR_FIELDS_H
#define DESCRIPTOR_FIELDS_H

#include <type_traits>
#include <stdexcept>
#include <optional>
#include "core/meta.h"
#include "core/vector.h"
#include "utilities/aliases.h"

#include "core/platform/column.h"

namespace olb {

/// Forward declaration of Expr type used for code generation
class Expr;

namespace descriptors {

// *INDENT-OFF*

/// \defgroup descriptor
//@{

/// Base of a field whose size is defined by [C,U_1,...,U_N]^T * [1,V_1,...V_N]
template <unsigned C, unsigned... U>
struct FIELD_BASE {
  FIELD_BASE() = default;

  /// Return value type of field
  /**
   * Most fields are stored using the same value type as the T type
   * parameter of their associated lattice. However this template
   * offers the possibility of declaring a different value type per
   * field. See TYPED_FIELD_BASE.
   **/
  template <typename T>
  using value_type = T;

  template <typename T>
  using column_type = AbstractColumn<T>;

  /// Get size of field for parameter vector V (strict)
  template <std::size_t... V>
  static constexpr std::size_t size(std::index_sequence<V...>)
  {
    static_assert(sizeof...(U) == sizeof...(V), "V size fits U params");
    return ((U*V) + ... + C);
  }

  /// Get size of field for parameter vector V (zero-padds or shortens V as required)
  template <std::size_t... V>
  static constexpr std::size_t size()
  {
    if constexpr (sizeof...(U) < sizeof...(V)) {
      return size(meta::take_n_sequence<sizeof...(U)>(std::index_sequence<V...>()));
    } else if constexpr (sizeof...(U) > sizeof...(V)) {
      using namespace meta;
      static_assert(is_zero_sequence(drop_n_sequence<sizeof...(V)>(std::index_sequence<U...>())),
                    "Dropped U entries are zero");
      return size(std::index_sequence<V...>() + zero_sequence<sizeof...(U) - sizeof...(V)>());
    } else {
      return size(std::index_sequence<V...>());
    }
    __builtin_unreachable();
  }

  // Initial value used for allocation of field data in FieldArrayD (optional)
  /**
   * Return value must be a correctly sized and typed olb::Vector.
   **/
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, DESCRIPTOR::template size<FIELD_BASE<C,U...>>()>{};
  }

  template <typename T, typename DESCRIPTOR, typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return true;
  }

  static constexpr bool isSerializable() {
    return true;
  }
};

struct SCALAR_FIELD     : public FIELD_BASE<1> { };
struct SPATIAL_FIELD    : public FIELD_BASE<0,1> { };
struct POPULATION_FIELD : public FIELD_BASE<0,0,1> { };

struct FIELD_BASE_CUSTOM_SIZE {
  template <typename T>
  using value_type = T;

  template <typename T>
  using column_type = AbstractColumn<T>;

  static constexpr bool isSerializable() {
    return true;
  }
};

struct NEIGHBOR_FIELD : public FIELD_BASE_CUSTOM_SIZE {
  template <unsigned D, unsigned...>
  static constexpr unsigned size() {
           if constexpr (D == 1) {
      return 2;
    } else if constexpr (D == 2) {
      return 8;
    } else if constexpr (D == 3) {
      return 26;
    }
  }

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, DESCRIPTOR::template size<NEIGHBOR_FIELD>()>{};
  }
};


/// Base of a descriptor field of scalar TYPE with dimensions A*B + B*Q + C
template <typename TYPE, unsigned C, unsigned... U>
struct TYPED_FIELD_BASE : public FIELD_BASE<C,U...> {
  template <typename T>
  // TODO Find better solution for this hacky workaround to instantiate the entire lattice with Expr
  using value_type = std::conditional_t<std::is_same_v<T,Expr>,
                                        std::conditional_t<std::is_pointer_v<TYPE> || std::is_enum_v<TYPE>,
                                                           TYPE,
                                                           Expr>,
                                        TYPE>;

  template <typename T>
  using column_type = AbstractColumn<value_type<T>>;

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, DESCRIPTOR::template size<FIELD_BASE<C,U...>>()>{};
  }

  static constexpr bool isSerializable() {
    return !std::is_pointer_v<TYPE>;
  }
};

/// Base of a field of scalar TYPE<BASE_TYPE> with dimensions same as FIELD_BASE
template <template<typename> typename TYPE, unsigned C, unsigned... U>
struct TEMPLATE_FIELD_BASE : public FIELD_BASE<C,U...> {
  template <typename T>
  using value_type = TYPE<T>;

  template <typename T>
  using column_type = AbstractColumn<value_type<T>>;

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, DESCRIPTOR::template size<FIELD_BASE<C,U...>>()>{};
  }

  static constexpr bool isSerializable() {
    return false;
  }
};

/// Base of a field of pointers to the individual columns of FIELD
template <typename FIELD>
struct POINTER_FIELD_BASE {
  template <typename T>
  using value_type = std::add_pointer_t<typename FIELD::template value_type<T>>;

  template <typename T>
  using column_type = AbstractColumn<value_type<T>>;

  template <unsigned D, unsigned Q>
  static constexpr unsigned size() {
    return FIELD::template size<D,Q>();
  }

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, DESCRIPTOR::template size<POINTER_FIELD_BASE>()>{};
  }

  static constexpr bool isSerializable() {
    return false;
  }
};

/// Base of a implicitly propagatable descriptor field
struct PROPAGATABLE_FIELD_BASE : public POPULATION_FIELD {
  template <typename T>
  using value_type = T;

  template <typename T>
  using column_type = AbstractCyclicColumn<T>;
};

/// Evaluates to true iff FIELD is propagatable (e.g. POPULATION)
template <typename FIELD>
using is_propagatable_field = typename std::is_base_of<PROPAGATABLE_FIELD_BASE, FIELD>::type;

/// Base of a tensor-valued descriptor field
struct TENSOR : public FIELD_BASE_CUSTOM_SIZE {
  template <unsigned D, unsigned...>
  static constexpr unsigned size() {
    return (D * (D+1)) / 2; // see `TensorVal` in `core/util.h`
  }

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, DESCRIPTOR::template size<TENSOR>()>{};
  }
};

/// Base of degrees-of-rotational-freedom-dimensionalized fields
struct ROTATION_FIELD_BASE : public FIELD_BASE_CUSTOM_SIZE {
  template <unsigned D, unsigned...>
  static constexpr unsigned size()
  {
    if constexpr (D == 2) {
      return 1;
    } else {
      return D;
    }
  }

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, DESCRIPTOR::template size<ROTATION_FIELD_BASE>()>{};
  }

};

/// Base of a descriptor field of pointer type
template <typename TYPE>
struct OBJECT_POINTER_FIELD_BASE : public TYPED_FIELD_BASE<std::add_pointer_t<TYPE>,1> {
  static_assert(std::is_object_v<TYPE>, "TYPE must be object");
};

/// \defgroup descriptor_fields Set of common descriptor fields
/// \ingroup descriptor
//@{

struct CELL_ID      : public TYPED_FIELD_BASE<std::size_t,1> { };
struct MATERIAL     : public TYPED_FIELD_BASE<int,        1> { };
struct LATTICE_TIME : public TYPED_FIELD_BASE<std::size_t,1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value >= 0;
  }
};

struct POPULATION : public PROPAGATABLE_FIELD_BASE { };

struct STATISTIC_GENERATED : public TYPED_FIELD_BASE<int,1> { };
struct STATISTIC           : public FIELD_BASE<2> { };

// Field types need to be distinct (i.e. not aliases)
// (Field size parametrized: Cs + Ds*D + Qs*Q)  Cs Ds Qs
struct VELOCITY             : public FIELD_BASE<0,  1, 0> {
};
struct EXTERNAL_VELOCITY    : public FIELD_BASE<0,  1, 0> { };
struct VELOCITY2            : public FIELD_BASE<0,  1, 0> { };
struct VX                   : public FIELD_BASE<1,  0, 0> { };
struct VY                   : public FIELD_BASE<1,  0, 0> { };
struct VZ                   : public FIELD_BASE<1,  0, 0> { };
struct FX                   : public FIELD_BASE<1,  0, 0> { };
struct FY                   : public FIELD_BASE<1,  0, 0> { };
struct FZ                   : public FIELD_BASE<1,  0, 0> { };
struct AVERAGE_VELOCITY     : public FIELD_BASE<0,  1, 0> { };
struct AVERAGE_DENSITY      : public FIELD_BASE<1,  0, 0> { };
struct AVERAGE_TKE          : public FIELD_BASE<0,  1, 0> { };
struct SOURCE               : public FIELD_BASE<1,  0, 0> { };
struct U_TAU                : public FIELD_BASE<1,  0, 0> { };
struct PRESSCORR            : public FIELD_BASE<1,  0, 0> { };
struct FORCE                : public FIELD_BASE<0,  1, 0> { };
struct EXTERNAL_FORCE       : public FIELD_BASE<0,  1, 0> { };
struct NABLARHO             : public FIELD_BASE<0,  1, 0> { };
struct TAU_EFF              : public FIELD_BASE<1,  0, 0> { };
struct RHO                  : public FIELD_BASE<1,  0, 0> { };
struct VISCOSITY            : public FIELD_BASE<1,  0, 0> { };
struct GAMMA                : public FIELD_BASE<1,  0, 0> { };
struct CUTOFF_KIN_ENERGY    : public FIELD_BASE<1,  0, 0> { };
struct CUTOFF_HEAT_FLUX     : public FIELD_BASE<1,  0, 0> { };
struct CHEM_POTENTIAL       : public FIELD_BASE<1,  0, 0> { };
struct ADDEND               : public FIELD_BASE<1,  0, 0> { };
struct V6                   : public FIELD_BASE<6,  0, 0> { };
struct V12                  : public FIELD_BASE<12, 0, 0> { };
struct OMEGA                : public FIELD_BASE<1,  0, 0> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value > 0;
  }
};
struct TAU                  : public FIELD_BASE<1,  0, 0> { };
struct INTERFACE_WIDTH      : public FIELD_BASE<1,  0, 0> { };
struct MAGIC                : public FIELD_BASE<1,  0, 0> { };
struct G                    : public FIELD_BASE<0,  1, 0> { };
struct EPSILON              : public FIELD_BASE<1,  0, 0> { };
struct BODY_FORCE           : public FIELD_BASE<0,  1, 0> { };
struct K                    : public FIELD_BASE<1,  0, 0> { };
struct NU                   : public FIELD_BASE<1,  0, 0> { };
struct VELOCITY_NUMERATOR   : public FIELD_BASE<0,  1, 0> { };
struct VELOCITY_DENOMINATOR : public FIELD_BASE<1,  0, 0> { };
struct ZETA                 : public FIELD_BASE<0,  0, 1> { };
struct LOCAL_DRAG           : public FIELD_BASE<0,  1, 0> { };
struct VELOCITY_SOLID       : public FIELD_BASE<0,  1, 0> { };
struct COORDINATE           : public FIELD_BASE<0,  1, 0> { };
struct F                    : public FIELD_BASE<0,  0, 1> { };
struct STRAINRATE           : public FIELD_BASE<6,  0, 0> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<STRAINRATE>()>(T(0));
  }
};
struct DENSITY           : public FIELD_BASE<1, 0, 0> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<DENSITY>()>(T(1.));
  }
};
struct NEIGHBOR             : public FIELD_BASE<1,  0, 0>{
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,1>(0.);
  }
};
struct DJDF                 : public FIELD_BASE<0,  0, 1> { };
struct DJDALPHA             : public FIELD_BASE<0,  1, 0> { };
struct AV_SHEAR             : public FIELD_BASE<1,  0, 0> { };
struct SHEAR_RATE_MAGNITUDE : public FIELD_BASE<1,  0, 0> { };
struct TAU_W                : public FIELD_BASE<1,  0, 0> { };
struct SCALAR               : public FIELD_BASE<1,  0, 0> { };
struct SMAGO_CONST          : public FIELD_BASE<1,  0, 0> { };
struct EFFECTIVE_OMEGA      : public FIELD_BASE<1,  0, 0> { };
struct VELO_GRAD            : public FIELD_BASE<0,  3, 0> { };
struct FIL_RHO              : public FIELD_BASE<1,  0, 0> { };
struct LOCAL_FIL_VEL_X      : public FIELD_BASE<1,  0, 0> { };
struct LOCAL_FIL_VEL_Y      : public FIELD_BASE<1,  0, 0> { };
struct LOCAL_FIL_VEL_Z      : public FIELD_BASE<1,  0, 0> { };
struct LOCAL_AV_DISS        : public FIELD_BASE<1,  0, 0> { };
struct LOCAL_AV_TKE         : public FIELD_BASE<1,  0, 0> { };
struct LOCAL_SIGMA_ADM      : public FIELD_BASE<1,  0, 0> { };
struct LOCAL_NU_EDDY        : public FIELD_BASE<1,  0, 0> { };
struct FILTERED_VEL_GRAD    : public FIELD_BASE<0,  3, 0> { };
struct ERROR_COVARIANCE     : public FIELD_BASE<1,  0, 0> { };
struct VARIANCE             : public FIELD_BASE<1,  0, 0> { };
struct TAU_SGS              : public FIELD_BASE<1,  0, 0> { };
struct FILTERED_POPULATION  : public FIELD_BASE<0,  0, 1> { };
struct INDICATE             : public FIELD_BASE<1,  0, 0> { };
struct BIOGAS_INSTANT       : public FIELD_BASE<1,  0, 0> { };
struct BIOGAS_CUMULATIVE    : public FIELD_BASE<1,  0, 0> { };
struct METHANE_INSTANT      : public FIELD_BASE<1,  0, 0> { };
struct METHANE_CUMULATIVE   : public FIELD_BASE<1,  0, 0> { };
struct CO2_INSTANT          : public FIELD_BASE<1,  0, 0> { };
struct CO2_CUMULATIVE       : public FIELD_BASE<1,  0, 0> { };
struct TEMPERATURE          : public FIELD_BASE<1,  0, 0> { };
struct SIGMA                : public FIELD_BASE<1,  0, 0> { };
struct PSI                  : public FIELD_BASE<1,  0, 0> { };
struct THETA                : public FIELD_BASE<1,  0, 0> { };
struct NORMGRADPSI          : public FIELD_BASE<1,  0, 0> { };
struct PSI0                 : public FIELD_BASE<1,  0, 0> { };
struct INTERPHASE_NORMAL    : public FIELD_BASE<0,  1, 0> { };
struct MASS                 : public FIELD_BASE<1,  0, 0> { };
struct CELL_TYPE            : public FIELD_BASE<1,  0, 0> { };
struct BOUNDARY             : public FIELD_BASE<1,  0, 0> { };
struct SOURCE_OLD           : public FIELD_BASE<1,  0, 0> { };
struct TOP                  : public FIELD_BASE<1,  0, 0> { };
struct BOTTOM               : public FIELD_BASE<1,  0, 0> { };
struct OLD_PHIU             : public FIELD_BASE<0,  1, 0> { };
struct CONTACT_DETECTION    : public TYPED_FIELD_BASE<std::size_t, 1, 0, 0> { };
struct POROSITY             : public FIELD_BASE<1,  0, 0> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<POROSITY>()>{1};
  }

  template <typename T, typename DESCRIPTOR, typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value >= 0 && value[0]<=1;
  }
};
struct SAMPLING_DISTANCE : public descriptors::FIELD_BASE<1,0,0> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<SAMPLING_DISTANCE>()>(3);
  }
};
struct POROSITY2            : public FIELD_BASE<1,  0, 0> { };
struct EUL2LAGR             : public FIELD_BASE<1,  0, 0> { };
struct LOCATION             : public FIELD_BASE<0,  1, 0> { };
struct Y1                   : public FIELD_BASE<0,  1, 0> { };
struct Y2                   : public FIELD_BASE<0,  1, 0> { };
struct VORTICITY            : public FIELD_BASE<1,  0, 0> { };

//TODO: This expression has been removed on master lately. As no obvious equivalent could be found immediately,
//      it is added back in to enable some functionality on feature/unifiedParticleFramework.
//      In case a proper equivalent exists, this expression can be removed for good!
template <
  typename WANTED_FIELD,
  typename CURRENT_FIELD,
  typename... FIELDS,
  // WANTED_FIELD equals the head of our field list, terminate recursion
  std::enable_if_t<std::is_same<WANTED_FIELD,CURRENT_FIELD>::value, int> = 0
  >
constexpr unsigned getIndexInFieldList()
{
  return 0;
}

template <
  typename WANTED_FIELD,
  typename CURRENT_FIELD,
  typename... FIELDS,
  // WANTED_FIELD doesn't equal the head of our field list
  std::enable_if_t<!std::is_same<WANTED_FIELD,CURRENT_FIELD>::value, int> = 0
  >
constexpr unsigned getIndexInFieldList()
{
  // Break compilation when WANTED_FIELD is not provided by list of fields
  static_assert(sizeof...(FIELDS) > 0, "Field not found.");

  return 1 + getIndexInFieldList<WANTED_FIELD,FIELDS...>();
}

//@}

//@}

}

namespace fields {

struct PHYS_R : public descriptors::FIELD_BASE<0,1> { };
struct BLOCK_LOWER : public descriptors::FIELD_BASE<0,1> { };
struct BLOCK_UPPER : public descriptors::FIELD_BASE<0,1> { };

namespace moments {

struct U : public descriptors::FIELD_BASE<0,1> { };

}

}

}

// *INDENT-ON*

#endif
