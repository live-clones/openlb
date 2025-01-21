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

#ifndef FSI_PRIMITIVE_F_H
#define FSI_PRIMITIVE_F_H

#include "fsi/fields.h"

#include "utilities/geometricOperations.h"

namespace olb {

struct CuboidPorosityF {
  /// Data fields in element store
  using data = meta::list<
    fields::fsi::ELEMENT_LOWER,
    fields::fsi::ELEMENT_UPPER,
    fields::fsi::ELEMENT_ROTATION,
    fields::fsi::ELEMENT_U_TRANSLATION,
    fields::fsi::ELEMENT_U_ROTATION
  >;

  /// Parameter fields required for computation
  using parameters = meta::list<
    fields::array_of<fields::fsi::ELEMENT_LOWER>,
    fields::array_of<fields::fsi::ELEMENT_UPPER>,
    fields::array_of<fields::fsi::ELEMENT_PIVOT>,
    fields::array_of<fields::fsi::ELEMENT_ROTATION>,
    fields::array_of<fields::fsi::ELEMENT_U_TRANSLATION>,
    fields::array_of<fields::fsi::ELEMENT_U_ROTATION>,
    // General non-FSI specific parameters
    fields::converter::PHYS_VELOCITY,
    fields::converter::PHYS_DELTA_X
  >;

  template <typename PARAMETERS, typename PHYS_R>
  int tag(PARAMETERS& params,
          PHYS_R& physR,
          unsigned iElement) any_platform {
    return params.template get<fields::array_of<fields::fsi::ELEMENT_TAG>>()[iElement];
  }

  template <typename PARAMETERS, typename PHYS_R>
  bool isInterior(PARAMETERS& params,
                  PHYS_R& physR,
                  unsigned iElement) any_platform {
    return false;
  }

  template <typename PARAMETERS, typename PHYS_R>
  auto compute(PARAMETERS& params,
               PHYS_R& physR,
               unsigned iElement) any_platform {
    using V = typename PARAMETERS::value_t;
    using DESCRIPTOR = typename PARAMETERS::descriptor_t;

    const auto lowerBounds = params.template get<fields::array_of<fields::fsi::ELEMENT_LOWER>>();
    const auto upperBounds = params.template get<fields::array_of<fields::fsi::ELEMENT_UPPER>>();
    const auto pivots = params.template get<fields::array_of<fields::fsi::ELEMENT_PIVOT>>();
    const auto rotations = params.template get<fields::array_of<fields::fsi::ELEMENT_ROTATION>>();

    const Vector<V,DESCRIPTOR::d> elementLowerR{lowerBounds, iElement};
    const Vector<V,DESCRIPTOR::d> elementUpperR{upperBounds, iElement};
    const Vector<V,DESCRIPTOR::d> elementPivotR{pivots, iElement};
    const FieldD<V,DESCRIPTOR,fields::fsi::ELEMENT_ROTATION> elementRotation{rotations, iElement};

    const Vector<V,DESCRIPTOR::d> shiftedPhysR = physR - elementPivotR;
    const auto rotatedPhysR = util::matrixVectorProduct(util::calculateRotationMatrix<V,DESCRIPTOR::d>(-1*elementRotation),
                                                        shiftedPhysR);
    const Vector<V,DESCRIPTOR::d> shiftedLowerR = elementLowerR - elementPivotR;
    const Vector<V,DESCRIPTOR::d> shiftedUpperR = elementUpperR - elementPivotR;

    const auto physDeltaX = params.template get<fields::converter::PHYS_DELTA_X>();

    if (shiftedLowerR-physDeltaX < rotatedPhysR && rotatedPhysR < shiftedUpperR+physDeltaX) {
      const auto centerR = V{0.5}*(elementUpperR - elementLowerR);
      const V eps = V{0.5}*physDeltaX;
      const V d = sdf::box(rotatedPhysR - centerR, centerR-0.5*eps);
      return util::min(util::max(d / eps, 0), 1);
    } else {
      return V{1};
    }
  }

  template <typename PARAMETERS, typename PHYS_R>
  auto computeU(PARAMETERS& params,
                PHYS_R& physR,
                unsigned iElement) any_platform {
    using V = typename PARAMETERS::value_t;
    using DESCRIPTOR = typename PARAMETERS::descriptor_t;

    const auto pivots = params.template get<fields::array_of<fields::fsi::ELEMENT_PIVOT>>();
    const auto translationUs = params.template get<fields::array_of<fields::fsi::ELEMENT_U_TRANSLATION>>();
    const auto rotationUs = params.template get<fields::array_of<fields::fsi::ELEMENT_U_ROTATION>>();

    const Vector<V,DESCRIPTOR::d> elementPivotR{pivots, iElement};
    const Vector<V,DESCRIPTOR::d> shiftedPhysR = physR - elementPivotR;

    Vector<V,DESCRIPTOR::d> u{translationUs, iElement};
    if constexpr (DESCRIPTOR::d == 2) {
      u[0] -= rotationUs[iElement] * shiftedPhysR[1];
      u[1] += rotationUs[iElement] * shiftedPhysR[0];
    } else {
      u += crossProduct(Vector<V,DESCRIPTOR::d>{rotationUs, iElement},
                        shiftedPhysR);
    }
    u *= params.template get<fields::converter::PHYS_VELOCITY>();
    return u;
  }
};

struct SpherePorosityF {
  struct RADIUS : public descriptors::FIELD_BASE<1> { };

  /// Data fields in element store
  using data = meta::list<
    RADIUS,
    fields::fsi::ELEMENT_U_TRANSLATION
  >;

  /// Parameter fields required for computation
  using parameters = meta::list<
    fields::array_of<RADIUS>,
    fields::array_of<fields::fsi::ELEMENT_PIVOT>,
    fields::array_of<fields::fsi::ELEMENT_U_TRANSLATION>,
    // General non-FSI specific parameters
    fields::converter::PHYS_VELOCITY,
    fields::converter::PHYS_DELTA_X
  >;

  template <typename PARAMETERS, typename PHYS_R>
  int tag(PARAMETERS& params,
          PHYS_R& physR,
          unsigned iElement) any_platform {
    return params.template get<fields::array_of<fields::fsi::ELEMENT_TAG>>()[iElement];
  }

  template <typename PARAMETERS, typename PHYS_R>
  bool isInterior(PARAMETERS& params,
                  PHYS_R& physR,
                  unsigned iElement) any_platform {
    return false;
  }

  template <typename PARAMETERS, typename PHYS_R>
  auto compute(PARAMETERS& params,
               PHYS_R& physR,
               unsigned iElement) any_platform {
    using V = typename PARAMETERS::value_t;
    using DESCRIPTOR = typename PARAMETERS::descriptor_t;

    const auto pivots = params.template get<fields::array_of<fields::fsi::ELEMENT_PIVOT>>();
    const auto radii  = params.template get<fields::array_of<RADIUS>>();

    const Vector<V,DESCRIPTOR::d> elementPivot{pivots, iElement};
    const V elementRadius = radii[iElement];

    const auto physDeltaX = params.template get<fields::converter::PHYS_DELTA_X>();

    const V d = sdf::sphere(physR - elementPivot, elementRadius);
    const V eps = V{0.5}*physDeltaX;
    if (d < eps) {
      return util::max(d / eps, 0);
    } else {
      return V{1};
    }
  }

  template <typename PARAMETERS, typename PHYS_R>
  auto computeU(PARAMETERS& params,
                PHYS_R& physR,
                unsigned iElement) any_platform {
    using V = typename PARAMETERS::value_t;
    using DESCRIPTOR = typename PARAMETERS::descriptor_t;

    const auto translationUs = params.template get<fields::array_of<fields::fsi::ELEMENT_U_TRANSLATION>>();

    Vector<V,DESCRIPTOR::d> u{translationUs, iElement};
    u *= params.template get<fields::converter::PHYS_VELOCITY>();
    return u;
  }
};

struct LinePorosityF {
  struct WIDTH : public descriptors::FIELD_BASE<1> { };

  struct POINT_A : public descriptors::FIELD_BASE<0,1> { };
  struct POINT_B : public descriptors::FIELD_BASE<0,1> { };

  struct POINT_A_U : public descriptors::FIELD_BASE<0,1> { };
  struct POINT_B_U : public descriptors::FIELD_BASE<0,1> { };

  struct POINT_A_F : public descriptors::FIELD_BASE<0,1> { };
  struct POINT_B_F : public descriptors::FIELD_BASE<0,1> { };

  /// Data fields in element store
  using data = meta::list<
    WIDTH,
    POINT_A,
    POINT_B,
    POINT_A_U,
    POINT_B_U,
    POINT_A_F,
    POINT_B_F,
    fields::fsi::ELEMENT_ROTATION,
    fields::fsi::ELEMENT_U_TRANSLATION,
    fields::fsi::ELEMENT_U_ROTATION
  >;

  /// Parameter fields required for computation
  using parameters = meta::list<
    fields::array_of<WIDTH>,
    fields::array_of<POINT_A>,
    fields::array_of<POINT_B>,
    fields::array_of<POINT_A_U>,
    fields::array_of<POINT_B_U>,
    fields::array_of<fields::fsi::ELEMENT_PIVOT>,
    // General non-FSI specific parameters
    fields::converter::PHYS_VELOCITY,
    fields::converter::PHYS_DELTA_X
  >;

  template <typename PARAMETERS, typename PHYS_R>
  int tag(PARAMETERS& params,
          PHYS_R& physR,
          unsigned iElement) any_platform {
    return params.template get<fields::array_of<fields::fsi::ELEMENT_TAG>>()[iElement];
  }

  template <typename PARAMETERS, typename PHYS_R>
  bool isInterior(PARAMETERS& params,
                  PHYS_R& physR,
                  unsigned iElement) any_platform {
    return false;
  }

  template <typename PARAMETERS, typename PHYS_R>
  auto compute(PARAMETERS& params,
               PHYS_R& physR,
               unsigned iElement) any_platform {
    using V = typename PARAMETERS::value_t;
    using DESCRIPTOR = typename PARAMETERS::descriptor_t;

    const auto as = params.template get<fields::array_of<POINT_A>>();
    const auto bs = params.template get<fields::array_of<POINT_B>>();
    const auto ws = params.template get<fields::array_of<WIDTH>>();

    const Vector<V,DESCRIPTOR::d> elementA{as, iElement};
    const Vector<V,DESCRIPTOR::d> elementB{bs, iElement};
    const V elementW = ws[iElement];

    const auto physDeltaX = params.template get<fields::converter::PHYS_DELTA_X>();

    const V d = sdf::line(physR, elementA, elementB, elementW);
    const V eps = V{0.5}*physDeltaX;
    if (d < eps) {
      return util::max(d / eps, 0);
    } else {
      return V{1};
    }
  }

  template <typename PARAMETERS, typename PHYS_R>
  auto computeU(PARAMETERS& params,
                PHYS_R& physR,
                unsigned iElement) any_platform {
    using V = typename PARAMETERS::value_t;
    using DESCRIPTOR = typename PARAMETERS::descriptor_t;

    const auto as = params.template get<fields::array_of<POINT_A>>();
    const auto bs = params.template get<fields::array_of<POINT_B>>();
    const auto aUs = params.template get<fields::array_of<POINT_A_U>>();
    const auto bUs = params.template get<fields::array_of<POINT_B_U>>();

    const Vector<V,DESCRIPTOR::d> a{as, iElement};
    const Vector<V,DESCRIPTOR::d> b{bs, iElement};

    const auto ba = b - a;
    const auto baNorm = norm(ba);
    const auto baHat = ba / baNorm;
    const auto pa = physR - a;
    const auto frac = sdf::clamp((pa * baHat) / baNorm, 0, 1);

    const Vector<V,DESCRIPTOR::d> aU{aUs, iElement};
    const Vector<V,DESCRIPTOR::d> bU{bUs, iElement};

    Vector<V,DESCRIPTOR::d> u = (1-frac)*bU + frac*aU;
    u *= params.template get<fields::converter::PHYS_VELOCITY>();
    return u;
  }
};

}

#endif
