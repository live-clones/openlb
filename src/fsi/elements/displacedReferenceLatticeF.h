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

#ifndef FSI_DISPLACED_REFERENCE_LATTICE_F_H
#define FSI_DISPLACED_REFERENCE_LATTICE_F_H

#include "fsi/fields.h"

#include "utilities/geometricOperations.h"

namespace olb {

namespace fields::fsi {

struct ELEMENT_REFERENCE_DISPLACEMENT : public descriptors::FIELD_BASE<0,1> { };
struct ELEMENT_REFERENCE_ORIGIN : public descriptors::TYPED_FIELD_BASE<int,0,1> { };
struct ELEMENT_REFERENCE_MODIFY : public descriptors::TYPED_FIELD_BASE<int,0,1> { };

}

struct DisplacedReferenceLatticePorosityF {
  /// Data fields in element store
  using data = meta::list<
    fields::fsi::ELEMENT_LOWER,
    fields::fsi::ELEMENT_UPPER,
    fields::fsi::ELEMENT_U_TRANSLATION,
    fields::fsi::ELEMENT_REFERENCE_DELTA_X,
    fields::fsi::ELEMENT_REFERENCE_PROJECTION,
    fields::fsi::ELEMENT_REFERENCE_POROSITIES,
    fields::fsi::ELEMENT_REFERENCE_DISPLACEMENT,
    fields::fsi::ELEMENT_REFERENCE_ORIGIN,
    fields::fsi::ELEMENT_REFERENCE_MODIFY
  >;

  /// Parameter fields required for computation
  using parameters = meta::list<
    fields::array_of<fields::fsi::ELEMENT_LOWER>,
    fields::array_of<fields::fsi::ELEMENT_UPPER>,
    fields::array_of<fields::fsi::ELEMENT_PIVOT>,
    fields::array_of<fields::fsi::ELEMENT_U_TRANSLATION>,
    fields::array_of<fields::fsi::ELEMENT_REFERENCE_DELTA_X>,
    fields::array_of<fields::fsi::ELEMENT_REFERENCE_PROJECTION>,
    fields::array_of<fields::fsi::ELEMENT_REFERENCE_POROSITIES>,
    fields::array_of<fields::fsi::ELEMENT_REFERENCE_DISPLACEMENT>,
    fields::array_of<fields::fsi::ELEMENT_REFERENCE_ORIGIN>,
    fields::array_of<fields::fsi::ELEMENT_REFERENCE_MODIFY>,
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
    const auto displacements = params.template get<fields::array_of<fields::fsi::ELEMENT_REFERENCE_DISPLACEMENT>>();

    const Vector<V,DESCRIPTOR::d> elementLowerR{lowerBounds, iElement};
    const Vector<V,DESCRIPTOR::d> elementUpperR{upperBounds, iElement};
    const Vector<V,DESCRIPTOR::d> elementDisplacementR{displacements, iElement};

    if (elementLowerR < physR && physR < elementUpperR) {
      const auto refDeltaX = params.template get<fields::array_of<fields::fsi::ELEMENT_REFERENCE_DELTA_X>>();
      const auto refProjection = params.template get<fields::array_of<fields::fsi::ELEMENT_REFERENCE_PROJECTION>>();
      const auto refLattice = params.template get<fields::array_of<fields::fsi::ELEMENT_REFERENCE_POROSITIES>>();

      const auto refOrigin = params.template get<fields::array_of<fields::fsi::ELEMENT_REFERENCE_ORIGIN>>();
      const auto refModify = params.template get<fields::array_of<fields::fsi::ELEMENT_REFERENCE_MODIFY>>();

      const Vector<std::size_t,DESCRIPTOR::d> projection{refProjection, iElement};
      Vector<int,DESCRIPTOR::d> discreteR = round((physR - elementLowerR + elementDisplacementR) / refDeltaX[iElement]);
      discreteR *= Vector<int,DESCRIPTOR::d>{refModify, iElement};
      discreteR += Vector<int,DESCRIPTOR::d>{refOrigin, iElement};

      auto porosity = refLattice[iElement][projection * discreteR];
      for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
        porosity += refLattice[iElement][projection * (discreteR + descriptors::c<DESCRIPTOR>(iPop))];
      }
      porosity /= DESCRIPTOR::q;
      return porosity;
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
    // TODO genericize
    u *= (physR[0] > V{0.084});
    return u;
  }
};

}

#endif
