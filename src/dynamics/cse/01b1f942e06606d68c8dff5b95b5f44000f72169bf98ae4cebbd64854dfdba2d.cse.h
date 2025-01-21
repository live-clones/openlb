/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021-24 Adrian Kummerlaender, Shota Ito
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

/*  ========================================================
 *  ==  WARNING: This is an automatically generated file, ==
 *  ==                  do not modify.                    ==
 *  ========================================================
 */

#pragma once


namespace olb {

namespace dynamics {

template <typename T, typename... FIELDS>
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::FreeEnergyInletOutletDensity, momenta::FreeEnergyInletOutletMomentum<1, -1>, momenta::RegularizedBoundaryStress<1, -1>, momenta::DefineSeparately>, equilibria::FreeEnergy, collision::FreeEnergyInletOutlet<1, -1>, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = x9 + V{-1};
auto x11 = cell[0]*x10;
auto x12 = cell.template getFieldComponent<descriptors::FORCE>(1)*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x13 = cell.template getFieldComponent<descriptors::FORCE>(0)*(V{1.5}*x12 + V{-1}) + V{1};
auto x14 = V{0.444444}*x13*x9;
auto x15 = cell[2]*x10;
auto x16 = cell[3]*x10;
auto x17 = cell[4]*x10;
auto x18 = cell[5]*x10;
auto x19 = cell[6]*x10;
auto x20 = V{0.111111}*x9;
auto x21 = x13*x20;
auto x22 = cell.template getFieldComponent<descriptors::FORCE>(0)*(-V{3}*cell.template getFieldComponent<descriptors::FORCE>(1) + V{3}*x12 + V{1});
auto x23 = x22 + V{-1};
auto x24 = cell.template getFieldComponent<descriptors::FORCE>(0)*x9;
auto x25 = V{0.166666666666667}*cell.template getFieldComponent<descriptors::FORCE>(0) + V{0.166666666666667}*x11 + V{0.166666666666667}*x15 + V{0.166666666666667}*x16 + V{0.166666666666667}*x17 + V{0.166666666666667}*x18 + V{0.166666666666667}*x19 + x21 - V{0.0277778333333333}*x23*x9 - V{0.0277777777777778}*x24 + V{-0.166666666666667};
auto x26 = V{0.111111111111111}*x24;
auto x27 = x21 + x26;
auto x28 = V{0.0277777777777778}*x24;
auto x29 = V{1} - x22;
cell[0] = V{0.555555555555556}*cell.template getFieldComponent<descriptors::FORCE>(0)*x9 - x11 - x14;
cell[1] = x25;
cell[2] = -x15 - x27;
cell[3] = -x16 - x28 - V{0.027778}*x29*x9;
cell[4] = -x17 - x20*x29 - x26;
cell[5] = -x18 + V{0.027778}*x23*x9 - x28;
cell[6] = -x19 - x27;
cell[7] = x25;
cell[8] = V{0.666666666666667}*cell.template getFieldComponent<descriptors::FORCE>(0) + V{0.666666666666667}*x11 + x14 + V{0.666666666666667}*x15 + V{0.666666666666667}*x16 + V{0.666666666666667}*x17 + V{0.666666666666667}*x18 + V{0.666666666666667}*x19 - V{0.111111333333333}*x23*x9 - V{0.111111111111111}*x24 + V{-0.666666666666667};
return { cell.template getFieldComponent<descriptors::FORCE>(0), x12 };
}
};

}

}
