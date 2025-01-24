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
struct CSE<SourcedAdvectionDiffusionBGKdynamics<T, descriptors::D3Q19<FIELDS...> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = cell.template getFieldComponent<descriptors::SOURCE>(0)*(V{0.5}*x19 + V{-1});
auto x21 = V{0.5}*cell.template getFieldComponent<descriptors::SOURCE>(0) + cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x22 = V{1}*x19 + V{-1};
auto x23 = V{0.0555555555555556}*x20;
auto x24 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x25 = x24 + V{-1};
auto x26 = x21 + V{1};
auto x27 = V{0.055556}*x19;
auto x28 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x29 = x28 + V{-1};
auto x30 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x31 = V{0.0277777777777778}*x20;
auto x32 = V{0.027778}*x19;
auto x33 = -x24;
auto x34 = x28 + V{1};
auto x35 = x30 + V{1};
auto x36 = -x28;
auto x37 = x24 + V{1};
auto x38 = -x30;
auto x0 = -cell[0]*(x19 + V{-1}) + V{0.333333}*x19*x21 - V{0.333333333333333}*x20;
auto x1 = -cell[1]*x22 - x23 - x27*(x25*x26 + V{1});
auto x2 = -cell[2]*x22 - x23 - x27*(x26*x29 + V{1});
auto x3 = -cell[3]*x22 - x23 - x27*(x26*(x30 + V{-1}) + V{1});
auto x4 = -cell[4]*x22 - x31 - x32*(x26*(x25 + x28) + V{1});
auto x5 = -cell[5]*x22 + V{0.027778}*x19*(x26*(x33 + x34) + V{-1}) - x31;
auto x6 = -cell[6]*x22 - x31 - x32*(x26*(x25 + x30) + V{1});
auto x7 = -cell[7]*x22 + V{0.027778}*x19*(x26*(x33 + x35) + V{-1}) - x31;
auto x8 = -cell[8]*x22 - x31 - x32*(x26*(x29 + x30) + V{1});
auto x9 = -cell[9]*x22 + V{0.027778}*x19*(x26*(x35 + x36) + V{-1}) - x31;
auto x10 = -cell[10]*x22 + V{0.055556}*x19*(x26*x37 + V{-1}) - x23;
auto x11 = -cell[11]*x22 + V{0.055556}*x19*(x26*x34 + V{-1}) - x23;
auto x12 = -cell[12]*x22 + V{0.055556}*x19*(x26*x35 + V{-1}) - x23;
auto x13 = -cell[13]*x22 + V{0.027778}*x19*(x26*(x28 + x37) + V{-1}) - x31;
auto x14 = -cell[14]*x22 + V{0.027778}*x19*(x26*(x36 + x37) + V{-1}) - x31;
auto x15 = -cell[15]*x22 + V{0.027778}*x19*(x26*(x30 + x37) + V{-1}) - x31;
auto x16 = -cell[16]*x22 + V{0.027778}*x19*(x26*(x37 + x38) + V{-1}) - x31;
auto x17 = -cell[17]*x22 + V{0.027778}*x19*(x26*(x30 + x34) + V{-1}) - x31;
auto x18 = -cell[18]*x22 + V{0.027778}*x19*(x26*(x34 + x38) + V{-1}) - x31;
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
cell[9] = x9;
cell[10] = x10;
cell[11] = x11;
cell[12] = x12;
cell[13] = x13;
cell[14] = x14;
cell[15] = x15;
cell[16] = x16;
cell[17] = x17;
cell[18] = x18;
return { x26, cell.template getFieldComponent<descriptors::VELOCITY>(0)*cell.template getFieldComponent<descriptors::VELOCITY>(0) + cell.template getFieldComponent<descriptors::VELOCITY>(1)*cell.template getFieldComponent<descriptors::VELOCITY>(1) + cell.template getFieldComponent<descriptors::VELOCITY>(2)*cell.template getFieldComponent<descriptors::VELOCITY>(2) };
}
};

}

}
