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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::ConstRhoBGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = parameters.template get<statistics::AVERAGE_RHO>();
auto x11 = x9 + V{-1};
auto x12 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x13 = x12 + V{1};
auto x14 = V{1} / ((x13)*(x13));
auto x15 = V{1.5}*x14;
auto x16 = cell[1] - cell[5];
auto x17 = -cell[4] + cell[8];
auto x18 = -cell[3] + cell[7] + x16 + x17;
auto x19 = x18*x18;
auto x20 = x15*x19;
auto x21 = cell[2] - cell[6];
auto x22 = cell[3] - cell[7] + x16 + x21;
auto x23 = x22*x22;
auto x24 = x15*x23;
auto x25 = x24 + V{-1};
auto x26 = x20 + x25;
auto x27 = V{1} / (x13);
auto x28 = -x27*(x10 + V{-1}) + V{1};
auto x29 = x12 + V{1};
auto x30 = x28*x29;
auto x31 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778};
auto x32 = V{4.5}*x14;
auto x33 = V{2}*cell[1] - V{2}*cell[5] + x17 + x21;
auto x34 = x32*(x33*x33);
auto x35 = V{3}*cell[7];
auto x36 = V{3}*cell[3];
auto x37 = V{3}*cell[1] - V{3}*cell[5];
auto x38 = x27*(-V{3}*cell[4] + V{3}*cell[8] + x35 - x36 + x37);
auto x39 = -x24;
auto x40 = x38 + x39;
auto x41 = x27*(V{3}*cell[2] - V{3}*cell[6] - x35 + x36 + x37);
auto x42 = V{1} - x20;
auto x43 = x41 + x42;
auto x44 = x34 + x40 + x43;
auto x45 = V{0.111111111111111}*cell[0] + V{0.111111111111111}*cell[1] + V{0.111111111111111}*cell[2] + V{0.111111111111111}*cell[3] + V{0.111111111111111}*cell[4] + V{0.111111111111111}*cell[5] + V{0.111111111111111}*cell[6] + V{0.111111111111111}*cell[7] + V{0.111111111111111}*cell[8] + V{0.111111111111111};
auto x46 = V{3}*x14;
auto x47 = x23*x46 + x42;
auto x48 = x41 + x47;
auto x49 = V{2}*cell[3] + cell[4] - V{2}*cell[7] - cell[8] + x21;
auto x50 = x32*(x49*x49);
auto x51 = -x38 + x39 + x43 + x50;
auto x52 = x19*x46;
auto x53 = x25 + x38 - x52;
auto x54 = x26 - x34 + x38 + x41;
auto x55 = -x41;
auto x56 = x47 + x55;
auto x57 = x40 + x42 + x50 + x55;
auto x58 = x40 + x52 + V{1};
auto x0 = -x11*(cell[0] + x26*(V{0.444444444444444}*cell[0] + V{0.444444444444444}*cell[1] + V{0.444444444444444}*cell[2] + V{0.444444444444444}*cell[3] + V{0.444444444444444}*cell[4] + V{0.444444444444444}*cell[5] + V{0.444444444444444}*cell[6] + V{0.444444444444444}*cell[7] + V{0.444444444444444}*cell[8] + V{0.444444444444444}) + V{0.444444444444444}) - V{0.444444444444444}*x26*x30 + V{-0.444444444444444};
auto x1 = -x11*(cell[1] - x31*x44 + V{0.0277777777777778}) + V{0.0277777777777778}*x28*x29*x44 + V{-0.0277777777777778};
auto x2 = -x11*(cell[2] - x45*x48 + V{0.111111111111111}) + V{0.111111111111111}*x28*x29*x48 + V{-0.111111111111111};
auto x3 = -x11*(cell[3] - x31*x51 + V{0.0277777777777778}) + V{0.0277777777777778}*x28*x29*x51 + V{-0.0277777777777778};
auto x4 = -x11*(cell[4] + x45*x53 + V{0.111111111111111}) - V{0.111111111111111}*x30*x53 + V{-0.111111111111111};
auto x5 = -x11*(cell[5] + x31*x54 + V{0.0277777777777778}) - V{0.0277777777777778}*x30*x54 + V{-0.0277777777777778};
auto x6 = -x11*(cell[6] - x45*x56 + V{0.111111111111111}) + V{0.111111111111111}*x28*x29*x56 + V{-0.111111111111111};
auto x7 = -x11*(cell[7] - x31*x57 + V{0.0277777777777778}) + V{0.0277777777777778}*x28*x29*x57 + V{-0.0277777777777778};
auto x8 = -x11*(cell[8] - x45*x58 + V{0.111111111111111}) + V{0.111111111111111}*x28*x29*x58 + V{-0.111111111111111};
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { -x10 + x12 + V{2}, V{1}*x14*(x19 + x23) };
}
};

}

}
