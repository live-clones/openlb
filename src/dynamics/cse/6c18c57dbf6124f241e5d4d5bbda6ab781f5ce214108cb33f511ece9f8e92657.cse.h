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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::fsi::HLBM>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = x9 + V{-1};
auto x11 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x12 = x11 + V{1};
auto x13 = x11 + V{1};
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
auto x28 = V{1}*cell[7];
auto x29 = V{1}*cell[3];
auto x30 = V{1}*cell[1] - V{1}*cell[5];
auto x31 = x27*(-V{1}*cell[4] + V{1}*cell[8] + x28 - x29 + x30);
auto x32 = cell.template getFieldComponent<descriptors::POROSITY>(0) + V{-1};
auto x33 = -x32;
auto x34 = -x31;
auto x35 = cell.template getFieldComponent<descriptors::VELOCITY>(1) + x34;
auto x36 = x33*x35;
auto x37 = x31 + x36;
auto x38 = x27*(V{1}*cell[2] - V{1}*cell[6] - x28 + x29 + x30);
auto x39 = cell.template getFieldComponent<descriptors::VELOCITY>(0) + x38;
auto x40 = x33*x39;
auto x41 = -x38 + x40;
auto x42 = V{-1} + V{1.5}*(x37*x37) + V{1.5}*(x41*x41);
auto x43 = V{0.0277777777777778}*x9;
auto x44 = V{4.5}*x14;
auto x45 = V{2}*cell[1] - V{2}*cell[5] + x17 + x21;
auto x46 = x44*(x45*x45);
auto x47 = V{3}*cell[7];
auto x48 = V{3}*cell[3];
auto x49 = V{3}*cell[1] - V{3}*cell[5];
auto x50 = x27*(-V{3}*cell[4] + V{3}*cell[8] + x47 - x48 + x49);
auto x51 = -x24;
auto x52 = x50 + x51;
auto x53 = x27*(V{3}*cell[2] - V{3}*cell[6] - x47 + x48 + x49);
auto x54 = V{1} - x20;
auto x55 = x53 + x54;
auto x56 = x46 + x52 + x55;
auto x57 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778};
auto x58 = x32*x39;
auto x59 = V{3}*x58;
auto x60 = x18*x27;
auto x61 = x32*x35;
auto x62 = x22*x27;
auto x63 = x58 + x62;
auto x64 = x60 - x61 + x63;
auto x65 = x50 + x53;
auto x66 = x38 + x58;
auto x67 = V{1.5}*(x66*x66);
auto x68 = -x34 - x61;
auto x69 = V{1.5}*(x68*x68);
auto x70 = V{3}*x61;
auto x71 = -x67 - x69 - x70 + V{1};
auto x72 = V{0.111111111111111}*cell[0] + V{0.111111111111111}*cell[1] + V{0.111111111111111}*cell[2] + V{0.111111111111111}*cell[3] + V{0.111111111111111}*cell[4] + V{0.111111111111111}*cell[5] + V{0.111111111111111}*cell[6] + V{0.111111111111111}*cell[7] + V{0.111111111111111}*cell[8] + V{0.111111111111111};
auto x73 = V{3}*x14;
auto x74 = x23*x73 + x54;
auto x75 = x53 + x74;
auto x76 = x40 - x62;
auto x77 = -x76;
auto x78 = -x53;
auto x79 = V{3}*x40;
auto x80 = x42 + x78 + x79;
auto x81 = -x50;
auto x82 = V{2}*cell[3] + cell[4] - V{2}*cell[7] - cell[8] + x21;
auto x83 = x44*(x82*x82);
auto x84 = x51 + x55 + x81 + x83;
auto x85 = x36 + x60;
auto x86 = -x76 - x85;
auto x87 = V{3}*x36 + x50;
auto x88 = V{0.111111111111111}*x9;
auto x89 = x19*x73;
auto x90 = x25 + x50 - x89;
auto x91 = x42 + x87;
auto x92 = x26 - x46 + x65;
auto x93 = x40 - x62 - x85;
auto x94 = x74 + x78;
auto x95 = x53 + x59 + x67 + x69 + V{-1};
auto x96 = x52 + x54 + x78 + x83;
auto x97 = -x60 + x61;
auto x98 = -x63 - x97;
auto x99 = x52 + x89 + V{1};
auto x100 = -x97;
auto x0 = -cell[0]*x10 + V{0.444444444444444}*x12*x26 - x42*(V{0.444444444444444}*cell[0] + V{0.444444444444444}*cell[1] + V{0.444444444444444}*cell[2] + V{0.444444444444444}*cell[3] + V{0.444444444444444}*cell[4] + V{0.444444444444444}*cell[5] + V{0.444444444444444}*cell[6] + V{0.444444444444444}*cell[7] + V{0.444444444444444}*cell[8] + V{0.444444444444444}) - V{0.444444444444444}*x9*(x12*x26 + V{1});
auto x1 = -cell[1]*x10 + x43*(x12*x56 + V{-1}) - x56*x57 + x57*(x59 + x65 + x71 + V{4.5}*(x64*x64));
auto x2 = -(cell[2]*x10 + x72*x75 + x72*(x80 - V{4.5}*x77*x77) - V{0.111111111111111}*x9*(x12*x75 + V{-1}));
auto x3 = -(cell[3]*x10 + x57*x84 + x57*(x80 + x87 - V{4.5}*x86*x86) - V{0.0277777777777778}*x9*(x12*x84 + V{-1}));
auto x4 = -(cell[4]*x10 - V{0.111111111111111}*x12*x90 + x72*(x91 - V{4.5}*x85*x85) + x88*(x12*x90 + V{1}));
auto x5 = -(cell[5]*x10 - V{0.0277777777777778}*x12*x92 + x43*(x12*x92 + V{1}) + x57*(x53 - x79 + x91 - V{4.5}*x93*x93));
auto x6 = -cell[6]*x10 - x72*x94 - x72*(x95 - V{4.5}*x63*x63) + x88*(x12*x94 + V{-1});
auto x7 = -cell[7]*x10 + x43*(x12*x96 + V{-1}) - x57*x96 - x57*(x70 + x81 + x95 - V{4.5}*x98*x98);
auto x8 = -cell[8]*x10 - x72*x99 + x72*(x50 + x71 + V{4.5}*(x100*x100)) + x88*(x12*x99 + V{-1});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x13, V{1}*x14*(x19 + x23) };
}
};

}

}
