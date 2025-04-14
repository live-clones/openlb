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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x21 = cell[15] + cell[17];
auto x22 = cell[12] + x21;
auto x23 = cell[11] + cell[18];
auto x24 = cell[10] + cell[14] + cell[16];
auto x25 = cell[2] + cell[8] + cell[9];
auto x26 = cell[13] + cell[3];
auto x27 = cell[0] + cell[1] + cell[4] + cell[5] + cell[6] + cell[7] + x22 + x23 + x24 + x25 + x26;
auto x28 = x27 + V{1};
auto x29 = V{1} / (x28);
auto x30 = x27 + V{1};
auto x31 = V{1} / ((x28)*(x28));
auto x32 = V{1}*x31;
auto x33 = x30*x32;
auto x34 = -cell[18];
auto x35 = -cell[3];
auto x36 = -cell[8];
auto x37 = cell[9] + x36;
auto x38 = x34 + x35 + x37;
auto x39 = -cell[6];
auto x40 = cell[7] + x39;
auto x41 = -cell[16] + x40;
auto x42 = x22 + x38 + x41;
auto x43 = x42*x42;
auto x44 = V{0.333333}*cell[0];
auto x45 = V{0.333333}*cell[10] - V{0.666667}*cell[17] - V{0.666667}*cell[18] + V{0.333333}*cell[1] - V{0.666667}*cell[8] - V{0.666667}*cell[9] + x44;
auto x46 = V{0.333333}*cell[11] - V{0.666667}*cell[15] - V{0.666667}*cell[16] + V{0.333333}*cell[2] - V{0.666667}*cell[6] - V{0.666667}*cell[7];
auto x47 = -V{0.666667}*cell[12] + V{0.333333}*cell[13] + V{0.333333}*cell[14] - V{0.666667}*cell[3] + V{0.333333}*cell[4] + V{0.333333}*cell[5] + x33*x43 + x45 + x46;
auto x48 = cell[13] + cell[17];
auto x49 = -cell[2];
auto x50 = -cell[9];
auto x51 = x23 + x36 + x49 + x50;
auto x52 = -cell[4];
auto x53 = cell[5] + x52;
auto x54 = -cell[14] + x53;
auto x55 = x48 + x51 + x54;
auto x56 = x55*x55;
auto x57 = V{0.333333}*cell[12] - V{0.666667}*cell[13] - V{0.666667}*cell[14] + V{0.333333}*cell[3] - V{0.666667}*cell[4] - V{0.666667}*cell[5];
auto x58 = -V{0.666667}*cell[11] + V{0.333333}*cell[15] + V{0.333333}*cell[16] - V{0.666667}*cell[2] + V{0.333333}*cell[6] + V{0.333333}*cell[7] + x33*x56 + x45 + x57;
auto x59 = cell[13] + cell[15];
auto x60 = -cell[1];
auto x61 = -cell[7];
auto x62 = x39 + x60 + x61;
auto x63 = -cell[5] + x52;
auto x64 = x24 + x59 + x62 + x63;
auto x65 = x64*x64;
auto x66 = -V{0.666667}*cell[10] + V{0.333333}*cell[17] + V{0.333333}*cell[18] - V{0.666667}*cell[1] + V{0.333333}*cell[8] + V{0.333333}*cell[9] + x33*x65 + x44 + x46 + x57;
auto x67 = x29*x64;
auto x68 = -cell[13] + cell[14] + x53 + x55*x67;
auto x69 = -cell[15] + cell[16];
auto x70 = x40 + x42*x67 + x69;
auto x71 = -cell[17];
auto x72 = cell[18] + x71;
auto x73 = x30*x31*x42*x55 + x37 + x72;
auto x74 = V{1} / (V{3.00000046417339}*util::sqrt(x29*(x20*x20)*util::sqrt(V{0.5}*(x47*x47) + V{0.5}*(x58*x58) + V{0.5}*(x66*x66) + x68*x68 + x70*x70 + x73*x73) + V{0.0277777691819762}/((x19)*(x19))) + V{0.5}/x19);
auto x75 = V{1.5}*x31;
auto x76 = x65*x75;
auto x77 = x56*x75;
auto x78 = x43*x75;
auto x79 = x77 + x78 + V{-1};
auto x80 = x76 + x79;
auto x81 = V{1} - x74;
auto x82 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x83 = V{3}*cell[14];
auto x84 = V{3}*cell[16];
auto x85 = V{3}*cell[5];
auto x86 = V{3}*cell[7];
auto x87 = V{3}*cell[13] - V{3}*cell[4];
auto x88 = V{3}*cell[15] - V{3}*cell[6];
auto x89 = x29*(V{3}*cell[10] - V{3}*cell[1] + x83 + x84 - x85 - x86 + x87 + x88);
auto x90 = V{3}*x31;
auto x91 = x65*x90;
auto x92 = V{3}*cell[18];
auto x93 = V{3}*cell[9];
auto x94 = V{3}*cell[17] - V{3}*cell[8];
auto x95 = x29*(V{3}*cell[11] - V{3}*cell[2] - x83 + x85 + x87 + x92 - x93 + x94);
auto x96 = x56*x90;
auto x97 = x76 + V{-1};
auto x98 = x29*(V{3}*cell[12] - V{3}*cell[3] - x84 + x86 + x88 - x92 + x93 + x94);
auto x99 = x43*x90;
auto x100 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x101 = V{4.5}*x31;
auto x102 = cell[10] + cell[16] + x62;
auto x103 = V{2}*cell[13] - V{2}*cell[4] + x102 + x21 + x51;
auto x104 = x101*(x103*x103);
auto x105 = x80 + x89;
auto x106 = -x95;
auto x107 = -cell[11] + V{2}*cell[14] + cell[15] - V{2}*cell[5] + x102 + x25 + x34 + x71;
auto x108 = -x107;
auto x109 = cell[10] + cell[14] + x60 + x63;
auto x110 = cell[12] + V{2}*cell[15] - V{2}*cell[6] + x109 + x38 + x48;
auto x111 = x101*(x110*x110);
auto x112 = -x98;
auto x113 = -cell[12] + x26;
auto x114 = V{2}*cell[16] - V{2}*cell[7] + cell[8] + x109 + x113 + x50 + x72;
auto x115 = -x114;
auto x116 = cell[11] + x49 + x54;
auto x117 = cell[12] + V{2}*cell[17] - V{2}*cell[8] + x116 + x35 + x41 + x59;
auto x118 = x101*(x117*x117);
auto x119 = x80 + x95;
auto x120 = V{2}*cell[18] + cell[6] - V{2}*cell[9] + x113 + x116 + x61 + x69;
auto x121 = -x120;
auto x122 = -x77;
auto x123 = V{1} - x78;
auto x124 = x122 + x123;
auto x125 = x124 + x89;
auto x126 = -x76;
auto x127 = x126 + x95;
auto x128 = x126 + x98;
auto x129 = -x89;
auto x130 = x80 + x98;
auto x0 = V{1}*cell[0]*x81 - x74*(x80*(V{0.333333333333333}*cell[0] + V{0.333333333333333}*cell[10] + V{0.333333333333333}*cell[11] + V{0.333333333333333}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] + V{0.333333333333333}) + V{0.333333333333333});
auto x1 = V{1}*cell[1]*x81 - x74*(x82*(x79 + x89 - x91) + V{0.0555555555555556});
auto x2 = V{1}*cell[2]*x81 - x74*(x82*(x78 + x95 - x96 + x97) + V{0.0555555555555556});
auto x3 = V{1}*cell[3]*x81 - x74*(x82*(x77 + x97 + x98 - x99) + V{0.0555555555555556});
auto x4 = V{1}*cell[4]*x81 - x74*(x100*(-x104 + x105 + x95) + V{0.0277777777777778});
auto x5 = V{1}*cell[5]*x81 - x74*(x100*(-x101*x108*x108 + x105 + x106) + V{0.0277777777777778});
auto x6 = V{1}*cell[6]*x81 - x74*(x100*(x105 - x111 + x98) + V{0.0277777777777778});
auto x7 = V{1}*cell[7]*x81 - x74*(x100*(-x101*x115*x115 + x105 + x112) + V{0.0277777777777778});
auto x8 = V{1}*cell[8]*x81 - x74*(x100*(-x118 + x119 + x98) + V{0.0277777777777778});
auto x9 = V{1}*cell[9]*x81 - x74*(x100*(-x101*x121*x121 + x112 + x119) + V{0.0277777777777778});
auto x10 = V{1}*cell[10]*x81 + x74*(x82*(x125 + x91) + V{-0.0555555555555556});
auto x11 = V{1}*cell[11]*x81 + x74*(x82*(x123 + x127 + x96) + V{-0.0555555555555556});
auto x12 = V{1}*cell[12]*x81 + x74*(x82*(x122 + x128 + x99 + V{1}) + V{-0.0555555555555556});
auto x13 = V{1}*cell[13]*x81 + x74*(x100*(x104 + x125 + x127) + V{-0.0277777777777778});
auto x14 = V{1}*cell[14]*x81 - x74*(x100*(-x101*x107*x107 + x119 + x129) + V{0.0277777777777778});
auto x15 = V{1}*cell[15]*x81 + x74*(x100*(x111 + x125 + x128) + V{-0.0277777777777778});
auto x16 = V{1}*cell[16]*x81 - x74*(x100*(-x101*x114*x114 + x129 + x130) + V{0.0277777777777778});
auto x17 = V{1}*cell[17]*x81 + x74*(x100*(x118 + x124 + x127 + x98) + V{-0.0277777777777778});
auto x18 = V{1}*cell[18]*x81 - x74*(x100*(-x101*x120*x120 + x106 + x130) + V{0.0277777777777778});
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
return { x28, x32*(x43 + x56 + x65) };
}
};

}

}
