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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::Kupershtokh>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x21 = cell.template getFieldComponent<descriptors::FORCE>(2);
auto x20 = cell.template getFieldComponent<descriptors::FORCE>(1);
auto x19 = cell.template getFieldComponent<descriptors::FORCE>(0);
auto x22 = parameters.template get<descriptors::OMEGA>();
auto x23 = x22 + V{-1};
auto x24 = cell[10] + cell[14];
auto x25 = cell[12] + cell[7];
auto x26 = cell[0] + cell[11] + cell[13] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[8] + cell[9] + x24 + x25;
auto x27 = x26 + V{1};
auto x28 = x26 + V{1};
auto x29 = V{1} / ((x28)*(x28));
auto x30 = V{1.5}*x29;
auto x31 = cell[13] - cell[4];
auto x32 = cell[15] - cell[6];
auto x33 = x31 + x32;
auto x34 = -cell[1];
auto x35 = cell[16] - cell[7];
auto x36 = x34 + x35;
auto x37 = -cell[5] + x24;
auto x38 = x33 + x36 + x37;
auto x39 = x38*x38;
auto x40 = x30*x39;
auto x41 = cell[17] - cell[8];
auto x42 = x31 + x41;
auto x43 = cell[18] - cell[9];
auto x44 = -cell[2];
auto x45 = cell[11] - cell[14] + cell[5] + x44;
auto x46 = x42 + x43 + x45;
auto x47 = x46*x46;
auto x48 = x30*x47;
auto x49 = x32 + x41;
auto x50 = -cell[3];
auto x51 = -cell[18] + cell[9];
auto x52 = x50 + x51;
auto x53 = -cell[16] + x25;
auto x54 = x49 + x52 + x53;
auto x55 = x54*x54;
auto x56 = x30*x55;
auto x57 = x48 + x56 + V{-1};
auto x58 = x40 + x57;
auto x59 = V{1} / (x28);
auto x60 = V{1}*cell[14];
auto x61 = V{1}*cell[16];
auto x62 = V{1}*cell[5];
auto x63 = V{1}*cell[7];
auto x64 = V{1}*cell[13] - V{1}*cell[4];
auto x65 = V{1}*cell[15] - V{1}*cell[6];
auto x66 = x19 + x59*(V{1}*cell[10] - V{1}*cell[1] + x60 + x61 - x62 - x63 + x64 + x65);
auto x67 = V{1.5}*(x66*x66);
auto x68 = V{1}*cell[18];
auto x69 = V{1}*cell[9];
auto x70 = V{1}*cell[17] - V{1}*cell[8];
auto x71 = x20 + x59*(V{1}*cell[11] - V{1}*cell[2] - x60 + x62 + x64 + x68 - x69 + x70);
auto x72 = V{1.5}*(x71*x71);
auto x73 = x21 + x59*(V{1}*cell[12] - V{1}*cell[3] - x61 + x63 + x65 - x68 + x69 + x70);
auto x74 = V{1.5}*(x73*x73);
auto x75 = x67 + x72 + x74 + V{-1};
auto x76 = V{0.0555555555555556}*x22;
auto x77 = V{3}*cell[14];
auto x78 = V{3}*cell[16];
auto x79 = V{3}*cell[5];
auto x80 = V{3}*cell[7];
auto x81 = V{3}*cell[13] - V{3}*cell[4];
auto x82 = V{3}*cell[15] - V{3}*cell[6];
auto x83 = x59*(V{3}*cell[10] - V{3}*cell[1] + x77 + x78 - x79 - x80 + x81 + x82);
auto x84 = V{3}*x29;
auto x85 = x39*x84;
auto x86 = x57 + x83 - x85;
auto x87 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x88 = x19 + x38*x59;
auto x89 = V{4.5}*(x88*x88);
auto x90 = V{3}*x19;
auto x91 = x83 + x90;
auto x92 = x75 + x91;
auto x93 = V{3}*cell[18];
auto x94 = V{3}*cell[9];
auto x95 = V{3}*cell[17] - V{3}*cell[8];
auto x96 = x59*(V{3}*cell[11] - V{3}*cell[2] - x77 + x79 + x81 + x93 - x94 + x95);
auto x97 = x47*x84;
auto x98 = x40 + V{-1};
auto x99 = x56 + x96 - x97 + x98;
auto x100 = x46*x59;
auto x101 = x100 + x20;
auto x102 = V{4.5}*(x101*x101);
auto x103 = V{3}*x20;
auto x104 = x103 + x96;
auto x105 = x104 + x75;
auto x106 = x59*(V{3}*cell[12] - V{3}*cell[3] - x78 + x80 + x82 - x93 + x94 + x95);
auto x107 = x55*x84;
auto x108 = x106 - x107 + x48 + x98;
auto x109 = x21 + x54*x59;
auto x110 = V{4.5}*(x109*x109);
auto x111 = V{3}*x21;
auto x112 = x106 + x111;
auto x113 = x112 + x75;
auto x114 = V{0.0277777777777778}*x22;
auto x115 = V{4.5}*x29;
auto x116 = cell[10] + x36;
auto x117 = cell[11] + V{2}*cell[13] - V{2}*cell[4] + x116 + x43 + x44 + x49;
auto x118 = x115*(x117*x117);
auto x119 = x58 + x83;
auto x120 = -x118 + x119 + x96;
auto x121 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x122 = x101 + x88;
auto x123 = V{4.5}*(x122*x122);
auto x124 = -x96;
auto x125 = -cell[17] + cell[8];
auto x126 = -cell[11] + V{2}*cell[14] + cell[2] - V{2}*cell[5] + x116 + x125 + x32 + x51;
auto x127 = -x126;
auto x128 = -x115*x127*x127 + x119 + x124;
auto x129 = -x100 - x20 + x88;
auto x130 = -x129;
auto x131 = -x103 + x124;
auto x132 = x34 + x37;
auto x133 = cell[12] + V{2}*cell[15] - V{2}*cell[6] + x132 + x42 + x52;
auto x134 = x115*(x133*x133);
auto x135 = x106 + x119 - x134;
auto x136 = x109 + x88;
auto x137 = V{4.5}*(x136*x136);
auto x138 = -x106;
auto x139 = -cell[12] + cell[3] + x31;
auto x140 = V{2}*cell[16] - V{2}*cell[7] + x125 + x132 + x139 + x43;
auto x141 = -x140;
auto x142 = -x115*x141*x141 + x119 + x138;
auto x143 = -x21 - x54*x59;
auto x144 = x143 + x88;
auto x145 = -x144;
auto x146 = -x111 + x138;
auto x147 = V{2}*cell[17] - V{2}*cell[8] + x33 + x45 + x50 + x53;
auto x148 = x115*(x147*x147);
auto x149 = x58 + x96;
auto x150 = x106 - x148 + x149;
auto x151 = x101 + x109;
auto x152 = V{4.5}*(x151*x151);
auto x153 = -cell[15] + V{2}*cell[18] + cell[6] - V{2}*cell[9] + x139 + x35 + x45;
auto x154 = -x153;
auto x155 = -x115*x154*x154 + x138 + x149;
auto x156 = x101 + x143;
auto x157 = -x156;
auto x158 = -x48;
auto x159 = V{1} - x56;
auto x160 = x158 + x159;
auto x161 = x160 + x83;
auto x162 = x161 + x85;
auto x163 = -x67 - x72 - x74;
auto x164 = x163 + V{1};
auto x165 = x164 + x91;
auto x166 = -x40;
auto x167 = x166 + x96;
auto x168 = x159 + x167 + x97;
auto x169 = x106 + V{1};
auto x170 = x107 + x158 + x166 + x169;
auto x171 = x111 + x163 + x169;
auto x172 = x118 + x161 + x167;
auto x173 = -x83;
auto x174 = -x115*x126*x126 + x149 + x173;
auto x175 = x173 - x90;
auto x176 = x106 + x134 + x161 + x166;
auto x177 = x106 + x58;
auto x178 = -x115*x140*x140 + x173 + x177;
auto x179 = x106 + x148 + x160 + x167;
auto x180 = -x115*x153*x153 + x124 + x177;
auto x0 = -cell[0]*x23 - V{0.333333333333333}*x22*(x27*x58 + V{1}) + V{0.333333333333333}*x27*x58 - x75*(V{0.333333333333333}*cell[0] + V{0.333333333333333}*cell[10] + V{0.333333333333333}*cell[11] + V{0.333333333333333}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] + V{0.333333333333333});
auto x1 = -cell[1]*x23 + V{0.0555555555555556}*x27*x86 - x76*(x27*x86 + V{1}) - x87*(-x89 + x92);
auto x2 = -cell[2]*x23 + V{0.0555555555555556}*x27*x99 - x76*(x27*x99 + V{1}) - x87*(-x102 + x105);
auto x3 = -cell[3]*x23 + V{0.0555555555555556}*x108*x27 - x76*(x108*x27 + V{1}) - x87*(-x110 + x113);
auto x4 = -cell[4]*x23 - x114*(x120*x27 + V{1}) + V{0.0277777777777778}*x120*x27 - x121*(x104 - x123 + x92);
auto x5 = -(cell[5]*x23 + x114*(x128*x27 + V{1}) + x121*(x131 + x92 - V{4.5}*x130*x130) - V{0.0277777777777778}*x128*x27);
auto x6 = -cell[6]*x23 - x114*(x135*x27 + V{1}) - x121*(x112 - x137 + x92) + V{0.0277777777777778}*x135*x27;
auto x7 = -(cell[7]*x23 + x114*(x142*x27 + V{1}) + x121*(x146 + x92 - V{4.5}*x145*x145) - V{0.0277777777777778}*x142*x27);
auto x8 = -cell[8]*x23 - x114*(x150*x27 + V{1}) - x121*(x105 + x112 - x152) + V{0.0277777777777778}*x150*x27;
auto x9 = -(cell[9]*x23 + x114*(x155*x27 + V{1}) + x121*(x105 + x146 - V{4.5}*x157*x157) - V{0.0277777777777778}*x155*x27);
auto x10 = -cell[10]*x23 - x162*x87 + x76*(x162*x27 + V{-1}) + x87*(x165 + x89);
auto x11 = -cell[11]*x23 - x168*x87 + x76*(x168*x27 + V{-1}) + x87*(x102 + x104 + x164);
auto x12 = -cell[12]*x23 - x170*x87 + x76*(x170*x27 + V{-1}) + x87*(x110 + x171);
auto x13 = -cell[13]*x23 + x114*(x172*x27 + V{-1}) - x121*x172 + x121*(x104 + x123 + x165);
auto x14 = -(cell[14]*x23 + x114*(x174*x27 + V{1}) + x121*(x105 + x175 - V{4.5}*x129*x129) - V{0.0277777777777778}*x174*x27);
auto x15 = -cell[15]*x23 + x114*(x176*x27 + V{-1}) - x121*x176 + x121*(x137 + x171 + x91);
auto x16 = -(cell[16]*x23 + x114*(x178*x27 + V{1}) + x121*(x113 + x175 - V{4.5}*x144*x144) - V{0.0277777777777778}*x178*x27);
auto x17 = -cell[17]*x23 + x114*(x179*x27 + V{-1}) - x121*x179 + x121*(x104 + x152 + x171);
auto x18 = -(cell[18]*x23 + x114*(x180*x27 + V{1}) + x121*(x113 + x131 - V{4.5}*x156*x156) - V{0.0277777777777778}*x180*x27);
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
return { x28, V{1}*x29*(x39 + x47 + x55) };
}
};

}

}
