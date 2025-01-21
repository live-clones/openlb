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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::Guo<momenta::Forced> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x22 = x21 + V{1};
auto x23 = x21 + V{1};
auto x24 = V{1} / (x23);
auto x25 = V{1}*cell[14];
auto x26 = V{1}*cell[16];
auto x27 = V{1}*cell[5];
auto x28 = V{1}*cell[7];
auto x29 = V{1}*cell[13] - V{1}*cell[4];
auto x30 = V{1}*cell[15] - V{1}*cell[6];
auto x31 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(0) + x24*(V{1}*cell[10] - V{1}*cell[1] + x25 + x26 - x27 - x28 + x29 + x30);
auto x32 = x31*x31;
auto x33 = V{1.5}*x32;
auto x34 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x35 = V{1}*cell[18];
auto x36 = V{1}*cell[9];
auto x37 = V{1}*cell[17] - V{1}*cell[8];
auto x38 = x24*(V{1}*cell[11] - V{1}*cell[2] - x25 + x27 + x29 + x35 - x36 + x37);
auto x39 = x34 + x38;
auto x40 = x39*x39;
auto x41 = V{1.5}*x40;
auto x42 = V{1}*cell[12] - V{1}*cell[3] - x26 + x28 + x30 - x35 + x36 + x37;
auto x43 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(2) + x24*x42;
auto x44 = x43*x43;
auto x45 = V{1.5}*x44;
auto x46 = x33 + x41 + x45 + V{-1};
auto x47 = V{0.5}*x19 + V{-1};
auto x48 = cell.template getFieldComponent<descriptors::FORCE>(0)*x31;
auto x49 = cell.template getFieldComponent<descriptors::FORCE>(1)*x39;
auto x50 = cell.template getFieldComponent<descriptors::FORCE>(2)*x43;
auto x51 = V{0.0555555555555556}*x19;
auto x52 = V{4.5}*cell[14];
auto x53 = V{4.5}*cell[16];
auto x54 = V{4.5}*cell[5];
auto x55 = V{4.5}*cell[7];
auto x56 = V{4.5}*cell[13] - V{4.5}*cell[4];
auto x57 = V{4.5}*cell[15] - V{4.5}*cell[6];
auto x58 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(0) + x24*(V{4.5}*cell[10] - V{4.5}*cell[1] + x52 + x53 - x54 - x55 + x56 + x57);
auto x59 = x31*x58;
auto x60 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x61 = V{3}*cell[14];
auto x62 = V{3}*cell[16];
auto x63 = V{3}*cell[5];
auto x64 = V{3}*cell[7];
auto x65 = V{3}*cell[13] - V{3}*cell[4];
auto x66 = V{3}*cell[15] - V{3}*cell[6];
auto x67 = x24*(V{3}*cell[10] - V{3}*cell[1] + x61 + x62 - x63 - x64 + x65 + x66);
auto x68 = x60 + x67;
auto x69 = x46 + x68;
auto x70 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x71 = V{6}*cell[14];
auto x72 = V{6}*cell[16];
auto x73 = V{6}*cell[5];
auto x74 = V{6}*cell[7];
auto x75 = V{6}*cell[13] - V{6}*cell[4];
auto x76 = V{6}*cell[15] - V{6}*cell[6];
auto x77 = x24*(V{6}*cell[10] - V{6}*cell[1] + x71 + x72 - x73 - x74 + x75 + x76);
auto x78 = x70 + x77;
auto x79 = x78 + V{-3};
auto x80 = V{0.055556}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x81 = V{0.166667}*x49;
auto x82 = V{0.166667}*x50;
auto x83 = x81 + x82;
auto x84 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x85 = V{4.5}*cell[18];
auto x86 = V{4.5}*cell[9];
auto x87 = V{4.5}*cell[17] - V{4.5}*cell[8];
auto x88 = x24*(V{4.5}*cell[11] - V{4.5}*cell[2] - x52 + x54 + x56 + x85 - x86 + x87);
auto x89 = x84 + x88;
auto x90 = x39*x89;
auto x91 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x92 = V{3}*cell[18];
auto x93 = V{3}*cell[9];
auto x94 = V{3}*cell[17] - V{3}*cell[8];
auto x95 = x24*(V{3}*cell[11] - V{3}*cell[2] - x61 + x63 + x65 + x92 - x93 + x94);
auto x96 = x91 + x95;
auto x97 = x46 + x96;
auto x98 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x99 = V{6}*cell[18];
auto x100 = V{6}*cell[9];
auto x101 = V{6}*cell[17] - V{6}*cell[8];
auto x102 = x24*(V{6}*cell[11] - V{6}*cell[2] - x100 + x101 - x71 + x73 + x75 + x99);
auto x103 = x102 + x98;
auto x104 = x103 + V{-3};
auto x105 = V{0.055556}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x106 = V{0.166667}*x48;
auto x107 = x106 + x82;
auto x108 = V{4.5}*cell[12] - V{4.5}*cell[3] - x53 + x55 + x57 - x85 + x86 + x87;
auto x109 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(2) + x108*x24;
auto x110 = x109*x43;
auto x111 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x112 = x24*(V{3}*cell[12] - V{3}*cell[3] - x62 + x64 + x66 - x92 + x93 + x94);
auto x113 = x111 + x112;
auto x114 = x113 + x46;
auto x115 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x116 = x24*(V{6}*cell[12] - V{6}*cell[3] + x100 + x101 - x72 + x74 + x76 - x99);
auto x117 = x115 + x116;
auto x118 = x117 + V{-3};
auto x119 = V{0.055556}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x120 = x106 + x81;
auto x121 = V{0.0277777777777778}*x19;
auto x122 = (x31 + x39)*(x58 + x89);
auto x123 = V{9}*cell[18];
auto x124 = V{9}*cell[5];
auto x125 = V{9}*cell[14];
auto x126 = V{9}*cell[9];
auto x127 = V{9}*cell[13] - V{9}*cell[4];
auto x128 = V{9}*cell[17] - V{9}*cell[8];
auto x129 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(1) + x24*(V{9}*cell[11] - V{9}*cell[2] + x123 + x124 - x125 - x126 + x127 + x128);
auto x130 = V{0.027778}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x131 = V{9}*cell[16];
auto x132 = V{9}*cell[7];
auto x133 = V{9}*cell[15] - V{9}*cell[6];
auto x134 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(0) + x24*(V{9}*cell[10] - V{9}*cell[1] - x124 + x125 + x127 + x131 - x132 + x133);
auto x135 = V{0.027778}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x136 = V{0.083333}*x50;
auto x137 = -x136;
auto x138 = x23*x47;
auto x139 = V{1}*x138;
auto x140 = x31 - x34 - x38;
auto x141 = x58 - x84 - x88;
auto x142 = -x91 - x95;
auto x143 = V{0.250002}*cell[14];
auto x144 = V{0.250002}*cell[16];
auto x145 = V{0.250002}*cell[5];
auto x146 = V{0.250002}*cell[7];
auto x147 = V{0.250002}*cell[13] - V{0.250002}*cell[4];
auto x148 = V{0.250002}*cell[15] - V{0.250002}*cell[6];
auto x149 = -V{0.125001}*cell.template getFieldComponent<descriptors::FORCE>(0) - x24*(V{0.250002}*cell[10] - V{0.250002}*cell[1] + x143 + x144 - x145 - x146 + x147 + x148) + V{0.083334};
auto x150 = V{0.166668}*cell[18];
auto x151 = V{0.166668}*cell[5];
auto x152 = V{0.166668}*cell[14];
auto x153 = V{0.166668}*cell[9];
auto x154 = V{0.166668}*cell[17] - V{0.166668}*cell[8];
auto x155 = V{0.166668}*cell[13] - V{0.166668}*cell[4];
auto x156 = V{0.083334}*cell.template getFieldComponent<descriptors::FORCE>(1) + x24*(V{0.166668}*cell[11] - V{0.166668}*cell[2] + x150 + x151 - x152 - x153 + x154 + x155);
auto x157 = x129 + V{3};
auto x158 = -x70 - x77;
auto x159 = (x109 + x58)*(x31 + x43);
auto x160 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(2) + x24*(V{9}*cell[12] - V{9}*cell[3] - x123 + x126 + x128 - x131 + x132 + x133);
auto x161 = V{0.027778}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x162 = V{0.083333}*x49;
auto x163 = -x162;
auto x164 = -V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(2) - x24*x42;
auto x165 = x164 + x31;
auto x166 = -V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(2) - x108*x24;
auto x167 = x166 + x58;
auto x168 = -x111 - x112;
auto x169 = V{0.166668}*cell[7];
auto x170 = V{0.166668}*cell[16];
auto x171 = V{0.166668}*cell[15] - V{0.166668}*cell[6];
auto x172 = V{0.083334}*cell.template getFieldComponent<descriptors::FORCE>(2) + x24*(V{0.166668}*cell[12] - V{0.166668}*cell[3] - x150 + x153 + x154 + x169 - x170 + x171);
auto x173 = x160 + V{3};
auto x174 = (x109 + x89)*(x39 + x43);
auto x175 = V{0.083333}*x48;
auto x176 = -x175;
auto x177 = x164 + x39;
auto x178 = x166 + x89;
auto x179 = V{0.250002}*cell[18];
auto x180 = V{0.250002}*cell[9];
auto x181 = V{0.250002}*cell[17] - V{0.250002}*cell[8];
auto x182 = -V{0.125001}*cell.template getFieldComponent<descriptors::FORCE>(1) - x24*(V{0.250002}*cell[11] - V{0.250002}*cell[2] - x143 + x145 + x147 + x179 - x180 + x181) + V{0.083334};
auto x183 = -x102 - x98;
auto x184 = -x33 - x41 - x45 + V{1};
auto x185 = x184 + x68;
auto x186 = x78 + V{3};
auto x187 = x184 + x96;
auto x188 = x103 + V{3};
auto x189 = x117 + V{3};
auto x190 = -x60 - x67;
auto x191 = V{0.083334}*cell.template getFieldComponent<descriptors::FORCE>(0) + x24*(V{0.166668}*cell[10] - V{0.166668}*cell[1] - x151 + x152 + x155 - x169 + x170 + x171);
auto x192 = x134 + V{3};
auto x193 = -V{0.125001}*cell.template getFieldComponent<descriptors::FORCE>(2) - x24*(V{0.250002}*cell[12] - V{0.250002}*cell[3] - x144 + x146 + x148 - x179 + x180 + x181) + V{0.083334};
auto x194 = -x115 - x116;
auto x0 = -cell[0]*x20 - V{0.333333333333333}*x19*(x22*x46 + V{1}) + V{1}*x23*x47*(x48 + x49 + x50);
auto x1 = -cell[1]*x20 + x23*x47*(-x79*x80 + x83) - x51*(x22*(-x59 + x69) + V{1});
auto x2 = -cell[2]*x20 + x23*x47*(-x104*x105 + x107) - x51*(x22*(-x90 + x97) + V{1});
auto x3 = -cell[3]*x20 + x23*x47*(-x118*x119 + x120) - x51*(x22*(-x110 + x114) + V{1});
auto x4 = -cell[4]*x20 - x121*(x22*(-x122 + x69 + x96) + V{1}) - x139*(x130*(x129 + x79) + x135*(x104 + x134) + x137);
auto x5 = -cell[5]*x20 - x121*(x22*(-x140*x141 + x142 + x69) + V{1}) + V{1}*x23*x47*(-cell.template getFieldComponent<descriptors::FORCE>(1)*(x149 + x156) + x130*(x157 + x158) + x136);
auto x6 = -cell[6]*x20 - x121*(x22*(x113 - x159 + x69) + V{1}) - x139*(x130*(x160 + x79) + x161*(x118 + x134) + x163);
auto x7 = -cell[7]*x20 - x121*(x22*(-x165*x167 + x168 + x69) + V{1}) + V{1}*x23*x47*(-cell.template getFieldComponent<descriptors::FORCE>(2)*(x149 + x172) + x130*(x158 + x173) + x162);
auto x8 = -cell[8]*x20 - x121*(x22*(x113 - x174 + x97) + V{1}) - x139*(x135*(x104 + x160) + x161*(x118 + x129) + x176);
auto x9 = -cell[9]*x20 - x121*(x22*(x168 - x177*x178 + x97) + V{1}) + V{1}*x23*x47*(-cell.template getFieldComponent<descriptors::FORCE>(2)*(x172 + x182) + x135*(x173 + x183) + x175);
auto x10 = -cell[10]*x20 + x138*(-x186*x80 + x83) + x51*(x22*(x185 + x59) + V{-1});
auto x11 = -cell[11]*x20 + x138*(-x105*x188 + x107) + x51*(x22*(x187 + x90) + V{-1});
auto x12 = -cell[12]*x20 + x138*(-x119*x189 + x120) + x51*(x22*(x110 + x113 + x184) + V{-1});
auto x13 = -cell[13]*x20 - x139*(x130*(x129 + x186) + x135*(x134 + x188) + x137) + V{0.0277777777777778}*x19*(x22*(x122 + x185 + x96) + V{-1});
auto x14 = -cell[14]*x20 - x121*(x22*(-x140*x141 + x190 + x97) + V{1}) + V{1}*x23*x47*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x182 + x191) + x135*(x183 + x192) + x136);
auto x15 = -cell[15]*x20 - x139*(x130*(x160 + x186) + x161*(x134 + x189) + x163) + V{0.0277777777777778}*x19*(x22*(x113 + x159 + x185) + V{-1});
auto x16 = -cell[16]*x20 - x121*(x22*(x114 - x165*x167 + x190) + V{1}) + V{1}*x23*x47*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x191 + x193) + x161*(x192 + x194) + x162);
auto x17 = -cell[17]*x20 - x139*(x135*(x160 + x188) + x161*(x129 + x189) + x176) + V{0.0277777777777778}*x19*(x22*(x113 + x174 + x187) + V{-1});
auto x18 = -cell[18]*x20 - x121*(x22*(x114 + x142 - x177*x178) + V{1}) + V{1}*x23*x47*(-cell.template getFieldComponent<descriptors::FORCE>(1)*(x156 + x193) + x161*(x157 + x194) + x175);
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
return { x23, x32 + x40 + x44 };
}
};

}

}
