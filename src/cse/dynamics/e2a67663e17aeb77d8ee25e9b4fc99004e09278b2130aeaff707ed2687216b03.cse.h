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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, forcing::Guo<momenta::Forced> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x21 = cell.template getFieldComponent<descriptors::FORCE>(2);
auto x19 = cell.template getFieldComponent<descriptors::FORCE>(0);
auto x22 = parameters.template get<descriptors::OMEGA>();
auto x23 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x20 = cell.template getFieldComponent<descriptors::FORCE>(1);
auto x24 = cell[14] + cell[5];
auto x25 = cell[18] + cell[9];
auto x26 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[15] + cell[16] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[7] + cell[8] + x24 + x25 + V{1};
auto x27 = V{1} / (x26);
auto x28 = V{0.333333333333333}*cell[13];
auto x29 = V{0.333333333333333}*cell[14];
auto x30 = V{0.333333333333333}*cell[4];
auto x31 = V{0.333333333333333}*cell[5];
auto x32 = V{1}*cell[12];
auto x33 = V{1}*cell[15];
auto x34 = V{1}*cell[17];
auto x35 = V{1}*cell[16];
auto x36 = V{1}*cell[18];
auto x37 = V{1}*cell[7];
auto x38 = V{1}*cell[6];
auto x39 = -x38;
auto x40 = x37 + x39;
auto x41 = V{1}*cell[9];
auto x42 = V{1}*cell[8];
auto x43 = -x42;
auto x44 = x41 + x43;
auto x45 = -V{1}*cell[3] + x32 + x33 + x34 - x35 - x36 + x40 + x44;
auto x46 = V{0.5}*x21 + x27*x45;
auto x47 = x46*x46;
auto x48 = V{0.333333333333333}*cell[0];
auto x49 = V{0.333333333333333}*cell[10];
auto x50 = V{0.333333333333333}*cell[1];
auto x51 = -V{0.666666666666667}*cell[17] - V{0.666666666666667}*cell[18] - V{0.666666666666667}*cell[8] - V{0.666666666666667}*cell[9] + x48 + x49 + x50;
auto x52 = V{0.333333333333333}*cell[11];
auto x53 = V{0.333333333333333}*cell[2];
auto x54 = -V{0.666666666666667}*cell[15] - V{0.666666666666667}*cell[16] - V{0.666666666666667}*cell[6] - V{0.666666666666667}*cell[7] + x52 + x53;
auto x55 = -V{0.666666666666667}*cell[12] - V{0.666666666666667}*cell[3] + x26*x47 + x28 + x29 + x30 + x31 + x51 + x54;
auto x56 = V{0.333333333333333}*cell[15];
auto x57 = V{0.333333333333333}*cell[16];
auto x58 = V{0.333333333333333}*cell[6];
auto x59 = V{0.333333333333333}*cell[7];
auto x60 = V{0.5}*x20;
auto x61 = V{1}*cell[11];
auto x62 = V{1}*cell[13];
auto x63 = V{1}*cell[14];
auto x64 = V{1}*cell[5];
auto x65 = V{1}*cell[4];
auto x66 = -x65;
auto x67 = x64 + x66;
auto x68 = x27*(-V{1}*cell[2] + x34 + x36 - x41 + x43 + x61 + x62 - x63 + x67);
auto x69 = x60 + x68;
auto x70 = x69*x69;
auto x71 = V{0.333333333333333}*cell[12];
auto x72 = V{0.333333333333333}*cell[3];
auto x73 = -V{0.666666666666667}*cell[13] - V{0.666666666666667}*cell[14] - V{0.666666666666667}*cell[4] - V{0.666666666666667}*cell[5] + x71 + x72;
auto x74 = -V{0.666666666666667}*cell[11] - V{0.666666666666667}*cell[2] + x26*x70 + x51 + x56 + x57 + x58 + x59 + x73;
auto x75 = V{0.333333333333333}*cell[17];
auto x76 = V{0.333333333333333}*cell[18];
auto x77 = V{0.333333333333333}*cell[8];
auto x78 = V{0.333333333333333}*cell[9];
auto x79 = V{1}*cell[10];
auto x80 = V{0.5}*x19 + x27*(-V{1}*cell[1] + x33 + x35 - x37 + x39 + x62 + x63 - x64 + x66 + x79);
auto x81 = x80*x80;
auto x82 = -V{0.666666666666667}*cell[10] - V{0.666666666666667}*cell[1] + x26*x81 + x48 + x54 + x73 + x75 + x76 + x77 + x78;
auto x83 = x26*x80;
auto x84 = -x33 + x35 + x40 + x46*x83;
auto x85 = x69*x83;
auto x86 = x26*x46*x69;
auto x87 = V{1} / (V{3.00000046417339}*util::sqrt(x27*(x23*x23)*util::sqrt((-cell[13] - cell[4] + x24 + x85)*(-x62 + x63 + x67 + x85) + (-cell[17] - cell[8] + x25 + x86)*(-x34 + x36 + x44 + x86) + V{0.5}*(x55*x55) + V{0.5}*(x74*x74) + V{0.5}*(x82*x82) + V{1}*(x84*x84)) + V{0.0277777691819762}/((x22)*(x22))) + V{0.5}/x22);
auto x88 = V{1.5}*x81;
auto x89 = V{1.5}*x70;
auto x90 = V{1.5}*x47;
auto x91 = x88 + x89 + x90 + V{-1};
auto x92 = V{1} - x87;
auto x93 = V{1.5}*x19;
auto x94 = V{3}*cell[14];
auto x95 = V{3}*cell[16];
auto x96 = V{3}*cell[5];
auto x97 = V{3}*cell[7];
auto x98 = V{3}*cell[13] - V{3}*cell[4];
auto x99 = V{3}*cell[15] - V{3}*cell[6];
auto x100 = x27*(V{3}*cell[10] - V{3}*cell[1] + x94 + x95 - x96 - x97 + x98 + x99);
auto x101 = x100 + x93;
auto x102 = x101*x19;
auto x103 = V{1.5}*x20;
auto x104 = V{3}*cell[18];
auto x105 = V{3}*cell[9];
auto x106 = V{3}*cell[17] - V{3}*cell[8];
auto x107 = x27*(V{3}*cell[11] - V{3}*cell[2] + x104 - x105 + x106 - x94 + x96 + x98);
auto x108 = x103 + x107;
auto x109 = x108*x20;
auto x110 = V{1.5}*x21;
auto x111 = x27*(V{3}*cell[12] - V{3}*cell[3] - x104 + x105 + x106 - x95 + x97 + x99);
auto x112 = x110 + x111;
auto x113 = x112*x21;
auto x114 = x109 + x113;
auto x115 = V{1} - V{0.5}*x87;
auto x116 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x117 = V{4.5}*cell[14];
auto x118 = V{4.5}*cell[16];
auto x119 = V{4.5}*cell[5];
auto x120 = V{4.5}*cell[7];
auto x121 = V{4.5}*cell[13] - V{4.5}*cell[4];
auto x122 = V{4.5}*cell[15] - V{4.5}*cell[6];
auto x123 = V{2.25}*x19 + x27*(V{4.5}*cell[10] - V{4.5}*cell[1] + x117 + x118 - x119 - x120 + x121 + x122);
auto x124 = x123*x80;
auto x125 = x101 + x91;
auto x126 = V{3}*x19;
auto x127 = V{6}*cell[14];
auto x128 = V{6}*cell[16];
auto x129 = V{6}*cell[5];
auto x130 = V{6}*cell[7];
auto x131 = V{6}*cell[13] - V{6}*cell[4];
auto x132 = V{6}*cell[15] - V{6}*cell[6];
auto x133 = x27*(V{6}*cell[10] - V{6}*cell[1] + x127 + x128 - x129 - x130 + x131 + x132);
auto x134 = x126 + x133;
auto x135 = x134 + V{-3};
auto x136 = V{0.0555555555555556}*x26;
auto x137 = V{2.25}*x20;
auto x138 = V{4.5}*cell[18];
auto x139 = V{4.5}*cell[9];
auto x140 = V{4.5}*cell[17] - V{4.5}*cell[8];
auto x141 = x27*(V{4.5}*cell[11] - V{4.5}*cell[2] - x117 + x119 + x121 + x138 - x139 + x140);
auto x142 = x137 + x141;
auto x143 = x142*x69;
auto x144 = x108 + x91;
auto x145 = V{3}*x20;
auto x146 = V{6}*cell[18];
auto x147 = V{6}*cell[9];
auto x148 = V{6}*cell[17] - V{6}*cell[8];
auto x149 = x27*(V{6}*cell[11] - V{6}*cell[2] - x127 + x129 + x131 + x146 - x147 + x148);
auto x150 = x145 + x149;
auto x151 = x150 + V{-3};
auto x152 = x102 + x113;
auto x153 = V{4.5}*cell[12] - V{4.5}*cell[3] - x118 + x120 + x122 - x138 + x139 + x140;
auto x154 = x153*x27 + V{2.25}*x21;
auto x155 = x154*x46;
auto x156 = x112 + x91;
auto x157 = V{3}*x21;
auto x158 = x27*(V{6}*cell[12] - V{6}*cell[3] - x128 + x130 + x132 - x146 + x147 + x148);
auto x159 = x157 + x158;
auto x160 = x159 + V{-3};
auto x161 = x102 + x109;
auto x162 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x163 = (x123 + x142)*(x69 + x80);
auto x164 = V{4.5}*x20;
auto x165 = V{9}*cell[18];
auto x166 = V{9}*cell[5];
auto x167 = V{9}*cell[14];
auto x168 = V{9}*cell[9];
auto x169 = V{9}*cell[13] - V{9}*cell[4];
auto x170 = V{9}*cell[17] - V{9}*cell[8];
auto x171 = x27*(V{9}*cell[11] - V{9}*cell[2] + x165 + x166 - x167 - x168 + x169 + x170);
auto x172 = x164 + x171;
auto x173 = V{4.5}*x19;
auto x174 = V{9}*cell[16];
auto x175 = V{9}*cell[7];
auto x176 = V{9}*cell[15] - V{9}*cell[6];
auto x177 = x27*(V{9}*cell[10] - V{9}*cell[1] - x166 + x167 + x169 + x174 - x175 + x176);
auto x178 = x173 + x177;
auto x179 = -x113;
auto x180 = V{0.0277777777777778}*x26;
auto x181 = -x60 - x68 + x80;
auto x182 = x123 - x137 - x141;
auto x183 = -x103 - x107;
auto x184 = x172 + V{3};
auto x185 = -x126 - x133;
auto x186 = x150 + V{3};
auto x187 = -x173 - x177;
auto x188 = (x123 + x154)*(x46 + x80);
auto x189 = V{4.5}*x21;
auto x190 = x27*(V{9}*cell[12] - V{9}*cell[3] - x165 + x168 + x170 - x174 + x175 + x176);
auto x191 = x189 + x190;
auto x192 = -x109;
auto x193 = -V{0.5}*x21 - x27*x45;
auto x194 = x193 + x80;
auto x195 = -x153*x27 - V{2.25}*x21;
auto x196 = x123 + x195;
auto x197 = -x110 - x111;
auto x198 = x191 + V{3};
auto x199 = x159 + V{3};
auto x200 = (x142 + x154)*(x46 + x69);
auto x201 = -x102;
auto x202 = x193 + x69;
auto x203 = x142 + x195;
auto x204 = -x145 - x149;
auto x205 = -x164 - x171;
auto x206 = -x88 - x89 - x90 + V{1};
auto x207 = x101 + x206;
auto x208 = x134 + V{3};
auto x209 = x108 + x206;
auto x210 = -x100 - x93;
auto x211 = x178 + V{3};
auto x212 = -x157 - x158;
auto x213 = -x189 - x190;
auto x0 = V{1}*cell[0]*x92 - V{0.333333333333333}*x115*x26*(x102 + x114) - x87*(x91*(x28 + x29 + x30 + x31 + x48 + x49 + x50 + x52 + x53 + x56 + x57 + x58 + x59 + x71 + x72 + x75 + x76 + x77 + x78 + V{0.333333333333333}) + V{0.333333333333333});
auto x1 = V{1}*cell[1]*x92 - x115*x136*(x114 - x135*x19) - x87*(x116*(-x124 + x125) + V{0.0555555555555556});
auto x2 = V{1}*cell[2]*x92 - x115*x136*(-x151*x20 + x152) - x87*(x116*(-x143 + x144) + V{0.0555555555555556});
auto x3 = V{1}*cell[3]*x92 - x115*x136*(-x160*x21 + x161) - x87*(x116*(-x155 + x156) + V{0.0555555555555556});
auto x4 = x115*x180*(x179 + x19*(x135 + x172) + x20*(x151 + x178)) + x65*x92 - x87*(x162*(x108 + x125 - x163) + V{0.0277777777777778});
auto x5 = V{1}*cell[5]*x92 - x115*x180*(x113 + x19*(x184 + x185) - x20*(x186 + x187)) - x87*(x162*(x125 - x181*x182 + x183) + V{0.0277777777777778});
auto x6 = x115*x180*(x19*(x135 + x191) + x192 + x21*(x160 + x178)) + x38*x92 - x87*(x162*(x112 + x125 - x188) + V{0.0277777777777778});
auto x7 = V{1}*cell[7]*x92 - x115*x180*(x109 + x19*(x185 + x198) - x21*(x187 + x199)) - x87*(x162*(x125 - x194*x196 + x197) + V{0.0277777777777778});
auto x8 = x115*x180*(x20*(x151 + x191) + x201 + x21*(x160 + x172)) + x42*x92 - x87*(x162*(x112 + x144 - x200) + V{0.0277777777777778});
auto x9 = V{1}*cell[9]*x92 - x115*x180*(x102 + x20*(x198 + x204) - x21*(x199 + x205)) - x87*(x162*(x144 + x197 - x202*x203) + V{0.0277777777777778});
auto x10 = -x115*x136*(x114 - x19*x208) + x79*x92 + x87*(x116*(x124 + x207) + V{-0.0555555555555556});
auto x11 = -x115*x136*(x152 - x186*x20) + x61*x92 + x87*(x116*(x143 + x209) + V{-0.0555555555555556});
auto x12 = -x115*x136*(x161 - x199*x21) + x32*x92 + x87*(x116*(x112 + x155 + x206) + V{-0.0555555555555556});
auto x13 = x115*x180*(x179 + x19*(x172 + x208) + x20*(x178 + x186)) + x62*x92 + x87*(x162*(x108 + x163 + x207) + V{-0.0277777777777778});
auto x14 = V{1}*cell[14]*x92 - x115*x180*(x113 - x19*(x205 + x208) + x20*(x204 + x211)) - x87*(x162*(x144 - x181*x182 + x210) + V{0.0277777777777778});
auto x15 = x115*x180*(x19*(x191 + x208) + x192 + x21*(x178 + x199)) + x33*x92 + x87*(x162*(x112 + x188 + x207) + V{-0.0277777777777778});
auto x16 = V{1}*cell[16]*x92 - x115*x180*(x109 - x19*(x208 + x213) + x21*(x211 + x212)) - x87*(x162*(x156 - x194*x196 + x210) + V{0.0277777777777778});
auto x17 = x115*x180*(x20*(x186 + x191) + x201 + x21*(x172 + x199)) + x34*x92 + x87*(x162*(x112 + x200 + x209) + V{-0.0277777777777778});
auto x18 = V{1}*cell[18]*x92 - x115*x180*(x102 - x20*(x186 + x213) + x21*(x184 + x212)) - x87*(x162*(x156 + x183 - x202*x203) + V{0.0277777777777778});
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
return { x26, x47 + x70 + x81 };
}
};

}

}
