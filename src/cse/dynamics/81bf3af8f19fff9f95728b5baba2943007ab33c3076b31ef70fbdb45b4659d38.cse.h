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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::Incompressible, collision::BGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = cell[13] - cell[4];
auto x21 = cell[15] - cell[6];
auto x22 = cell[10] + cell[14] + cell[16];
auto x23 = -cell[1] - cell[5] - cell[7] + x20 + x21 + x22;
auto x24 = x23*x23;
auto x25 = cell[17] - cell[8];
auto x26 = cell[11] + cell[18] + cell[5];
auto x27 = -cell[14] - cell[2] - cell[9] + x20 + x25 + x26;
auto x28 = x27*x27;
auto x29 = cell[12] + cell[7] + cell[9];
auto x30 = -cell[16] - cell[18] - cell[3] + x21 + x25 + x29;
auto x31 = x30*x30;
auto x32 = x19 + V{-1};
auto x33 = V{0.111111111111111}*cell[10];
auto x34 = V{0.0555555555555556}*cell[0];
auto x35 = V{0.0555555555555556}*cell[12];
auto x36 = V{0.0555555555555556}*cell[3];
auto x37 = V{0.111111111111111}*cell[13];
auto x38 = -V{0.0833333333333333}*x31;
auto x39 = V{0.222222222222222}*cell[4] + x34 + x35 + x36 - x37 + x38;
auto x40 = V{0.0555555555555556}*cell[11];
auto x41 = V{0.0555555555555556}*cell[2];
auto x42 = V{0.111111111111111}*cell[15];
auto x43 = -V{0.0833333333333333}*x28;
auto x44 = V{0.222222222222222}*cell[6] + x40 + x41 - x42 + x43;
auto x45 = V{0.111111111111111}*cell[14];
auto x46 = V{0.222222222222222}*cell[5] - x45;
auto x47 = V{0.111111111111111}*cell[16];
auto x48 = V{0.222222222222222}*cell[7] - x47;
auto x49 = V{0.0555555555555556}*cell[17];
auto x50 = V{0.0555555555555556}*cell[18];
auto x51 = V{0.0555555555555556}*cell[8];
auto x52 = V{0.0555555555555556}*cell[9];
auto x53 = V{0.166666666666667}*x24 + x49 + x50 + x51 + x52;
auto x54 = V{0.111111111111111}*cell[11];
auto x55 = V{0.111111111111111}*cell[5];
auto x56 = V{0.222222222222222}*cell[14] - x55;
auto x57 = V{0.0555555555555556}*cell[10];
auto x58 = V{0.0555555555555556}*cell[1];
auto x59 = V{0.111111111111111}*cell[17];
auto x60 = -V{0.0833333333333333}*x24;
auto x61 = V{0.222222222222222}*cell[8] + x57 + x58 - x59 + x60;
auto x62 = V{0.111111111111111}*cell[18];
auto x63 = V{0.222222222222222}*cell[9] - x62;
auto x64 = V{0.0555555555555556}*cell[15];
auto x65 = V{0.0555555555555556}*cell[16];
auto x66 = V{0.0555555555555556}*cell[6];
auto x67 = V{0.0555555555555556}*cell[7];
auto x68 = V{0.166666666666667}*x28 + x64 + x65 + x66 + x67;
auto x69 = V{0.111111111111111}*cell[12];
auto x70 = V{0.111111111111111}*cell[7];
auto x71 = V{0.222222222222222}*cell[16] + x34 - x70;
auto x72 = V{0.111111111111111}*cell[9];
auto x73 = V{0.222222222222222}*cell[18] - x72;
auto x74 = V{0.0555555555555556}*cell[13];
auto x75 = V{0.0555555555555556}*cell[14];
auto x76 = V{0.0555555555555556}*cell[4];
auto x77 = V{0.0555555555555556}*cell[5];
auto x78 = V{0.166666666666667}*x31 + x74 + x75 + x76 + x77;
auto x79 = V{0.0277777777777778}*cell[0];
auto x80 = -V{0.0416666666666667}*x24;
auto x81 = -V{0.0416666666666667}*x28;
auto x82 = -V{0.0416666666666667}*x31;
auto x83 = -x50 + x72 + x79 + x80 + x81 + x82;
auto x84 = -x65 + x70;
auto x85 = V{0.111111111111111}*cell[1];
auto x86 = -x57;
auto x87 = V{0.111111111111111}*cell[8];
auto x88 = -x49 + x87;
auto x89 = x85 + x86 + x88;
auto x90 = V{0.111111111111111}*cell[2];
auto x91 = -x40;
auto x92 = V{0.111111111111111}*cell[6];
auto x93 = -x64 + x92;
auto x94 = x90 + x91 + x93;
auto x95 = V{0.0277777777777778}*cell[12];
auto x96 = V{0.0277777777777778}*cell[3];
auto x97 = V{1}*cell[10];
auto x98 = V{1}*cell[17];
auto x99 = -V{1}*cell[1];
auto x100 = V{1}*cell[8];
auto x101 = -x100 + x97 + x98 + x99;
auto x102 = V{1}*cell[18];
auto x103 = V{1}*cell[9];
auto x104 = x102 - x103;
auto x105 = V{1}*cell[11];
auto x106 = V{1}*cell[15];
auto x107 = V{1}*cell[2];
auto x108 = -x107;
auto x109 = V{1}*cell[6];
auto x110 = -x109;
auto x111 = x105 + x106 + x108 + x110;
auto x112 = V{1}*cell[16];
auto x113 = V{1}*cell[7];
auto x114 = x112 - x113;
auto x115 = V{4.5}*cell[10];
auto x116 = V{4.5}*cell[17];
auto x117 = -V{4.5}*cell[1];
auto x118 = V{4.5}*cell[8];
auto x119 = x115 + x116 + x117 - x118;
auto x120 = V{4.5}*cell[18];
auto x121 = V{4.5}*cell[9];
auto x122 = x120 - x121;
auto x123 = V{4.5}*cell[11];
auto x124 = V{4.5}*cell[15];
auto x125 = V{4.5}*cell[2];
auto x126 = -x125;
auto x127 = V{4.5}*cell[6];
auto x128 = -x127;
auto x129 = x123 + x124 + x126 + x128;
auto x130 = V{4.5}*cell[16];
auto x131 = V{4.5}*cell[7];
auto x132 = x130 - x131;
auto x133 = V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[5] + x95 + x96 + V{0.0277777777777778}*(V{2}*cell[13] - V{2}*cell[4] + x101 + x104 + x111 + x114)*(V{9}*cell[13] - V{9}*cell[4] + x119 + x122 + x129 + x132);
auto x134 = -x51 + x59;
auto x135 = x134 + x85 + x86;
auto x136 = -x102 + x103;
auto x137 = x100 + x97 - x98 + x99;
auto x138 = -x120 + x121;
auto x139 = x115 - x116 + x117 + x118;
auto x140 = V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[4] + x95 + x96 + V{0.0277777777777778}*(V{2}*cell[14] - V{2}*cell[5] - x105 + x106 + x107 + x110 + x114 + x136 + x137)*(V{9}*cell[14] - V{9}*cell[5] - x123 + x124 + x125 + x128 + x132 + x138 + x139);
auto x141 = -x52 + x62;
auto x142 = -x41;
auto x143 = x79 + x80 + x81 + x82;
auto x144 = x143 + x84;
auto x145 = x142 + x144 + x54;
auto x146 = x141 + x143;
auto x147 = x55 - x75;
auto x148 = V{0.111111111111111}*cell[3];
auto x149 = -x35;
auto x150 = V{0.111111111111111}*cell[4];
auto x151 = x150 - x74;
auto x152 = x148 + x149 + x151;
auto x153 = V{0.0277777777777778}*cell[11];
auto x154 = V{0.0277777777777778}*cell[2];
auto x155 = V{1}*cell[14];
auto x156 = V{1}*cell[5];
auto x157 = V{1}*cell[13] - V{1}*cell[4];
auto x158 = x155 - x156 + x157;
auto x159 = V{1}*cell[12];
auto x160 = V{1}*cell[3];
auto x161 = x159 - x160;
auto x162 = V{4.5}*cell[14];
auto x163 = V{4.5}*cell[5];
auto x164 = V{4.5}*cell[13] - V{4.5}*cell[4];
auto x165 = x162 - x163 + x164;
auto x166 = V{4.5}*cell[12];
auto x167 = V{4.5}*cell[3];
auto x168 = x166 - x167;
auto x169 = V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[7] + x153 + x154 + V{0.0277777777777778}*(V{2}*cell[15] - V{2}*cell[6] + x101 + x136 + x158 + x161)*(V{9}*cell[15] - V{9}*cell[6] + x119 + x138 + x165 + x168);
auto x170 = -x36;
auto x171 = x151 + x170 + x69;
auto x172 = -x159 + x160;
auto x173 = -x166 + x167;
auto x174 = V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[6] + x153 + x154 + V{0.0277777777777778}*(V{2}*cell[16] - V{2}*cell[7] + x104 + x137 + x158 + x172)*(V{9}*cell[16] - V{9}*cell[7] + x122 + x139 + x165 + x173);
auto x175 = V{0.0277777777777778}*cell[10];
auto x176 = V{0.0277777777777778}*cell[1];
auto x177 = x45 - x77;
auto x178 = x175 + x176 + x177;
auto x179 = -x155 + x156 + x157;
auto x180 = -x162 + x163 + x164;
auto x181 = V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778}*(V{2}*cell[17] - V{2}*cell[8] + x111 - x112 + x113 + x161 + x179)*(V{9}*cell[17] - V{9}*cell[8] + x129 - x130 + x131 + x168 + x180);
auto x182 = x47 - x67;
auto x183 = x143 + x182;
auto x184 = V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*(V{2}*cell[18] - V{2}*cell[9] + x105 - x106 + x108 + x109 + x114 + x172 + x179)*(V{9}*cell[18] - V{9}*cell[9] + x123 - x124 + x126 + x127 + x132 + x173 + x180);
auto x185 = x42 - x66;
auto x186 = x185 + x90 + x91;
auto x187 = V{0.222222222222222}*cell[13] - x150 + x35 + x36 + x38;
auto x188 = V{0.222222222222222}*cell[15] + x40 + x41 + x43 - x92;
auto x189 = V{0.222222222222222}*cell[17] + x34 + x57 + x58 + x60 - x87;
auto x190 = x142 + x54;
auto x191 = -x58;
auto x192 = x182 + x191 + x33;
auto x193 = x37 - x76;
auto x194 = x170 + x193 + x69;
auto x195 = x177 + x191 + x33;
auto x196 = x148 + x149 + x193;
auto x197 = x147 + x175 + x176;
auto x198 = cell[0] + cell[13] + cell[15] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + x22 + x26 + x29 + V{1};
auto x0 = -cell[0]*x32 + x19*(V{0.333333333333333}*cell[0] + V{0.333333333333333}*cell[10] + V{0.333333333333333}*cell[11] + V{0.333333333333333}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] - V{0.5}*x24 - V{0.5}*x28 - V{0.5}*x31);
auto x1 = -cell[1]*x32 + x19*(V{0.222222222222222}*cell[1] - x33 + x39 + x44 + x46 + x48 + x53);
auto x2 = -cell[2]*x32 + x19*(V{0.222222222222222}*cell[2] + x39 - x54 + x56 + x61 + x63 + x68);
auto x3 = -cell[3]*x32 + x19*(V{0.222222222222222}*cell[3] + x44 + x61 - x69 + x71 + x73 + x78);
auto x4 = -cell[4]*x32 + x19*(-V{0.138888888888889}*cell[13] + V{0.194444444444444}*cell[4] + x133 + x83 + x84 + x89 + x94);
auto x5 = -cell[5]*x32 + x19*(-V{0.138888888888889}*cell[14] + V{0.194444444444444}*cell[5] + x135 + x140 + x141 + x145 + x93);
auto x6 = -cell[6]*x32 + x19*(-V{0.138888888888889}*cell[15] + V{0.194444444444444}*cell[6] + x146 + x147 + x152 + x169 + x89);
auto x7 = -cell[7]*x32 + x19*(-V{0.138888888888889}*cell[16] + V{0.194444444444444}*cell[7] + x135 + x147 + x171 + x174 + x83);
auto x8 = -cell[8]*x32 + x19*(-V{0.138888888888889}*cell[17] + V{0.194444444444444}*cell[8] + x152 + x178 + x181 + x183 + x94);
auto x9 = -cell[9]*x32 + x19*(-V{0.138888888888889}*cell[18] + V{0.194444444444444}*cell[9] + x144 + x171 + x178 + x184 + x186);
auto x10 = -cell[10]*x32 + x19*(V{0.222222222222222}*cell[10] + x187 + x188 + x53 + x56 + x71 - x85);
auto x11 = -cell[11]*x32 + x19*(V{0.222222222222222}*cell[11] + x187 + x189 + x46 + x68 + x73 - x90);
auto x12 = -cell[12]*x32 + x19*(V{0.222222222222222}*cell[12] - x148 + x188 + x189 + x48 + x63 + x78);
auto x13 = -cell[13]*x32 + x19*(V{0.194444444444444}*cell[13] - V{0.138888888888889}*cell[4] + x133 + x134 + x146 + x185 + x190 + x192);
auto x14 = -cell[14]*x32 + x19*(V{0.194444444444444}*cell[14] - V{0.138888888888889}*cell[5] + x140 + x186 + x192 + x83 + x88);
auto x15 = -cell[15]*x32 + x19*(V{0.194444444444444}*cell[15] - V{0.138888888888889}*cell[6] + x134 + x169 + x194 + x195 + x83);
auto x16 = -cell[16]*x32 + x19*(V{0.194444444444444}*cell[16] - V{0.138888888888889}*cell[7] + x146 + x174 + x195 + x196 + x88);
auto x17 = -cell[17]*x32 + x19*(V{0.194444444444444}*cell[17] - V{0.138888888888889}*cell[8] + x145 + x181 + x185 + x194 + x197);
auto x18 = -cell[18]*x32 + x19*(V{0.194444444444444}*cell[18] - V{0.138888888888889}*cell[9] + x183 + x184 + x190 + x196 + x197 + x93);
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
return { x198, V{1}*(x24 + x28 + x31)/((x198)*(x198)) };
}
};

}

}
