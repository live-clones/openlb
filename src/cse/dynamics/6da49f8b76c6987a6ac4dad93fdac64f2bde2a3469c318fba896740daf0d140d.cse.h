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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::TRT, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = parameters.template get<collision::TRT::MAGIC>();
auto x21 = cell[10] + cell[14];
auto x22 = cell[12] + cell[7];
auto x23 = cell[0] + cell[11] + cell[13] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[8] + cell[9] + x21 + x22 + V{1};
auto x24 = V{1} / ((x23)*(x23));
auto x25 = V{1.5}*x24;
auto x26 = cell[13] - cell[4];
auto x27 = cell[15] - cell[6];
auto x28 = x26 + x27;
auto x29 = -cell[1];
auto x30 = cell[16] - cell[7];
auto x31 = x29 + x30;
auto x32 = -cell[5] + x21;
auto x33 = x28 + x31 + x32;
auto x34 = x33*x33;
auto x35 = x25*x34;
auto x36 = cell[17] - cell[8];
auto x37 = x26 + x36;
auto x38 = cell[18] - cell[9];
auto x39 = -cell[14];
auto x40 = -cell[2];
auto x41 = cell[11] + cell[5] + x39 + x40;
auto x42 = x37 + x38 + x41;
auto x43 = x42*x42;
auto x44 = x25*x43;
auto x45 = x27 + x36;
auto x46 = -cell[3];
auto x47 = -cell[18];
auto x48 = cell[9] + x47;
auto x49 = x46 + x48;
auto x50 = -cell[16];
auto x51 = x22 + x50;
auto x52 = x45 + x49 + x51;
auto x53 = x52*x52;
auto x54 = x25*x53;
auto x55 = x44 + x54 + V{-1};
auto x56 = x35 + x55;
auto x57 = V{1} / (-x20/(V{0.5} - V{1}/x19) + V{0.5});
auto x58 = V{0.5}*cell[1];
auto x59 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x60 = V{1} / (x23);
auto x61 = V{3}*cell[14];
auto x62 = V{3}*cell[16];
auto x63 = V{3}*cell[5];
auto x64 = V{3}*cell[7];
auto x65 = V{3}*cell[13] - V{3}*cell[4];
auto x66 = V{3}*cell[15] - V{3}*cell[6];
auto x67 = x60*(V{3}*cell[10] - V{3}*cell[1] + x61 + x62 - x63 - x64 + x65 + x66);
auto x68 = -V{3}*x24*x34 + x55 + x67;
auto x69 = V{3}*x24;
auto x70 = -x44;
auto x71 = V{1} - x54;
auto x72 = x70 + x71;
auto x73 = x67 + x72;
auto x74 = x59*(x34*x69 + x73);
auto x75 = V{0.5}*cell[10] - x74;
auto x76 = x58 + x59*x68;
auto x77 = x19*(x75 + x76 + V{0.0555555555555556});
auto x78 = V{0.5}*cell[2];
auto x79 = V{3}*cell[18];
auto x80 = V{3}*cell[9];
auto x81 = V{3}*cell[17] - V{3}*cell[8];
auto x82 = x60*(V{3}*cell[11] - V{3}*cell[2] - x61 + x63 + x65 + x79 - x80 + x81);
auto x83 = x35 + V{-1};
auto x84 = -V{3}*x24*x43 + x54 + x82 + x83;
auto x85 = -x35;
auto x86 = x82 + x85;
auto x87 = x59*(x43*x69 + x71 + x86);
auto x88 = V{0.5}*cell[11] - x87;
auto x89 = x59*x84 + x78;
auto x90 = x19*(x88 + x89 + V{0.0555555555555556});
auto x91 = V{0.5}*cell[3];
auto x92 = x60*(V{3}*cell[12] - V{3}*cell[3] - x62 + x64 + x66 - x79 + x80 + x81);
auto x93 = -V{3}*x24*x53 + x44 + x83 + x92;
auto x94 = x85 + x92;
auto x95 = x59*(x53*x69 + x70 + x94 + V{1});
auto x96 = V{0.5}*cell[12] - x95;
auto x97 = x59*x93 + x91;
auto x98 = x19*(x96 + x97 + V{0.0555555555555556});
auto x99 = V{0.5}*cell[4];
auto x100 = V{0.0138888888888889}*cell[0] + V{0.0138888888888889}*cell[10] + V{0.0138888888888889}*cell[11] + V{0.0138888888888889}*cell[12] + V{0.0138888888888889}*cell[13] + V{0.0138888888888889}*cell[14] + V{0.0138888888888889}*cell[15] + V{0.0138888888888889}*cell[16] + V{0.0138888888888889}*cell[17] + V{0.0138888888888889}*cell[18] + V{0.0138888888888889}*cell[1] + V{0.0138888888888889}*cell[2] + V{0.0138888888888889}*cell[3] + V{0.0138888888888889}*cell[4] + V{0.0138888888888889}*cell[5] + V{0.0138888888888889}*cell[6] + V{0.0138888888888889}*cell[7] + V{0.0138888888888889}*cell[8] + V{0.0138888888888889}*cell[9] + V{0.0138888888888889};
auto x101 = cell[10] + x31;
auto x102 = cell[11] + V{2}*cell[13] - V{2}*cell[4] + x101 + x38 + x40 + x45;
auto x103 = x102*x102;
auto x104 = x56 + x67;
auto x105 = -V{4.5}*x103*x24 + x104 + x82;
auto x106 = V{4.5}*x24;
auto x107 = x100*(x103*x106 + x73 + x86);
auto x108 = V{0.5}*cell[13] - x107;
auto x109 = x100*x105 + x99;
auto x110 = x19*(x108 + x109 + V{0.0277777777777778});
auto x111 = -x82;
auto x112 = -cell[11];
auto x113 = -cell[17];
auto x114 = cell[8] + x113;
auto x115 = V{2}*cell[14] + cell[2] - V{2}*cell[5] + x101 + x112 + x114 + x27 + x48;
auto x116 = -x115;
auto x117 = x104 + x111 - V{4.5}*x24*x116*x116;
auto x118 = x56 - x67;
auto x119 = x118 - V{4.5}*x24*x115*x115 + x82;
auto x120 = V{0.5}*cell[14];
auto x121 = V{0.5}*cell[5];
auto x122 = x120 - x121;
auto x123 = x100*x117;
auto x124 = x100*x119;
auto x125 = x19*(x120 + x121 + x123 + x124 + V{0.0277777777777778});
auto x126 = V{0.5}*cell[6];
auto x127 = x29 + x32;
auto x128 = cell[12] + V{2}*cell[15] - V{2}*cell[6] + x127 + x37 + x49;
auto x129 = x128*x128;
auto x130 = x104 - V{4.5}*x129*x24 + x92;
auto x131 = x100*(x106*x129 + x73 + x94);
auto x132 = V{0.5}*cell[15] - x131;
auto x133 = x100*x130 + x126;
auto x134 = x19*(x132 + x133 + V{0.0277777777777778});
auto x135 = -x92;
auto x136 = -cell[12];
auto x137 = cell[3] + x136 + x26;
auto x138 = V{2}*cell[16] - V{2}*cell[7] + x114 + x127 + x137 + x38;
auto x139 = -x138;
auto x140 = x104 + x135 - V{4.5}*x24*x139*x139;
auto x141 = x118 - V{4.5}*x24*x138*x138 + x92;
auto x142 = V{0.5}*cell[16];
auto x143 = V{0.5}*cell[7];
auto x144 = x142 - x143;
auto x145 = x100*x140;
auto x146 = x100*x141;
auto x147 = x19*(x142 + x143 + x145 + x146 + V{0.0277777777777778});
auto x148 = V{0.5}*cell[8];
auto x149 = V{2}*cell[17] - V{2}*cell[8] + x28 + x41 + x46 + x51;
auto x150 = x149*x149;
auto x151 = x56 + x82;
auto x152 = -V{4.5}*x150*x24 + x151 + x92;
auto x153 = x100*(x106*x150 + x72 + x86 + x92);
auto x154 = V{0.5}*cell[17] - x153;
auto x155 = x100*x152 + x148;
auto x156 = x19*(x154 + x155 + V{0.0277777777777778});
auto x157 = -cell[15];
auto x158 = V{2}*cell[18] + cell[6] - V{2}*cell[9] + x137 + x157 + x30 + x41;
auto x159 = -x158;
auto x160 = x135 + x151 - V{4.5}*x24*x159*x159;
auto x161 = x111 - V{4.5}*x24*x158*x158 + x56 + x92;
auto x162 = V{0.5}*cell[18];
auto x163 = V{0.5}*cell[9];
auto x164 = x162 - x163;
auto x165 = x100*x160;
auto x166 = x100*x161;
auto x167 = x19*(x162 + x163 + x165 + x166 + V{0.0277777777777778});
auto x0 = cell[0] - x19*(V{1}*cell[0] + x56*(V{0.333333333333333}*cell[0] + V{0.333333333333333}*cell[10] + V{0.333333333333333}*cell[11] + V{0.333333333333333}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] + V{0.333333333333333}) + V{0.333333333333333});
auto x1 = cell[1] + x57*(-x58 - x59*x68 + x75) - x77;
auto x2 = cell[2] + x57*(-x59*x84 - x78 + x88) - x90;
auto x3 = cell[3] + x57*(-x59*x93 - x91 + x96) - x98;
auto x4 = cell[4] - x110 + x57*(-x100*x105 + x108 - x99);
auto x5 = cell[5] - x125 + x57*(-x100*x117 + x100*x119 + x122);
auto x6 = cell[6] - x134 + x57*(-x100*x130 - x126 + x132);
auto x7 = cell[7] - x147 + x57*(-x100*x140 + x100*x141 + x144);
auto x8 = cell[8] - x156 + x57*(-x100*x152 - x148 + x154);
auto x9 = cell[9] - x167 + x57*(-x100*x160 + x100*x161 + x164);
auto x10 = cell[10] - x57*(V{0.5}*cell[10] - x74 - x76) - x77;
auto x11 = -x112 - x57*(V{0.5}*cell[11] - x87 - x89) - x90;
auto x12 = -x136 - x57*(V{0.5}*cell[12] - x95 - x97) - x98;
auto x13 = cell[13] - x110 - x57*(V{0.5}*cell[13] - x107 - x109);
auto x14 = -x125 - x39 - x57*(x122 - x123 + x124);
auto x15 = -x134 - x157 - x57*(V{0.5}*cell[15] - x131 - x133);
auto x16 = -x147 - x50 - x57*(x144 - x145 + x146);
auto x17 = -x113 - x156 - x57*(V{0.5}*cell[17] - x153 - x155);
auto x18 = -x167 - x47 - x57*(x164 - x165 + x166);
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
return { x23, V{1}*x24*(x34 + x43 + x53) };
}
};

}

}
