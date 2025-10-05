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
struct CSE<Dual<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::PlainGuo>, T, descriptors::D3Q19<FIELDS...> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{1}*cell[0];
auto x22 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x23 = cell.template getFieldComponent<opti::F>(11) + cell.template getFieldComponent<opti::F>(13) + cell.template getFieldComponent<opti::F>(18);
auto x24 = cell.template getFieldComponent<opti::F>(12) + cell.template getFieldComponent<opti::F>(15) + cell.template getFieldComponent<opti::F>(9);
auto x25 = cell.template getFieldComponent<opti::F>(1) + cell.template getFieldComponent<opti::F>(4) + cell.template getFieldComponent<opti::F>(6);
auto x26 = cell.template getFieldComponent<opti::F>(0) + cell.template getFieldComponent<opti::F>(10) + cell.template getFieldComponent<opti::F>(14) + cell.template getFieldComponent<opti::F>(16) + cell.template getFieldComponent<opti::F>(17) + cell.template getFieldComponent<opti::F>(2) + cell.template getFieldComponent<opti::F>(3) + cell.template getFieldComponent<opti::F>(5) + cell.template getFieldComponent<opti::F>(7) + cell.template getFieldComponent<opti::F>(8) + x23 + x24 + x25;
auto x27 = x26 + V{1};
auto x28 = V{1} / (x27);
auto x29 = V{6}*cell.template getFieldComponent<opti::F>(4);
auto x30 = V{6}*cell.template getFieldComponent<opti::F>(6);
auto x31 = -V{6}*cell.template getFieldComponent<opti::F>(14) + V{6}*cell.template getFieldComponent<opti::F>(5);
auto x32 = -V{6}*cell.template getFieldComponent<opti::F>(16) + V{6}*cell.template getFieldComponent<opti::F>(7);
auto x33 = x28*(-V{6}*cell.template getFieldComponent<opti::F>(1) + V{6}*cell.template getFieldComponent<opti::F>(10) + V{6}*cell.template getFieldComponent<opti::F>(13) + V{6}*cell.template getFieldComponent<opti::F>(15) - x29 - x30 - x31 - x32);
auto x34 = x22 + x33;
auto x35 = x34 + V{-3};
auto x36 = V{9}*cell.template getFieldComponent<opti::F>(18);
auto x37 = V{9}*cell.template getFieldComponent<opti::F>(4);
auto x38 = V{9}*cell.template getFieldComponent<opti::F>(9);
auto x39 = -V{9}*cell.template getFieldComponent<opti::F>(14) + V{9}*cell.template getFieldComponent<opti::F>(5);
auto x40 = V{9}*cell.template getFieldComponent<opti::F>(17) - V{9}*cell.template getFieldComponent<opti::F>(8);
auto x41 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(1) + x28*(V{9}*cell.template getFieldComponent<opti::F>(11) + V{9}*cell.template getFieldComponent<opti::F>(13) - V{9}*cell.template getFieldComponent<opti::F>(2) + x36 - x37 - x38 + x39 + x40);
auto x42 = V{0.027778}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x43 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x44 = V{6}*cell.template getFieldComponent<opti::F>(18);
auto x45 = V{6}*cell.template getFieldComponent<opti::F>(9);
auto x46 = V{6}*cell.template getFieldComponent<opti::F>(17) - V{6}*cell.template getFieldComponent<opti::F>(8);
auto x47 = x28*(V{6}*cell.template getFieldComponent<opti::F>(11) + V{6}*cell.template getFieldComponent<opti::F>(13) - V{6}*cell.template getFieldComponent<opti::F>(2) - x29 + x31 + x44 - x45 + x46);
auto x48 = x43 + x47;
auto x49 = x48 + V{-3};
auto x50 = V{9}*cell.template getFieldComponent<opti::F>(6);
auto x51 = -V{9}*cell.template getFieldComponent<opti::F>(16) + V{9}*cell.template getFieldComponent<opti::F>(7);
auto x52 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(0) + x28*(-V{9}*cell.template getFieldComponent<opti::F>(1) + V{9}*cell.template getFieldComponent<opti::F>(10) + V{9}*cell.template getFieldComponent<opti::F>(13) + V{9}*cell.template getFieldComponent<opti::F>(15) - x37 - x39 - x50 - x51);
auto x53 = V{0.027778}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x54 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x55 = V{1}*cell.template getFieldComponent<opti::F>(9);
auto x56 = V{1}*cell.template getFieldComponent<opti::F>(18);
auto x57 = V{1}*cell.template getFieldComponent<opti::F>(6);
auto x58 = V{1}*cell.template getFieldComponent<opti::F>(17) - V{1}*cell.template getFieldComponent<opti::F>(8);
auto x59 = -V{1}*cell.template getFieldComponent<opti::F>(16) + V{1}*cell.template getFieldComponent<opti::F>(7);
auto x60 = x28*(V{1}*cell.template getFieldComponent<opti::F>(12) + V{1}*cell.template getFieldComponent<opti::F>(15) - V{1}*cell.template getFieldComponent<opti::F>(3) + x55 - x56 - x57 + x58 + x59);
auto x61 = x54 + x60;
auto x62 = cell.template getFieldComponent<descriptors::FORCE>(2)*x61;
auto x63 = V{0.083333}*x62;
auto x64 = -x63;
auto x65 = V{0.250002}*cell.template getFieldComponent<opti::F>(4);
auto x66 = V{0.250002}*cell.template getFieldComponent<opti::F>(6);
auto x67 = -V{0.250002}*cell.template getFieldComponent<opti::F>(14) + V{0.250002}*cell.template getFieldComponent<opti::F>(5);
auto x68 = -V{0.250002}*cell.template getFieldComponent<opti::F>(16) + V{0.250002}*cell.template getFieldComponent<opti::F>(7);
auto x69 = -V{0.125001}*cell.template getFieldComponent<descriptors::FORCE>(0) - x28*(-V{0.250002}*cell.template getFieldComponent<opti::F>(1) + V{0.250002}*cell.template getFieldComponent<opti::F>(10) + V{0.250002}*cell.template getFieldComponent<opti::F>(13) + V{0.250002}*cell.template getFieldComponent<opti::F>(15) - x65 - x66 - x67 - x68) + V{0.083334};
auto x70 = V{0.166668}*cell.template getFieldComponent<opti::F>(18);
auto x71 = V{0.166668}*cell.template getFieldComponent<opti::F>(4);
auto x72 = V{0.166668}*cell.template getFieldComponent<opti::F>(9);
auto x73 = V{0.166668}*cell.template getFieldComponent<opti::F>(17) - V{0.166668}*cell.template getFieldComponent<opti::F>(8);
auto x74 = -V{0.166668}*cell.template getFieldComponent<opti::F>(14) + V{0.166668}*cell.template getFieldComponent<opti::F>(5);
auto x75 = V{0.083334}*cell.template getFieldComponent<descriptors::FORCE>(1) + x28*(V{0.166668}*cell.template getFieldComponent<opti::F>(11) + V{0.166668}*cell.template getFieldComponent<opti::F>(13) - V{0.166668}*cell.template getFieldComponent<opti::F>(2) + x70 - x71 - x72 + x73 + x74);
auto x76 = x41 + V{3};
auto x77 = -x22 - x33;
auto x78 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(2) + x28*(V{9}*cell.template getFieldComponent<opti::F>(12) + V{9}*cell.template getFieldComponent<opti::F>(15) - V{9}*cell.template getFieldComponent<opti::F>(3) - x36 + x38 + x40 - x50 + x51);
auto x79 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x80 = x28*(V{6}*cell.template getFieldComponent<opti::F>(12) + V{6}*cell.template getFieldComponent<opti::F>(15) - V{6}*cell.template getFieldComponent<opti::F>(3) - x30 + x32 - x44 + x45 + x46);
auto x81 = x79 + x80;
auto x82 = x81 + V{-3};
auto x83 = V{0.027778}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x84 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x85 = V{1}*cell.template getFieldComponent<opti::F>(4);
auto x86 = -V{1}*cell.template getFieldComponent<opti::F>(14) + V{1}*cell.template getFieldComponent<opti::F>(5);
auto x87 = x28*(V{1}*cell.template getFieldComponent<opti::F>(11) + V{1}*cell.template getFieldComponent<opti::F>(13) - V{1}*cell.template getFieldComponent<opti::F>(2) - x55 + x56 + x58 - x85 + x86);
auto x88 = x84 + x87;
auto x89 = cell.template getFieldComponent<descriptors::FORCE>(1)*x88;
auto x90 = V{0.083333}*x89;
auto x91 = -x90;
auto x92 = V{0.166668}*cell.template getFieldComponent<opti::F>(6);
auto x93 = -V{0.166668}*cell.template getFieldComponent<opti::F>(16) + V{0.166668}*cell.template getFieldComponent<opti::F>(7);
auto x94 = V{0.083334}*cell.template getFieldComponent<descriptors::FORCE>(2) + x28*(V{0.166668}*cell.template getFieldComponent<opti::F>(12) + V{0.166668}*cell.template getFieldComponent<opti::F>(15) - V{0.166668}*cell.template getFieldComponent<opti::F>(3) - x70 + x72 + x73 - x92 + x93);
auto x95 = x78 + V{3};
auto x96 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x97 = V{1}*cell.template getFieldComponent<opti::F>(1) - V{1}*cell.template getFieldComponent<opti::F>(10) - V{1}*cell.template getFieldComponent<opti::F>(13) - V{1}*cell.template getFieldComponent<opti::F>(15) + x57 + x59 + x85 + x86;
auto x98 = -x28*x97 + x96;
auto x99 = cell.template getFieldComponent<descriptors::FORCE>(0)*x98;
auto x100 = V{0.083333}*x99;
auto x101 = -x100;
auto x102 = V{0.250002}*cell.template getFieldComponent<opti::F>(18);
auto x103 = V{0.250002}*cell.template getFieldComponent<opti::F>(9);
auto x104 = V{0.250002}*cell.template getFieldComponent<opti::F>(17) - V{0.250002}*cell.template getFieldComponent<opti::F>(8);
auto x105 = -V{0.125001}*cell.template getFieldComponent<descriptors::FORCE>(1) - x28*(V{0.250002}*cell.template getFieldComponent<opti::F>(11) + V{0.250002}*cell.template getFieldComponent<opti::F>(13) - V{0.250002}*cell.template getFieldComponent<opti::F>(2) + x102 - x103 + x104 - x65 + x67) + V{0.083334};
auto x106 = -x43 - x47;
auto x107 = x34 + V{3};
auto x108 = x48 + V{3};
auto x109 = V{0.083334}*cell.template getFieldComponent<descriptors::FORCE>(0) + x28*(-V{0.166668}*cell.template getFieldComponent<opti::F>(1) + V{0.166668}*cell.template getFieldComponent<opti::F>(10) + V{0.166668}*cell.template getFieldComponent<opti::F>(13) + V{0.166668}*cell.template getFieldComponent<opti::F>(15) - x71 - x74 - x92 - x93);
auto x110 = x52 + V{3};
auto x111 = x81 + V{3};
auto x112 = -V{0.125001}*cell.template getFieldComponent<descriptors::FORCE>(2) - x28*(V{0.250002}*cell.template getFieldComponent<opti::F>(12) + V{0.250002}*cell.template getFieldComponent<opti::F>(15) - V{0.250002}*cell.template getFieldComponent<opti::F>(3) - x102 + x103 + x104 - x66 + x68) + V{0.083334};
auto x113 = -x79 - x80;
auto x114 = V{0.055556}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x115 = V{0.166667}*x89;
auto x116 = V{0.166667}*x62;
auto x117 = x115 + x116;
auto x118 = V{0.055556}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x119 = V{0.166667}*x99;
auto x120 = x116 + x119;
auto x121 = V{0.055556}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x122 = x115 + x119;
auto x123 = V{0.5}*x19 + V{-1};
auto x124 = V{1}*x123;
auto x125 = V{0.0277777777777778}*x19;
auto x126 = x88 + x98;
auto x127 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x128 = V{4.5}*cell.template getFieldComponent<opti::F>(13);
auto x129 = V{4.5}*cell.template getFieldComponent<opti::F>(18);
auto x130 = V{4.5}*cell.template getFieldComponent<opti::F>(4);
auto x131 = V{4.5}*cell.template getFieldComponent<opti::F>(9);
auto x132 = -V{4.5}*cell.template getFieldComponent<opti::F>(14) + V{4.5}*cell.template getFieldComponent<opti::F>(5);
auto x133 = V{4.5}*cell.template getFieldComponent<opti::F>(17) - V{4.5}*cell.template getFieldComponent<opti::F>(8);
auto x134 = x28*(V{4.5}*cell.template getFieldComponent<opti::F>(11) - V{4.5}*cell.template getFieldComponent<opti::F>(2) + x128 + x129 - x130 - x131 + x132 + x133);
auto x135 = x127 + x134;
auto x136 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x137 = V{4.5}*cell.template getFieldComponent<opti::F>(6);
auto x138 = V{4.5}*cell.template getFieldComponent<opti::F>(15);
auto x139 = -V{4.5}*cell.template getFieldComponent<opti::F>(16) + V{4.5}*cell.template getFieldComponent<opti::F>(7);
auto x140 = V{4.5}*cell.template getFieldComponent<opti::F>(1) - V{4.5}*cell.template getFieldComponent<opti::F>(10) - x128 + x130 + x132 + x137 - x138 + x139;
auto x141 = x136 - x140*x28;
auto x142 = x135 + x141;
auto x143 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x144 = V{3}*cell.template getFieldComponent<opti::F>(13);
auto x145 = V{3}*cell.template getFieldComponent<opti::F>(18);
auto x146 = V{3}*cell.template getFieldComponent<opti::F>(4);
auto x147 = V{3}*cell.template getFieldComponent<opti::F>(9);
auto x148 = -V{3}*cell.template getFieldComponent<opti::F>(14) + V{3}*cell.template getFieldComponent<opti::F>(5);
auto x149 = V{3}*cell.template getFieldComponent<opti::F>(17) - V{3}*cell.template getFieldComponent<opti::F>(8);
auto x150 = x28*(V{3}*cell.template getFieldComponent<opti::F>(11) - V{3}*cell.template getFieldComponent<opti::F>(2) + x144 + x145 - x146 - x147 + x148 + x149);
auto x151 = x143 + x150;
auto x152 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x153 = V{3}*cell.template getFieldComponent<opti::F>(6);
auto x154 = V{3}*cell.template getFieldComponent<opti::F>(15);
auto x155 = -V{3}*cell.template getFieldComponent<opti::F>(16) + V{3}*cell.template getFieldComponent<opti::F>(7);
auto x156 = V{3}*cell.template getFieldComponent<opti::F>(1) - V{3}*cell.template getFieldComponent<opti::F>(10) - x144 + x146 + x148 + x153 - x154 + x155;
auto x157 = -x156*x28;
auto x158 = V{1.5}*(x88*x88);
auto x159 = V{1.5}*(x61*x61);
auto x160 = x158 + x159 + V{-1};
auto x161 = x160 + V{1.5}*(x98*x98);
auto x162 = x152 + x157 + x161;
auto x163 = -x84 - x87 + x98;
auto x164 = -x127 - x134 + x141;
auto x165 = -x143 - x150;
auto x166 = x61 + x98;
auto x167 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x168 = x28*(V{4.5}*cell.template getFieldComponent<opti::F>(12) - V{4.5}*cell.template getFieldComponent<opti::F>(3) - x129 + x131 + x133 - x137 + x138 + x139);
auto x169 = x167 + x168;
auto x170 = x141 + x169;
auto x171 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x172 = x28*(V{3}*cell.template getFieldComponent<opti::F>(12) - V{3}*cell.template getFieldComponent<opti::F>(3) - x145 + x147 + x149 - x153 + x154 + x155);
auto x173 = x171 + x172;
auto x174 = -x54 - x60;
auto x175 = x174 + x98;
auto x176 = -x167 - x168;
auto x177 = x141 + x176;
auto x178 = -x171 - x172;
auto x179 = x61 + x88;
auto x180 = x135 + x169;
auto x181 = x179*x180;
auto x182 = x151 + x161;
auto x183 = x174 + x88;
auto x184 = x135 + x176;
auto x185 = -x152;
auto x186 = -x157 + x185;
auto x187 = x161 + x173;
auto x188 = V{0.0555555555555556}*x19;
auto x189 = x135*x88;
auto x190 = x169*x61;
auto x191 = x28*x97;
auto x192 = -x191 + x96;
auto x193 = x192 + x88;
auto x194 = x140*x28;
auto x195 = x136 - x194;
auto x196 = x135 + x195;
auto x197 = V{1.5}*(x192*x192);
auto x198 = -x158 - x159 - x197 + V{1};
auto x199 = x151 + x198;
auto x200 = x156*x28;
auto x201 = x152 - x200;
auto x202 = x192 + x61;
auto x203 = x169 + x195;
auto x204 = x173 + x198;
auto x205 = -cell.template getFieldComponent<opti::F>(14) + cell.template getFieldComponent<opti::F>(5);
auto x206 = -cell.template getFieldComponent<opti::F>(16) + cell.template getFieldComponent<opti::F>(7);
auto x207 = -cell.template getFieldComponent<opti::F>(10) - cell.template getFieldComponent<opti::F>(13) - cell.template getFieldComponent<opti::F>(15) + x205 + x206 + x25;
auto x208 = V{0.250002}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x209 = V{0.250002}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x210 = V{0.166668}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x211 = V{0.333336}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x212 = V{0.166667}*cell[12] + V{0.166667}*cell[3];
auto x213 = V{0.166667}*cell[11] + V{0.166667}*cell[2];
auto x214 = x124*(-V{1}*cell.template getFieldComponent<descriptors::FORCE>(0)*(V{0.083333}*cell[17] + V{0.083333}*cell[18] + V{0.083333}*cell[8] + V{0.083333}*cell[9] + x21 + x212 + x213) + cell[10]*x211 + cell[13]*x208 + cell[13]*x210 - cell[14]*x208 + cell[14]*x210 + cell[15]*x209 + cell[15]*x210 - cell[16]*x209 + cell[16]*x210 + cell[1]*x211 + cell[4]*x208 + cell[4]*x210 - cell[5]*x208 + cell[5]*x210 + cell[6]*x209 + cell[6]*x210 - cell[7]*x209 + cell[7]*x210);
auto x215 = cell.template getFieldComponent<opti::F>(17) - cell.template getFieldComponent<opti::F>(8);
auto x216 = -cell.template getFieldComponent<opti::F>(2) - cell.template getFieldComponent<opti::F>(4) - cell.template getFieldComponent<opti::F>(9) + x205 + x215 + x23;
auto x217 = V{0.250002}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x218 = V{0.166668}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x219 = V{0.333336}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x220 = V{0.166667}*cell[10] + V{0.166667}*cell[1] + x21;
auto x221 = -V{1}*cell.template getFieldComponent<descriptors::FORCE>(1)*(V{0.083333}*cell[15] + V{0.083333}*cell[16] + V{0.083333}*cell[6] + V{0.083333}*cell[7] + x212 + x220) + cell[11]*x219 + cell[13]*x217 + cell[13]*x218 - cell[14]*x217 + cell[14]*x218 + cell[17]*x209 + cell[17]*x218 - cell[18]*x209 + cell[18]*x218 + cell[2]*x219 + cell[4]*x217 + cell[4]*x218 - cell[5]*x217 + cell[5]*x218 + cell[8]*x209 + cell[8]*x218 - cell[9]*x209 + cell[9]*x218;
auto x222 = -cell.template getFieldComponent<opti::F>(18) - cell.template getFieldComponent<opti::F>(3) - cell.template getFieldComponent<opti::F>(6) + x206 + x215 + x24;
auto x223 = V{0.166668}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x224 = V{0.333336}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x225 = -V{1}*cell.template getFieldComponent<descriptors::FORCE>(2)*(V{0.083333}*cell[13] + V{0.083333}*cell[14] + V{0.083333}*cell[4] + V{0.083333}*cell[5] + x213 + x220) + cell[12]*x224 + cell[15]*x217 + cell[15]*x223 - cell[16]*x217 + cell[16]*x223 + cell[17]*x208 + cell[17]*x223 - cell[18]*x208 + cell[18]*x223 + cell[3]*x224 + cell[6]*x217 + cell[6]*x223 - cell[7]*x217 + cell[7]*x223 + cell[8]*x208 + cell[8]*x223 - cell[9]*x208 + cell[9]*x223;
auto x226 = V{0.0833333333333333}*cell[14];
auto x227 = V{0.0833333333333333}*cell[9];
auto x228 = V{0.0833333333333333}*cell[18];
auto x229 = V{0.0833333333333333}*cell[5];
auto x230 = V{0.125}*x183;
auto x231 = cell[18]*x230;
auto x232 = cell[9]*x230;
auto x233 = V{0.25}*x88;
auto x234 = V{0.0277777777777778}*x184;
auto x235 = cell[18]*x234;
auto x236 = cell[9]*x234;
auto x237 = V{0.0555555555555556}*x135;
auto x238 = V{0.125}*cell[14]*x163;
auto x239 = x191 - x96;
auto x240 = V{0.125}*cell[5]*(-x239 - x88);
auto x241 = V{1}*cell[0] + V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[12] + V{0.0833333333333333}*cell[13] + V{0.0833333333333333}*cell[14] + V{0.0833333333333333}*cell[15] + V{0.0833333333333333}*cell[16] + V{0.0833333333333333}*cell[17] + V{0.0833333333333333}*cell[18] + V{0.166666666666667}*cell[1] + V{0.166666666666667}*cell[2] + V{0.166666666666667}*cell[3] + V{0.0833333333333333}*cell[4] + V{0.0833333333333333}*cell[5] + V{0.0833333333333333}*cell[6] + V{0.0833333333333333}*cell[7] + V{0.0833333333333333}*cell[8] + V{0.0833333333333333}*cell[9];
auto x242 = V{0.0277777777777778}*cell[14]*x164;
auto x243 = -x136 + x194;
auto x244 = V{0.0277777777777778}*cell[5]*(-x135 - x243);
auto x245 = V{0.125}*cell[13]*x126 + V{0.0277777777777778}*cell[13]*x142 - V{0.0833333333333333}*cell[13] + V{0.125}*cell[4]*x193 + V{0.0277777777777778}*cell[4]*x196 + V{0.0833333333333333}*cell[4];
auto x246 = V{0.125}*x179;
auto x247 = V{0.0277777777777778}*x180;
auto x248 = cell[17]*x246 + cell[17]*x247 - V{0.0833333333333333}*cell[17] + cell[8]*x246 + cell[8]*x247 + V{0.0833333333333333}*cell[8];
auto x249 = cell[11]*x233 + cell[11]*x237 - V{0.166666666666667}*cell[11] + cell[2]*x233 + cell[2]*x237 + V{0.166666666666667}*cell[2] + x226 + x227 - x228 - x229 + x231 + x232 + x235 + x236 - x238 - x240 - x241*x88 - x242 - x244 + x245 + x248;
auto x250 = V{1} / ((x27)*(x27));
auto x251 = x26 + V{1};
auto x252 = V{1}*x19*x251;
auto x253 = x250*x252;
auto x254 = V{0.0833333333333333}*cell[16];
auto x255 = V{0.0833333333333333}*cell[7];
auto x256 = V{0.25}*x61;
auto x257 = V{0.0555555555555556}*x169;
auto x258 = V{0.125}*cell[16]*x175;
auto x259 = V{0.125}*cell[7]*(-x239 - x61);
auto x260 = V{0.0277777777777778}*cell[16]*x177;
auto x261 = V{0.0277777777777778}*cell[7]*(-x169 - x243);
auto x262 = V{0.125}*cell[15]*x166 + V{0.0277777777777778}*cell[15]*x170 - V{0.0833333333333333}*cell[15] + V{0.125}*cell[6]*x202 + V{0.0277777777777778}*cell[6]*x203 + V{0.0833333333333333}*cell[6];
auto x263 = cell[12]*x256 + cell[12]*x257 - V{0.166666666666667}*cell[12] + cell[3]*x256 + cell[3]*x257 + V{0.166666666666667}*cell[3] - x227 + x228 - x231 - x232 - x235 - x236 - x241*x61 + x248 + x254 - x255 - x258 - x259 - x260 - x261 + x262;
auto x264 = V{0.0555555555555556}*cell[10]*x141 + V{0.25}*cell[10]*x98 - V{0.166666666666667}*cell[10] + V{0.25}*cell[1]*x192 + V{0.0555555555555556}*cell[1]*x195 + V{0.166666666666667}*cell[1] - x192*x241 - x226 + x229 + x238 + x240 + x242 + x244 + x245 - x254 + x255 + x258 + x259 + x260 + x261 + x262;
auto x265 = V{0.333333333333333}*cell[0]*x161*x19 + cell[10]*x188*(-x141*x98 + x162) + cell[11]*x188*(x182 - x189) + cell[12]*x188*(x187 - x190) + cell[13]*x125*(-x126*x142 + x151 + x162) + cell[14]*x125*(x162 - x163*x164 + x165) + cell[15]*x125*(x162 - x166*x170 + x173) + cell[16]*x125*(x162 - x175*x177 + x178) + cell[17]*x125*(x173 - x181 + x182) + cell[18]*x125*(x178 + x182 - x183*x184) - V{0.0555555555555556}*cell[1]*x19*(-x160 - x185 + x192*x195 - x197 - x200) - V{0.0555555555555556}*cell[2]*x19*(x189 + x199) - V{0.0555555555555556}*cell[3]*x19*(x190 + x204) - V{0.0277777777777778}*cell[4]*x19*(x193*x196 + x199 + x201) + cell[5]*x125*(-x163*x164 + x182 + x186) - V{0.0277777777777778}*cell[6]*x19*(x201 + x202*x203 + x204) + cell[7]*x125*(-x175*x177 + x186 + x187) - V{0.0277777777777778}*cell[8]*x19*(x173 + x181 + x199) + cell[9]*x125*(x165 - x183*x184 + x187) - V{1}*x123*x216*x221*x28 - V{1}*x123*x222*x225*x28 + x124*(-cell[0]*(x62 + x89 + x99) - cell[10]*(-x114*x35 + x117) - cell[11]*(-x118*x49 + x120) - cell[12]*(-x121*x82 + x122) + cell[13]*(x42*(x35 + x41) + x53*(x49 + x52) + x64) + cell[14]*(cell.template getFieldComponent<descriptors::FORCE>(1)*(x69 + x75) - x42*(x76 + x77) - x63) + cell[15]*(x42*(x35 + x78) + x83*(x52 + x82) + x91) + cell[16]*(cell.template getFieldComponent<descriptors::FORCE>(2)*(x69 + x94) - x42*(x77 + x95) - x90) + cell[17]*(x101 + x53*(x49 + x78) + x83*(x41 + x82)) + cell[18]*(cell.template getFieldComponent<descriptors::FORCE>(2)*(x105 + x94) - x100 - x53*(x106 + x95)) - cell[1]*(-x107*x114 + x117) - cell[2]*(-x108*x118 + x120) - cell[3]*(-x111*x121 + x122) + cell[4]*(x42*(x107 + x41) + x53*(x108 + x52) + x64) + cell[5]*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x105 + x109) - x53*(x106 + x110) - x63) + cell[6]*(x42*(x107 + x78) + x83*(x111 + x52) + x91) + cell[7]*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x109 + x112) - x83*(x110 + x113) - x90) + cell[8]*(x101 + x53*(x108 + x78) + x83*(x111 + x41)) + cell[9]*(cell.template getFieldComponent<descriptors::FORCE>(1)*(x112 + x75) - x100 - x83*(x113 + x76))) - V{1}*x19*x207*x250*x251*x264 + x207*x214*x28 + x216*x249*x253 + x222*x253*x263;
auto x266 = -cell.template getFieldComponent<opti::DJDF>(0) + x20*x21 + x265;
auto x267 = V{1}*x20;
auto x268 = x252*x28;
auto x269 = x264*x268;
auto x270 = x214 + x265 - x269;
auto x271 = -cell.template getFieldComponent<opti::DJDF>(10) + cell[1]*x267 + x270;
auto x272 = x124*x221;
auto x273 = x249*x268;
auto x274 = x272 - x273;
auto x275 = x265 + x274;
auto x276 = -cell.template getFieldComponent<opti::DJDF>(11) + cell[2]*x267 + x275;
auto x277 = x124*x225;
auto x278 = x263*x268;
auto x279 = x277 - x278;
auto x280 = x265 + x279;
auto x281 = -cell.template getFieldComponent<opti::DJDF>(12) + cell[3]*x267 + x280;
auto x282 = -cell.template getFieldComponent<opti::DJDF>(13) + cell[4]*x267 + x270 + x274;
auto x283 = -x272 + x273;
auto x284 = -cell.template getFieldComponent<opti::DJDF>(14) + cell[5]*x267 + x270 + x283;
auto x285 = -cell.template getFieldComponent<opti::DJDF>(15) + cell[6]*x267 + x270 + x279;
auto x286 = -x277 + x278;
auto x287 = -cell.template getFieldComponent<opti::DJDF>(16) + cell[7]*x267 + x270 + x286;
auto x288 = -cell.template getFieldComponent<opti::DJDF>(17) + cell[8]*x267 + x275 + x279;
auto x289 = -cell.template getFieldComponent<opti::DJDF>(18) + cell[9]*x267 + x275 + x286;
auto x290 = -x214 + x269;
auto x291 = x265 + x290;
auto x292 = -cell.template getFieldComponent<opti::DJDF>(1) + cell[10]*x267 + x291;
auto x293 = x265 + x283;
auto x294 = -cell.template getFieldComponent<opti::DJDF>(2) + cell[11]*x267 + x293;
auto x295 = -cell.template getFieldComponent<opti::DJDF>(3) + cell[12]*x267 + x265 + x286;
auto x296 = -cell.template getFieldComponent<opti::DJDF>(4) + cell[13]*x267 + x283 + x291;
auto x297 = -cell.template getFieldComponent<opti::DJDF>(5) + cell[14]*x267 + x275 + x290;
auto x298 = -cell.template getFieldComponent<opti::DJDF>(6) + cell[15]*x267 + x286 + x291;
auto x299 = -cell.template getFieldComponent<opti::DJDF>(7) + cell[16]*x267 + x280 + x290;
auto x300 = -cell.template getFieldComponent<opti::DJDF>(8) + cell[17]*x267 + x286 + x293;
auto x301 = -cell.template getFieldComponent<opti::DJDF>(9) + cell[18]*x267 + x280 + x283;
cell[0] = -x266;
cell[1] = -x271;
cell[2] = -x276;
cell[3] = -x281;
cell[4] = -x282;
cell[5] = -x284;
cell[6] = -x285;
cell[7] = -x287;
cell[8] = -x288;
cell[9] = -x289;
cell[10] = -x292;
cell[11] = -x294;
cell[12] = -x295;
cell[13] = -x296;
cell[14] = -x297;
cell[15] = -x298;
cell[16] = -x299;
cell[17] = -x300;
cell[18] = -x301;
return { V{1} - x266, x266*x266 + x271*x271 + x276*x276 + x281*x281 + x282*x282 + x284*x284 + x285*x285 + x287*x287 + x288*x288 + x289*x289 + x292*x292 + x294*x294 + x295*x295 + x296*x296 + x297*x297 + x298*x298 + x299*x299 + x300*x300 + x301*x301 };
}
};

}

}
