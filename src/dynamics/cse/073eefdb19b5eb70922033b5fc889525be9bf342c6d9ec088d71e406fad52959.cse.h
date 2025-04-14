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
struct CSE<dynamics::Tuple<T, descriptors::D3Q27<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, forcing::Guo<momenta::ForcedWithStress> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x27 = parameters.template get<descriptors::OMEGA>();
auto x28 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x29 = cell[13] + cell[1] + cell[5] + cell[7];
auto x30 = cell[18] + cell[25] + cell[2] + cell[9];
auto x31 = cell[10] + cell[20] + cell[22] + cell[24] + cell[3];
auto x32 = cell[0] + cell[11] + cell[12] + cell[14] + cell[15] + cell[16] + cell[17] + cell[19] + cell[21] + cell[23] + cell[26] + cell[4] + cell[6] + cell[8] + x29 + x30 + x31;
auto x33 = x32 + V{1};
auto x34 = V{1} / (x33);
auto x35 = x32 + V{1};
auto x36 = x35/((x33)*(x33));
auto x37 = -cell[23];
auto x38 = cell[10] + cell[11] - cell[17] - cell[24] + cell[4] + x37;
auto x39 = cell[12] - cell[19] - cell[25] + cell[6];
auto x40 = -cell[14] - cell[18] - cell[20] - cell[26] + x29 + x38 + x39;
auto x41 = x34*x35;
auto x42 = V{0.666667}*cell[10];
auto x43 = V{0.666667}*cell[11];
auto x44 = V{0.666667}*cell[12];
auto x45 = V{0.666667}*cell[13];
auto x46 = V{0.666667}*cell[23];
auto x47 = V{0.666667}*cell[24];
auto x48 = V{0.666667}*cell[25];
auto x49 = V{0.666667}*cell[26];
auto x50 = -V{0.333333}*cell[0];
auto x51 = -V{0.333333}*cell[16] + V{0.666667}*cell[17] + V{0.666667}*cell[18] - V{0.333333}*cell[3] + V{0.666667}*cell[4] + V{0.666667}*cell[5] + x42 + x43 + x44 + x45 + x46 + x47 + x48 + x49 + x50;
auto x52 = -V{0.333333}*cell[15] + V{0.666667}*cell[19] + V{0.666667}*cell[20] - V{0.333333}*cell[2] + V{0.666667}*cell[6] + V{0.666667}*cell[7];
auto x53 = -cell.template getFieldComponent<descriptors::FORCE>(0)*x40*x41 + V{0.666667}*cell[14] + V{0.666667}*cell[1] - V{0.333333}*cell[21] - V{0.333333}*cell[22] - V{0.333333}*cell[8] - V{0.333333}*cell[9] - x36*x40*x40 + x51 + x52;
auto x54 = -cell[13] - cell[21] + cell[26] + cell[8];
auto x55 = -cell[12] - cell[15] - cell[22] - cell[5] + x30 + x38 + x54;
auto x56 = -V{0.333333}*cell[14] - V{0.333333}*cell[1] + V{0.666667}*cell[21] + V{0.666667}*cell[22] + V{0.666667}*cell[8] + V{0.666667}*cell[9];
auto x57 = -cell.template getFieldComponent<descriptors::FORCE>(1)*x41*x55 + V{0.666667}*cell[15] - V{0.333333}*cell[19] - V{0.333333}*cell[20] + V{0.666667}*cell[2] - V{0.333333}*cell[6] - V{0.333333}*cell[7] - x36*x55*x55 + x51 + x56;
auto x58 = -cell[11] - cell[16] - cell[7] - cell[9] + x31 + x37 + x39 + x54;
auto x59 = -cell.template getFieldComponent<descriptors::FORCE>(2)*x41*x58 + V{0.666667}*cell[16] - V{0.333333}*cell[17] - V{0.333333}*cell[18] + V{0.666667}*cell[3] - V{0.333333}*cell[4] - V{0.333333}*cell[5] - x36*x58*x58 + x42 + x43 + x44 + x45 + x46 + x47 + x48 + x49 + x50 + x52 + x56;
auto x60 = V{1}*cell[18];
auto x61 = V{1}*cell[11];
auto x62 = -x61;
auto x63 = V{1}*cell[4];
auto x64 = x41*(cell.template getFieldComponent<descriptors::FORCE>(0)*x55 + cell.template getFieldComponent<descriptors::FORCE>(1)*x40);
auto x65 = V{1}*x55;
auto x66 = x34*x40;
auto x67 = V{1}*cell[12];
auto x68 = V{1}*cell[25];
auto x69 = V{1}*cell[10];
auto x70 = -x69;
auto x71 = -V{1}*cell[23];
auto x72 = x67 + x68 + x70 + x71;
auto x73 = V{1}*cell[13];
auto x74 = V{1}*cell[26];
auto x75 = x73 + x74;
auto x76 = V{1}*cell[5];
auto x77 = -V{1}*cell[17];
auto x78 = V{1}*cell[24];
auto x79 = -x78;
auto x80 = x76 + x77 + x79;
auto x81 = V{2}*cell[13];
auto x82 = V{2}*cell[26];
auto x83 = V{2}*cell[11];
auto x84 = V{2}*cell[24];
auto x85 = V{2}*x55;
auto x86 = -V{2}*cell[10] + V{2}*cell[12] - V{2}*cell[23] + V{2}*cell[25];
auto x87 = V{1}*cell[22];
auto x88 = -x74;
auto x89 = V{1}*cell[8];
auto x90 = x41*(cell.template getFieldComponent<descriptors::FORCE>(1)*x58 + cell.template getFieldComponent<descriptors::FORCE>(2)*x55);
auto x91 = x36*x58;
auto x92 = x61 + x78;
auto x93 = V{1}*cell[9];
auto x94 = -x73;
auto x95 = -V{1}*cell[21];
auto x96 = x93 + x94 + x95;
auto x97 = V{1}*cell[20];
auto x98 = V{1}*cell[6];
auto x99 = -x67 + x71;
auto x100 = V{1}*cell[7];
auto x101 = -V{1}*cell[19];
auto x102 = -x68;
auto x103 = x100 + x101 + x102;
auto x104 = -x103 - V{0.5}*x41*(cell.template getFieldComponent<descriptors::FORCE>(0)*x58 + cell.template getFieldComponent<descriptors::FORCE>(2)*x40) - V{1}*x58*x66 - x70 - x75 - x92 - x97 + x98 - x99;
auto x105 = V{1} / (V{3.00000046417339}*util::sqrt(x34*(x28*x28)*util::sqrt(V{0.5}*(-x60 - x62 + x63 - V{0.5}*x64 - x65*x66 - x72 - x75 - x80)*(V{2}*cell[17] - V{2}*cell[18] + V{2}*cell[4] - V{2}*cell[5] - V{1}*x64 - x66*x85 - x81 - x82 + x83 + x84 - x86) + V{0.5}*(-x65*x91 - x72 - x87 - x88 + x89 - V{0.5}*x90 - x92 - x96)*(V{2}*cell[21] - V{2}*cell[22] + V{2}*cell[8] - V{2}*cell[9] + x81 + x82 - x83 - x84 - x85*x91 - x86 - V{1}*x90) + x104*x104 + V{0.5}*(x53*x53) + V{0.5}*(x57*x57) + V{0.5}*(x59*x59)) + V{0.0277777691819762}/((x27)*(x27))) + V{0.5}/x27);
auto x106 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x107 = x61 + x63 + x69;
auto x108 = x67 + x71 + x98;
auto x109 = -V{1}*cell[14] + V{1}*cell[1] + x103 + x107 + x108 - x60 + x73 + x80 + x88 - x97;
auto x110 = -x109*x34;
auto x111 = x106 + x110;
auto x112 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x113 = x74 + x89;
auto x114 = -V{1}*cell[15] + V{1}*cell[2] + x107 + x113 + x60 + x68 - x76 + x77 + x79 - x87 + x96 + x99;
auto x115 = -x114*x34;
auto x116 = x112 + x115;
auto x117 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x118 = -V{1}*cell[16] + V{1}*cell[3] - x100 + x101 + x102 + x108 + x113 + x62 + x69 + x78 + x87 - x93 + x94 + x95 + x97;
auto x119 = -x118*x34;
auto x120 = x117 + x119;
auto x121 = V{-1} + V{1.5}*(x111*x111) + V{1.5}*(x116*x116) + V{1.5}*(x120*x120);
auto x122 = V{1} - x105;
auto x123 = x106 - x109*x34;
auto x124 = cell.template getFieldComponent<descriptors::FORCE>(0)*x123;
auto x125 = x114*x34;
auto x126 = x112 - x125;
auto x127 = cell.template getFieldComponent<descriptors::FORCE>(1)*x126;
auto x128 = x118*x34;
auto x129 = x117 - x128;
auto x130 = cell.template getFieldComponent<descriptors::FORCE>(2)*x129;
auto x131 = V{1} - V{0.5}*x105;
auto x132 = V{0.0740740740740741}*cell[0] + V{0.0740740740740741}*cell[10] + V{0.0740740740740741}*cell[11] + V{0.0740740740740741}*cell[12] + V{0.0740740740740741}*cell[13] + V{0.0740740740740741}*cell[14] + V{0.0740740740740741}*cell[15] + V{0.0740740740740741}*cell[16] + V{0.0740740740740741}*cell[17] + V{0.0740740740740741}*cell[18] + V{0.0740740740740741}*cell[19] + V{0.0740740740740741}*cell[1] + V{0.0740740740740741}*cell[20] + V{0.0740740740740741}*cell[21] + V{0.0740740740740741}*cell[22] + V{0.0740740740740741}*cell[23] + V{0.0740740740740741}*cell[24] + V{0.0740740740740741}*cell[25] + V{0.0740740740740741}*cell[26] + V{0.0740740740740741}*cell[2] + V{0.0740740740740741}*cell[3] + V{0.0740740740740741}*cell[4] + V{0.0740740740740741}*cell[5] + V{0.0740740740740741}*cell[6] + V{0.0740740740740741}*cell[7] + V{0.0740740740740741}*cell[8] + V{0.0740740740740741}*cell[9] + V{0.0740740740740741};
auto x133 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x134 = V{4.5}*cell[13];
auto x135 = V{4.5}*cell[5];
auto x136 = V{4.5}*cell[7];
auto x137 = V{4.5}*cell[18];
auto x138 = V{4.5}*cell[20];
auto x139 = V{4.5}*cell[26];
auto x140 = V{4.5}*cell[10];
auto x141 = V{4.5}*cell[11];
auto x142 = -V{4.5}*cell[23];
auto x143 = V{4.5}*cell[24];
auto x144 = -V{4.5}*cell[17] + V{4.5}*cell[4] + x140 + x141 + x142 - x143;
auto x145 = V{4.5}*cell[12];
auto x146 = V{4.5}*cell[25];
auto x147 = -V{4.5}*cell[19] + V{4.5}*cell[6] + x145 - x146;
auto x148 = -V{4.5}*cell[14] + V{4.5}*cell[1] + x134 + x135 + x136 - x137 - x138 - x139 + x144 + x147;
auto x149 = -x148*x34;
auto x150 = x133 + x149;
auto x151 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x152 = V{3}*cell[13];
auto x153 = V{3}*cell[5];
auto x154 = V{3}*cell[7];
auto x155 = V{3}*cell[18];
auto x156 = V{3}*cell[20];
auto x157 = V{3}*cell[26];
auto x158 = V{3}*cell[10];
auto x159 = V{3}*cell[11];
auto x160 = -V{3}*cell[23];
auto x161 = V{3}*cell[24];
auto x162 = -V{3}*cell[17] + V{3}*cell[4] + x158 + x159 + x160 - x161;
auto x163 = V{3}*cell[12];
auto x164 = V{3}*cell[25];
auto x165 = -V{3}*cell[19] + V{3}*cell[6] + x163 - x164;
auto x166 = -V{3}*cell[14] + V{3}*cell[1] + x152 + x153 + x154 - x155 - x156 - x157 + x162 + x165;
auto x167 = -x166*x34;
auto x168 = x121 + x151 + x167;
auto x169 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x170 = V{6}*cell[13];
auto x171 = V{6}*cell[5];
auto x172 = V{6}*cell[7];
auto x173 = V{6}*cell[18];
auto x174 = V{6}*cell[20];
auto x175 = V{6}*cell[26];
auto x176 = V{6}*cell[10];
auto x177 = V{6}*cell[11];
auto x178 = -V{6}*cell[23];
auto x179 = V{6}*cell[24];
auto x180 = -V{6}*cell[17] + V{6}*cell[4] + x176 + x177 + x178 - x179;
auto x181 = V{6}*cell[12];
auto x182 = V{6}*cell[25];
auto x183 = -V{6}*cell[19] + V{6}*cell[6] + x181 - x182;
auto x184 = x34*(-V{6}*cell[14] + V{6}*cell[1] + x170 + x171 + x172 - x173 - x174 - x175 + x180 + x183);
auto x185 = -x169 + x184 + V{3};
auto x186 = V{0.074074}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x187 = V{0.222222}*x127;
auto x188 = V{0.222222}*x130;
auto x189 = x187 + x188;
auto x190 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x191 = V{4.5}*cell[9];
auto x192 = V{4.5}*cell[22];
auto x193 = -V{4.5}*cell[21] + V{4.5}*cell[8] - x134 + x139;
auto x194 = -V{4.5}*cell[15] + V{4.5}*cell[2] - x135 + x137 + x144 - x145 + x146 + x191 - x192 + x193;
auto x195 = -x194*x34;
auto x196 = x190 + x195;
auto x197 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x198 = V{3}*cell[9];
auto x199 = V{3}*cell[22];
auto x200 = -V{3}*cell[21] + V{3}*cell[8] - x152 + x157;
auto x201 = -V{3}*cell[15] + V{3}*cell[2] - x153 + x155 + x162 - x163 + x164 + x198 - x199 + x200;
auto x202 = -x201*x34;
auto x203 = x197 + x202;
auto x204 = x121 + x203;
auto x205 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x206 = V{6}*cell[9];
auto x207 = V{6}*cell[22];
auto x208 = -V{6}*cell[21] + V{6}*cell[8] - x170 + x175;
auto x209 = x34*(-V{6}*cell[15] + V{6}*cell[2] - x171 + x173 + x180 - x181 + x182 + x206 - x207 + x208);
auto x210 = -x205 + x209 + V{3};
auto x211 = V{0.074074}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x212 = V{0.222222}*x124;
auto x213 = x188 + x212;
auto x214 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x215 = -V{4.5}*cell[16] + V{4.5}*cell[3] - x136 + x138 + x140 - x141 + x142 + x143 + x147 - x191 + x192 + x193;
auto x216 = -x215*x34;
auto x217 = x214 + x216;
auto x218 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x219 = -V{3}*cell[16] + V{3}*cell[3] - x154 + x156 + x158 - x159 + x160 + x161 + x165 - x198 + x199 + x200;
auto x220 = -x219*x34;
auto x221 = x218 + x220;
auto x222 = x121 + x221;
auto x223 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x224 = x34*(-V{6}*cell[16] + V{6}*cell[3] - x172 + x174 + x176 - x177 + x178 + x179 + x183 - x206 + x207 + x208);
auto x225 = -x223 + x224 + V{3};
auto x226 = V{0.074074}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x227 = x187 + x212;
auto x228 = V{0.0185185185185185}*cell[0] + V{0.0185185185185185}*cell[10] + V{0.0185185185185185}*cell[11] + V{0.0185185185185185}*cell[12] + V{0.0185185185185185}*cell[13] + V{0.0185185185185185}*cell[14] + V{0.0185185185185185}*cell[15] + V{0.0185185185185185}*cell[16] + V{0.0185185185185185}*cell[17] + V{0.0185185185185185}*cell[18] + V{0.0185185185185185}*cell[19] + V{0.0185185185185185}*cell[1] + V{0.0185185185185185}*cell[20] + V{0.0185185185185185}*cell[21] + V{0.0185185185185185}*cell[22] + V{0.0185185185185185}*cell[23] + V{0.0185185185185185}*cell[24] + V{0.0185185185185185}*cell[25] + V{0.0185185185185185}*cell[26] + V{0.0185185185185185}*cell[2] + V{0.0185185185185185}*cell[3] + V{0.0185185185185185}*cell[4] + V{0.0185185185185185}*cell[5] + V{0.0185185185185185}*cell[6] + V{0.0185185185185185}*cell[7] + V{0.0185185185185185}*cell[8] + V{0.0185185185185185}*cell[9] + V{0.0185185185185185};
auto x229 = x111 + x116;
auto x230 = x150 + x196;
auto x231 = x168 + x203;
auto x232 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x233 = V{9}*cell[18];
auto x234 = V{9}*cell[25];
auto x235 = V{9}*cell[9];
auto x236 = V{9}*cell[12];
auto x237 = V{9}*cell[22];
auto x238 = V{9}*cell[5];
auto x239 = V{9}*cell[10];
auto x240 = V{9}*cell[11];
auto x241 = -V{9}*cell[23];
auto x242 = V{9}*cell[24];
auto x243 = -V{9}*cell[17] + V{9}*cell[4] + x239 + x240 + x241 - x242;
auto x244 = V{9}*cell[26];
auto x245 = V{9}*cell[13];
auto x246 = -V{9}*cell[21] + V{9}*cell[8] + x244 - x245;
auto x247 = x34*(-V{9}*cell[15] + V{9}*cell[2] + x233 + x234 + x235 - x236 - x237 - x238 + x243 + x246);
auto x248 = -x232 + x247;
auto x249 = x185 + x248;
auto x250 = V{0.018519}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x251 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x252 = V{9}*cell[7];
auto x253 = V{9}*cell[20];
auto x254 = -V{9}*cell[19] + V{9}*cell[6] - x234 + x236;
auto x255 = x34*(-V{9}*cell[14] + V{9}*cell[1] - x233 + x238 + x243 - x244 + x245 + x252 - x253 + x254);
auto x256 = -x251 + x255;
auto x257 = x210 + x256;
auto x258 = V{0.018519}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x259 = V{0.055556}*x130;
auto x260 = V{1}*x33;
auto x261 = -x112;
auto x262 = x111 - x115 + x261;
auto x263 = -x190;
auto x264 = x150 - x195 + x263;
auto x265 = -x197;
auto x266 = -x202 + x265;
auto x267 = x168 + x266;
auto x268 = V{0.166671}*cell[13];
auto x269 = V{0.166671}*cell[5];
auto x270 = V{0.166671}*cell[7];
auto x271 = V{0.166671}*cell[18];
auto x272 = V{0.166671}*cell[20];
auto x273 = V{0.166671}*cell[26];
auto x274 = V{0.166671}*cell[10];
auto x275 = V{0.166671}*cell[11];
auto x276 = -V{0.166671}*cell[23];
auto x277 = V{0.166671}*cell[24];
auto x278 = -V{0.166671}*cell[17] + V{0.166671}*cell[4] + x274 + x275 + x276 - x277;
auto x279 = V{0.166671}*cell[12];
auto x280 = V{0.166671}*cell[25];
auto x281 = -V{0.166671}*cell[19] + V{0.166671}*cell[6] + x279 - x280;
auto x282 = -V{0.0833355}*cell.template getFieldComponent<descriptors::FORCE>(0) + x34*(-V{0.166671}*cell[14] + V{0.166671}*cell[1] + x268 + x269 + x270 - x271 - x272 - x273 + x278 + x281) + V{0.055557};
auto x283 = V{0.111114}*cell[18];
auto x284 = V{0.111114}*cell[25];
auto x285 = V{0.111114}*cell[9];
auto x286 = V{0.111114}*cell[12];
auto x287 = V{0.111114}*cell[22];
auto x288 = V{0.111114}*cell[5];
auto x289 = V{0.111114}*cell[10];
auto x290 = V{0.111114}*cell[26];
auto x291 = V{0.111114}*cell[13];
auto x292 = -V{0.111114}*cell[23];
auto x293 = -V{0.111114}*cell[21] + V{0.111114}*cell[8] + x289 + x290 - x291 + x292;
auto x294 = V{0.111114}*cell[11];
auto x295 = V{0.111114}*cell[24];
auto x296 = -V{0.111114}*cell[17] + V{0.111114}*cell[4] + x294 - x295;
auto x297 = V{0.055557}*cell.template getFieldComponent<descriptors::FORCE>(1) - x34*(-V{0.111114}*cell[15] + V{0.111114}*cell[2] + x283 + x284 + x285 - x286 - x287 - x288 + x293 + x296);
auto x298 = x232 - x247;
auto x299 = x185 + x298;
auto x300 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x301 = x34*(-V{9}*cell[16] + V{9}*cell[3] - x235 + x237 + x239 - x240 + x241 + x242 + x246 - x252 + x253 + x254);
auto x302 = -x300 + x301;
auto x303 = x225 + x256;
auto x304 = V{0.018519}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x305 = V{0.055556}*x127;
auto x306 = -x117;
auto x307 = -x119 + x306;
auto x308 = x111 + x307;
auto x309 = -x214;
auto x310 = -x216 + x309;
auto x311 = x150 + x310;
auto x312 = -x218;
auto x313 = -x220 + x312;
auto x314 = V{0.111114}*cell[20];
auto x315 = V{0.111114}*cell[7];
auto x316 = -V{0.111114}*cell[19] + V{0.111114}*cell[6] - x284 + x286;
auto x317 = V{0.055557}*cell.template getFieldComponent<descriptors::FORCE>(2) - x34*(-V{0.111114}*cell[16] + V{0.111114}*cell[3] - x285 + x287 + x293 - x294 + x295 + x314 - x315 + x316);
auto x318 = x300 - x301;
auto x319 = x116 + x120;
auto x320 = x196 + x217;
auto x321 = x204 + x221;
auto x322 = x210 + x302;
auto x323 = x225 + x248;
auto x324 = V{0.055556}*x124;
auto x325 = x116 + x307;
auto x326 = x196 + x310;
auto x327 = V{0.166671}*cell[9];
auto x328 = V{0.166671}*cell[22];
auto x329 = -V{0.166671}*cell[21] + V{0.166671}*cell[8] - x268 + x273;
auto x330 = -V{0.0833355}*cell.template getFieldComponent<descriptors::FORCE>(1) + x34*(-V{0.166671}*cell[15] + V{0.166671}*cell[2] - x269 + x271 + x278 - x279 + x280 + x327 - x328 + x329) + V{0.055557};
auto x331 = x210 + x318;
auto x332 = V{0.00462962962962963}*cell[0] + V{0.00462962962962963}*cell[10] + V{0.00462962962962963}*cell[11] + V{0.00462962962962963}*cell[12] + V{0.00462962962962963}*cell[13] + V{0.00462962962962963}*cell[14] + V{0.00462962962962963}*cell[15] + V{0.00462962962962963}*cell[16] + V{0.00462962962962963}*cell[17] + V{0.00462962962962963}*cell[18] + V{0.00462962962962963}*cell[19] + V{0.00462962962962963}*cell[1] + V{0.00462962962962963}*cell[20] + V{0.00462962962962963}*cell[21] + V{0.00462962962962963}*cell[22] + V{0.00462962962962963}*cell[23] + V{0.00462962962962963}*cell[24] + V{0.00462962962962963}*cell[25] + V{0.00462962962962963}*cell[26] + V{0.00462962962962963}*cell[2] + V{0.00462962962962963}*cell[3] + V{0.00462962962962963}*cell[4] + V{0.00462962962962963}*cell[5] + V{0.00462962962962963}*cell[6] + V{0.00462962962962963}*cell[7] + V{0.00462962962962963}*cell[8] + V{0.00462962962962963}*cell[9] + V{0.00462962962962963};
auto x333 = V{0.00463}*x33;
auto x334 = x223 - x224 + V{3};
auto x335 = x205 - x209 + V{3};
auto x336 = x166*x34;
auto x337 = x128 + x306;
auto x338 = x125 + x261;
auto x339 = x133 - x148*x34;
auto x340 = x215*x34;
auto x341 = x309 + x340;
auto x342 = x194*x34;
auto x343 = x263 + x342;
auto x344 = x123*x123;
auto x345 = x126*x126;
auto x346 = x129*x129;
auto x347 = V{1.5}*x344 + V{1.5}*x345 + V{1.5}*x346 + V{-1};
auto x348 = x201*x34;
auto x349 = x265 + x348;
auto x350 = x347 + x349;
auto x351 = x219*x34;
auto x352 = x312 + x351;
auto x353 = x350 + x352;
auto x354 = x318 + x335;
auto x355 = x298 + x334;
auto x356 = -x151;
auto x357 = x336 + x347 + x356;
auto x358 = x169 - x184 + V{3};
auto x359 = x190 - x342;
auto x360 = x214 - x340;
auto x361 = x123 + x126;
auto x362 = x339 + x359;
auto x363 = x349 + x357;
auto x364 = x298 + x358;
auto x365 = x251 - x255;
auto x366 = x335 + x365;
auto x367 = -x167 + x356;
auto x368 = V{0.055557}*cell.template getFieldComponent<descriptors::FORCE>(0) - x34*(-V{0.111114}*cell[14] + V{0.111114}*cell[1] - x283 + x288 + x289 - x290 + x291 + x292 + x296 - x314 + x315 + x316);
auto x369 = x123 + x129;
auto x370 = x339 + x360;
auto x371 = x352 + x357;
auto x372 = x318 + x358;
auto x373 = x334 + x365;
auto x374 = -V{0.0833355}*cell.template getFieldComponent<descriptors::FORCE>(2) + x34*(-V{0.166671}*cell[16] + V{0.166671}*cell[3] - x270 + x272 + x274 - x275 + x276 + x277 + x281 - x327 + x328 + x329) + V{0.055557};
auto x375 = x225 + x365;
auto x0 = V{1}*cell[0]*x122 - x105*(x121*(V{0.296296296296296}*cell[0] + V{0.296296296296296}*cell[10] + V{0.296296296296296}*cell[11] + V{0.296296296296296}*cell[12] + V{0.296296296296296}*cell[13] + V{0.296296296296296}*cell[14] + V{0.296296296296296}*cell[15] + V{0.296296296296296}*cell[16] + V{0.296296296296296}*cell[17] + V{0.296296296296296}*cell[18] + V{0.296296296296296}*cell[19] + V{0.296296296296296}*cell[1] + V{0.296296296296296}*cell[20] + V{0.296296296296296}*cell[21] + V{0.296296296296296}*cell[22] + V{0.296296296296296}*cell[23] + V{0.296296296296296}*cell[24] + V{0.296296296296296}*cell[25] + V{0.296296296296296}*cell[26] + V{0.296296296296296}*cell[2] + V{0.296296296296296}*cell[3] + V{0.296296296296296}*cell[4] + V{0.296296296296296}*cell[5] + V{0.296296296296296}*cell[6] + V{0.296296296296296}*cell[7] + V{0.296296296296296}*cell[8] + V{0.296296296296296}*cell[9] + V{0.296296296296296}) + V{0.296296296296296}) - V{0.888889}*x131*x33*(x124 + x127 + x130);
auto x1 = V{1}*cell[1]*x122 - x105*(x132*(-x111*x150 + x168) + V{0.0740740740740741}) - x131*x33*(x185*x186 + x189);
auto x2 = V{1}*cell[2]*x122 - x105*(x132*(-x116*x196 + x204) + V{0.0740740740740741}) - x131*x33*(x210*x211 + x213);
auto x3 = V{1}*cell[3]*x122 - x105*(x132*(-x120*x217 + x222) + V{0.0740740740740741}) - x131*x33*(x225*x226 + x227);
auto x4 = V{1}*cell[4]*x122 - x105*(x228*(-x229*x230 + x231) + V{0.0185185185185185}) - x131*x260*(x249*x250 + x257*x258 + x259);
auto x5 = V{1}*cell[5]*x122 - x105*(x228*(-x262*x264 + x267) + V{0.0185185185185185}) - x131*x260*(-cell.template getFieldComponent<descriptors::FORCE>(1)*(x282 + x297) + x250*x299 + x259);
auto x6 = V{1}*cell[6]*x122 - x105*(x228*(x168 + x221 - (x111 + x120)*(x150 + x217)) + V{0.0185185185185185}) - x131*x260*(x250*(x185 + x302) + x303*x304 + x305);
auto x7 = V{1}*cell[7]*x122 - x105*(x228*(x168 - x308*x311 + x313) + V{0.0185185185185185}) - x131*x260*(-cell.template getFieldComponent<descriptors::FORCE>(2)*(x282 + x317) + x250*(x185 + x318) + x305);
auto x8 = V{1}*cell[8]*x122 - x105*(x228*(-x319*x320 + x321) + V{0.0185185185185185}) - x131*x260*(x258*x322 + x304*x323 + x324);
auto x9 = V{1}*cell[9]*x122 - x105*(x228*(x204 + x313 - x325*x326) + V{0.0185185185185185}) - x131*x260*(-cell.template getFieldComponent<descriptors::FORCE>(2)*(x317 + x330) + x258*x331 + x324);
auto x10 = V{1}*cell[10]*x122 - x105*(x332*(x221 + x231 - (x120 + x229)*(x217 + x230)) + V{0.00462962962962963}) - x131*x333*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x249 + x302) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x257 + x302) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x248 + x303));
auto x11 = -x105*(x332*(x231 + x313 - (-x229 - x307)*(-x230 - x310)) + V{0.00462962962962963}) + x122*x61 + x131*x333*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x249 + x318) - cell.template getFieldComponent<descriptors::FORCE>(1)*(x257 + x318) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x248 + x256 + x334));
auto x12 = -x105*(x332*(x221 + x267 - (-x120 - x262)*(-x217 - x264)) + V{0.00462962962962963}) + x122*x67 + x131*x333*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x299 + x302) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x256 + x302 + x335) - cell.template getFieldComponent<descriptors::FORCE>(2)*(x298 + x303));
auto x13 = -x105*(x332*(x151 - x336 + x353 - (-x123 - x337 - x338)*(-x339 - x341 - x343)) + V{0.00462962962962963}) + x122*x73 + x131*x333*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x299 + x318) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x256 + x354) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x256 + x355));
auto x14 = V{1}*cell[14]*x122 - x105*(x132*(-x123*x339 + x357) + V{0.0740740740740741}) - x131*x33*(-x186*x358 + x189);
auto x15 = V{1}*cell[15]*x122 - x105*(x132*(-x126*x359 + x350) + V{0.0740740740740741}) - x131*x33*(-x211*x335 + x213);
auto x16 = V{1}*cell[16]*x122 - x105*(x132*(-x129*x360 + x347 + x352) + V{0.0740740740740741}) - x131*x33*(-x226*x334 + x227);
auto x17 = V{1}*cell[17]*x122 - x105*(x228*(-x361*x362 + x363) + V{0.0185185185185185}) + x131*x260*(x250*x364 + x258*x366 - x259);
auto x18 = V{1}*cell[18]*x122 - x105*(x228*(x204 - x262*x264 + x367) + V{0.0185185185185185}) - x131*x260*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x330 + x368) + x258*(x210 + x365) + x259);
auto x19 = V{1}*cell[19]*x122 - x105*(x228*(-x369*x370 + x371) + V{0.0185185185185185}) + x131*x260*(x250*x372 + x304*x373 - x305);
auto x20 = V{1}*cell[20]*x122 - x105*(x228*(x222 - x308*x311 + x367) + V{0.0185185185185185}) - x131*x260*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x368 + x374) + x304*x375 + x305);
auto x21 = V{1}*cell[21]*x122 - x105*(x228*(x353 - (x126 + x129)*(x359 + x360)) + V{0.0185185185185185}) + x131*x260*(x258*x354 + x304*x355 - x324);
auto x22 = V{1}*cell[22]*x122 - x105*(x228*(x222 + x266 - x325*x326) + V{0.0185185185185185}) - x131*x260*(-cell.template getFieldComponent<descriptors::FORCE>(1)*(x297 + x374) + x304*(x225 + x298) + x324);
auto x23 = V{1}*cell[23]*x122 - x105*(x332*(x352 + x363 - (x129 + x361)*(x360 + x362)) + V{0.00462962962962963}) + x131*x333*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x318 + x364) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x318 + x366) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x298 + x373));
auto x24 = -x105*(x332*(x218 - x351 + x363 - (x337 + x361)*(x341 + x362)) + V{0.00462962962962963}) + x122*x78 + x131*x333*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x302 + x364) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x302 + x366) - cell.template getFieldComponent<descriptors::FORCE>(2)*(x298 + x375));
auto x25 = -x105*(x332*(x197 - x348 + x371 - (x338 + x369)*(x343 + x370)) + V{0.00462962962962963}) + x122*x68 + x131*x333*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x248 + x372) - cell.template getFieldComponent<descriptors::FORCE>(1)*(x331 + x365) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x248 + x373));
auto x26 = -x105*(x332*(x321 + x367 - (x106 + x110 - x319)*(x133 + x149 - x320)) + V{0.00462962962962963}) + x122*x74 + x131*x333*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(-x248 - x302 - x358) + cell.template getFieldComponent<descriptors::FORCE>(1)*(-x322 - x365) - cell.template getFieldComponent<descriptors::FORCE>(2)*(x323 + x365));
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
cell[19] = x19;
cell[20] = x20;
cell[21] = x21;
cell[22] = x22;
cell[23] = x23;
cell[24] = x24;
cell[25] = x25;
cell[26] = x26;
return { x33, x344 + x345 + x346 };
}
};

}

}
