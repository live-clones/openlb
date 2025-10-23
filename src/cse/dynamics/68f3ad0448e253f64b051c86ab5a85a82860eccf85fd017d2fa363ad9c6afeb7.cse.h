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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::ThirdOrderRLB, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell[13] + cell[15];
auto x22 = cell[10] + x21;
auto x23 = cell[11] + cell[5];
auto x24 = cell[12] + cell[7] + cell[9];
auto x25 = cell[2] + cell[8];
auto x26 = cell[18] + cell[3];
auto x27 = cell[16] + cell[6];
auto x28 = cell[0] + cell[14] + cell[17] + cell[1] + cell[4] + x22 + x23 + x24 + x25 + x26 + x27;
auto x29 = x28 + V{1};
auto x30 = V{1} / (x29);
auto x31 = V{0.5}*x30;
auto x32 = -cell[1];
auto x33 = -cell[7];
auto x34 = -cell[6];
auto x35 = cell[16] + x34;
auto x36 = x32 + x33 + x35;
auto x37 = -cell[4];
auto x38 = cell[14] + x37;
auto x39 = -cell[5] + x38;
auto x40 = x22 + x36 + x39;
auto x41 = x40*x40;
auto x42 = cell[13] + cell[17];
auto x43 = -cell[2];
auto x44 = -cell[9];
auto x45 = -cell[8];
auto x46 = cell[18] + x45;
auto x47 = x43 + x44 + x46;
auto x48 = -cell[14] + x23 + x37;
auto x49 = x42 + x47 + x48;
auto x50 = x49*x49;
auto x51 = cell[15] + cell[17];
auto x52 = -cell[18];
auto x53 = -cell[3];
auto x54 = x45 + x52 + x53;
auto x55 = -cell[16] + x34;
auto x56 = x24 + x51 + x54 + x55;
auto x57 = x56*x56;
auto x58 = V{0.333333333333333}*cell[0];
auto x59 = V{0.333333333333333}*cell[10];
auto x60 = V{0.333333333333333}*cell[11];
auto x61 = V{0.333333333333333}*cell[12];
auto x62 = V{0.333333333333333}*cell[13];
auto x63 = V{0.333333333333333}*cell[14];
auto x64 = V{0.333333333333333}*cell[15];
auto x65 = V{0.333333333333333}*cell[16];
auto x66 = V{0.333333333333333}*cell[17];
auto x67 = V{0.333333333333333}*cell[18];
auto x68 = V{0.333333333333333}*cell[1];
auto x69 = V{0.333333333333333}*cell[2];
auto x70 = V{0.333333333333333}*cell[3];
auto x71 = V{0.333333333333333}*cell[4];
auto x72 = V{0.333333333333333}*cell[5];
auto x73 = V{0.333333333333333}*cell[6];
auto x74 = V{0.333333333333333}*cell[7];
auto x75 = V{0.333333333333333}*cell[8];
auto x76 = V{0.333333333333333}*cell[9];
auto x77 = V{1} / ((x29)*(x29));
auto x78 = V{1.5}*x77;
auto x79 = x57*x78;
auto x80 = x41*x78;
auto x81 = x50*x78;
auto x82 = x80 + x81 + V{-1};
auto x83 = x79 + x82;
auto x84 = x31*x40;
auto x85 = V{0.666666666666667}*cell[12];
auto x86 = V{0.666666666666667}*cell[3];
auto x87 = V{1}*x30;
auto x88 = x57*x87;
auto x89 = V{0.666666666666667}*cell[17];
auto x90 = V{0.666666666666667}*cell[18];
auto x91 = V{0.666666666666667}*cell[8];
auto x92 = V{0.666666666666667}*cell[9];
auto x93 = x58 + x59 + x68 - x89 - x90 - x91 - x92;
auto x94 = V{0.666666666666667}*cell[15];
auto x95 = V{0.666666666666667}*cell[16];
auto x96 = V{0.666666666666667}*cell[6];
auto x97 = V{0.666666666666667}*cell[7];
auto x98 = x60 + x69 - x94 - x95 - x96 - x97;
auto x99 = x62 + x63 + x71 + x72 - x85 - x86 + x88 + x93 + x98;
auto x100 = V{0.666666666666667}*cell[11];
auto x101 = V{0.666666666666667}*cell[2];
auto x102 = x50*x87;
auto x103 = V{0.666666666666667}*cell[13];
auto x104 = V{0.666666666666667}*cell[14];
auto x105 = V{0.666666666666667}*cell[4];
auto x106 = V{0.666666666666667}*cell[5];
auto x107 = -x103 - x104 - x105 - x106 + x61 + x70;
auto x108 = -x100 - x101 + x102 + x107 + x64 + x65 + x73 + x74 + x93;
auto x109 = x30*x40;
auto x110 = x109*x49;
auto x111 = -cell[13] + cell[5] + x110 + x38;
auto x112 = x111*x87;
auto x113 = -cell[15];
auto x114 = x109*x56;
auto x115 = cell[7] + x113 + x114 + x35;
auto x116 = x56*x87;
auto x117 = V{0.166666666666667}*x30;
auto x118 = V{6.93889390390723e-18}*cell[0];
auto x119 = V{0.0833333333333333}*cell[12];
auto x120 = V{0.0833333333333333}*cell[3];
auto x121 = V{0.0833333333333333}*x30;
auto x122 = -V{0.0833333333333333}*cell[13] - V{0.0833333333333333}*cell[14] - V{0.0833333333333333}*cell[4] - V{0.0833333333333333}*cell[5] + x118 + x119 + x120 - x121*x57;
auto x123 = V{0.0833333333333333}*cell[11];
auto x124 = V{0.0833333333333333}*cell[2];
auto x125 = -V{0.0833333333333333}*cell[15] - V{0.0833333333333333}*cell[16] - V{0.0833333333333333}*cell[6] - V{0.0833333333333333}*cell[7] - x121*x50 + x123 + x124;
auto x126 = -V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] - V{0.166666666666667}*cell[1] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x117*x41 + x122 + x125;
auto x127 = x28 + V{1};
auto x128 = V{3}*x77;
auto x129 = x128*x41;
auto x130 = -x79;
auto x131 = V{1} - x81;
auto x132 = util::pow(x29, -3);
auto x133 = V{6.000003}*x132;
auto x134 = x133*x40;
auto x135 = V{2.999997}*x132;
auto x136 = x135*x49;
auto x137 = V{3}*cell[14];
auto x138 = V{3}*cell[16];
auto x139 = V{3}*cell[5];
auto x140 = V{3}*cell[7];
auto x141 = V{3}*cell[13] - V{3}*cell[4];
auto x142 = V{3}*cell[15] - V{3}*cell[6];
auto x143 = x30*(V{3}*cell[10] - V{3}*cell[1] + x137 + x138 - x139 - x140 + x141 + x142);
auto x144 = -x143;
auto x145 = x134*x50;
auto x146 = x136*x41;
auto x147 = x144 + x145 + x146;
auto x148 = x134*x57 + x136*x57 + x147;
auto x149 = x31*x49;
auto x150 = V{0.666666666666667}*cell[10];
auto x151 = V{0.666666666666667}*cell[1];
auto x152 = x41*x87;
auto x153 = x107 - x150 - x151 + x152 + x58 + x66 + x67 + x75 + x76 + x98;
auto x154 = x30*x49;
auto x155 = x154*x56;
auto x156 = -cell[17];
auto x157 = cell[9] + x156;
auto x158 = x155 + x157 + x46;
auto x159 = V{0.0833333333333333}*cell[10];
auto x160 = V{0.0833333333333333}*cell[1];
auto x161 = -V{0.0833333333333333}*cell[17] - V{0.0833333333333333}*cell[18] - V{0.0833333333333333}*cell[8] - V{0.0833333333333333}*cell[9] - x121*x41 + x159 + x160;
auto x162 = -V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] - V{0.166666666666667}*cell[2] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] + x117*x50 + x122 + x161;
auto x163 = -x80;
auto x164 = x128*x50;
auto x165 = x133*x49;
auto x166 = x135*x40;
auto x167 = V{3}*cell[18];
auto x168 = V{3}*cell[9];
auto x169 = V{3}*cell[17] - V{3}*cell[8];
auto x170 = x30*(V{3}*cell[11] - V{3}*cell[2] - x137 + x139 + x141 + x167 - x168 + x169);
auto x171 = -x170;
auto x172 = x165*x41;
auto x173 = x166*x50;
auto x174 = x171 + x172 + x173;
auto x175 = x165*x57 + x166*x57 + x174;
auto x176 = x31*x56;
auto x177 = x40*x87;
auto x178 = x49*x87;
auto x179 = -V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] - V{0.166666666666667}*cell[3] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] + x117*x57 + x118 + x125 + x161;
auto x180 = x128*x57;
auto x181 = x30*(V{3}*cell[12] - V{3}*cell[3] - x138 + x140 + x142 - x167 + x168 + x169);
auto x182 = -x181;
auto x183 = V{9}*x132;
auto x184 = x183*x56;
auto x185 = x182 + x184*x41 + x184*x50;
auto x186 = x131 + x163;
auto x187 = -x119;
auto x188 = -x120;
auto x189 = V{0.0416666666666667}*x30;
auto x190 = x189*x57;
auto x191 = V{0.25000025}*x109;
auto x192 = x108*x191;
auto x193 = V{0.5000005}*x111;
auto x194 = x154*x193;
auto x195 = V{2.49999999985601e-07}*x109;
auto x196 = x195*x99;
auto x197 = x30*x56;
auto x198 = V{4.99999999971202e-07}*x197;
auto x199 = x115*x198;
auto x200 = -V{0.0416666666666667}*cell[0];
auto x201 = V{0.0833333333333333}*x30;
auto x202 = V{0.0416666666666667}*cell[10] + V{6.93889390390723e-18}*cell[17] + V{6.93889390390723e-18}*cell[18] + V{0.0416666666666667}*cell[1] + V{6.93889390390723e-18}*cell[8] + V{6.93889390390723e-18}*cell[9] + x200 - x201*x41;
auto x203 = V{0.0416666666666667}*cell[11] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[16] + V{0.0416666666666667}*cell[2] + V{6.93889390390723e-18}*cell[6] + V{6.93889390390723e-18}*cell[7] - x201*x50;
auto x204 = x187 + x188 + x190 + x192 + x194 - x196 - x199 + x202 + x203;
auto x205 = V{0.25000025}*x154;
auto x206 = x153*x205;
auto x207 = x109*x193;
auto x208 = V{2.49999999985601e-07}*x154;
auto x209 = x208*x99;
auto x210 = x158*x198;
auto x211 = x206 + x207 - x209 - x210;
auto x212 = V{0.25}*x110;
auto x213 = V{0.375}*cell[13] - V{0.125}*cell[14] + V{0.375}*cell[4] - V{0.125}*cell[5] - x212;
auto x214 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x215 = V{4.5}*x77;
auto x216 = cell[10] + x36;
auto x217 = cell[11] + V{2}*cell[13] - V{2}*cell[4] + x216 + x47 + x51;
auto x218 = x215*(x217*x217);
auto x219 = x143 + x83;
auto x220 = V{18}*x132;
auto x221 = x170 - x183*x40*x57 - x183*x49*x57 + x220*x40*x50 + x220*x41*x49;
auto x222 = -x206 - x207 + x209 + x210;
auto x223 = -V{0.125}*cell[13] + V{0.375}*cell[14] - V{0.125}*cell[4] + V{0.375}*cell[5] + x212;
auto x224 = -cell[11] + V{2}*cell[14] + cell[15] - V{2}*cell[5] + x157 + x216 + x25 + x52;
auto x225 = -x224;
auto x226 = V{3.000006}*x132;
auto x227 = x226*x49*x57;
auto x228 = V{6.000012}*x132;
auto x229 = x228*x40*x50;
auto x230 = x226*x40*x57;
auto x231 = x228*x41*x49;
auto x232 = -x123;
auto x233 = -x124;
auto x234 = x189*x50;
auto x235 = x191*x99;
auto x236 = V{0.5000005}*x115;
auto x237 = x197*x236;
auto x238 = x108*x195;
auto x239 = V{4.99999999971202e-07}*x154;
auto x240 = x111*x239;
auto x241 = V{0.0416666666666667}*cell[12] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[14] + V{0.0416666666666667}*cell[3] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[5] - x201*x57;
auto x242 = x202 + x232 + x233 + x234 + x235 + x237 - x238 - x240 + x241;
auto x243 = V{0.25000025}*x197;
auto x244 = x153*x243;
auto x245 = x109*x236;
auto x246 = V{2.49999999985601e-07}*x197;
auto x247 = x108*x246;
auto x248 = x158*x239;
auto x249 = x244 + x245 - x247 - x248;
auto x250 = V{0.25}*x114;
auto x251 = V{0.375}*cell[15] - V{0.125}*cell[16] + V{0.375}*cell[6] - V{0.125}*cell[7] - x250;
auto x252 = cell[10] + x32 + x39;
auto x253 = cell[12] + V{2}*cell[15] - V{2}*cell[6] + cell[9] + x252 + x42 + x54;
auto x254 = x215*(x253*x253);
auto x255 = x132*x56;
auto x256 = V{9.000009}*x255;
auto x257 = x256*x41;
auto x258 = V{8.99999999948164e-06}*x255;
auto x259 = x258*x50;
auto x260 = x181 + x257 - x259;
auto x261 = V{5.999994}*x132;
auto x262 = x261*x49*x57;
auto x263 = V{12.000006}*x132;
auto x264 = x263*x40*x57;
auto x265 = -x145;
auto x266 = -x146;
auto x267 = x219 + x262 + x264 + x265 + x266;
auto x268 = -x244 - x245 + x247 + x248;
auto x269 = -V{0.125}*cell[15] + V{0.375}*cell[16] - V{0.125}*cell[6] + V{0.375}*cell[7] + x250;
auto x270 = -cell[12] + cell[13];
auto x271 = V{2}*cell[16] - V{2}*cell[7] + cell[8] + x156 + x252 + x26 + x270 + x44;
auto x272 = -x271;
auto x273 = -x159;
auto x274 = -x160;
auto x275 = x189*x41;
auto x276 = x205*x99;
auto x277 = V{0.5000005}*x158;
auto x278 = x197*x277;
auto x279 = x153*x208;
auto x280 = V{4.99999999971202e-07}*x109;
auto x281 = x111*x280;
auto x282 = x200 + x203 + x241 + x273 + x274 + x275 + x276 + x278 - x279 - x281;
auto x283 = x108*x243;
auto x284 = x154*x277;
auto x285 = x153*x246;
auto x286 = x115*x280;
auto x287 = x283 + x284 - x285 - x286;
auto x288 = V{0.25}*x155;
auto x289 = V{0.375}*cell[17] - V{0.125}*cell[18] + V{0.375}*cell[8] - V{0.125}*cell[9] - x288;
auto x290 = x43 + x48;
auto x291 = cell[12] + V{2}*cell[17] + cell[7] - V{2}*cell[8] + x21 + x290 + x53 + x55;
auto x292 = x215*(x291*x291);
auto x293 = x256*x50;
auto x294 = x258*x41;
auto x295 = x181 + x293 - x294;
auto x296 = x261*x40*x57;
auto x297 = x263*x49*x57;
auto x298 = -x172;
auto x299 = -x173;
auto x300 = x170 + x83;
auto x301 = x296 + x297 + x298 + x299 + x300;
auto x302 = -x283 - x284 + x285 + x286;
auto x303 = -V{0.125}*cell[17] + V{0.375}*cell[18] - V{0.125}*cell[8] + V{0.375}*cell[9] + x288;
auto x304 = V{2}*cell[18] + cell[3] - V{2}*cell[9] + x113 + x27 + x270 + x290 + x33;
auto x305 = -x304;
auto x306 = -x58;
auto x307 = x306 - x59 - x68 + x89 + x90 + x91 + x92;
auto x308 = x103 + x104 + x105 + x106 - x61 - x70;
auto x309 = x100 + x101 - x102 + x307 + x308 - x64 - x65 - x73 - x74;
auto x310 = -x60 - x69 + x94 + x95 + x96 + x97;
auto x311 = x307 + x310 - x62 - x63 - x71 - x72 + x85 + x86 - x88;
auto x312 = -x111;
auto x313 = -x115;
auto x314 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x315 = x79 + V{-1};
auto x316 = x150 + x151 - x152 + x306 + x308 + x310 - x66 - x67 - x75 - x76;
auto x317 = -x158;
auto x318 = x187 + x188 + x190 - x192 - x194 + x196 + x199 + x202 + x203;
auto x319 = x130 + x186;
auto x320 = x143 + x319;
auto x321 = x202 + x232 + x233 + x234 - x235 - x237 + x238 + x240 + x241;
auto x322 = x200 + x203 + x241 + x273 + x274 + x275 - x276 - x278 + x279 + x281;
auto x0 = x20*(-V{0.5}*cell[0] + V{4.16333634234434e-17}*cell[10] + V{4.16333634234434e-17}*cell[11] + V{4.16333634234434e-17}*cell[12] + V{0.5}*cell[13] + V{0.5}*cell[14] + V{0.5}*cell[15] + V{0.5}*cell[16] + V{0.5}*cell[17] + V{0.5}*cell[18] + V{4.16333634234434e-17}*cell[1] + V{4.16333634234434e-17}*cell[2] + V{4.16333634234434e-17}*cell[3] + V{0.5}*cell[4] + V{0.5}*cell[5] + V{0.5}*cell[6] + V{0.5}*cell[7] + V{0.5}*cell[8] + V{0.5}*cell[9] - x31*x41 - x31*x50 - x31*x57) - x83*(x58 + x59 + x60 + x61 + x62 + x63 + x64 + x65 + x66 + x67 + x68 + x69 + x70 + x71 + x72 + x73 + x74 + x75 + x76 + V{0.333333333333333}) + V{-0.333333333333333};
auto x1 = V{0.0555555555555556}*x127*(x129 + x130 + x131 + x148) - x20*(-x108*x84 - x112*x49 - x115*x116 - x126 - x84*x99) + V{-0.0555555555555556};
auto x2 = V{0.0555555555555556}*x127*(x130 + x163 + x164 + x175 + V{1}) - x20*(-x112*x40 - x116*x158 - x149*x153 - x149*x99 - x162) + V{-0.0555555555555556};
auto x3 = V{0.0555555555555556}*x127*(x180 + x185 + x186) - x20*(-x108*x176 - x115*x177 - x153*x176 - x158*x178 - x179) + V{-0.0555555555555556};
auto x4 = -x20*(x204 + x211 + x213) - x214*(-x218 + x219 + x221) + V{-0.0277777777777778};
auto x5 = -(x20*(x204 + x222 + x223) + x214*(x171 - x215*x225*x225 + x219 + x227 + x229 - x230 - x231) + V{0.0277777777777778});
auto x6 = -x20*(x242 + x249 + x251) - x214*(-x254 + x260 + x267) + V{-0.0277777777777778};
auto x7 = -(x20*(x242 + x268 + x269) + x214*(x182 - x215*x272*x272 - x257 + x259 + x267) + V{0.0277777777777778});
auto x8 = -x20*(x282 + x287 + x289) - x214*(-x292 + x295 + x301) + V{-0.0277777777777778};
auto x9 = -(x20*(x282 + x302 + x303) + x214*(x182 - x215*x305*x305 - x293 + x294 + x301) + V{0.0277777777777778});
auto x10 = x20*(x116*x313 + x126 + x178*x312 + x309*x84 + x311*x84) - x314*(-x129 + x148 + x315 + x81) + V{-0.0555555555555556};
auto x11 = x20*(x116*x317 + x149*x311 + x149*x316 + x162 + x177*x312) - x314*(-x164 + x175 + x315 + x80) + V{-0.0555555555555556};
auto x12 = x20*(x176*x309 + x176*x316 + x177*x313 + x178*x317 + x179) - x314*(-x180 + x185 + x82) + V{-0.0555555555555556};
auto x13 = V{0.0277777777777778}*x127*(x218 + x221 + x320) - x20*(x213 + x222 + x318) + V{-0.0277777777777778};
auto x14 = -(x20*(x211 + x223 + x318) + x214*(x144 - x215*x224*x224 - x227 - x229 + x230 + x231 + x300) + V{0.0277777777777778});
auto x15 = V{0.0277777777777778}*x127*(x254 + x260 + x262 + x264 + x265 + x266 + x320) - x20*(x251 + x268 + x321) + V{-0.0277777777777778};
auto x16 = -(x20*(x249 + x269 + x321) + x214*(x147 - x215*x271*x271 + x260 - x262 - x264 + x83) + V{0.0277777777777778});
auto x17 = V{0.0277777777777778}*x127*(x170 + x292 + x295 + x296 + x297 + x298 + x299 + x319) - x20*(x289 + x302 + x322) + V{-0.0277777777777778};
auto x18 = -(x20*(x287 + x303 + x322) + x214*(x174 - x215*x304*x304 + x295 - x296 - x297 + x83) + V{0.0277777777777778});
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
return { x29, V{1}*x77*(x41 + x50 + x57) };
}
};

}

}
