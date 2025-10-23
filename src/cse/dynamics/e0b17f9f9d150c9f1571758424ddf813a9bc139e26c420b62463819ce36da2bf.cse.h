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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::ParameterFromCell<collision::LES::SMAGORINSKY, collision::SmagorinskyEffectiveOmega<collision::RLBThirdOrder> > >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x24 = cell.template getFieldComponent<descriptors::WMPOROSITY>(0);
auto x21 = cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x22 = cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x25 = cell.template getFieldComponent<descriptors::WMVELOCITY>(0);
auto x19 = cell.template getFieldComponent<collision::LES::SMAGORINSKY>(0);
auto x26 = cell.template getFieldComponent<descriptors::WMVELOCITY>(1);
auto x28 = parameters.template get<descriptors::OMEGA>();
auto x23 = cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x27 = cell.template getFieldComponent<descriptors::WMVELOCITY>(2);
auto x20 = cell.template getFieldComponent<descriptors::POROSITY>(0);
auto x29 = V{0.5}/x28;
auto x30 = V{0.0277777691819762}/((x28)*(x28));
auto x31 = cell[12] + cell[17] + cell[7];
auto x32 = cell[11] + cell[13] + cell[18];
auto x33 = cell[10] + cell[14] + cell[16];
auto x34 = cell[0] + cell[15] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[8] + cell[9] + x31 + x32 + x33;
auto x35 = x34 + V{1};
auto x36 = V{1} / (x35);
auto x37 = x19*x19;
auto x38 = x20 + V{-1};
auto x39 = -x38;
auto x40 = -cell[8];
auto x41 = cell[9] + x40;
auto x42 = cell[15] - cell[6];
auto x43 = V{1}*x20*x36;
auto x44 = x43*(-cell[16] - cell[18] - cell[3] + x31 + x41 + x42);
auto x45 = x23*x39 + x44;
auto x46 = x24*x45;
auto x47 = x24 + V{-1};
auto x48 = -x47;
auto x49 = x27*x48;
auto x50 = x46 + x49;
auto x51 = x50*x50;
auto x52 = x35*x51;
auto x53 = V{0.333333333333333}*cell[13];
auto x54 = V{0.333333333333333}*cell[14];
auto x55 = V{0.333333333333333}*cell[4];
auto x56 = V{0.333333333333333}*cell[5];
auto x57 = V{0.666666666666667}*cell[12];
auto x58 = V{0.666666666666667}*cell[3];
auto x59 = V{0.333333333333333}*cell[0];
auto x60 = V{0.333333333333333}*cell[10];
auto x61 = V{0.333333333333333}*cell[1];
auto x62 = V{0.666666666666667}*cell[17];
auto x63 = V{0.666666666666667}*cell[18];
auto x64 = V{0.666666666666667}*cell[8];
auto x65 = V{0.666666666666667}*cell[9];
auto x66 = x59 + x60 + x61 - x62 - x63 - x64 - x65;
auto x67 = V{0.333333333333333}*cell[11];
auto x68 = V{0.333333333333333}*cell[2];
auto x69 = V{0.666666666666667}*cell[15];
auto x70 = V{0.666666666666667}*cell[16];
auto x71 = V{0.666666666666667}*cell[6];
auto x72 = V{0.666666666666667}*cell[7];
auto x73 = x67 + x68 - x69 - x70 - x71 - x72;
auto x74 = x53 + x54 + x55 + x56 - x57 - x58 + x66 + x73;
auto x75 = x52 + x74;
auto x76 = -cell[4];
auto x77 = cell[5] + x76;
auto x78 = x43*(-cell[14] + cell[17] - cell[2] - cell[9] + x32 + x40 + x77);
auto x79 = x22*x39 + x78;
auto x80 = x24*x79;
auto x81 = x26*x48;
auto x82 = x80 + x81;
auto x83 = x82*x82;
auto x84 = x35*x83;
auto x85 = V{0.333333333333333}*cell[15];
auto x86 = V{0.333333333333333}*cell[16];
auto x87 = V{0.333333333333333}*cell[6];
auto x88 = V{0.333333333333333}*cell[7];
auto x89 = V{0.666666666666667}*cell[11];
auto x90 = V{0.666666666666667}*cell[2];
auto x91 = V{0.333333333333333}*cell[12];
auto x92 = V{0.333333333333333}*cell[3];
auto x93 = V{0.666666666666667}*cell[13];
auto x94 = V{0.666666666666667}*cell[14];
auto x95 = V{0.666666666666667}*cell[4];
auto x96 = V{0.666666666666667}*cell[5];
auto x97 = x91 + x92 - x93 - x94 - x95 - x96;
auto x98 = x66 + x85 + x86 + x87 + x88 - x89 - x90 + x97;
auto x99 = x84 + x98;
auto x100 = x43*(cell[13] - cell[1] - cell[5] - cell[7] + x33 + x42 + x76);
auto x101 = x100 + x21*x39;
auto x102 = x101*x24;
auto x103 = x25*x48;
auto x104 = x102 + x103;
auto x105 = x104*x104;
auto x106 = x105*x35;
auto x107 = V{0.333333333333333}*cell[17];
auto x108 = V{0.333333333333333}*cell[18];
auto x109 = V{0.333333333333333}*cell[8];
auto x110 = V{0.333333333333333}*cell[9];
auto x111 = V{0.666666666666667}*cell[10];
auto x112 = V{0.666666666666667}*cell[1];
auto x113 = x107 + x108 + x109 + x110 - x111 - x112 + x59 + x73 + x97;
auto x114 = x106 + x113;
auto x115 = x104*x35;
auto x116 = -V{1}*cell[15] + V{1}*cell[16] - V{1}*cell[6] + V{1}*cell[7];
auto x117 = x115*x50 + x116;
auto x118 = x115*x82;
auto x119 = -cell[13] + cell[14] + x77;
auto x120 = -V{1}*cell[13] + V{1}*cell[14] - V{1}*cell[4] + V{1}*cell[5];
auto x121 = x118 + x120;
auto x122 = x35*x50*x82;
auto x123 = -cell[17] + cell[18] + x41;
auto x124 = -V{1}*cell[17] + V{1}*cell[18] - V{1}*cell[8] + V{1}*cell[9];
auto x125 = x122 + x124;
auto x126 = V{1} - V{1} / (x29 + V{3.00000046417339}*util::sqrt(x30 + x36*x37*util::sqrt(x121*(x118 + x119) + x125*(x122 + x123) + V{0.5}*(x114*x114) + V{1}*(x117*x117) + V{0.5}*(x75*x75) + V{0.5}*(x99*x99))));
auto x127 = V{0.5}*cell[13];
auto x128 = V{0.5}*cell[14];
auto x129 = V{0.5}*cell[15];
auto x130 = V{0.5}*cell[16];
auto x131 = V{0.5}*cell[17];
auto x132 = V{0.5}*cell[18];
auto x133 = V{0.5}*cell[4];
auto x134 = V{0.5}*cell[5];
auto x135 = V{0.5}*cell[6];
auto x136 = V{0.5}*cell[7];
auto x137 = V{0.5}*cell[8];
auto x138 = V{0.5}*cell[9];
auto x139 = V{0.5}*cell[0];
auto x140 = V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + x127 + x128 + x129 + x130 + x131 + x132 + x133 + x134 + x135 + x136 + x137 + x138 + x139 + V{0.5};
auto x141 = V{1.5}*x51;
auto x142 = V{1.5}*x105;
auto x143 = V{1.5}*x83;
auto x144 = x142 + x143 + V{-1};
auto x145 = x141 + x144;
auto x146 = -x23*x38 + x44;
auto x147 = x146*x24;
auto x148 = x27*x47;
auto x149 = x147 - x148;
auto x150 = x149*x149;
auto x151 = x150*x35 + x74;
auto x152 = -x22*x38 + x78;
auto x153 = x152*x24;
auto x154 = x26*x47;
auto x155 = x153 - x154;
auto x156 = x155*x155;
auto x157 = x156*x35 + x98;
auto x158 = x100 - x21*x38;
auto x159 = x158*x24;
auto x160 = x25*x47;
auto x161 = x159 - x160;
auto x162 = x161*x161;
auto x163 = x113 + x162*x35;
auto x164 = x161*x35;
auto x165 = x149*x164;
auto x166 = x116 + x165;
auto x167 = x155*x164;
auto x168 = x120 + x167;
auto x169 = x149*x155*x35;
auto x170 = x124 + x169;
auto x171 = V{1} - V{1} / (x29 + V{3.00000046417339}*util::sqrt(x30 + x36*x37*util::sqrt(x168*(x119 + x167) + x170*(x123 + x169) + V{0.5}*(x151*x151) + V{0.5}*(x157*x157) + V{0.5}*(x163*x163) + V{1}*(x166*x166))));
auto x172 = V{0.5}*x158*x24 - V{0.5}*x25*x47;
auto x173 = V{1}*x121;
auto x174 = V{1}*x149;
auto x175 = V{0.166666666666667}*cell[10];
auto x176 = V{0.166666666666667}*cell[11];
auto x177 = V{0.166666666666667}*cell[12];
auto x178 = V{0.166666666666667}*cell[13];
auto x179 = V{0.166666666666667}*cell[14];
auto x180 = V{0.166666666666667}*cell[15];
auto x181 = V{0.166666666666667}*cell[16];
auto x182 = V{0.166666666666667}*cell[17];
auto x183 = V{0.166666666666667}*cell[18];
auto x184 = V{0.166666666666667}*cell[1];
auto x185 = V{0.166666666666667}*cell[2];
auto x186 = V{0.166666666666667}*cell[3];
auto x187 = V{0.166666666666667}*cell[4];
auto x188 = V{0.166666666666667}*cell[5];
auto x189 = V{0.166666666666667}*cell[6];
auto x190 = V{0.166666666666667}*cell[7];
auto x191 = V{0.166666666666667}*cell[8];
auto x192 = V{0.166666666666667}*cell[9];
auto x193 = V{0.166666666666667}*cell[0] + x175 + x176 + x177 + x178 + x179 + x180 + x181 + x182 + x183 + x184 + x185 + x186 + x187 + x188 + x189 + x190 + x191 + x192 + V{0.166666666666667};
auto x194 = V{0.0833333333333333}*cell[0] + V{0.0833333333333333}*cell[10] + V{0.0833333333333333}*cell[11] + V{0.0833333333333333}*cell[12] + V{0.0833333333333333}*cell[13] + V{0.0833333333333333}*cell[14] + V{0.0833333333333333}*cell[15] + V{0.0833333333333333}*cell[16] + V{0.0833333333333333}*cell[17] + V{0.0833333333333333}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0833333333333333}*cell[2] + V{0.0833333333333333}*cell[3] + V{0.0833333333333333}*cell[4] + V{0.0833333333333333}*cell[5] + V{0.0833333333333333}*cell[6] + V{0.0833333333333333}*cell[7] + V{0.0833333333333333}*cell[8] + V{0.0833333333333333}*cell[9] + V{0.0833333333333333};
auto x195 = V{6.93889390390723e-18}*cell[0];
auto x196 = V{0.0833333333333333}*cell[12];
auto x197 = V{0.0833333333333333}*cell[3];
auto x198 = V{0.0833333333333333}*cell[13];
auto x199 = V{0.0833333333333333}*cell[14];
auto x200 = V{0.0833333333333333}*cell[4];
auto x201 = V{0.0833333333333333}*cell[5];
auto x202 = x195 + x196 + x197 - x198 - x199 - x200 - x201;
auto x203 = -x150*x194 + x202;
auto x204 = V{0.0833333333333333}*cell[11];
auto x205 = V{0.0833333333333333}*cell[2];
auto x206 = V{0.0833333333333333}*cell[15];
auto x207 = V{0.0833333333333333}*cell[16];
auto x208 = V{0.0833333333333333}*cell[6];
auto x209 = V{0.0833333333333333}*cell[7];
auto x210 = x204 + x205 - x206 - x207 - x208 - x209;
auto x211 = -x156*x194 + x210;
auto x212 = -x175 + x182 + x183 - x184 + x191 + x192;
auto x213 = x34 + V{1};
auto x214 = V{1.5}*x156;
auto x215 = V{6.000003}*x158*x24 - V{6.000003}*x25*x47;
auto x216 = V{2.999997}*x152*x24 - V{2.999997}*x26*x47;
auto x217 = V{1.5}*x150;
auto x218 = V{1} - x217;
auto x219 = -V{3}*x159 + V{3}*x160;
auto x220 = x156*x215 + x162*x216 + x219;
auto x221 = V{0.5}*x152*x24 - V{0.5}*x26*x47;
auto x222 = V{0.0833333333333333}*cell[10];
auto x223 = V{0.0833333333333333}*cell[1];
auto x224 = V{0.0833333333333333}*cell[17];
auto x225 = V{0.0833333333333333}*cell[18];
auto x226 = V{0.0833333333333333}*cell[8];
auto x227 = V{0.0833333333333333}*cell[9];
auto x228 = x222 + x223 - x224 - x225 - x226 - x227;
auto x229 = -x162*x194 + x228;
auto x230 = -x176 + x180 + x181 - x185 + x189 + x190;
auto x231 = V{1.5}*x162;
auto x232 = V{6.000003}*x152*x24 - V{6.000003}*x26*x47;
auto x233 = V{2.999997}*x158*x24 - V{2.999997}*x25*x47;
auto x234 = -V{3}*x153 + V{3}*x154;
auto x235 = x156*x233 + x162*x232 + x234;
auto x236 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x237 = V{9}*x24*x45 + V{9}*x27*x48;
auto x238 = x105*x237;
auto x239 = x237*x83;
auto x240 = x144 - V{3}*x51;
auto x241 = V{3}*x46;
auto x242 = V{3}*x49;
auto x243 = x241 + x242;
auto x244 = V{0.5}*x146*x24 - V{0.5}*x27*x47;
auto x245 = x117*x161;
auto x246 = x125*x155;
auto x247 = -x177 + x178 + x179 - x186 + x187 + x188 + x195;
auto x248 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x249 = x104 + x82;
auto x250 = V{3}*x102;
auto x251 = V{3}*x103;
auto x252 = x145 + x250 + x251;
auto x253 = V{3}*x80;
auto x254 = V{3}*x81;
auto x255 = x253 + x254;
auto x256 = -x196;
auto x257 = -x197;
auto x258 = V{0.25000025}*x158*x24 - V{0.25000025}*x25*x47;
auto x259 = x258*x99;
auto x260 = V{0.5000005}*x121;
auto x261 = x155*x260;
auto x262 = V{0.0416666666666667}*cell[0] + V{0.0416666666666667}*cell[10] + V{0.0416666666666667}*cell[11] + V{0.0416666666666667}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0416666666666667}*cell[1] + V{0.0416666666666667}*cell[2] + V{0.0416666666666667}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + V{0.0416666666666667};
auto x263 = x150*x262;
auto x264 = V{2.49999999985601e-07}*x158*x24 - V{2.49999999985601e-07}*x25*x47;
auto x265 = x264*x75;
auto x266 = V{4.99999999971202e-07}*x149;
auto x267 = x117*x266;
auto x268 = -V{0.0416666666666667}*cell[0];
auto x269 = V{0.0833333333333333}*cell[0] + x196 + x197 + x198 + x199 + x200 + x201 + x204 + x205 + x206 + x207 + x208 + x209 + x222 + x223 + x224 + x225 + x226 + x227 + V{0.0833333333333333};
auto x270 = V{0.0416666666666667}*cell[10] + V{6.93889390390723e-18}*cell[17] + V{6.93889390390723e-18}*cell[18] + V{0.0416666666666667}*cell[1] + V{6.93889390390723e-18}*cell[8] + V{6.93889390390723e-18}*cell[9] - x162*x269 + x268;
auto x271 = V{0.0416666666666667}*cell[11] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[16] + V{0.0416666666666667}*cell[2] + V{6.93889390390723e-18}*cell[6] + V{6.93889390390723e-18}*cell[7] - x156*x269;
auto x272 = x256 + x257 + x259 + x261 + x263 - x265 - x267 + x270 + x271;
auto x273 = V{0.25000025}*x152*x24 - V{0.25000025}*x26*x47;
auto x274 = x114*x273;
auto x275 = x161*x260;
auto x276 = V{2.49999999985601e-07}*x152*x24 - V{2.49999999985601e-07}*x26*x47;
auto x277 = x276*x75;
auto x278 = x125*x266;
auto x279 = x274 + x275 - x277 - x278;
auto x280 = V{0.25}*x167;
auto x281 = V{0.375}*cell[13] - V{0.125}*cell[14] + V{0.375}*cell[4] - V{0.125}*cell[5] - x280;
auto x282 = x104 - x80 - x81;
auto x283 = -x282;
auto x284 = x51*(V{3.000006}*x24*x79 + V{3.000006}*x26*x48);
auto x285 = x83*(V{6.000012}*x101*x24 + V{6.000012}*x25*x48);
auto x286 = x51*(V{3.000006}*x101*x24 + V{3.000006}*x25*x48);
auto x287 = x105*(V{6.000012}*x24*x79 + V{6.000012}*x26*x48);
auto x288 = -x253 - x254;
auto x289 = -V{0.125}*cell[13] + V{0.375}*cell[14] - V{0.125}*cell[4] + V{0.375}*cell[5] + x280;
auto x290 = x104 + x50;
auto x291 = V{9.000009}*x24*x45 + V{9.000009}*x27*x48;
auto x292 = x105*x291;
auto x293 = V{8.99999999948164e-06}*x24*x45 + V{8.99999999948164e-06}*x27*x48;
auto x294 = x293*x83;
auto x295 = x243 + x292 - x294;
auto x296 = x51*(V{5.999994}*x24*x79 + V{5.999994}*x26*x48);
auto x297 = x51*(V{12.000006}*x101*x24 + V{12.000006}*x25*x48);
auto x298 = V{6.000003}*x101*x24 + V{6.000003}*x25*x48;
auto x299 = x298*x83;
auto x300 = V{2.999997}*x24*x79 + V{2.999997}*x26*x48;
auto x301 = x105*x300;
auto x302 = x252 + x296 + x297 - x299 - x301;
auto x303 = -x204;
auto x304 = -x205;
auto x305 = x258*x75;
auto x306 = V{0.5000005}*x149;
auto x307 = x117*x306;
auto x308 = x156*x262;
auto x309 = x264*x99;
auto x310 = V{4.99999999971202e-07}*x121;
auto x311 = x155*x310;
auto x312 = V{0.0416666666666667}*cell[12] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[14] + V{0.0416666666666667}*cell[3] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[5] - x150*x269;
auto x313 = x270 + x303 + x304 + x305 + x307 + x308 - x309 - x311 + x312;
auto x314 = V{0.25000025}*x146*x24 - V{0.25000025}*x27*x47;
auto x315 = x114*x314;
auto x316 = V{0.5000005}*x245;
auto x317 = V{2.49999999985601e-07}*x146*x24 - V{2.49999999985601e-07}*x27*x47;
auto x318 = x317*x99;
auto x319 = V{4.99999999971202e-07}*x246;
auto x320 = x315 + x316 - x318 - x319;
auto x321 = V{0.25}*x165;
auto x322 = V{0.375}*cell[15] - V{0.125}*cell[16] + V{0.375}*cell[6] - V{0.125}*cell[7] - x321;
auto x323 = -x24*x45 - x27*x48;
auto x324 = x104 + x323;
auto x325 = -x324;
auto x326 = -x241 - x242;
auto x327 = -V{0.125}*cell[15] + V{0.375}*cell[16] - V{0.125}*cell[6] + V{0.375}*cell[7] + x321;
auto x328 = x50 + x82;
auto x329 = x291*x83;
auto x330 = x105*x293;
auto x331 = x243 + x329 - x330;
auto x332 = x51*(V{5.999994}*x101*x24 + V{5.999994}*x25*x48);
auto x333 = x51*(V{12.000006}*x24*x79 + V{12.000006}*x26*x48);
auto x334 = V{6.000003}*x24*x79 + V{6.000003}*x26*x48;
auto x335 = x105*x334;
auto x336 = V{2.999997}*x101*x24 + V{2.999997}*x25*x48;
auto x337 = x336*x83;
auto x338 = x145 + x255;
auto x339 = x332 + x333 - x335 - x337 + x338;
auto x340 = -x222;
auto x341 = -x223;
auto x342 = x273*x75;
auto x343 = x125*x306;
auto x344 = x162*x262;
auto x345 = x114*x276;
auto x346 = x161*x310;
auto x347 = x268 + x271 + x312 + x340 + x341 + x342 + x343 + x344 - x345 - x346;
auto x348 = x314*x99;
auto x349 = V{0.5000005}*x246;
auto x350 = x114*x317;
auto x351 = V{4.99999999971202e-07}*x245;
auto x352 = x348 + x349 - x350 - x351;
auto x353 = V{0.25}*x169;
auto x354 = V{0.375}*cell[17] - V{0.125}*cell[18] + V{0.375}*cell[8] - V{0.125}*cell[9] - x353;
auto x355 = x323 + x82;
auto x356 = -x355;
auto x357 = -V{0.125}*cell[17] + V{0.375}*cell[18] - V{0.125}*cell[8] + V{0.375}*cell[9] + x353;
auto x358 = V{0.5}*x101*x24 + V{0.5}*x25*x48;
auto x359 = -x59;
auto x360 = x359 - x60 - x61 + x62 + x63 + x64 + x65;
auto x361 = -x91 - x92 + x93 + x94 + x95 + x96;
auto x362 = x360 + x361 - x84 - x85 - x86 - x87 - x88 + x89 + x90;
auto x363 = -x67 - x68 + x69 + x70 + x71 + x72;
auto x364 = x360 + x363 - x52 - x53 - x54 - x55 - x56 + x57 + x58;
auto x365 = -V{1}*x121;
auto x366 = -x117;
auto x367 = V{1}*x50;
auto x368 = -x194*x51 + x202;
auto x369 = -x194*x83 + x210;
auto x370 = x141 + V{-1};
auto x371 = -x250 - x251;
auto x372 = x299 + x301 + x371;
auto x373 = V{0.5}*x24*x79 + V{0.5}*x26*x48;
auto x374 = -x106 - x107 - x108 - x109 - x110 + x111 + x112 + x359 + x361 + x363;
auto x375 = -x125;
auto x376 = -x105*x194 + x228;
auto x377 = x288 + x335 + x337;
auto x378 = V{0.5}*x24*x45 + V{0.5}*x27*x48;
auto x379 = x155 + x161;
auto x380 = x214 + x217 + x231 + V{-1};
auto x381 = V{0.5000005}*x168;
auto x382 = x256 + x257 + x263 + x270 + x271;
auto x383 = x149 + x161;
auto x384 = V{8.99999999948164e-06}*x146*x24 - V{8.99999999948164e-06}*x27*x47;
auto x385 = V{9.000009}*x146*x24 - V{9.000009}*x27*x47;
auto x386 = -V{3}*x147 + V{3}*x148 + x380;
auto x387 = V{4.99999999971202e-07}*x155;
auto x388 = x161*x166;
auto x389 = x270 + x303 + x304 + x308 + x312;
auto x390 = x149 + x155;
auto x391 = x268 + x271 + x312 + x340 + x341 + x344;
auto x0 = -V{1}*x126*(V{4.16333634234434e-17}*cell[10] + V{4.16333634234434e-17}*cell[11] + V{4.16333634234434e-17}*cell[12] + V{4.16333634234434e-17}*cell[1] + V{4.16333634234434e-17}*cell[2] + V{4.16333634234434e-17}*cell[3] - x105*x140 + x127 + x128 + x129 + x130 + x131 + x132 + x133 + x134 + x135 + x136 + x137 + x138 - x139 - x140*x51 - x140*x83) - x145*(x107 + x108 + x109 + x110 + x53 + x54 + x55 + x56 + x59 + x60 + x61 + x67 + x68 + x85 + x86 + x87 + x88 + x91 + x92 + V{0.333333333333333}) + V{-0.333333333333333};
auto x1 = -V{1}*x171*(x117*x174 + x155*x173 + x162*x193 + x172*x75 + x172*x99 + x203 + x211 + x212) + V{0.0555555555555556}*x213*(x150*x215 + x150*x216 + V{3}*x162 - x214 + x218 + x220) + V{-0.0555555555555556};
auto x2 = -V{1}*x171*(x114*x221 + x125*x174 + x156*x193 + x161*x173 + x203 + x221*x75 + x229 + x230) + V{0.0555555555555556}*x213*(x150*x232 + x150*x233 + V{3}*x156 + x218 - x231 + x235) + V{-0.0555555555555556};
auto x3 = V{1}*x126*(-x114*x244 - x150*x193 - x211 - x229 - x244*x99 - V{1}*x245 - V{1}*x246 - x247) - x236*(-x238 - x239 + x240 + x243) + V{-0.0555555555555556};
auto x4 = -(-V{1}*x126*(x272 + x279 + x281) + x248*(x105*(V{18}*x24*x79 + V{18}*x26*x48) + x252 + x255 - x51*(V{9}*x101*x24 + V{9}*x25*x48) - x51*(V{9}*x24*x79 + V{9}*x26*x48) + x83*(V{18}*x101*x24 + V{18}*x25*x48) - V{4.5}*x249*x249) + V{0.0277777777777778});
auto x5 = -(-V{1}*x126*(x272 - x274 - x275 + x277 + x278 + x289) + x248*(x252 + x284 + x285 - x286 - x287 + x288 - V{4.5}*x283*x283) + V{0.0277777777777778});
auto x6 = -(-V{1}*x126*(x313 + x320 + x322) + x248*(x295 + x302 - V{4.5}*x290*x290) + V{0.0277777777777778});
auto x7 = -(-V{1}*x126*(x313 - x315 - x316 + x318 + x319 + x327) + x248*(-x292 + x294 + x302 + x326 - V{4.5}*x325*x325) + V{0.0277777777777778});
auto x8 = -(-V{1}*x126*(x347 + x352 + x354) + x248*(x331 + x339 - V{4.5}*x328*x328) + V{0.0277777777777778});
auto x9 = -(-V{1}*x126*(x347 - x348 - x349 + x350 + x351 + x357) + x248*(x326 - x329 + x330 + x339 - V{4.5}*x356*x356) + V{0.0277777777777778});
auto x10 = -V{1}*x126*(x105*x193 + x212 + x358*x362 + x358*x364 + x365*x82 + x366*x367 + x368 + x369) - x236*(-V{3}*x105 + x143 + x298*x51 + x300*x51 + x370 + x372) + V{-0.0555555555555556};
auto x11 = -V{1}*x126*(x104*x365 + x193*x83 + x230 + x364*x373 + x367*x375 + x368 + x373*x374 + x376) - x236*(x142 + x334*x51 + x336*x51 + x370 + x377 - V{3}*x83) + V{-0.0555555555555556};
auto x12 = -V{1}*x126*(V{1}*x104*x366 + x193*x51 + x247 + x362*x378 + x369 + x374*x378 + V{1}*x375*x82 + x376) - x236*(x238 + x239 + x240 + x326) + V{-0.0555555555555556};
auto x13 = -(-V{1}*x171*(x151*x264 + x151*x276 - x155*x381 - x157*x258 - x161*x381 - x163*x273 + x166*x266 + x170*x266 + x281 + x382) + x248*(x150*(V{9}*x152*x24 - V{9}*x26*x47) + x150*(V{9}*x158*x24 - V{9}*x25*x47) - x156*(V{18}*x158*x24 - V{18}*x25*x47) - x162*(V{18}*x152*x24 - V{18}*x26*x47) + x219 + x234 + x380 - V{4.5}*x379*x379) + V{0.0277777777777778});
auto x14 = -(-V{1}*x126*(-x259 - x261 + x265 + x267 + x279 + x289 + x382) + x248*(-x284 - x285 + x286 + x287 + x338 + x371 - V{4.5}*x282*x282) + V{0.0277777777777778});
auto x15 = -(-V{1}*x171*(-x151*x258 + x157*x264 + x157*x317 - x163*x314 - x166*x306 + x168*x387 + x170*x387 + x322 - V{0.5000005}*x388 + x389) + x248*(-x150*(V{5.999994}*x152*x24 - V{5.999994}*x26*x47) - x150*(V{12.000006}*x158*x24 - V{12.000006}*x25*x47) + x156*x384 - x162*x385 + x220 + x386 - V{4.5}*x383*x383) + V{0.0277777777777778});
auto x16 = -(-V{1}*x126*(-x305 - x307 + x309 + x311 + x320 + x327 + x389) + x248*(x145 + x295 - x296 - x297 + x372 - V{4.5}*x324*x324) + V{0.0277777777777778});
auto x17 = -(-V{1}*x171*(-x151*x273 - V{0.5000005}*x155*x170 - x157*x314 + V{4.99999999971202e-07}*x161*x168 + x163*x276 + x163*x317 - x170*x306 + x354 + V{4.99999999971202e-07}*x388 + x391) + x248*(-x150*(V{12.000006}*x152*x24 - V{12.000006}*x26*x47) - x150*(V{5.999994}*x158*x24 - V{5.999994}*x25*x47) - x156*x385 + x162*x384 + x235 + x386 - V{4.5}*x390*x390) + V{0.0277777777777778});
auto x18 = -(-V{1}*x126*(-x342 - x343 + x345 + x346 + x352 + x357 + x391) + x248*(x145 + x331 - x332 - x333 + x377 - V{4.5}*x355*x355) + V{0.0277777777777778});
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
return { x35, x150 + x156 + x162 };
}
};

}

}
