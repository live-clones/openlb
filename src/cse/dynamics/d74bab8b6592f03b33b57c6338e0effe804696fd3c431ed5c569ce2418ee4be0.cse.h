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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::DualPorousBGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = cell.template getFieldComponent<opti::F>(1) + cell.template getFieldComponent<opti::F>(4);
auto x21 = cell.template getFieldComponent<opti::F>(12) + cell.template getFieldComponent<opti::F>(15);
auto x22 = cell.template getFieldComponent<opti::F>(11) + cell.template getFieldComponent<opti::F>(18);
auto x23 = cell.template getFieldComponent<opti::F>(0) + cell.template getFieldComponent<opti::F>(10) + cell.template getFieldComponent<opti::F>(13) + cell.template getFieldComponent<opti::F>(14) + cell.template getFieldComponent<opti::F>(16) + cell.template getFieldComponent<opti::F>(17) + cell.template getFieldComponent<opti::F>(2) + cell.template getFieldComponent<opti::F>(3) + cell.template getFieldComponent<opti::F>(5) + cell.template getFieldComponent<opti::F>(6) + cell.template getFieldComponent<opti::F>(7) + cell.template getFieldComponent<opti::F>(8) + cell.template getFieldComponent<opti::F>(9) + x20 + x21 + x22;
auto x24 = x23 + V{1};
auto x25 = V{1} / (x24);
auto x26 = x25*(x23 + V{1});
auto x27 = V{0.0277777777777778}*cell[13];
auto x28 = cell.template getFieldComponent<descriptors::POROSITY>(0)*x25;
auto x29 = V{3}*x28;
auto x30 = cell.template getFieldComponent<opti::F>(17) - cell.template getFieldComponent<opti::F>(8);
auto x31 = -cell.template getFieldComponent<opti::F>(9);
auto x32 = -cell.template getFieldComponent<opti::F>(14) + cell.template getFieldComponent<opti::F>(5);
auto x33 = cell.template getFieldComponent<opti::F>(18) + x31 + x32;
auto x34 = -cell.template getFieldComponent<opti::F>(2);
auto x35 = cell.template getFieldComponent<opti::F>(11) + cell.template getFieldComponent<opti::F>(13) - cell.template getFieldComponent<opti::F>(4) + x34;
auto x36 = x30 + x33 + x35;
auto x37 = x28*x36;
auto x38 = x37 + V{1};
auto x39 = x36*x38;
auto x40 = -cell.template getFieldComponent<opti::F>(16) + cell.template getFieldComponent<opti::F>(7);
auto x41 = x32 + x40;
auto x42 = -cell.template getFieldComponent<opti::F>(10);
auto x43 = -cell.template getFieldComponent<opti::F>(15) + cell.template getFieldComponent<opti::F>(6);
auto x44 = x42 + x43;
auto x45 = -cell.template getFieldComponent<opti::F>(13) + x20;
auto x46 = x41 + x44 + x45;
auto x47 = -x46;
auto x48 = x28*x47;
auto x49 = x48 + V{1};
auto x50 = x47*x49;
auto x51 = x30 + x40;
auto x52 = -cell.template getFieldComponent<opti::F>(18) + cell.template getFieldComponent<opti::F>(9);
auto x53 = -cell.template getFieldComponent<opti::F>(3);
auto x54 = -cell.template getFieldComponent<opti::F>(6) + x21 + x53;
auto x55 = x51 + x52 + x54;
auto x56 = x55*x55;
auto x57 = x28*x56;
auto x58 = x50 + x57;
auto x59 = V{3}*cell.template getFieldComponent<opti::F>(18);
auto x60 = V{3}*cell.template getFieldComponent<opti::F>(4);
auto x61 = V{3}*cell.template getFieldComponent<opti::F>(9);
auto x62 = -V{3}*cell.template getFieldComponent<opti::F>(14) + V{3}*cell.template getFieldComponent<opti::F>(5);
auto x63 = V{3}*cell.template getFieldComponent<opti::F>(17) - V{3}*cell.template getFieldComponent<opti::F>(8);
auto x64 = x25*(V{3}*cell.template getFieldComponent<opti::F>(11) + V{3}*cell.template getFieldComponent<opti::F>(13) - V{3}*cell.template getFieldComponent<opti::F>(2) + x59 - x60 - x61 + x62 + x63);
auto x65 = V{1} / ((x24)*(x24));
auto x66 = V{4.5}*x65;
auto x67 = cell.template getFieldComponent<opti::F>(1) + x44;
auto x68 = -cell.template getFieldComponent<opti::F>(17) + cell.template getFieldComponent<opti::F>(8);
auto x69 = cell.template getFieldComponent<opti::F>(11) + V{2}*cell.template getFieldComponent<opti::F>(13) - cell.template getFieldComponent<opti::F>(2) - V{2}*cell.template getFieldComponent<opti::F>(4) - x40 - x52 - x67 - x68;
auto x70 = -x66*x69*x69;
auto x71 = V{3}*cell.template getFieldComponent<opti::F>(6);
auto x72 = -V{3}*cell.template getFieldComponent<opti::F>(16) + V{3}*cell.template getFieldComponent<opti::F>(7);
auto x73 = V{3}*cell.template getFieldComponent<opti::F>(1) - V{3}*cell.template getFieldComponent<opti::F>(10) - V{3}*cell.template getFieldComponent<opti::F>(13) - V{3}*cell.template getFieldComponent<opti::F>(15) + x60 + x62 + x71 + x72;
auto x74 = -x25*x73;
auto x75 = V{1.5}*x65;
auto x76 = x47*x47;
auto x77 = x75*x76;
auto x78 = x36*x36;
auto x79 = x75*x78;
auto x80 = x56*x75;
auto x81 = x79 + x80 + V{-1};
auto x82 = x77 + x81;
auto x83 = x74 + x82;
auto x84 = x64 + x70 + x83;
auto x85 = V{0.0277777777777778}*cell[14];
auto x86 = x37 + V{-1};
auto x87 = -x86;
auto x88 = -x36*x87;
auto x89 = -x64;
auto x90 = -V{2}*cell.template getFieldComponent<opti::F>(14) + V{2}*cell.template getFieldComponent<opti::F>(5) + x22 + x31 + x34 + x51 + x67;
auto x91 = -x66*x90*x90 + x83 + x89;
auto x92 = V{0.0277777777777778}*cell[15];
auto x93 = x28*x55;
auto x94 = x93 + V{1};
auto x95 = x55*x94;
auto x96 = x28*x78;
auto x97 = x50 + x96;
auto x98 = x25*(V{3}*cell.template getFieldComponent<opti::F>(12) + V{3}*cell.template getFieldComponent<opti::F>(15) - V{3}*cell.template getFieldComponent<opti::F>(3) - x59 + x61 + x63 - x71 + x72);
auto x99 = x42 + x45;
auto x100 = -cell.template getFieldComponent<opti::F>(12) + cell.template getFieldComponent<opti::F>(3);
auto x101 = V{2}*cell.template getFieldComponent<opti::F>(15) - V{2}*cell.template getFieldComponent<opti::F>(6) - x100 - x33 - x68 - x99;
auto x102 = -x66*x101*x101;
auto x103 = x102 + x83 + x98;
auto x104 = V{0.0277777777777778}*cell[16];
auto x105 = x93 + V{-1};
auto x106 = -x105;
auto x107 = -x106*x55;
auto x108 = -x98;
auto x109 = cell.template getFieldComponent<opti::F>(12) - V{2}*cell.template getFieldComponent<opti::F>(16) + V{2}*cell.template getFieldComponent<opti::F>(7) + x30 + x32 + x52 + x53 + x99;
auto x110 = x108 - x66*x109*x109 + x83;
auto x111 = V{0.0277777777777778}*cell[17];
auto x112 = x28*x76;
auto x113 = x112 + x39;
auto x114 = V{2}*cell.template getFieldComponent<opti::F>(17) - V{2}*cell.template getFieldComponent<opti::F>(8) + x35 + x41 + x54;
auto x115 = x66*(x114*x114);
auto x116 = x64 + x82;
auto x117 = -x115 + x116 + x98;
auto x118 = V{0.0277777777777778}*cell[18];
auto x119 = cell.template getFieldComponent<opti::F>(16) + V{2}*cell.template getFieldComponent<opti::F>(18) - cell.template getFieldComponent<opti::F>(7) - V{2}*cell.template getFieldComponent<opti::F>(9) + x100 + x32 + x35 + x43;
auto x120 = -x119;
auto x121 = x108 + x116 - x66*x120*x120;
auto x122 = V{0.0277777777777778}*cell[5];
auto x123 = x48 + V{-1};
auto x124 = -x123;
auto x125 = -x124*x47;
auto x126 = x39 + x57;
auto x127 = -x74;
auto x128 = -x90;
auto x129 = x116 + x127 - x66*x128*x128;
auto x130 = V{0.0277777777777778}*cell[7];
auto x131 = x95 + x96;
auto x132 = -x109;
auto x133 = x82 + x98;
auto x134 = x127 + x133 - x66*x132*x132;
auto x135 = V{0.0277777777777778}*cell[9];
auto x136 = x112 + x95;
auto x137 = x133 - x66*x119*x119 + x89;
auto x138 = V{0.0555555555555556}*cell[10];
auto x139 = x46*x49;
auto x140 = x57 + x96;
auto x141 = V{3}*x65;
auto x142 = -x141*x76 + x74 + x81;
auto x143 = V{0.0555555555555556}*cell[11];
auto x144 = x46*x46;
auto x145 = x144*x28;
auto x146 = x141*x78;
auto x147 = x77 + V{-1};
auto x148 = -x146 + x147 + x64 + x80;
auto x149 = V{0.0555555555555556}*cell[12];
auto x150 = x141*x56;
auto x151 = x147 - x150 + x79 + x98;
auto x152 = V{0.333333333333333}*cell[0];
auto x153 = cell.template getFieldComponent<descriptors::POROSITY>(0)*cell.template getFieldComponent<descriptors::POROSITY>(0);
auto x154 = V{0.0277777777777778}*cell[4];
auto x155 = x124*x47;
auto x156 = x36*x87;
auto x157 = -x57;
auto x158 = x156 + x157;
auto x159 = x144*x75;
auto x160 = x25*x73 + x81;
auto x161 = x159 + x160;
auto x162 = -x161 - x70 - x89;
auto x163 = V{0.0277777777777778}*cell[6];
auto x164 = x106*x55;
auto x165 = -x96;
auto x166 = x164 + x165;
auto x167 = -x102 - x108 - x161;
auto x168 = V{0.0277777777777778}*cell[8];
auto x169 = -x112;
auto x170 = x156 + x164;
auto x171 = V{1} - x159;
auto x172 = x171 + x64 - x80;
auto x173 = -x79 + x98;
auto x174 = x115 + x172 + x173;
auto x175 = V{0.0555555555555556}*cell[1];
auto x176 = x123*x46;
auto x177 = V{3}*x144*x65 - x160;
auto x178 = V{0.0555555555555556}*cell[2];
auto x179 = x36*x86;
auto x180 = x146 + x172;
auto x181 = V{0.0555555555555556}*cell[3];
auto x182 = x105*x55;
auto x183 = x150 + x171 + x173;
auto x184 = cell.template getFieldComponent<opti::DJDF>(0) + cell[0] - x19*(cell[0] + x26*(x103*x92*(x29*(x95 + x97) + V{1}) + x104*x110*(-x29*(-x107 - x97) + V{1}) + x111*x117*(x29*(x113 + x95) + V{1}) + x118*x121*(-x29*(-x107 - x113) + V{1}) + x122*x129*(-x29*(-x125 - x126) + V{1}) + x130*x134*(-x29*(-x125 - x131) + V{1}) + x135*x137*(-x29*(-x136 - x88) + V{1}) + x138*x142*(x29*(-x139 + x140) + V{1}) + x143*x148*(x29*(x126 + x145) + V{1}) + x149*x151*(x29*(x131 + x145) + V{1}) + x152*x82*(x141*x153*(x144 + x56 + x78) + V{1}) - x154*x162*(-x29*(x155 + x158) + V{1}) - x163*x167*(-x29*(x155 + x166) + V{1}) - x168*x174*(-x29*(x169 + x170) + V{1}) - x175*x177*(x29*(x140 - x176) + V{1}) - x178*x180*(x29*(x145 + x179 + x57) + V{1}) - x181*x183*(x29*(x145 + x182 + x96) + V{1}) + x27*x84*(x29*(x39 + x58) + V{1}) + x85*x91*(-x29*(-x58 - x88) + V{1})));
auto x185 = V{3}*cell.template getFieldComponent<descriptors::POROSITY>(0);
auto x186 = x25*x39;
auto x187 = x25*x47;
auto x188 = x187 + V{-1};
auto x189 = x188*x49;
auto x190 = cell.template getFieldComponent<descriptors::POROSITY>(0)*x65;
auto x191 = x190*x56;
auto x192 = x189 + x191;
auto x193 = x179*x25;
auto x194 = x25*x95;
auto x195 = x190*x78;
auto x196 = x189 + x195;
auto x197 = x182*x25;
auto x198 = -x188;
auto x199 = -cell.template getFieldComponent<descriptors::POROSITY>(0)*x198*x47;
auto x200 = x199 + x95;
auto x201 = -x39;
auto x202 = x198*x47;
auto x203 = cell.template getFieldComponent<descriptors::POROSITY>(0)*x202;
auto x204 = x164 + x203;
auto x205 = x123*x188;
auto x206 = x191 + x205;
auto x207 = x195 + x205;
auto x208 = -x95;
auto x209 = x203 + x208;
auto x210 = V{3}*x25;
auto x211 = x188*x46;
auto x212 = x25*x78;
auto x213 = x25*x56;
auto x214 = x212 + x213;
auto x215 = cell.template getFieldComponent<opti::DJDF>(10) + cell[1] - x19*(cell[1] + x26*(x103*x92*(x185*(x194 + x196) + V{1}) + x104*x110*(x185*(x196 + x197) + V{1}) + x111*x117*(-x29*(-x200 - x39) + V{1}) + x118*x121*(-x29*(x201 + x204) + V{1}) + x122*x129*(x185*(x186 + x206) + V{1}) + x130*x134*(x185*(x194 + x207) + V{1}) + x135*x137*(-x29*(x156 + x209) + V{1}) + x138*x142*(x185*(x192 + x195) + V{1}) + x143*x148*(-x29*(-x126 - x199) + V{1}) + x149*x151*(-x29*(-x131 - x199) + V{1}) + x152*x82*(x153*x210*(-x211 + x214) + V{1}) - x154*x162*(x185*(x193 + x206) + V{1}) - x163*x167*(x185*(x197 + x207) + V{1}) - x168*x174*(-x29*(x170 + x203) + V{1}) - x175*x177*(x185*(x195 + x206) + V{1}) - x178*x180*(-x29*(x158 + x203) + V{1}) - x181*x183*(-x29*(x166 + x203) + V{1}) + x27*x84*(x185*(x186 + x192) + V{1}) + x85*x91*(x185*(x192 + x193) + V{1})));
auto x216 = x25*x36;
auto x217 = x216 + V{-1};
auto x218 = x217*x38;
auto x219 = -x139*x25;
auto x220 = x191 + x219;
auto x221 = x217*x86;
auto x222 = -x217;
auto x223 = -cell.template getFieldComponent<descriptors::POROSITY>(0)*x222*x36;
auto x224 = x223 + x50;
auto x225 = -x50;
auto x226 = x222*x36;
auto x227 = cell.template getFieldComponent<descriptors::POROSITY>(0)*x226;
auto x228 = x164 + x227;
auto x229 = x144*x190;
auto x230 = x218 + x229;
auto x231 = -x176*x25;
auto x232 = x191 + x231;
auto x233 = x155 + x227;
auto x234 = x221 + x229;
auto x235 = x217*x36;
auto x236 = x144*x25;
auto x237 = x213 + x236;
auto x238 = cell.template getFieldComponent<opti::DJDF>(11) + cell[2] - x19*(cell[2] + x26*(x103*x92*(-x29*(-x224 - x95) + V{1}) + x104*x110*(-x29*(x225 + x228) + V{1}) + x111*x117*(x185*(x194 + x230) + V{1}) + x118*x121*(x185*(x197 + x230) + V{1}) + x122*x129*(x185*(x218 + x232) + V{1}) + x130*x134*(-x29*(x208 + x233) + V{1}) + x135*x137*(x185*(x194 + x234) + V{1}) + x138*x142*(-x29*(-x223 - x58) + V{1}) + x143*x148*(x185*(x191 + x230) + V{1}) + x149*x151*(-x29*(-x136 - x223) + V{1}) + x152*x82*(x153*x210*(x235 + x237) + V{1}) - x154*x162*(x185*(x221 + x232) + V{1}) - x163*x167*(-x29*(x155 + x228) + V{1}) - x168*x174*(x185*(x197 + x234) + V{1}) - x175*x177*(-x29*(x157 + x233) + V{1}) - x178*x180*(x185*(x191 + x234) + V{1}) - x181*x183*(-x29*(x169 + x228) + V{1}) + x27*x84*(x185*(x218 + x220) + V{1}) + x85*x91*(x185*(x220 + x221) + V{1})));
auto x239 = x25*x55;
auto x240 = x239 + V{-1};
auto x241 = -x240;
auto x242 = -cell.template getFieldComponent<descriptors::POROSITY>(0)*x241*x55;
auto x243 = x39 + x50;
auto x244 = x241*x55;
auto x245 = cell.template getFieldComponent<descriptors::POROSITY>(0)*x244;
auto x246 = x156 + x245;
auto x247 = x240*x94;
auto x248 = x195 + x219;
auto x249 = x105*x240;
auto x250 = x186 + x229;
auto x251 = x155 + x245;
auto x252 = x195 + x247;
auto x253 = x193 + x229;
auto x254 = x240*x55;
auto x255 = x212 + x236;
auto x256 = x195 + x249;
auto x257 = cell.template getFieldComponent<opti::DJDF>(12) + cell[3] - x19*(cell[3] + x26*(x103*x92*(x185*(x247 + x248) + V{1}) + x104*x110*(x185*(x248 + x249) + V{1}) + x111*x117*(x185*(x247 + x250) + V{1}) + x118*x121*(x185*(x249 + x250) + V{1}) + x122*x129*(-x29*(x201 + x251) + V{1}) + x130*x134*(x185*(x231 + x252) + V{1}) + x135*x137*(x185*(x247 + x253) + V{1}) + x138*x142*(-x29*(-x242 - x97) + V{1}) + x143*x148*(-x29*(-x113 - x242) + V{1}) + x149*x151*(x185*(x229 + x252) + V{1}) + x152*x82*(x153*x210*(x254 + x255) + V{1}) - x154*x162*(-x29*(x155 + x246) + V{1}) - x163*x167*(x185*(x231 + x256) + V{1}) - x168*x174*(x185*(x249 + x253) + V{1}) - x175*x177*(-x29*(x165 + x251) + V{1}) - x178*x180*(-x29*(x169 + x246) + V{1}) - x181*x183*(x185*(x229 + x256) + V{1}) + x27*x84*(-x29*(-x242 - x243) + V{1}) + x85*x91*(-x29*(x225 + x246) + V{1})));
auto x258 = -x194;
auto x259 = x198*x49;
auto x260 = x226*x28;
auto x261 = x259 + x260;
auto x262 = x164*x25;
auto x263 = x222*x38;
auto x264 = x202*x28;
auto x265 = x263 + x264;
auto x266 = -cell.template getFieldComponent<descriptors::POROSITY>(0)*x222*x25*x36;
auto x267 = x124*x198;
auto x268 = x194 + x267;
auto x269 = x222*x87;
auto x270 = -cell.template getFieldComponent<descriptors::POROSITY>(0)*x198*x25*x47;
auto x271 = x194 + x270;
auto x272 = x235*x28;
auto x273 = -x211*x28;
auto x274 = x191 + x273;
auto x275 = -x267;
auto x276 = x260 + x262;
auto x277 = -x269;
auto x278 = x262 + x264;
auto x279 = cell.template getFieldComponent<opti::DJDF>(13) + cell[4] - x19*(cell[4] + x26*(x103*x92*(-x185*(x258 + x261) + V{1}) + x104*x110*(-x185*(x261 + x262) + V{1}) + x111*x117*(-x185*(x258 + x265) + V{1}) + x118*x121*(-x185*(x262 + x265) + V{1}) + x122*x129*(x185*(x206 + x218) + V{1}) + x130*x134*(-x185*(-x266 - x268) + V{1}) + x135*x137*(-x185*(-x269 - x271) + V{1}) + x138*x142*(x185*(x192 + x272) + V{1}) + x143*x148*(x185*(x218 + x274) + V{1}) + x149*x151*(-x29*(x209 + x227) + V{1}) + x152*x82*(-x153*x210*(x202 - x213 + x226) + V{1}) - x154*x162*(x185*(x206 + x221) + V{1}) - x163*x167*(-x185*(x275 + x276) + V{1}) - x168*x174*(-x185*(x277 + x278) + V{1}) - x175*x177*(x185*(x206 + x272) + V{1}) - x178*x180*(x185*(x221 + x274) + V{1}) - x181*x183*(-x29*(x204 + x227) + V{1}) + x27*x84*(x185*(x192 + x218) + V{1}) + x85*x91*(x185*(x192 + x221) + V{1})));
auto x280 = x216 + V{1};
auto x281 = x280*x38;
auto x282 = x280*x86;
auto x283 = -x198*x49;
auto x284 = x280*x36;
auto x285 = x28*x284;
auto x286 = x194 + x285;
auto x287 = -x285;
auto x288 = x262 + x287;
auto x289 = -x281;
auto x290 = x280*x87;
auto x291 = cell.template getFieldComponent<descriptors::POROSITY>(0)*x284;
auto x292 = -x198*x47;
auto x293 = x213 + x284;
auto x294 = -x106*x25*x55;
auto x295 = x285 + x294;
auto x296 = -x291;
auto x297 = cell.template getFieldComponent<opti::DJDF>(14) + cell[5] - x19*(cell[5] + x26*(x103*x92*(-x185*(-x283 - x286) + V{1}) + x104*x110*(-x185*(x259 + x288) + V{1}) + x111*x117*(-x185*(-x271 - x281) + V{1}) + x118*x121*(-x185*(x278 + x289) + V{1}) + x122*x129*(x185*(x206 + x281) + V{1}) + x130*x134*(x185*(x268 + x285) + V{1}) + x135*x137*(-x185*(x258 + x264 + x290) + V{1}) + x138*x142*(x185*(x192 + x285) + V{1}) + x143*x148*(x185*(x274 + x281) + V{1}) + x149*x151*(-x29*(-x200 - x291) + V{1}) + x152*x82*(-x153*x210*(-x292 - x293) + V{1}) - x154*x162*(x185*(x206 + x282) + V{1}) - x163*x167*(-x185*(-x267 - x295) + V{1}) - x168*x174*(-x185*(x278 + x290) + V{1}) - x175*x177*(x185*(x206 + x285) + V{1}) - x178*x180*(x185*(x274 + x282) + V{1}) - x181*x183*(-x29*(x204 + x296) + V{1}) + x27*x84*(x185*(x192 + x281) + V{1}) + x85*x91*(x185*(x192 + x282) + V{1})));
auto x298 = -x186;
auto x299 = x244*x28;
auto x300 = x259 + x299;
auto x301 = x156*x25;
auto x302 = x241*x94;
auto x303 = x264 + x302;
auto x304 = x106*x241;
auto x305 = x186 + x270;
auto x306 = -cell.template getFieldComponent<descriptors::POROSITY>(0)*x241*x25*x55;
auto x307 = x186 + x267;
auto x308 = x254*x28;
auto x309 = x299 + x301;
auto x310 = -x304;
auto x311 = x264 + x301;
auto x312 = cell.template getFieldComponent<opti::DJDF>(15) + cell[6] - x19*(cell[6] + x26*(x103*x92*(x185*(x196 + x247) + V{1}) + x104*x110*(x185*(x196 + x249) + V{1}) + x111*x117*(-x185*(x298 + x303) + V{1}) + x118*x121*(-x185*(-x304 - x305) + V{1}) + x122*x129*(-x185*(-x306 - x307) + V{1}) + x130*x134*(x185*(x207 + x247) + V{1}) + x135*x137*(-x185*(x301 + x303) + V{1}) + x138*x142*(x185*(x196 + x308) + V{1}) + x143*x148*(-x29*(x201 + x203 + x245) + V{1}) + x149*x151*(x185*(x252 + x273) + V{1}) + x152*x82*(-x153*x210*(x202 - x212 + x244) + V{1}) - x154*x162*(-x185*(x275 + x309) + V{1}) - x163*x167*(x185*(x207 + x249) + V{1}) - x168*x174*(-x185*(x310 + x311) + V{1}) - x175*x177*(x185*(x207 + x308) + V{1}) - x178*x180*(-x29*(x203 + x246) + V{1}) - x181*x183*(x185*(x256 + x273) + V{1}) + x27*x84*(-x185*(x298 + x300) + V{1}) + x85*x91*(-x185*(x300 + x301) + V{1})));
auto x313 = x239 + V{1};
auto x314 = x313*x55;
auto x315 = x28*x314;
auto x316 = x186 + x315;
auto x317 = -x315;
auto x318 = x301 + x317;
auto x319 = x313*x94;
auto x320 = x105*x313;
auto x321 = x106*x313;
auto x322 = -x319;
auto x323 = cell.template getFieldComponent<descriptors::POROSITY>(0)*x314;
auto x324 = x323 + x39;
auto x325 = x195 + x273;
auto x326 = x212 + x314;
auto x327 = -x25*x36*x87;
auto x328 = x315 + x327;
auto x329 = -x323;
auto x330 = x156 + x329;
auto x331 = cell.template getFieldComponent<opti::DJDF>(16) + cell[7] - x19*(cell[7] + x26*(x103*x92*(x185*(x196 + x319) + V{1}) + x104*x110*(x185*(x196 + x320) + V{1}) + x111*x117*(-x185*(-x305 - x319) + V{1}) + x118*x121*(-x185*(x264 + x298 + x321) + V{1}) + x122*x129*(x185*(x307 + x315) + V{1}) + x130*x134*(x185*(x207 + x319) + V{1}) + x135*x137*(-x185*(x311 + x322) + V{1}) + x138*x142*(x185*(x196 + x315) + V{1}) + x143*x148*(-x29*(-x199 - x324) + V{1}) + x149*x151*(x185*(x319 + x325) + V{1}) + x152*x82*(-x153*x210*(-x292 - x326) + V{1}) - x154*x162*(-x185*(-x267 - x328) + V{1}) - x163*x167*(x185*(x207 + x320) + V{1}) - x168*x174*(-x185*(x311 + x321) + V{1}) - x175*x177*(x185*(x207 + x315) + V{1}) - x178*x180*(-x29*(x203 + x330) + V{1}) - x181*x183*(x185*(x320 + x325) + V{1}) + x27*x84*(-x185*(-x283 - x316) + V{1}) + x85*x91*(-x185*(x259 + x318) + V{1})));
auto x332 = x25*x50;
auto x333 = -x332;
auto x334 = x263 + x299;
auto x335 = x269 + x332;
auto x336 = x260 + x302;
auto x337 = x266 + x332;
auto x338 = x155*x25;
auto x339 = x229 + x272;
auto x340 = x25*x76;
auto x341 = x299 + x338;
auto x342 = x260 + x338;
auto x343 = cell.template getFieldComponent<opti::DJDF>(17) + cell[8] - x19*(cell[8] + x26*(x103*x92*(-x185*(x333 + x336) + V{1}) + x104*x110*(-x185*(-x304 - x337) + V{1}) + x111*x117*(x185*(x230 + x247) + V{1}) + x118*x121*(x185*(x230 + x249) + V{1}) + x122*x129*(-x185*(x334 + x338) + V{1}) + x130*x134*(-x185*(x336 + x338) + V{1}) + x135*x137*(x185*(x234 + x247) + V{1}) + x138*x142*(-x29*(x225 + x227 + x245) + V{1}) + x143*x148*(x185*(x230 + x308) + V{1}) + x149*x151*(x185*(x247 + x339) + V{1}) + x152*x82*(-x153*x210*(x226 + x244 - x340) + V{1}) - x154*x162*(-x185*(x277 + x341) + V{1}) - x163*x167*(-x185*(x310 + x342) + V{1}) - x168*x174*(x185*(x234 + x249) + V{1}) - x175*x177*(-x29*(x233 + x245) + V{1}) - x178*x180*(x185*(x234 + x308) + V{1}) - x181*x183*(x185*(x249 + x339) + V{1}) + x27*x84*(-x185*(x333 + x334) + V{1}) + x85*x91*(-x185*(-x306 - x335) + V{1})));
auto x344 = -x222*x38;
auto x345 = x315 + x332;
auto x346 = x317 + x338;
auto x347 = -x222*x36;
auto x348 = x314 + x340;
auto x349 = -x124*x25*x47;
auto x350 = x315 + x349;
auto x351 = cell.template getFieldComponent<opti::DJDF>(18) + cell[9] - x19*(cell[9] + x26*(x103*x92*(-x185*(-x319 - x337) + V{1}) + x104*x110*(-x185*(x260 + x321 + x333) + V{1}) + x111*x117*(x185*(x230 + x319) + V{1}) + x118*x121*(x185*(x230 + x320) + V{1}) + x122*x129*(-x185*(x263 + x346) + V{1}) + x130*x134*(-x185*(x322 + x342) + V{1}) + x135*x137*(x185*(x234 + x319) + V{1}) + x138*x142*(-x29*(-x224 - x323) + V{1}) + x143*x148*(x185*(x230 + x315) + V{1}) + x149*x151*(x185*(x319 + x339) + V{1}) + x152*x82*(-x153*x210*(-x347 - x348) + V{1}) - x154*x162*(-x185*(-x269 - x350) + V{1}) - x163*x167*(-x185*(x321 + x342) + V{1}) - x168*x174*(x185*(x234 + x320) + V{1}) - x175*x177*(-x29*(x233 + x329) + V{1}) - x178*x180*(x185*(x234 + x315) + V{1}) - x181*x183*(x185*(x320 + x339) + V{1}) + x27*x84*(-x185*(-x344 - x345) + V{1}) + x85*x91*(x185*(x315 + x335) + V{1})));
auto x352 = x187 + V{1};
auto x353 = x352*x49;
auto x354 = x191 + x353;
auto x355 = x195 + x353;
auto x356 = x352*x47;
auto x357 = cell.template getFieldComponent<descriptors::POROSITY>(0)*x356;
auto x358 = x357 + x39;
auto x359 = x123*x352;
auto x360 = x191 + x359;
auto x361 = x195 + x359;
auto x362 = x357 + x88;
auto x363 = x352*x46;
auto x364 = -x357;
auto x365 = x107 + x357;
auto x366 = cell.template getFieldComponent<opti::DJDF>(1) + cell[10] - x19*(cell[10] + x26*(x103*x92*(x185*(x194 + x355) + V{1}) + x104*x110*(x185*(x197 + x355) + V{1}) + x111*x117*(x29*(x358 + x95) + V{1}) + x118*x121*(-x29*(-x107 - x358) + V{1}) + x122*x129*(x185*(x186 + x360) + V{1}) + x130*x134*(x185*(x194 + x361) + V{1}) + x135*x137*(-x29*(-x362 - x95) + V{1}) + x138*x142*(x185*(x195 + x354) + V{1}) + x143*x148*(x29*(x126 + x357) + V{1}) + x149*x151*(x29*(x131 + x357) + V{1}) + x152*x82*(x153*x210*(x214 - x363) + V{1}) - x154*x162*(x185*(x193 + x360) + V{1}) - x163*x167*(x185*(x197 + x361) + V{1}) - x168*x174*(-x29*(x170 + x364) + V{1}) - x175*x177*(x185*(x195 + x360) + V{1}) - x178*x180*(-x29*(-x362 - x57) + V{1}) - x181*x183*(-x29*(-x365 - x96) + V{1}) + x27*x84*(x185*(x186 + x354) + V{1}) + x85*x91*(x185*(x193 + x354) + V{1})));
auto x367 = x291 + x50;
auto x368 = x229 + x281;
auto x369 = x125 + x291;
auto x370 = x229 + x282;
auto x371 = cell.template getFieldComponent<opti::DJDF>(2) + cell[11] - x19*(cell[11] + x26*(x103*x92*(x29*(x367 + x95) + V{1}) + x104*x110*(-x29*(-x107 - x367) + V{1}) + x111*x117*(x185*(x194 + x368) + V{1}) + x118*x121*(x185*(x197 + x368) + V{1}) + x122*x129*(x185*(x232 + x281) + V{1}) + x130*x134*(-x29*(-x369 - x95) + V{1}) + x135*x137*(x185*(x194 + x370) + V{1}) + x138*x142*(x29*(x291 + x58) + V{1}) + x143*x148*(x185*(x191 + x368) + V{1}) + x149*x151*(x29*(x136 + x291) + V{1}) + x152*x82*(x153*x210*(x237 + x284) + V{1}) - x154*x162*(x185*(x232 + x282) + V{1}) - x163*x167*(-x29*(x155 + x164 + x296) + V{1}) - x168*x174*(x185*(x197 + x370) + V{1}) - x175*x177*(-x29*(-x369 - x57) + V{1}) - x178*x180*(x185*(x191 + x370) + V{1}) - x181*x183*(-x29*(-x107 - x112 - x291) + V{1}) + x27*x84*(x185*(x220 + x281) + V{1}) + x85*x91*(x185*(x220 + x282) + V{1})));
auto x372 = x323 + x88;
auto x373 = x195 + x319;
auto x374 = x195 + x320;
auto x375 = cell.template getFieldComponent<opti::DJDF>(3) + cell[12] - x19*(cell[12] + x26*(x103*x92*(x185*(x248 + x319) + V{1}) + x104*x110*(x185*(x248 + x320) + V{1}) + x111*x117*(x185*(x250 + x319) + V{1}) + x118*x121*(x185*(x250 + x320) + V{1}) + x122*x129*(-x29*(-x125 - x324) + V{1}) + x130*x134*(x185*(x231 + x373) + V{1}) + x135*x137*(x185*(x253 + x319) + V{1}) + x138*x142*(x29*(x323 + x97) + V{1}) + x143*x148*(x29*(x113 + x323) + V{1}) + x149*x151*(x185*(x229 + x373) + V{1}) + x152*x82*(x153*x210*(x255 + x314) + V{1}) - x154*x162*(-x29*(x155 + x330) + V{1}) - x163*x167*(x185*(x231 + x374) + V{1}) - x168*x174*(x185*(x253 + x320) + V{1}) - x175*x177*(-x29*(-x125 - x323 - x96) + V{1}) - x178*x180*(-x29*(-x112 - x372) + V{1}) - x181*x183*(x185*(x229 + x374) + V{1}) + x27*x84*(x29*(x243 + x323) + V{1}) + x85*x91*(-x29*(-x372 - x50) + V{1})));
auto x376 = x28*x356;
auto x377 = x281 + x376;
auto x378 = -x124*x352;
auto x379 = -x280*x87;
auto x380 = x194 + x376;
auto x381 = -x28*x363;
auto x382 = x191 + x381;
auto x383 = x357 + x95;
auto x384 = x124*x352;
auto x385 = -x376;
auto x386 = x262 + x385;
auto x387 = cell.template getFieldComponent<opti::DJDF>(4) + cell[13] - x19*(cell[13] + x26*(x103*x92*(x185*(x286 + x353) + V{1}) + x104*x110*(-x185*(-x295 - x353) + V{1}) + x111*x117*(x185*(x194 + x377) + V{1}) + x118*x121*(-x185*(-x294 - x377) + V{1}) + x122*x129*(x185*(x281 + x360) + V{1}) + x130*x134*(-x185*(-x286 - x378) + V{1}) + x135*x137*(-x185*(-x379 - x380) + V{1}) + x138*x142*(x185*(x285 + x354) + V{1}) + x143*x148*(x185*(x281 + x382) + V{1}) + x149*x151*(x29*(x291 + x383) + V{1}) + x152*x82*(x153*x210*(x293 + x356) + V{1}) - x154*x162*(x185*(x282 + x360) + V{1}) - x163*x167*(-x185*(x288 + x384) + V{1}) - x168*x174*(-x185*(x290 + x386) + V{1}) - x175*x177*(x185*(x285 + x360) + V{1}) - x178*x180*(x185*(x282 + x382) + V{1}) - x181*x183*(-x29*(-x291 - x365) + V{1}) + x27*x84*(x185*(x281 + x354) + V{1}) + x85*x91*(x185*(x282 + x354) + V{1})));
auto x388 = -x353;
auto x389 = cell.template getFieldComponent<opti::DJDF>(5) + cell[14] - x19*(cell[14] + x26*(x103*x92*(-x185*(-x194 - x266 - x353) + V{1}) + x104*x110*(-x185*(x276 + x388) + V{1}) + x111*x117*(-x185*(-x344 - x380) + V{1}) + x118*x121*(-x185*(x263 + x386) + V{1}) + x122*x129*(x185*(x218 + x360) + V{1}) + x130*x134*(-x185*(x258 + x260 + x384) + V{1}) + x135*x137*(x185*(x269 + x380) + V{1}) + x138*x142*(x185*(x272 + x354) + V{1}) + x143*x148*(x185*(x218 + x382) + V{1}) + x149*x151*(-x29*(-x223 - x383) + V{1}) + x152*x82*(-x153*x210*(-x213 - x347 - x356) + V{1}) - x154*x162*(x185*(x221 + x360) + V{1}) - x163*x167*(-x185*(x276 + x384) + V{1}) - x168*x174*(-x185*(-x269 - x294 - x376) + V{1}) - x175*x177*(x185*(x272 + x360) + V{1}) - x178*x180*(x185*(x221 + x382) + V{1}) - x181*x183*(-x29*(x228 + x364) + V{1}) + x27*x84*(x185*(x218 + x354) + V{1}) + x85*x91*(x185*(x221 + x354) + V{1})));
auto x390 = x186 + x376;
auto x391 = -x106*x313;
auto x392 = x327 + x376;
auto x393 = x301 + x385;
auto x394 = cell.template getFieldComponent<opti::DJDF>(6) + cell[15] - x19*(cell[15] + x26*(x103*x92*(x185*(x319 + x355) + V{1}) + x104*x110*(x185*(x320 + x355) + V{1}) + x111*x117*(x185*(x319 + x390) + V{1}) + x118*x121*(-x185*(-x390 - x391) + V{1}) + x122*x129*(-x185*(-x316 - x378) + V{1}) + x130*x134*(x185*(x319 + x361) + V{1}) + x135*x137*(-x185*(-x319 - x392) + V{1}) + x138*x142*(x185*(x315 + x355) + V{1}) + x143*x148*(x29*(x324 + x357) + V{1}) + x149*x151*(x185*(x373 + x381) + V{1}) + x152*x82*(x153*x210*(x326 + x356) + V{1}) - x154*x162*(-x185*(x318 + x384) + V{1}) - x163*x167*(x185*(x320 + x361) + V{1}) - x168*x174*(-x185*(x321 + x393) + V{1}) - x175*x177*(x185*(x315 + x361) + V{1}) - x178*x180*(-x29*(-x323 - x362) + V{1}) - x181*x183*(x185*(x374 + x381) + V{1}) + x27*x84*(x185*(x316 + x353) + V{1}) + x85*x91*(-x185*(-x328 - x353) + V{1})));
auto x395 = -x241*x94;
auto x396 = -x241*x55;
auto x397 = cell.template getFieldComponent<opti::DJDF>(7) + cell[16] - x19*(cell[16] + x26*(x103*x92*(x185*(x252 + x353) + V{1}) + x104*x110*(x185*(x256 + x353) + V{1}) + x111*x117*(-x185*(-x390 - x395) + V{1}) + x118*x121*(x185*(x304 + x390) + V{1}) + x122*x129*(-x185*(x298 + x299 + x384) + V{1}) + x130*x134*(x185*(x252 + x359) + V{1}) + x135*x137*(-x185*(x302 + x393) + V{1}) + x138*x142*(x185*(x308 + x355) + V{1}) + x143*x148*(-x29*(-x242 - x358) + V{1}) + x149*x151*(x185*(x252 + x381) + V{1}) + x152*x82*(-x153*x210*(-x212 - x356 - x396) + V{1}) - x154*x162*(-x185*(x309 + x384) + V{1}) - x163*x167*(x185*(x256 + x359) + V{1}) - x168*x174*(-x185*(-x304 - x392) + V{1}) - x175*x177*(x185*(x308 + x361) + V{1}) - x178*x180*(-x29*(x246 + x364) + V{1}) - x181*x183*(x185*(x256 + x381) + V{1}) + x27*x84*(-x185*(-x186 - x306 - x353) + V{1}) + x85*x91*(-x185*(x309 + x388) + V{1})));
auto x398 = x285 + x332;
auto x399 = x285 + x319;
auto x400 = x287 + x338;
auto x401 = x229 + x285;
auto x402 = cell.template getFieldComponent<opti::DJDF>(8) + cell[17] - x19*(cell[17] + x26*(x103*x92*(x185*(x319 + x398) + V{1}) + x104*x110*(-x185*(-x391 - x398) + V{1}) + x111*x117*(x185*(x319 + x368) + V{1}) + x118*x121*(x185*(x320 + x368) + V{1}) + x122*x129*(-x185*(-x281 - x350) + V{1}) + x130*x134*(-x185*(-x349 - x399) + V{1}) + x135*x137*(x185*(x319 + x370) + V{1}) + x138*x142*(x29*(x323 + x367) + V{1}) + x143*x148*(x185*(x315 + x368) + V{1}) + x149*x151*(x185*(x229 + x399) + V{1}) + x152*x82*(x153*x210*(x284 + x348) + V{1}) - x154*x162*(-x185*(x290 + x346) + V{1}) - x163*x167*(-x185*(x321 + x400) + V{1}) - x168*x174*(x185*(x320 + x370) + V{1}) - x175*x177*(-x29*(-x323 - x369) + V{1}) - x178*x180*(x185*(x315 + x370) + V{1}) - x181*x183*(x185*(x320 + x401) + V{1}) + x27*x84*(x185*(x281 + x345) + V{1}) + x85*x91*(-x185*(-x345 - x379) + V{1})));
auto x403 = cell.template getFieldComponent<opti::DJDF>(9) + cell[18] - x19*(cell[18] + x26*(x103*x92*(-x185*(-x395 - x398) + V{1}) + x104*x110*(x185*(x304 + x398) + V{1}) + x111*x117*(x185*(x247 + x368) + V{1}) + x118*x121*(x185*(x249 + x368) + V{1}) + x122*x129*(-x185*(x289 + x341) + V{1}) + x130*x134*(-x185*(x302 + x400) + V{1}) + x135*x137*(x185*(x247 + x370) + V{1}) + x138*x142*(-x29*(-x242 - x367) + V{1}) + x143*x148*(x185*(x308 + x368) + V{1}) + x149*x151*(x185*(x247 + x401) + V{1}) + x152*x82*(-x153*x210*(-x284 - x340 - x396) + V{1}) - x154*x162*(-x185*(x290 + x341) + V{1}) - x163*x167*(-x185*(-x285 - x304 - x349) + V{1}) - x168*x174*(x185*(x249 + x370) + V{1}) - x175*x177*(-x29*(x251 + x296) + V{1}) - x178*x180*(x185*(x308 + x370) + V{1}) - x181*x183*(x185*(x249 + x401) + V{1}) + x27*x84*(-x185*(-x281 - x306 - x332) + V{1}) + x85*x91*(-x185*(x290 + x299 + x333) + V{1})));
cell[0] = x184;
cell[1] = x215;
cell[2] = x238;
cell[3] = x257;
cell[4] = x279;
cell[5] = x297;
cell[6] = x312;
cell[7] = x331;
cell[8] = x343;
cell[9] = x351;
cell[10] = x366;
cell[11] = x371;
cell[12] = x375;
cell[13] = x387;
cell[14] = x389;
cell[15] = x394;
cell[16] = x397;
cell[17] = x402;
cell[18] = x403;
return { x184 + V{1}, x184*x184 + x215*x215 + x238*x238 + x257*x257 + x279*x279 + x297*x297 + x312*x312 + x331*x331 + x343*x343 + x351*x351 + x366*x366 + x371*x371 + x375*x375 + x387*x387 + x389*x389 + x394*x394 + x397*x397 + x402*x402 + x403*x403 };
}
};

}

}
