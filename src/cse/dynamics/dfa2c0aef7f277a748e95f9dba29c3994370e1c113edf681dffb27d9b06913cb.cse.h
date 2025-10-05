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
struct CSE<Dual<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, T, descriptors::D3Q19<FIELDS...> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = V{1}*x19 + V{-1};
auto x21 = V{0.0277777777777778}*x19;
auto x22 = cell.template getFieldComponent<opti::F>(17) - cell.template getFieldComponent<opti::F>(8);
auto x23 = -cell.template getFieldComponent<opti::F>(9);
auto x24 = -cell.template getFieldComponent<opti::F>(14) + cell.template getFieldComponent<opti::F>(5);
auto x25 = cell.template getFieldComponent<opti::F>(18) + x23 + x24;
auto x26 = -cell.template getFieldComponent<opti::F>(2);
auto x27 = cell.template getFieldComponent<opti::F>(11) + cell.template getFieldComponent<opti::F>(13);
auto x28 = -cell.template getFieldComponent<opti::F>(4) + x26 + x27;
auto x29 = x22 + x25 + x28;
auto x30 = cell.template getFieldComponent<opti::F>(12) + cell.template getFieldComponent<opti::F>(15);
auto x31 = cell.template getFieldComponent<opti::F>(0) + cell.template getFieldComponent<opti::F>(1) + cell.template getFieldComponent<opti::F>(10) + cell.template getFieldComponent<opti::F>(14) + cell.template getFieldComponent<opti::F>(16) + cell.template getFieldComponent<opti::F>(17) + cell.template getFieldComponent<opti::F>(18) + cell.template getFieldComponent<opti::F>(2) + cell.template getFieldComponent<opti::F>(3) + cell.template getFieldComponent<opti::F>(4) + cell.template getFieldComponent<opti::F>(5) + cell.template getFieldComponent<opti::F>(6) + cell.template getFieldComponent<opti::F>(7) + cell.template getFieldComponent<opti::F>(8) + cell.template getFieldComponent<opti::F>(9) + x27 + x30;
auto x32 = x31 + V{1};
auto x33 = cell.template getFieldComponent<descriptors::POROSITY>(0)/x32;
auto x34 = V{3}*x33;
auto x35 = x29*x34;
auto x36 = V{1} / ((x32)*(x32));
auto x37 = V{4.5}*x36;
auto x38 = cell.template getFieldComponent<descriptors::POROSITY>(0)*cell.template getFieldComponent<descriptors::POROSITY>(0);
auto x39 = -cell.template getFieldComponent<opti::F>(16) + cell.template getFieldComponent<opti::F>(7);
auto x40 = -cell.template getFieldComponent<opti::F>(18) + cell.template getFieldComponent<opti::F>(9);
auto x41 = -cell.template getFieldComponent<opti::F>(17) + cell.template getFieldComponent<opti::F>(8);
auto x42 = -cell.template getFieldComponent<opti::F>(10);
auto x43 = -cell.template getFieldComponent<opti::F>(15) + cell.template getFieldComponent<opti::F>(6);
auto x44 = cell.template getFieldComponent<opti::F>(1) + x42 + x43;
auto x45 = cell.template getFieldComponent<opti::F>(11) + V{2}*cell.template getFieldComponent<opti::F>(13) - cell.template getFieldComponent<opti::F>(2) - V{2}*cell.template getFieldComponent<opti::F>(4) - x39 - x40 - x41 - x44;
auto x46 = -x37*x38*x45*x45;
auto x47 = x24 + x39;
auto x48 = cell.template getFieldComponent<opti::F>(1) - cell.template getFieldComponent<opti::F>(13) + cell.template getFieldComponent<opti::F>(4) + x42;
auto x49 = x43 + x47 + x48;
auto x50 = -x49;
auto x51 = x34*x50;
auto x52 = V{1.5}*x36;
auto x53 = x38*(x50*x50);
auto x54 = x52*x53;
auto x55 = x38*(x29*x29);
auto x56 = x52*x55;
auto x57 = x22 + x39;
auto x58 = -cell.template getFieldComponent<opti::F>(3);
auto x59 = -cell.template getFieldComponent<opti::F>(6) + x30 + x58;
auto x60 = x40 + x57 + x59;
auto x61 = x38*(x60*x60);
auto x62 = x52*x61;
auto x63 = x56 + x62 + V{-1};
auto x64 = x54 + x63;
auto x65 = x51 + x64;
auto x66 = -x35;
auto x67 = cell.template getFieldComponent<opti::F>(11) - V{2}*cell.template getFieldComponent<opti::F>(14) + cell.template getFieldComponent<opti::F>(18) + V{2}*cell.template getFieldComponent<opti::F>(5) + x23 + x26 + x44 + x57;
auto x68 = x34*x60;
auto x69 = -cell.template getFieldComponent<opti::F>(12) + cell.template getFieldComponent<opti::F>(3);
auto x70 = V{2}*cell.template getFieldComponent<opti::F>(15) - V{2}*cell.template getFieldComponent<opti::F>(6) - x25 - x41 - x48 - x69;
auto x71 = -x37*x38*x70*x70;
auto x72 = -x68;
auto x73 = cell.template getFieldComponent<opti::F>(12) - V{2}*cell.template getFieldComponent<opti::F>(16) + V{2}*cell.template getFieldComponent<opti::F>(7) + x22 + x24 + x40 + x48 + x58;
auto x74 = V{2}*cell.template getFieldComponent<opti::F>(17) - V{2}*cell.template getFieldComponent<opti::F>(8) + x28 + x47 + x59;
auto x75 = x37*x38*(x74*x74);
auto x76 = x35 + x64;
auto x77 = cell.template getFieldComponent<opti::F>(16) + V{2}*cell.template getFieldComponent<opti::F>(18) - cell.template getFieldComponent<opti::F>(7) - V{2}*cell.template getFieldComponent<opti::F>(9) + x24 + x28 + x43 + x69;
auto x78 = -x77;
auto x79 = -x51;
auto x80 = -x67;
auto x81 = -x73;
auto x82 = x64 + x68;
auto x83 = V{0.0555555555555556}*x19;
auto x84 = V{3}*x36;
auto x85 = x55*x84;
auto x86 = x54 + V{-1};
auto x87 = x61*x84;
auto x88 = x49*x49;
auto x89 = x38*x52*x88;
auto x90 = x34*x49 + x63;
auto x91 = x89 + x90;
auto x92 = V{1} - x89;
auto x93 = x35 - x62 + x92;
auto x94 = -x56 + x68;
auto x95 = V{0.0833333333333333}*cell[14];
auto x96 = V{0.0833333333333333}*cell[9];
auto x97 = V{0.0833333333333333}*cell[18];
auto x98 = V{0.0833333333333333}*cell[5];
auto x99 = V{0.25}*x33;
auto x100 = x77*x99;
auto x101 = cell[18]*x100;
auto x102 = cell[5]*x67*x99;
auto x103 = cell[9]*x100;
auto x104 = x29*x33;
auto x105 = V{0.5}*x104;
auto x106 = cell[14]*x80*x99;
auto x107 = V{1}*cell[0] + V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[12] + V{0.0833333333333333}*cell[13] + V{0.0833333333333333}*cell[14] + V{0.0833333333333333}*cell[15] + V{0.0833333333333333}*cell[16] + V{0.0833333333333333}*cell[17] + V{0.0833333333333333}*cell[18] + V{0.166666666666667}*cell[1] + V{0.166666666666667}*cell[2] + V{0.166666666666667}*cell[3] + V{0.0833333333333333}*cell[4] + V{0.0833333333333333}*cell[5] + V{0.0833333333333333}*cell[6] + V{0.0833333333333333}*cell[7] + V{0.0833333333333333}*cell[8] + V{0.0833333333333333}*cell[9];
auto x108 = x45*x99;
auto x109 = cell[13]*x108 - V{0.0833333333333333}*cell[13] + cell[4]*x108 + V{0.0833333333333333}*cell[4];
auto x110 = x74*x99;
auto x111 = cell[17]*x110 - V{0.0833333333333333}*cell[17] + cell[8]*x110 + V{0.0833333333333333}*cell[8];
auto x112 = cell[11]*x105 - V{0.166666666666667}*cell[11] + cell[2]*x105 + V{0.166666666666667}*cell[2] + x101 + x102 + x103 - x104*x107 - x106 + x109 + x111 + x95 + x96 - x97 - x98;
auto x113 = x31 + V{1};
auto x114 = V{1}*x113*x19;
auto x115 = cell.template getFieldComponent<descriptors::POROSITY>(0)*x114*x36;
auto x116 = V{0.0833333333333333}*cell[16];
auto x117 = V{0.0833333333333333}*cell[7];
auto x118 = cell[7]*x73*x99;
auto x119 = x33*x60;
auto x120 = V{0.5}*x119;
auto x121 = cell[16]*x81*x99;
auto x122 = x70*x99;
auto x123 = cell[15]*x122 - V{0.0833333333333333}*cell[15] + cell[6]*x122 + V{0.0833333333333333}*cell[6];
auto x124 = cell[12]*x120 - V{0.166666666666667}*cell[12] + cell[3]*x120 + V{0.166666666666667}*cell[3] - x101 - x103 - x107*x119 + x111 + x116 - x117 + x118 - x121 + x123 - x96 + x97;
auto x125 = x33*x49;
auto x126 = V{0.5}*x125;
auto x127 = -cell[10]*x126 - V{0.166666666666667}*cell[10] - cell[1]*x126 + V{0.166666666666667}*cell[1] - x102 + x106 + x107*x125 + x109 - x116 + x117 - x118 + x121 + x123 - x95 + x98;
auto x128 = -V{1}*cell.template getFieldComponent<descriptors::POROSITY>(0)*x113*x127*x19*x36*x49 + V{0.333333333333333}*cell[0]*x19*x64 + cell[10]*x83*(x51 - x53*x84 + x63) + cell[11]*x83*(x35 + x62 - x85 + x86) + cell[12]*x83*(x56 + x68 + x86 - x87) + cell[13]*x21*(x35 + x46 + x65) + cell[14]*x21*(-x37*x38*x67*x67 + x65 + x66) + cell[15]*x21*(x65 + x68 + x71) + cell[16]*x21*(-x37*x38*x73*x73 + x65 + x72) + cell[17]*x21*(x68 - x75 + x76) + cell[18]*x21*(-x37*x38*x78*x78 + x72 + x76) - V{0.0555555555555556}*cell[1]*x19*(V{3}*x36*x38*x88 - x90) - V{0.0555555555555556}*cell[2]*x19*(x85 + x93) - V{0.0555555555555556}*cell[3]*x19*(x87 + x92 + x94) - V{0.0277777777777778}*cell[4]*x19*(-x46 - x66 - x91) + cell[5]*x21*(-x37*x38*x80*x80 + x76 + x79) - V{0.0277777777777778}*cell[6]*x19*(-x71 - x72 - x91) + cell[7]*x21*(-x37*x38*x81*x81 + x79 + x82) - V{0.0277777777777778}*cell[8]*x19*(x75 + x93 + x94) + cell[9]*x21*(-x37*x38*x77*x77 + x66 + x82) + x112*x115*x29 + x115*x124*x60;
auto x129 = -cell.template getFieldComponent<opti::DJDF>(0) + cell[0]*x20 + x128;
auto x130 = x114*x33;
auto x131 = x127*x130;
auto x132 = x128 - x131;
auto x133 = -cell.template getFieldComponent<opti::DJDF>(10) + cell[1]*x20 + x132;
auto x134 = x112*x130;
auto x135 = -x134;
auto x136 = x128 + x135;
auto x137 = -cell.template getFieldComponent<opti::DJDF>(11) + cell[2]*x20 + x136;
auto x138 = x124*x130;
auto x139 = -x138;
auto x140 = x128 + x139;
auto x141 = -cell.template getFieldComponent<opti::DJDF>(12) + cell[3]*x20 + x140;
auto x142 = -cell.template getFieldComponent<opti::DJDF>(13) + cell[4]*x20 + x132 + x135;
auto x143 = -cell.template getFieldComponent<opti::DJDF>(14) + cell[5]*x20 + x132 + x134;
auto x144 = -cell.template getFieldComponent<opti::DJDF>(15) + cell[6]*x20 + x132 + x139;
auto x145 = -cell.template getFieldComponent<opti::DJDF>(16) + cell[7]*x20 + x132 + x138;
auto x146 = -cell.template getFieldComponent<opti::DJDF>(17) + cell[8]*x20 + x136 + x139;
auto x147 = -cell.template getFieldComponent<opti::DJDF>(18) + cell[9]*x20 + x136 + x138;
auto x148 = x128 + x131;
auto x149 = -cell.template getFieldComponent<opti::DJDF>(1) + cell[10]*x20 + x148;
auto x150 = x128 + x134;
auto x151 = -cell.template getFieldComponent<opti::DJDF>(2) + cell[11]*x20 + x150;
auto x152 = -cell.template getFieldComponent<opti::DJDF>(3) + cell[12]*x20 + x128 + x138;
auto x153 = -cell.template getFieldComponent<opti::DJDF>(4) + cell[13]*x20 + x134 + x148;
auto x154 = -cell.template getFieldComponent<opti::DJDF>(5) + cell[14]*x20 + x131 + x136;
auto x155 = -cell.template getFieldComponent<opti::DJDF>(6) + cell[15]*x20 + x138 + x148;
auto x156 = -cell.template getFieldComponent<opti::DJDF>(7) + cell[16]*x20 + x131 + x140;
auto x157 = -cell.template getFieldComponent<opti::DJDF>(8) + cell[17]*x20 + x138 + x150;
auto x158 = -cell.template getFieldComponent<opti::DJDF>(9) + cell[18]*x20 + x134 + x140;
cell[0] = -x129;
cell[1] = -x133;
cell[2] = -x137;
cell[3] = -x141;
cell[4] = -x142;
cell[5] = -x143;
cell[6] = -x144;
cell[7] = -x145;
cell[8] = -x146;
cell[9] = -x147;
cell[10] = -x149;
cell[11] = -x151;
cell[12] = -x152;
cell[13] = -x153;
cell[14] = -x154;
cell[15] = -x155;
cell[16] = -x156;
cell[17] = -x157;
cell[18] = -x158;
return { V{1} - x129, x129*x129 + x133*x133 + x137*x137 + x141*x141 + x142*x142 + x143*x143 + x144*x144 + x145*x145 + x146*x146 + x147*x147 + x149*x149 + x151*x151 + x152*x152 + x153*x153 + x154*x154 + x155*x155 + x156*x156 + x157*x157 + x158*x158 };
}
};

}

}
