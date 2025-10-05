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
struct CSE<Dual<ZouHeDynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, 1, 1>, T, descriptors::D3Q19<FIELDS...> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{1}*x20;
auto x22 = cell[0]*x21;
auto x23 = V{0.25}*x19 + V{-0.25};
auto x24 = cell[17] - cell[18];
auto x25 = x23*x24;
auto x26 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{1});
auto x27 = cell.template getFieldComponent<opti::F>(10) + cell.template getFieldComponent<opti::F>(16);
auto x28 = cell.template getFieldComponent<opti::F>(1) + cell.template getFieldComponent<opti::F>(6);
auto x29 = cell.template getFieldComponent<opti::F>(0) + cell.template getFieldComponent<opti::F>(15) + cell.template getFieldComponent<opti::F>(3) + cell.template getFieldComponent<opti::F>(7) + V{1};
auto x30 = x26*(V{2}*cell.template getFieldComponent<opti::F>(11) + cell.template getFieldComponent<opti::F>(12) + V{2}*cell.template getFieldComponent<opti::F>(13) + V{2}*cell.template getFieldComponent<opti::F>(17) + V{2}*cell.template getFieldComponent<opti::F>(18) + V{2}*cell.template getFieldComponent<opti::F>(5) + x27 + x28 + x29);
auto x31 = V{0.0277777777777778}*x30;
auto x32 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x33 = -V{4.5}*x32*x32;
auto x34 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x35 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x36 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x37 = V{1.5}*x36;
auto x38 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x39 = V{1.5}*x38;
auto x40 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x41 = V{1.5}*x40;
auto x42 = x39 + x41 + V{-1};
auto x43 = x37 + x42;
auto x44 = -x35 + x43;
auto x45 = x34 + x44;
auto x46 = x33 + x45;
auto x47 = -x34;
auto x48 = x35 + x43;
auto x49 = x33 + x47 + x48;
auto x50 = x27 - x31*x49;
auto x51 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x52 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x51;
auto x53 = -V{4.5}*x52*x52;
auto x54 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x55 = -x54;
auto x56 = x48 + x55;
auto x57 = x53 + x56;
auto x58 = x44 + x53 + x54;
auto x59 = cell.template getFieldComponent<opti::F>(12) + x31*x58;
auto x60 = -x31*x57 + x59;
auto x61 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x62 = x61*x61;
auto x63 = x34 + x48 - V{4.5}*x62;
auto x64 = x31*x63;
auto x65 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x66 = x65*x65;
auto x67 = x48 + x54 - V{4.5}*x66;
auto x68 = x31*x67;
auto x69 = V{0.0555555555555556}*x30;
auto x70 = x35 - V{3}*x36 + x42;
auto x71 = -x39;
auto x72 = V{1} - x41;
auto x73 = x71 + x72;
auto x74 = x35 + x73;
auto x75 = -x37;
auto x76 = x34 + x75;
auto x77 = V{4.5}*x62 + x74 + x76;
auto x78 = x31*x77;
auto x79 = x28 - x78;
auto x80 = x54 + x75;
auto x81 = V{4.5}*x66 + x74 + x80;
auto x82 = x31*x81;
auto x83 = V{3}*x36 + x74;
auto x84 = V{2}*cell.template getFieldComponent<opti::F>(11) + V{2}*cell.template getFieldComponent<opti::F>(13) + V{2}*cell.template getFieldComponent<opti::F>(17) + V{2}*cell.template getFieldComponent<opti::F>(18) + V{2}*cell.template getFieldComponent<opti::F>(5) + x29 - x69*x83 - x82;
auto x85 = -x64 - x68 - x69*x70 + x79 + x84;
auto x86 = V{1} / (x31*x46 + x50 + x60 + x85);
auto x87 = cell.template getFieldComponent<opti::F>(15) - cell.template getFieldComponent<opti::F>(6);
auto x88 = -cell.template getFieldComponent<opti::F>(16) + cell.template getFieldComponent<opti::F>(7);
auto x89 = -cell.template getFieldComponent<opti::F>(3) + x82 + x87 + x88;
auto x90 = x68 + x89;
auto x91 = x86*(x60 + x90);
auto x92 = V{1}*x25*x91;
auto x93 = cell[13] - cell[5];
auto x94 = -cell.template getFieldComponent<opti::F>(1) - cell.template getFieldComponent<opti::F>(7) + x78 + x87;
auto x95 = -x32;
auto x96 = x45 - V{4.5}*x95*x95;
auto x97 = x31*x96 + x50;
auto x98 = -x64 - x94 - x97;
auto x99 = -x20;
auto x100 = V{0.25}*x99;
auto x101 = V{0.75}*x99;
auto x102 = -x23;
auto x103 = -x31*x67;
auto x104 = -x31*x49;
auto x105 = -x31*x96;
auto x106 = x104 - x105 + x27;
auto x107 = -x52;
auto x108 = x56 - V{4.5}*x107*x107;
auto x109 = cell.template getFieldComponent<opti::F>(12) - x108*x31 + x31*x58;
auto x110 = -x31*x63;
auto x111 = x110 + x79;
auto x112 = V{1} / (x103 + x106 + x109 + x111 - x69*x70 + x84);
auto x113 = x106 - x110 + x94;
auto x114 = -x103 + x109 + x89;
auto x115 = V{1}*x102;
auto x116 = x115*x93;
auto x117 = x115*x24;
auto x118 = -x102*x112*(x113*x93 + x114*x24) + x112*x113*x116 + x112*x114*x117;
auto x119 = x81*(cell[17]*x101 + cell[18]*x100 + x118);
auto x120 = x77*(cell[13]*x101 + cell[5]*x100 + x118);
auto x121 = V{0.25}*x20;
auto x122 = cell[13]*x121;
auto x123 = V{0.75}*x20;
auto x124 = -x108*x31 + x59;
auto x125 = x124 + x85 + x97;
auto x126 = V{1} / (x125);
auto x127 = -x24;
auto x128 = x124 + x90;
auto x129 = x127*x128;
auto x130 = -x93;
auto x131 = -cell.template getFieldComponent<opti::F>(10) - cell.template getFieldComponent<opti::F>(15) - x104 + x105 + x111 + x88;
auto x132 = x130*x131;
auto x133 = -V{1}*x126*x127*x128*x23 + V{1}*x126*x130*x131*x23 + x126*x23*(x129 - x132);
auto x134 = x49*(-cell[5]*x123 - x122 - x133);
auto x135 = cell[17]*x121;
auto x136 = x108*(-cell[18]*x123 - x133 - x135);
auto x137 = cell[18]*x121;
auto x138 = x67*(-cell[17]*x123 - x133 - x137);
auto x139 = cell[5]*x121;
auto x140 = x63*(-cell[13]*x123 - x133 - x139);
auto x141 = V{1}*cell[11];
auto x142 = x83*(x118 + x141*x99);
auto x143 = x141*x20;
auto x144 = x70*(-x133 - x143);
auto x145 = x96*(cell[13]*x100 + cell[5]*x101 + x118);
auto x146 = x58*(cell[17]*x100 + cell[18]*x101 + x118);
auto x147 = x26*(cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x116 + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x117 + V{0.0277777777777778}*x119 + V{0.0277777777777778}*x120 + V{0.0277777777777778}*x134 + V{0.0277777777777778}*x136 + V{0.0277777777777778}*x138 + V{0.0277777777777778}*x140 + V{0.0555555555555556}*x142 + V{0.0555555555555556}*x144 - V{0.0277777777777778}*x145 - V{0.0277777777777778}*x146);
auto x148 = x19*x26;
auto x149 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x150 = V{4.5}*(x149*x149);
auto x151 = x34 + x43;
auto x152 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x51;
auto x153 = -x152;
auto x154 = V{3}*x38;
auto x155 = x37 + V{-1};
auto x156 = V{3}*x40;
auto x157 = V{0.333333333333333}*cell[0]*x43 + V{0.0555555555555556}*cell[10]*(-x154 + x155 + x34 + x41) + V{0.0555555555555556}*cell[11]*x70 + V{0.0555555555555556}*cell[12]*(x155 - x156 + x39 + x54) + V{0.0277777777777778}*cell[13]*x63 + V{0.0277777777777778}*cell[14]*x96 + V{0.0277777777777778}*cell[15]*(-x150 + x151 + x54) + V{0.0277777777777778}*cell[16]*(x151 + x55 - V{4.5}*x153*x153) + V{0.0277777777777778}*cell[17]*x67 + V{0.0277777777777778}*cell[18]*x108 - V{0.0555555555555556}*cell[1]*(x154 + x72 + x76) - V{0.0555555555555556}*cell[2]*x83 - V{0.0555555555555556}*cell[3]*(x156 + x71 + x80 + V{1}) - V{0.0277777777777778}*cell[4]*x77 + V{0.0277777777777778}*cell[5]*x49 - V{0.0277777777777778}*cell[6]*(x150 + x54 + x73 + x76) + V{0.0277777777777778}*cell[7]*(x43 + x47 + x54 - V{4.5}*x152*x152) - V{0.0277777777777778}*cell[8]*x81 + V{0.0277777777777778}*cell[9]*x58;
auto x158 = V{1}*x148*x157;
auto x159 = V{0.0138888888888889}*x30;
auto x160 = V{1} / ((V{0.5}*cell.template getFieldComponent<opti::F>(0) + V{0.5}*cell.template getFieldComponent<opti::F>(1) + V{0.5}*cell.template getFieldComponent<opti::F>(10) + cell.template getFieldComponent<opti::F>(11) + V{0.5}*cell.template getFieldComponent<opti::F>(12) + cell.template getFieldComponent<opti::F>(13) + V{0.5}*cell.template getFieldComponent<opti::F>(15) + V{0.5}*cell.template getFieldComponent<opti::F>(16) + cell.template getFieldComponent<opti::F>(17) + cell.template getFieldComponent<opti::F>(18) + V{0.5}*cell.template getFieldComponent<opti::F>(3) + cell.template getFieldComponent<opti::F>(5) + V{0.5}*cell.template getFieldComponent<opti::F>(6) + V{0.5}*cell.template getFieldComponent<opti::F>(7) + x159*x46 - x159*x49 - x159*x57 + x159*x58 - x159*x63 - x159*x67 - x159*x77 - x159*x81 - x31*x70 - x31*x83 + V{0.5})*(V{0.5}*cell.template getFieldComponent<opti::F>(0) + V{0.5}*cell.template getFieldComponent<opti::F>(1) + V{0.5}*cell.template getFieldComponent<opti::F>(10) + cell.template getFieldComponent<opti::F>(11) + V{0.5}*cell.template getFieldComponent<opti::F>(12) + cell.template getFieldComponent<opti::F>(13) + V{0.5}*cell.template getFieldComponent<opti::F>(15) + V{0.5}*cell.template getFieldComponent<opti::F>(16) + cell.template getFieldComponent<opti::F>(17) + cell.template getFieldComponent<opti::F>(18) + V{0.5}*cell.template getFieldComponent<opti::F>(3) + cell.template getFieldComponent<opti::F>(5) + V{0.5}*cell.template getFieldComponent<opti::F>(6) + V{0.5}*cell.template getFieldComponent<opti::F>(7) + x159*x46 - x159*x49 - x159*x57 + x159*x58 - x159*x63 - x159*x67 - x159*x77 - x159*x81 - x31*x70 - x31*x83 + V{0.5}));
auto x161 = -x129 + x132;
auto x162 = -V{0.25}*x125*x160*x161*x23 + x147 + x158;
auto x163 = x162 - V{1}*x23*x86*x93*x98 + x92;
auto x164 = cell.template getFieldComponent<opti::DJDF>(0) - x163 - x22;
auto x165 = x122 - x139 + x163;
auto x166 = cell.template getFieldComponent<opti::DJDF>(10) - cell[1]*x21 - x165;
auto x167 = V{2}*x102;
auto x168 = -V{0.5}*x125*x160*x161*x23 + V{2}*x148*x157 - V{2}*x23*x86*x93*x98 + V{2}*x25*x91 + x26*(cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x167*x93 + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x167*x24 + V{0.0555555555555556}*x119 + V{0.0555555555555556}*x120 + V{0.0555555555555556}*x134 + V{0.0555555555555556}*x136 + V{0.0555555555555556}*x138 + V{0.0555555555555556}*x140 + V{0.111111111111111}*x142 + V{0.111111111111111}*x144 - V{0.0555555555555556}*x145 - V{0.0555555555555556}*x146);
auto x169 = x135 - x137;
auto x170 = x163 + x169;
auto x171 = cell.template getFieldComponent<opti::DJDF>(12) - cell[3]*x21 - x170;
auto x172 = cell.template getFieldComponent<opti::DJDF>(15) - cell[6]*x21 - x165 - x169;
auto x173 = -x135 + x137;
auto x174 = cell.template getFieldComponent<opti::DJDF>(16) - cell[7]*x21 - x165 - x173;
auto x175 = -x122 + x139;
auto x176 = x163 + x175;
auto x177 = cell.template getFieldComponent<opti::DJDF>(1) - cell[10]*x21 - x176;
auto x178 = cell.template getFieldComponent<opti::DJDF>(3) - cell[12]*x21 - x163 - x173;
auto x179 = cell.template getFieldComponent<opti::DJDF>(6) - cell[15]*x21 - x173 - x176;
auto x180 = cell.template getFieldComponent<opti::DJDF>(7) - cell[16]*x21 - x170 - x175;
auto x181 = V{0.5}*x20;
auto x182 = x162 - x23*x86*x93*x98 + x25*x91;
auto x183 = V{0.5}*cell.template getFieldComponent<opti::DJDF>(11) - cell[11]*x181 - cell[2]*x181 - x182;
auto x184 = V{0.5}*cell.template getFieldComponent<opti::DJDF>(13) - cell[13]*x181 - cell[4]*x181 - x182;
auto x185 = V{0.5}*cell.template getFieldComponent<opti::DJDF>(17) - cell[17]*x181 - cell[8]*x181 - x182;
auto x186 = V{0.5}*cell.template getFieldComponent<opti::DJDF>(18) - cell[18]*x181 - cell[9]*x181 - x182;
auto x187 = V{0.5}*cell.template getFieldComponent<opti::DJDF>(5) - cell[14]*x181 - cell[5]*x181 - x182;
cell[0] = x164;
cell[1] = x166;
cell[2] = cell.template getFieldComponent<opti::DJDF>(11) - cell[2]*x21 - x143 - x168;
cell[3] = x171;
cell[4] = cell.template getFieldComponent<opti::DJDF>(13) - cell[13]*x21 - cell[4]*x21 - x168;
cell[5] = cell.template getFieldComponent<opti::DJDF>(14);
cell[6] = x172;
cell[7] = x174;
cell[8] = cell.template getFieldComponent<opti::DJDF>(17) - cell[17]*x21 - cell[8]*x21 - x168;
cell[9] = cell.template getFieldComponent<opti::DJDF>(18) - cell[18]*x21 - cell[9]*x21 - x168;
cell[10] = x177;
cell[11] = cell.template getFieldComponent<opti::DJDF>(2);
cell[12] = x178;
cell[13] = cell.template getFieldComponent<opti::DJDF>(4);
cell[14] = cell.template getFieldComponent<opti::DJDF>(5) - cell[14]*x21 - cell[5]*x21 - x168;
cell[15] = x179;
cell[16] = x180;
cell[17] = cell.template getFieldComponent<opti::DJDF>(8);
cell[18] = cell.template getFieldComponent<opti::DJDF>(9);
return { cell.template getFieldComponent<opti::DJDF>(0) + x125*x160*x161*(V{0.0625}*x19 + V{-0.0625}) - x147 - x158 - x22 + V{1}*x23*x86*x93*x98 - x92 + V{1}, cell.template getFieldComponent<opti::DJDF>(14)*cell.template getFieldComponent<opti::DJDF>(14) + cell.template getFieldComponent<opti::DJDF>(2)*cell.template getFieldComponent<opti::DJDF>(2) + cell.template getFieldComponent<opti::DJDF>(4)*cell.template getFieldComponent<opti::DJDF>(4) + cell.template getFieldComponent<opti::DJDF>(8)*cell.template getFieldComponent<opti::DJDF>(8) + cell.template getFieldComponent<opti::DJDF>(9)*cell.template getFieldComponent<opti::DJDF>(9) + x164*x164 + x166*x166 + x171*x171 + x172*x172 + x174*x174 + x177*x177 + x178*x178 + x179*x179 + x180*x180 + V{4}*(x183*x183) + V{4}*(x184*x184) + V{4}*(x185*x185) + V{4}*(x186*x186) + V{4}*(x187*x187) };
}
};

}

}
