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
struct CSE<Dual<ZouHeDynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, 1, -1>, T, descriptors::D3Q19<FIELDS...> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{1}*x20;
auto x22 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{-1};
auto x23 = V{1} / (x22);
auto x24 = cell.template getFieldComponent<opti::F>(0) + cell.template getFieldComponent<opti::F>(15) + cell.template getFieldComponent<opti::F>(16) + cell.template getFieldComponent<opti::F>(3) + cell.template getFieldComponent<opti::F>(6) + cell.template getFieldComponent<opti::F>(7) + V{1};
auto x25 = cell.template getFieldComponent<opti::F>(1) + cell.template getFieldComponent<opti::F>(10) + cell.template getFieldComponent<opti::F>(12);
auto x26 = V{2}*cell.template getFieldComponent<opti::F>(14) + V{2}*cell.template getFieldComponent<opti::F>(2) + V{2}*cell.template getFieldComponent<opti::F>(4) + V{2}*cell.template getFieldComponent<opti::F>(8) + V{2}*cell.template getFieldComponent<opti::F>(9) + x24 + x25;
auto x27 = x23*x26;
auto x28 = V{0.0277777777777778}*x27;
auto x29 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x30 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x31 = -x30;
auto x32 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x33 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x34 = V{1.5}*x33;
auto x35 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x36 = V{1.5}*x35;
auto x37 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x38 = V{1.5}*x37;
auto x39 = x36 + x38 + V{-1};
auto x40 = x34 + x39;
auto x41 = -x32 + x40;
auto x42 = x29 + x41 - V{4.5}*x31*x31;
auto x43 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x44 = V{4.5}*(x43*x43);
auto x45 = -x36;
auto x46 = V{1} - x38;
auto x47 = x45 + x46;
auto x48 = x32 + x47;
auto x49 = -x34;
auto x50 = x29 + x49;
auto x51 = x44 + x48 + x50;
auto x52 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x53 = V{4.5}*(x52*x52);
auto x54 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x55 = x49 + x54;
auto x56 = x48 + x53 + x55;
auto x57 = x28*x56;
auto x58 = x32 + x40;
auto x59 = x29 - x44 + x58;
auto x60 = -x53 + x54 + x58;
auto x61 = x28*x60;
auto x62 = -x29;
auto x63 = x58 + x62 - V{4.5}*x30*x30;
auto x64 = V{0.0555555555555556}*x27;
auto x65 = V{3}*x33;
auto x66 = x48 + x65;
auto x67 = x32 + x39 - x65;
auto x68 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x69 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x68;
auto x70 = x41 + x54 - V{4.5}*x69*x69;
auto x71 = -x54;
auto x72 = -x69;
auto x73 = x58 + x71 - V{4.5}*x72*x72;
auto x74 = cell.template getFieldComponent<opti::F>(3) + x28*x70 - x28*x73;
auto x75 = V{2}*cell.template getFieldComponent<opti::F>(14) + V{2}*cell.template getFieldComponent<opti::F>(2) + V{2}*cell.template getFieldComponent<opti::F>(4) + V{2}*cell.template getFieldComponent<opti::F>(8) + V{2}*cell.template getFieldComponent<opti::F>(9);
auto x76 = V{1} / (cell.template getFieldComponent<opti::F>(0) + cell.template getFieldComponent<opti::F>(15) + cell.template getFieldComponent<opti::F>(16) + cell.template getFieldComponent<opti::F>(6) + cell.template getFieldComponent<opti::F>(7) + x25 + x28*x42 - x28*x51 - x28*x59 - x28*x63 - x57 - x61 - x64*x66 - x64*x67 + x74 + x75 + V{1});
auto x77 = V{0.25}*x19;
auto x78 = x77 + V{-0.25};
auto x79 = -x78;
auto x80 = cell[8] - cell[9];
auto x81 = x79*x80;
auto x82 = V{0.25}*x81;
auto x83 = -cell.template getFieldComponent<opti::F>(15) + cell.template getFieldComponent<opti::F>(6);
auto x84 = cell.template getFieldComponent<opti::F>(16) - cell.template getFieldComponent<opti::F>(7);
auto x85 = cell.template getFieldComponent<opti::F>(12) - x57 - x61 - x74 - x83 - x84;
auto x86 = cell[14] - cell[4];
auto x87 = -x78*x86;
auto x88 = V{0.25}*x87;
auto x89 = -x26/x22;
auto x90 = V{0.0277777777777778}*x89;
auto x91 = -x59*x90;
auto x92 = x51*x90;
auto x93 = -cell.template getFieldComponent<opti::F>(16) + cell.template getFieldComponent<opti::F>(7);
auto x94 = -x42*x90;
auto x95 = -x63*x90;
auto x96 = cell.template getFieldComponent<opti::F>(1) + x94 - x95;
auto x97 = -cell.template getFieldComponent<opti::F>(10) + x83 + x91 - x92 + x93 + x96;
auto x98 = x82*x85 + x88*x97;
auto x99 = x76*x85;
auto x100 = x76*x97;
auto x101 = V{1}*x100*x87 - V{4}*x76*x98 + V{1}*x81*x99;
auto x102 = V{0.75}*x20;
auto x103 = V{0.25}*x20;
auto x104 = cell[14]*x103;
auto x105 = x101 + x104;
auto x106 = x51*(-cell[4]*x102 - x105);
auto x107 = cell[9]*x103;
auto x108 = x101 + x107;
auto x109 = x56*(-cell[8]*x102 - x108);
auto x110 = -x20;
auto x111 = V{0.25}*x110;
auto x112 = V{0.75}*x110;
auto x113 = -x70*x90;
auto x114 = V{0.0555555555555556}*x89;
auto x115 = -x73*x90;
auto x116 = cell.template getFieldComponent<opti::F>(10) - x91 + x92;
auto x117 = cell.template getFieldComponent<opti::F>(12) + x56*x90 + x60*x90;
auto x118 = V{1} / (x113 + x114*x66 + x114*x67 - x115 + x116 + x117 + x24 + x75 + x96);
auto x119 = cell.template getFieldComponent<opti::F>(15) - cell.template getFieldComponent<opti::F>(6);
auto x120 = -cell.template getFieldComponent<opti::F>(1) + x116 + x119 + x84 - x94 + x95;
auto x121 = -x80;
auto x122 = -cell.template getFieldComponent<opti::F>(3) - x113 + x115 + x117 + x119 + x93;
auto x123 = V{1}*x79;
auto x124 = x123*x86;
auto x125 = x118*x120*x124 + x118*x121*x122*x123 - x118*x79*(x120*x86 + x121*x122);
auto x126 = x59*(cell[14]*x111 + cell[4]*x112 + x125);
auto x127 = x63*(cell[14]*x112 + cell[4]*x111 + x125);
auto x128 = x73*(cell[8]*x111 + cell[9]*x112 + x125);
auto x129 = x60*(cell[8]*x112 + cell[9]*x111 + x125);
auto x130 = V{1}*cell[2];
auto x131 = x130*x20;
auto x132 = x66*(-x101 - x131);
auto x133 = x67*(x110*x130 + x125);
auto x134 = cell[4]*x103;
auto x135 = x101 + x134;
auto x136 = x42*(-cell[14]*x102 - x135);
auto x137 = cell[8]*x103;
auto x138 = x101 + x137;
auto x139 = x70*(-cell[9]*x102 - x138);
auto x140 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x121*x78;
auto x141 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x124 + V{0.0277777777777778}*x106 + V{0.0277777777777778}*x109 + V{0.0277777777777778}*x126 + V{0.0277777777777778}*x127 + V{0.0277777777777778}*x128 + V{0.0277777777777778}*x129 + V{0.0555555555555556}*x132 + V{0.0555555555555556}*x133 - V{0.0277777777777778}*x136 - V{0.0277777777777778}*x139 + V{1}*x140;
auto x142 = x19*x23;
auto x143 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x144 = V{4.5}*(x143*x143);
auto x145 = x29 + x40;
auto x146 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x68;
auto x147 = -x146;
auto x148 = V{3}*x35;
auto x149 = x34 + V{-1};
auto x150 = V{3}*x37;
auto x151 = -(V{0.333333333333333}*cell[0]*x40 + V{0.0555555555555556}*cell[10]*(-x148 + x149 + x29 + x38) + V{0.0555555555555556}*cell[11]*x67 + V{0.0555555555555556}*cell[12]*(x149 - x150 + x36 + x54) + V{0.0277777777777778}*cell[13]*x59 + V{0.0277777777777778}*cell[14]*x42 + V{0.0277777777777778}*cell[15]*(-x144 + x145 + x54) + V{0.0277777777777778}*cell[16]*(x145 + x71 - V{4.5}*x147*x147) + V{0.0277777777777778}*cell[17]*x60 + V{0.0277777777777778}*cell[18]*x73 - V{0.0555555555555556}*cell[1]*(x148 + x46 + x50) - V{0.0555555555555556}*cell[2]*x66 - V{0.0555555555555556}*cell[3]*(x150 + x45 + x55 + V{1}) - V{0.0277777777777778}*cell[4]*x51 + V{0.0277777777777778}*cell[5]*x63 - V{0.0277777777777778}*cell[6]*(x144 + x47 + x50 + x54) + V{0.0277777777777778}*cell[7]*(x40 + x54 + x62 - V{4.5}*x146*x146) - V{0.0277777777777778}*cell[8]*x56 + V{0.0277777777777778}*cell[9]*x70);
auto x152 = x141*x23 + V{1}*x142*x151;
auto x153 = -cell.template getFieldComponent<opti::DJDF>(0) + cell[0]*x21 + x101 + x152;
auto x154 = x105 - x134 + x152;
auto x155 = -x137;
auto x156 = x108 + x152 + x155;
auto x157 = V{2}*x100*x87 + V{2}*x142*x151 + x23*(-V{2}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x79*x86 + V{0.0555555555555556}*x106 + V{0.0555555555555556}*x109 + V{0.0555555555555556}*x126 + V{0.0555555555555556}*x127 + V{0.0555555555555556}*x128 + V{0.0555555555555556}*x129 + V{0.111111111111111}*x132 + V{0.111111111111111}*x133 - V{0.0555555555555556}*x136 - V{0.0555555555555556}*x139 + V{2}*x140) - V{4}*x76*(V{0.5}*x81*x85 + V{0.5}*x87*x97) + V{2}*x81*x99;
auto x158 = -x107;
auto x159 = x137 + x158;
auto x160 = -x104;
auto x161 = x135 + x152 + x160;
auto x162 = V{0.25}*x141*x23 + x151*x23*x77 - x76*x98;
auto x163 = x100*x88 + x162 + x82*x99;
auto x164 = -V{0.25}*cell.template getFieldComponent<opti::DJDF>(0) + cell[0]*x103 + x163;
auto x165 = V{0.0625}*x20;
auto x166 = cell[4]*x165;
auto x167 = cell[14]*x165;
auto x168 = x163 + x166 - x167;
auto x169 = -V{0.25}*cell.template getFieldComponent<opti::DJDF>(1) + cell[10]*x103 + x168;
auto x170 = x163 - x166 + x167;
auto x171 = -V{0.25}*cell.template getFieldComponent<opti::DJDF>(10) + cell[1]*x103 + x170;
auto x172 = cell[9]*x165;
auto x173 = cell[8]*x165;
auto x174 = x172 - x173;
auto x175 = -V{0.25}*cell.template getFieldComponent<opti::DJDF>(12) + cell[3]*x103 + x163 + x174;
auto x176 = -x172 + x173;
auto x177 = -V{0.25}*cell.template getFieldComponent<opti::DJDF>(3) + cell[12]*x103 + x163 + x176;
auto x178 = -V{0.25}*cell.template getFieldComponent<opti::DJDF>(15) + cell[6]*x103 + x170 + x174;
auto x179 = -V{0.25}*cell.template getFieldComponent<opti::DJDF>(16) + cell[7]*x103 + x170 + x176;
auto x180 = -V{0.25}*cell.template getFieldComponent<opti::DJDF>(6) + cell[15]*x103 + x168 + x176;
auto x181 = -V{0.25}*cell.template getFieldComponent<opti::DJDF>(7) + cell[16]*x103 + x168 + x174;
auto x182 = V{0.125}*x20;
auto x183 = V{0.25}*x100*x87 + x162 + V{0.25}*x81*x99;
auto x184 = -V{0.125}*cell.template getFieldComponent<opti::DJDF>(14) + cell[14]*x182 + cell[5]*x182 + x183;
auto x185 = -V{0.125}*cell.template getFieldComponent<opti::DJDF>(2) + cell[11]*x182 + cell[2]*x182 + x183;
auto x186 = -V{0.125}*cell.template getFieldComponent<opti::DJDF>(4) + cell[13]*x182 + cell[4]*x182 + x183;
auto x187 = -V{0.125}*cell.template getFieldComponent<opti::DJDF>(8) + cell[17]*x182 + cell[8]*x182 + x183;
auto x188 = -V{0.125}*cell.template getFieldComponent<opti::DJDF>(9) + cell[18]*x182 + cell[9]*x182 + x183;
cell[0] = -x153;
cell[1] = cell.template getFieldComponent<opti::DJDF>(10) - cell[1]*x21 - x154;
cell[2] = cell.template getFieldComponent<opti::DJDF>(11);
cell[3] = cell.template getFieldComponent<opti::DJDF>(12) - cell[3]*x21 - x156;
cell[4] = cell.template getFieldComponent<opti::DJDF>(13);
cell[5] = cell.template getFieldComponent<opti::DJDF>(14) - cell[14]*x21 - cell[5]*x21 - x157;
cell[6] = cell.template getFieldComponent<opti::DJDF>(15) - cell[6]*x21 - x107 - x154 - x155;
cell[7] = cell.template getFieldComponent<opti::DJDF>(16) - cell[7]*x21 - x154 - x159;
cell[8] = cell.template getFieldComponent<opti::DJDF>(17);
cell[9] = cell.template getFieldComponent<opti::DJDF>(18);
cell[10] = cell.template getFieldComponent<opti::DJDF>(1) - cell[10]*x21 - x161;
cell[11] = cell.template getFieldComponent<opti::DJDF>(2) - cell[11]*x21 - x131 - x157;
cell[12] = cell.template getFieldComponent<opti::DJDF>(3) - cell[12]*x21 - x138 - x152 - x158;
cell[13] = cell.template getFieldComponent<opti::DJDF>(4) - cell[13]*x21 - cell[4]*x21 - x157;
cell[14] = cell.template getFieldComponent<opti::DJDF>(5);
cell[15] = cell.template getFieldComponent<opti::DJDF>(6) - cell[15]*x21 - x159 - x161;
cell[16] = cell.template getFieldComponent<opti::DJDF>(7) - cell[16]*x21 - x134 - x156 - x160;
cell[17] = cell.template getFieldComponent<opti::DJDF>(8) - cell[17]*x21 - cell[8]*x21 - x157;
cell[18] = cell.template getFieldComponent<opti::DJDF>(9) - cell[18]*x21 - cell[9]*x21 - x157;
return { V{1} - x153, cell.template getFieldComponent<opti::DJDF>(11)*cell.template getFieldComponent<opti::DJDF>(11) + cell.template getFieldComponent<opti::DJDF>(13)*cell.template getFieldComponent<opti::DJDF>(13) + cell.template getFieldComponent<opti::DJDF>(17)*cell.template getFieldComponent<opti::DJDF>(17) + cell.template getFieldComponent<opti::DJDF>(18)*cell.template getFieldComponent<opti::DJDF>(18) + cell.template getFieldComponent<opti::DJDF>(5)*cell.template getFieldComponent<opti::DJDF>(5) + 16*(x164*x164) + 16*(x169*x169) + 16*(x171*x171) + 16*(x175*x175) + 16*(x177*x177) + 16*(x178*x178) + 16*(x179*x179) + 16*(x180*x180) + 16*(x181*x181) + V{64}*(x184*x184) + V{64}*(x185*x185) + V{64}*(x186*x186) + V{64}*(x187*x187) + V{64}*(x188*x188) };
}
};

}

}
