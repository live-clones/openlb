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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerCornerDensity3D<-1, 1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerCornerStress3D<-1, 1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{-1};
auto x22 = V{1} / (x21);
auto x23 = V{0.1666665}*cell[0];
auto x24 = V{0.1666665}*cell[12] + V{0.1666665}*cell[3] + V{0.333333}*cell[5] + x23 + V{0.1666665};
auto x25 = V{0.1666665}*cell[11] + V{0.1666665}*cell[2] + V{0.333333}*cell[7];
auto x26 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{1});
auto x27 = V{0.1666665}*cell[10] + V{0.333333}*cell[17] + V{0.1666665}*cell[1];
auto x28 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{1});
auto x29 = cell[0] + cell[10] + V{2}*cell[17] + cell[1] + V{1};
auto x30 = cell[12] + cell[3] + V{2}*cell[5];
auto x31 = V{2}*cell[11] + V{2}*cell[13] + cell[15] + cell[16] + V{2}*cell[18] + cell[6] + cell[7] + x29 + x30;
auto x32 = cell[11] + cell[2] + V{2}*cell[7];
auto x33 = V{2}*cell[12] + cell[13] + cell[14] + V{2}*cell[15] + cell[4] + cell[5] + V{2}*cell[9] + x29 + x32;
auto x34 = cell[0] + cell[17] + cell[18] + V{2}*cell[1] + V{2}*cell[4] + V{2}*cell[6] + cell[8] + cell[9] + x30 + x32 + V{1};
auto x35 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x36 = V{1.5}*x35;
auto x37 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x38 = V{1.5}*x37;
auto x39 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x40 = V{1.5}*x39;
auto x41 = x38 + x40 + V{-1};
auto x42 = x36 + x41;
auto x43 = x42*(-V{0.111111111111111}*x22*x34 + V{0.111111111111111}*x26*x31 + V{0.111111111111111}*x28*x33);
auto x44 = -1/x21;
auto x45 = V{0.166666666666667}*x26*x31 + V{0.166666666666667}*x28*x33 + V{0.166666666666667}*x34*x44;
auto x46 = -V{0.00925925925925926}*x22*x34 + V{0.00925925925925926}*x26*x31 + V{0.00925925925925926}*x28*x33;
auto x47 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x48 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x47;
auto x49 = -x48;
auto x50 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x51 = -x50;
auto x52 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x53 = x42 + x52;
auto x54 = x51 + x53;
auto x55 = x54 - V{4.5}*x49*x49;
auto x56 = V{5.47820073070147e-33}*x26*x28*x31*x33;
auto x57 = V{7.40148683083438e-17}*x34*x44*x56;
auto x58 = V{3.08148791101958e-33}*x26*x28*x31*x33;
auto x59 = V{5.55111512312578e-17}*x34*x44*x58;
auto x60 = -V{0.0185185185185185}*x22*x34 + V{0.0185185185185185}*x26*x31 + V{0.0185185185185185}*x28*x33;
auto x61 = V{3}*x37;
auto x62 = -x36;
auto x63 = V{1} - x40;
auto x64 = x62 + x63;
auto x65 = x52 + x64;
auto x66 = x61 + x65;
auto x67 = x60*x66;
auto x68 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x69 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x70 = V{4.5}*(x69*x69);
auto x71 = -x38;
auto x72 = x50 + x71;
auto x73 = x64 + x68 + x70 + x72;
auto x74 = x36 + V{-1};
auto x75 = x40 + x52 - x61 + x74;
auto x76 = x60*x75;
auto x77 = V{1.66533453693773e-16}*cell[10] + V{1.66533453693773e-16}*cell[1] + x60*x73 + x67 - x76;
auto x78 = -x37*x59 + x77;
auto x79 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x80 = -x79;
auto x81 = -x68;
auto x82 = x53 + x81;
auto x83 = x82 - V{4.5}*x80*x80;
auto x84 = V{1.36955018267537e-33}*x26*x28*x31*x33;
auto x85 = V{3.70074341541719e-17}*x34*x44*x84;
auto x86 = V{3}*x35;
auto x87 = x63 + x72 + x86;
auto x88 = x60*x87;
auto x89 = x41 + x50 - x86;
auto x90 = x60*x89;
auto x91 = V{1.11022302462516e-16}*cell[12] + V{3.33066907387547e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{3.33066907387547e-16}*cell[5] + x88 - x90;
auto x92 = -x35*x85 - x60*x83 + x91;
auto x93 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x94 = V{4.5}*(x93*x93);
auto x95 = x50 + x53 - x94;
auto x96 = -x46*x95;
auto x97 = x65 + x72 + x94;
auto x98 = x46*x97;
auto x99 = -x52;
auto x100 = -V{4.5}*x48*x48;
auto x101 = x42 + x50;
auto x102 = x100 + x101 + x99;
auto x103 = -x102*x46;
auto x104 = V{1.38777878078145e-17}*cell[0];
auto x105 = V{1.38777878078145e-17}*cell[12] + V{1.38777878078145e-17}*cell[3] + V{2.77555756156289e-17}*cell[5] + x104 + V{1.38777878078145e-17};
auto x106 = V{1.38777878078145e-17}*cell[11] + V{1.38777878078145e-17}*cell[2] + V{2.77555756156289e-17}*cell[7];
auto x107 = x22*(V{1.38777878078145e-17}*cell[17] + V{1.38777878078145e-17}*cell[18] + V{2.77555756156289e-17}*cell[1] + V{2.77555756156289e-17}*cell[4] + V{2.77555756156289e-17}*cell[6] + V{1.38777878078145e-17}*cell[8] + V{1.38777878078145e-17}*cell[9] + x105 + x106);
auto x108 = V{1.38777878078145e-17}*cell[10] + V{2.77555756156289e-17}*cell[17] + V{1.38777878078145e-17}*cell[1];
auto x109 = -x26*(V{2.77555756156289e-17}*cell[11] + V{2.77555756156289e-17}*cell[13] + V{1.38777878078145e-17}*cell[15] + V{1.38777878078145e-17}*cell[16] + V{2.77555756156289e-17}*cell[18] + V{1.38777878078145e-17}*cell[6] + V{1.38777878078145e-17}*cell[7] + x105 + x108);
auto x110 = -x28*(V{2.77555756156289e-17}*cell[12] + V{1.38777878078145e-17}*cell[13] + V{1.38777878078145e-17}*cell[14] + V{2.77555756156289e-17}*cell[15] + V{1.38777878078145e-17}*cell[4] + V{1.38777878078145e-17}*cell[5] + V{2.77555756156289e-17}*cell[9] + x104 + x106 + x108 + V{1.38777878078145e-17});
auto x111 = -x43;
auto x112 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x113 = V{4.5}*(x112*x112);
auto x114 = x68 + x71;
auto x115 = x113 + x114 + x65;
auto x116 = -V{0.037037037037037}*x22*x34 + V{0.037037037037037}*x26*x31 + V{0.037037037037037}*x28*x33;
auto x117 = V{3}*x39;
auto x118 = x114 + x117 + x62 + V{1};
auto x119 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x47;
auto x120 = -V{4.5}*x119*x119;
auto x121 = x101 + x120 + x81;
auto x122 = V{2.22044604925031e-16}*cell[11] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{4.44089209850063e-16}*cell[17] + V{2.22044604925031e-16}*cell[18] + V{2.22044604925031e-16}*cell[2] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[8] + V{2.22044604925031e-16}*cell[9] + x103 + x107 + x109 + x110 + x111 + x115*x60 + x116*x118 - x121*x60 + x96 + x98 + V{4.44089209850063e-16};
auto x123 = -x119;
auto x124 = x42 + x68;
auto x125 = x124 + x51;
auto x126 = x125 - V{4.5}*x123*x123;
auto x127 = x118*x60;
auto x128 = -x117 + x38 + x68 + x74;
auto x129 = x128*x60;
auto x130 = x107 + x109 + x110 + x111 + x127 - x129 + V{4.44089209850063e-16};
auto x131 = x130 - x55*x60;
auto x132 = -x113 + x53 + x68;
auto x133 = -x132*x46;
auto x134 = x115*x46;
auto x135 = -V{4.5}*x79*x79;
auto x136 = x124 + x135 + x99;
auto x137 = -x136*x46;
auto x138 = V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.44089209850063e-16}*cell[13] + V{4.9960036108132e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{2.22044604925031e-16}*cell[17] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.44089209850063e-16}*cell[4] + V{4.9960036108132e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + V{2.22044604925031e-16}*cell[8] + x116*x87 + x133 + x134 + x137 + x60*x97;
auto x139 = x101 + x68 - x70;
auto x140 = -x139*x46;
auto x141 = x46*x73;
auto x142 = -x121*x46;
auto x143 = V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] - x116*x75 - x132*x60 + x140 + x141 + x142 - x60*x95;
auto x144 = -x22*(-x126*x46 + x131 + x143 - x37*x57 - x39*x59 + x92) + x26*(x122 - x39*x57 - x46*x55 + x78 + x92) + x28*(-x126*x60 + x131 + x138 - x35*x57 - x39*x85 - x46*x83 + x78);
auto x145 = V{0.166666666666667}*cell[10];
auto x146 = V{0.166666666666667}*cell[1];
auto x147 = V{0.0555555555555556}*x26*x31 + V{0.0555555555555556}*x28*x33 + V{0.0555555555555556}*x34*x44;
auto x148 = V{0.0833333333333334}*cell[13];
auto x149 = V{0.0833333333333334}*cell[14];
auto x150 = V{0.0833333333333334}*cell[4];
auto x151 = V{0.0833333333333334}*cell[5];
auto x152 = V{3.46944695195361e-18}*cell[0];
auto x153 = V{3.46944695195361e-18}*cell[10] + V{6.93889390390723e-18}*cell[17] + V{3.46944695195361e-18}*cell[1] + x152 + V{3.46944695195361e-18};
auto x154 = V{3.46944695195361e-18}*cell[12] + V{3.46944695195361e-18}*cell[3] + V{6.93889390390723e-18}*cell[5];
auto x155 = V{6.93889390390723e-18}*cell[11] + V{6.93889390390723e-18}*cell[13] + V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[18] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + x153 + x154;
auto x156 = x155*x26;
auto x157 = V{3.46944695195361e-18}*cell[11] + V{3.46944695195361e-18}*cell[2] + V{6.93889390390723e-18}*cell[7];
auto x158 = V{6.93889390390723e-18}*cell[12] + V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{6.93889390390723e-18}*cell[15] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + V{6.93889390390723e-18}*cell[9] + x153 + x157;
auto x159 = x158*x28;
auto x160 = x22*(V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{6.93889390390723e-18}*cell[1] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[6] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x152 + x154 + x157 + V{3.46944695195361e-18});
auto x161 = -x160;
auto x162 = V{0.0277777777777778}*x26*x31 + V{0.0277777777777778}*x28*x33 + V{0.0277777777777778}*x34*x44;
auto x163 = x162*x35;
auto x164 = V{0.0833333333333333}*cell[12] + V{0.0833333333333333}*cell[3] - x148 - x149 - x150 - x151 + x156 + x159 + x161 - x163 + V{0.0555555555555555};
auto x165 = V{0.0833333333333334}*cell[15];
auto x166 = V{0.0833333333333334}*cell[16];
auto x167 = V{0.0833333333333334}*cell[6];
auto x168 = V{0.0833333333333334}*cell[7];
auto x169 = x162*x39;
auto x170 = V{0.0833333333333333}*cell[11] + V{0.0833333333333333}*cell[2] - x165 - x166 - x167 - x168 - x169;
auto x171 = V{0.166666666666667}*cell[11];
auto x172 = V{0.166666666666667}*cell[2];
auto x173 = V{0.0833333333333334}*cell[17];
auto x174 = V{0.0833333333333334}*cell[18];
auto x175 = V{0.0833333333333334}*cell[8];
auto x176 = V{0.0833333333333334}*cell[9];
auto x177 = x162*x37;
auto x178 = V{0.0833333333333333}*cell[10] + V{0.0833333333333333}*cell[1] - x173 - x174 - x175 - x176 - x177;
auto x179 = V{0.166666666666667}*cell[12];
auto x180 = V{0.166666666666667}*cell[3];
auto x181 = x26*x31 + x28*x33;
auto x182 = x181 - x22*x34;
auto x183 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x182;
auto x184 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x183;
auto x185 = -V{0.0138888888888889}*x22*x34 + V{0.0138888888888889}*x26*x31 + V{0.0138888888888889}*x28*x33;
auto x186 = V{0.013888875}*cell[0];
auto x187 = V{0.013888875}*cell[12] + V{0.013888875}*cell[3] + V{0.02777775}*cell[5] + x186 + V{0.013888875};
auto x188 = V{0.013888875}*cell[11] + V{0.013888875}*cell[2] + V{0.02777775}*cell[7];
auto x189 = x22*(V{0.013888875}*cell[17] + V{0.013888875}*cell[18] + V{0.02777775}*cell[1] + V{0.02777775}*cell[4] + V{0.02777775}*cell[6] + V{0.013888875}*cell[8] + V{0.013888875}*cell[9] + x187 + x188);
auto x190 = V{0.013888875}*cell[10] + V{0.02777775}*cell[17] + V{0.013888875}*cell[1];
auto x191 = -x26*(V{0.02777775}*cell[11] + V{0.02777775}*cell[13] + V{0.013888875}*cell[15] + V{0.013888875}*cell[16] + V{0.02777775}*cell[18] + V{0.013888875}*cell[6] + V{0.013888875}*cell[7] + x187 + x190);
auto x192 = -x28*(V{0.02777775}*cell[12] + V{0.013888875}*cell[13] + V{0.013888875}*cell[14] + V{0.02777775}*cell[15] + V{0.013888875}*cell[4] + V{0.013888875}*cell[5] + V{0.02777775}*cell[9] + x186 + x188 + x190 + V{0.013888875});
auto x193 = -V{0.0277777777777778}*x22*x34 + V{0.0277777777777778}*x26*x31 + V{0.0277777777777778}*x28*x33;
auto x194 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x189 + x191 + x192 - x193*x37 + V{0.0138888472222222};
auto x195 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] - x193*x39;
auto x196 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x185*x35 + x194 + x195;
auto x197 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] - x184 + x196;
auto x198 = x135 + x82;
auto x199 = -x198*x46;
auto x200 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] + x184 + x196;
auto x201 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x183;
auto x202 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x193*x35;
auto x203 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x185*x39 + x194 + x202;
auto x204 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] - x201 + x203;
auto x205 = x100 + x54;
auto x206 = -x205*x46;
auto x207 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] + x201 + x203;
auto x208 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x182;
auto x209 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x185*x37 + x189 + x191 + x192 + x195 + x202 + V{0.0138888472222222};
auto x210 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x208 + x209;
auto x211 = x120 + x125;
auto x212 = -x211*x46;
auto x213 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x208 + x209;
auto x214 = -V{7.40148683083438e-17}*x22*x34*x56;
auto x215 = -V{5.55111512312578e-17}*x22*x34*x58;
auto x216 = -x215*x37 + x77;
auto x217 = -V{3.70074341541719e-17}*x22*x34*x84;
auto x218 = -x198*x60 - x217*x35 + x91;
auto x219 = x130 - x205*x60;
auto x220 = -x22*(x143 + x212 - x214*x37 - x215*x39 + x218 + x219) + x26*(x122 + x206 - x214*x39 + x216 + x218) + x28*(x138 + x199 - x211*x60 - x214*x35 + x216 - x217*x39 + x219);
auto x221 = x181 + x34*x44;
auto x222 = -x155*x26;
auto x223 = -x158*x28;
auto x224 = -V{0.0833333333333333}*cell[12] - V{0.0833333333333333}*cell[3] + x148 + x149 + x150 + x151 + x160 + x163 + x222 + x223 + V{-0.0555555555555555};
auto x225 = -V{0.0833333333333333}*cell[11] - V{0.0833333333333333}*cell[2] + x165 + x166 + x167 + x168 + x169;
auto x226 = -V{0.0833333333333333}*cell[10] - V{0.0833333333333333}*cell[1] + x173 + x174 + x175 + x176 + x177;
auto x0 = -x19*(V{0.111111111111111}*x144*x42 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x22*(V{0.1666665}*cell[17] + V{0.1666665}*cell[18] + V{0.333333}*cell[1] + V{0.333333}*cell[4] + V{0.333333}*cell[6] + V{0.1666665}*cell[8] + V{0.1666665}*cell[9] + x24 + x25) - x26*(V{0.333333}*cell[11] + V{0.333333}*cell[13] + V{0.1666665}*cell[15] + V{0.1666665}*cell[16] + V{0.333333}*cell[18] + V{0.1666665}*cell[6] + V{0.1666665}*cell[7] + x24 + x27) - x28*(V{0.333333}*cell[12] + V{0.1666665}*cell[13] + V{0.1666665}*cell[14] + V{0.333333}*cell[15] + V{0.1666665}*cell[4] + V{0.1666665}*cell[5] + V{0.333333}*cell[9] + x23 + x25 + x27 + V{0.1666665}) - x35*x45 - x37*x45 - x39*x45 + x43 + V{0.833332833333333});
auto x1 = -x19*(V{0.0185185185185185}*x144*x75 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] - x145 - x146 + x147*x37 + x164 + x170 + x76);
auto x2 = -x19*(V{0.0185185185185185}*x128*x144 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] + x129 + x147*x39 + x164 - x171 - x172 + x178);
auto x3 = -x19*(V{0.0185185185185185}*x144*x89 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] + x147*x35 + x156 + x159 + x161 + x170 + x178 - x179 - x180 + x90 + V{0.0555555555555555});
auto x4 = -x19*(V{0.00925925925925926}*x132*x144 + V{0.0277777777777778}) - x20*(x133 + x197);
auto x5 = -x19*(V{0.00925925925925926}*x144*x83 + V{0.0277777777777778}) - x20*(x199 + x200);
auto x6 = -x19*(V{0.00925925925925926}*x144*x95 + V{0.0277777777777778}) - x20*(x204 + x96);
auto x7 = -x19*(V{0.00925925925925926}*x144*x55 + V{0.0277777777777778}) - x20*(x206 + x207);
auto x8 = -x19*(V{0.00925925925925926}*x139*x144 + V{0.0277777777777778}) - x20*(x140 + x210);
auto x9 = -x19*(V{0.00925925925925926}*x126*x144 + V{0.0277777777777778}) - x20*(x212 + x213);
auto x10 = x19*(V{0.0185185185185185}*x220*x66 + V{-0.0555555555555556}) + x20*(V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] - x145 - x146 + V{0.0555555555555556}*x221*x37 - x224 - x225 - x67);
auto x11 = x19*(V{0.0185185185185185}*x118*x220 + V{-0.0555555555555556}) + x20*(V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] - x127 - x171 - x172 + V{0.0555555555555556}*x221*x39 - x224 - x226);
auto x12 = x19*(V{0.0185185185185185}*x220*x87 + V{-0.0555555555555556}) + x20*(V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] - x160 - x179 - x180 + V{0.0555555555555556}*x221*x35 - x222 - x223 - x225 - x226 - x88 + V{0.0555555555555555});
auto x13 = x19*(V{0.00925925925925926}*x115*x220 + V{-0.0277777777777778}) - x20*(x134 + x197);
auto x14 = -x19*(V{0.00925925925925926}*x136*x144 + V{0.0277777777777778}) - x20*(x137 + x200);
auto x15 = x19*(V{0.00925925925925926}*x220*x97 + V{-0.0277777777777778}) - x20*(x204 + x98);
auto x16 = -x19*(V{0.00925925925925926}*x102*x144 + V{0.0277777777777778}) - x20*(x103 + x207);
auto x17 = x19*(V{0.00925925925925926}*x220*x73 + V{-0.0277777777777778}) - x20*(x141 + x210);
auto x18 = -x19*(V{0.00925925925925926}*x121*x144 + V{0.0277777777777778}) - x20*(x142 + x213);
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
return { V{0.333333333333333}*x220, x35 + x37 + x39 };
}
};

}

}
