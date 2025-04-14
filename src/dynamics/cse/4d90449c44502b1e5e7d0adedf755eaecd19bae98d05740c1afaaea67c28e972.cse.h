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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerCornerDensity3D<1, -1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerCornerStress3D<1, -1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{-1};
auto x22 = V{1} / (x21);
auto x23 = V{0.1666665}*cell[0];
auto x24 = V{0.1666665}*cell[12] + V{0.333333}*cell[14] + V{0.1666665}*cell[3] + x23 + V{0.1666665};
auto x25 = V{0.1666665}*cell[10] + V{0.1666665}*cell[1] + V{0.333333}*cell[9];
auto x26 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{1});
auto x27 = V{0.1666665}*cell[11] + V{0.333333}*cell[15] + V{0.1666665}*cell[2];
auto x28 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{1});
auto x29 = cell[0] + cell[11] + V{2}*cell[15] + cell[2] + V{1};
auto x30 = cell[12] + V{2}*cell[14] + cell[3];
auto x31 = V{2}*cell[10] + V{2}*cell[13] + V{2}*cell[16] + cell[17] + cell[18] + cell[8] + cell[9] + x29 + x30;
auto x32 = cell[10] + cell[1] + V{2}*cell[9];
auto x33 = V{2}*cell[12] + cell[13] + cell[14] + V{2}*cell[17] + cell[4] + cell[5] + V{2}*cell[7] + x29 + x32;
auto x34 = cell[0] + cell[15] + cell[16] + V{2}*cell[2] + V{2}*cell[4] + cell[6] + cell[7] + V{2}*cell[8] + x30 + x32 + V{1};
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
auto x46 = V{3.08148791101958e-33}*x26*x28*x31*x33;
auto x47 = V{5.55111512312578e-17}*x34*x44*x46;
auto x48 = -V{0.00925925925925926}*x22*x34 + V{0.00925925925925926}*x26*x31 + V{0.00925925925925926}*x28*x33;
auto x49 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x50 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x49;
auto x51 = -x50;
auto x52 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x53 = -x52;
auto x54 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x55 = x42 + x54;
auto x56 = x53 + x55;
auto x57 = x56 - V{4.5}*x51*x51;
auto x58 = V{1.36955018267537e-33}*x26*x28*x31*x33;
auto x59 = V{3.70074341541719e-17}*x34*x44*x58;
auto x60 = V{5.47820073070147e-33}*x26*x28*x31*x33;
auto x61 = V{7.40148683083438e-17}*x34*x44*x60;
auto x62 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x63 = V{4.5}*(x62*x62);
auto x64 = x52 + x55 - x63;
auto x65 = -x48*x64;
auto x66 = -x36;
auto x67 = V{1} - x38;
auto x68 = x66 + x67;
auto x69 = x54 + x68;
auto x70 = -x40;
auto x71 = x52 + x70;
auto x72 = x63 + x69 + x71;
auto x73 = x48*x72;
auto x74 = -x54;
auto x75 = -V{4.5}*x50*x50;
auto x76 = x42 + x52;
auto x77 = x74 + x75 + x76;
auto x78 = -x48*x77;
auto x79 = -V{0.0185185185185185}*x22*x34 + V{0.0185185185185185}*x26*x31 + V{0.0185185185185185}*x28*x33;
auto x80 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x81 = V{4.5}*(x80*x80);
auto x82 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x83 = x70 + x82;
auto x84 = x69 + x81 + x83;
auto x85 = -V{0.037037037037037}*x22*x34 + V{0.037037037037037}*x26*x31 + V{0.037037037037037}*x28*x33;
auto x86 = V{3}*x37;
auto x87 = x66 + x83 + x86 + V{1};
auto x88 = -x82;
auto x89 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x49;
auto x90 = -V{4.5}*x89*x89;
auto x91 = x76 + x88 + x90;
auto x92 = V{3}*x39;
auto x93 = x69 + x92;
auto x94 = x79*x93;
auto x95 = V{1.38777878078145e-17}*cell[0];
auto x96 = V{1.38777878078145e-17}*cell[12] + V{2.77555756156289e-17}*cell[14] + V{1.38777878078145e-17}*cell[3] + x95 + V{1.38777878078145e-17};
auto x97 = V{1.38777878078145e-17}*cell[10] + V{1.38777878078145e-17}*cell[1] + V{2.77555756156289e-17}*cell[9];
auto x98 = x22*(V{1.38777878078145e-17}*cell[15] + V{1.38777878078145e-17}*cell[16] + V{2.77555756156289e-17}*cell[2] + V{2.77555756156289e-17}*cell[4] + V{1.38777878078145e-17}*cell[6] + V{1.38777878078145e-17}*cell[7] + V{2.77555756156289e-17}*cell[8] + x96 + x97);
auto x99 = V{1.38777878078145e-17}*cell[11] + V{2.77555756156289e-17}*cell[15] + V{1.38777878078145e-17}*cell[2];
auto x100 = -x26*(V{2.77555756156289e-17}*cell[10] + V{2.77555756156289e-17}*cell[13] + V{2.77555756156289e-17}*cell[16] + V{1.38777878078145e-17}*cell[17] + V{1.38777878078145e-17}*cell[18] + V{1.38777878078145e-17}*cell[8] + V{1.38777878078145e-17}*cell[9] + x96 + x99);
auto x101 = -x28*(V{2.77555756156289e-17}*cell[12] + V{1.38777878078145e-17}*cell[13] + V{1.38777878078145e-17}*cell[14] + V{2.77555756156289e-17}*cell[17] + V{1.38777878078145e-17}*cell[4] + V{1.38777878078145e-17}*cell[5] + V{2.77555756156289e-17}*cell[7] + x95 + x97 + x99 + V{1.38777878078145e-17});
auto x102 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x103 = V{4.5}*(x102*x102);
auto x104 = x103 + x68 + x71 + x82;
auto x105 = x36 + V{-1};
auto x106 = x105 + x38 + x54 - x92;
auto x107 = x106*x79;
auto x108 = -x43;
auto x109 = x100 + x101 + x104*x79 - x107 + x108 + x94 + x98 + V{4.44089209850063e-16};
auto x110 = V{3}*x35;
auto x111 = x110 + x67 + x71;
auto x112 = x111*x79;
auto x113 = -x110 + x41 + x52;
auto x114 = x113*x79;
auto x115 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x116 = -V{4.5}*x115*x115;
auto x117 = x116 + x55 + x88;
auto x118 = V{3.33066907387547e-16}*cell[14] + V{3.33066907387547e-16}*cell[5] + x112 - x114 - x117*x79;
auto x119 = V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{3.33066907387547e-16}*cell[13] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + x109 + x118 + x65 + x73 + x78 + x79*x84 - x79*x91 + x85*x87;
auto x120 = -x115;
auto x121 = x42 + x82;
auto x122 = x121 + x74;
auto x123 = x122 - V{4.5}*x120*x120;
auto x124 = -x89;
auto x125 = x121 + x53;
auto x126 = x125 - V{4.5}*x124*x124;
auto x127 = x79*x87;
auto x128 = x105 + x40 + x82 - x86;
auto x129 = x128*x79;
auto x130 = V{1.66533453693773e-16}*cell[10] + V{2.22044604925031e-16}*cell[17] + V{1.66533453693773e-16}*cell[1] + V{2.22044604925031e-16}*cell[8] + x127 - x129;
auto x131 = x130 - x37*x47 - x57*x79;
auto x132 = x55 - x81 + x82;
auto x133 = -x132*x48;
auto x134 = x48*x84;
auto x135 = -x117*x48;
auto x136 = V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.44089209850063e-16}*cell[13] + V{4.9960036108132e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.44089209850063e-16}*cell[4] + V{4.9960036108132e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + x109 + x111*x85 + x133 + x134 + x135 + x72*x79;
auto x137 = V{1.04856185861083e-33}*x26*x28*x31*x33;
auto x138 = -x103 + x76 + x82;
auto x139 = -x138*x48;
auto x140 = x104*x48;
auto x141 = -x48*x91;
auto x142 = V{2.22044604925031e-16}*cell[11] + V{9.71445146547012e-17}*cell[12] + V{1.11022302462516e-16}*cell[13] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{2.22044604925031e-16}*cell[18] + V{2.22044604925031e-16}*cell[2] + V{9.71445146547012e-17}*cell[3] + V{1.11022302462516e-16}*cell[4] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{2.22044604925031e-16}*cell[9] + x100 + x101 - x106*x85 + x108 + x118 - x132*x79 + x139 + x140 + x141 - x64*x79 + x98 + V{4.44089209850063e-16};
auto x143 = -x22*(-x126*x48 + x131 - V{3.23815048849004e-17}*x137*x34*x35*x44 + x142 - x39*x61) + x26*(x119 - x35*x59 - x37*x61 - x39*x47 - x48*x57) + x28*(-x123*x48 - x126*x79 + x131 + x136 - x35*x61 - x39*x59);
auto x144 = V{0.166666666666667}*cell[10];
auto x145 = V{0.166666666666667}*cell[1];
auto x146 = V{0.0555555555555556}*x26*x31 + V{0.0555555555555556}*x28*x33 + V{0.0555555555555556}*x34*x44;
auto x147 = V{0.0833333333333334}*cell[13];
auto x148 = V{0.0833333333333334}*cell[14];
auto x149 = V{0.0833333333333334}*cell[4];
auto x150 = V{0.0833333333333334}*cell[5];
auto x151 = V{3.46944695195361e-18}*cell[0];
auto x152 = V{3.46944695195361e-18}*cell[11] + V{6.93889390390723e-18}*cell[15] + V{3.46944695195361e-18}*cell[2] + x151 + V{3.46944695195361e-18};
auto x153 = V{3.46944695195361e-18}*cell[12] + V{6.93889390390723e-18}*cell[14] + V{3.46944695195361e-18}*cell[3];
auto x154 = V{6.93889390390723e-18}*cell[10] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[16] + V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x152 + x153;
auto x155 = x154*x26;
auto x156 = V{3.46944695195361e-18}*cell[10] + V{3.46944695195361e-18}*cell[1] + V{6.93889390390723e-18}*cell[9];
auto x157 = V{6.93889390390723e-18}*cell[12] + V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{6.93889390390723e-18}*cell[17] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + V{6.93889390390723e-18}*cell[7] + x152 + x156;
auto x158 = x157*x28;
auto x159 = x22*(V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[2] + V{6.93889390390723e-18}*cell[4] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + V{6.93889390390723e-18}*cell[8] + x151 + x153 + x156 + V{3.46944695195361e-18});
auto x160 = -x159;
auto x161 = V{0.0277777777777778}*x26*x31 + V{0.0277777777777778}*x28*x33 + V{0.0277777777777778}*x34*x44;
auto x162 = x161*x35;
auto x163 = V{0.0833333333333333}*cell[12] + V{0.0833333333333333}*cell[3] - x147 - x148 - x149 - x150 + x155 + x158 + x160 - x162 + V{0.0555555555555555};
auto x164 = V{0.0833333333333334}*cell[15];
auto x165 = V{0.0833333333333334}*cell[16];
auto x166 = V{0.0833333333333334}*cell[6];
auto x167 = V{0.0833333333333334}*cell[7];
auto x168 = x161*x39;
auto x169 = V{0.0833333333333333}*cell[11] + V{0.0833333333333333}*cell[2] - x164 - x165 - x166 - x167 - x168;
auto x170 = V{0.166666666666667}*cell[11];
auto x171 = V{0.166666666666667}*cell[2];
auto x172 = V{0.0833333333333334}*cell[17];
auto x173 = V{0.0833333333333334}*cell[18];
auto x174 = V{0.0833333333333334}*cell[8];
auto x175 = V{0.0833333333333334}*cell[9];
auto x176 = x161*x37;
auto x177 = V{0.0833333333333333}*cell[10] + V{0.0833333333333333}*cell[1] - x172 - x173 - x174 - x175 - x176;
auto x178 = V{0.166666666666667}*cell[12];
auto x179 = V{0.166666666666667}*cell[3];
auto x180 = x26*x31 + x28*x33;
auto x181 = x180 - x22*x34;
auto x182 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x181;
auto x183 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x182;
auto x184 = -V{0.0138888888888889}*x22*x34 + V{0.0138888888888889}*x26*x31 + V{0.0138888888888889}*x28*x33;
auto x185 = V{0.013888875}*cell[0];
auto x186 = V{0.013888875}*cell[12] + V{0.02777775}*cell[14] + V{0.013888875}*cell[3] + x185 + V{0.013888875};
auto x187 = V{0.013888875}*cell[10] + V{0.013888875}*cell[1] + V{0.02777775}*cell[9];
auto x188 = x22*(V{0.013888875}*cell[15] + V{0.013888875}*cell[16] + V{0.02777775}*cell[2] + V{0.02777775}*cell[4] + V{0.013888875}*cell[6] + V{0.013888875}*cell[7] + V{0.02777775}*cell[8] + x186 + x187);
auto x189 = V{0.013888875}*cell[11] + V{0.02777775}*cell[15] + V{0.013888875}*cell[2];
auto x190 = -x26*(V{0.02777775}*cell[10] + V{0.02777775}*cell[13] + V{0.02777775}*cell[16] + V{0.013888875}*cell[17] + V{0.013888875}*cell[18] + V{0.013888875}*cell[8] + V{0.013888875}*cell[9] + x186 + x189);
auto x191 = -x28*(V{0.02777775}*cell[12] + V{0.013888875}*cell[13] + V{0.013888875}*cell[14] + V{0.02777775}*cell[17] + V{0.013888875}*cell[4] + V{0.013888875}*cell[5] + V{0.02777775}*cell[7] + x185 + x187 + x189 + V{0.013888875});
auto x192 = -V{0.0277777777777778}*x22*x34 + V{0.0277777777777778}*x26*x31 + V{0.0277777777777778}*x28*x33;
auto x193 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x188 + x190 + x191 - x192*x37 + V{0.0138888472222222};
auto x194 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] - x192*x39;
auto x195 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x184*x35 + x193 + x194;
auto x196 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] - x183 + x195;
auto x197 = -x48*(x116 + x122);
auto x198 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] + x183 + x195;
auto x199 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x182;
auto x200 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x192*x35;
auto x201 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x184*x39 + x193 + x200;
auto x202 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] - x199 + x201;
auto x203 = x125 + x90;
auto x204 = -x203*x48;
auto x205 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] + x199 + x201;
auto x206 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x181;
auto x207 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x184*x37 + x188 + x190 + x191 + x194 + x200 + V{0.0138888472222222};
auto x208 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x206 + x207;
auto x209 = x56 + x75;
auto x210 = -x209*x48;
auto x211 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x206 + x207;
auto x212 = -V{5.55111512312578e-17}*x22*x34*x46;
auto x213 = -V{3.70074341541719e-17}*x22*x34*x58;
auto x214 = -V{7.40148683083438e-17}*x22*x34*x60;
auto x215 = x130 - x209*x79 - x212*x37;
auto x216 = -x22*(V{3.23815048849004e-17}*x137*x22*x34*x35 + x142 + x204 - x214*x39 + x215) + x26*(x119 + x210 - x212*x39 - x213*x35 - x214*x37) + x28*(x136 + x197 - x203*x79 - x213*x39 - x214*x35 + x215);
auto x217 = x180 + x34*x44;
auto x218 = -x154*x26;
auto x219 = -x157*x28;
auto x220 = -V{0.0833333333333333}*cell[12] - V{0.0833333333333333}*cell[3] + x147 + x148 + x149 + x150 + x159 + x162 + x218 + x219 + V{-0.0555555555555555};
auto x221 = -V{0.0833333333333333}*cell[11] - V{0.0833333333333333}*cell[2] + x164 + x165 + x166 + x167 + x168;
auto x222 = -V{0.0833333333333333}*cell[10] - V{0.0833333333333333}*cell[1] + x172 + x173 + x174 + x175 + x176;
auto x0 = -x19*(V{0.111111111111111}*x143*x42 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x22*(V{0.1666665}*cell[15] + V{0.1666665}*cell[16] + V{0.333333}*cell[2] + V{0.333333}*cell[4] + V{0.1666665}*cell[6] + V{0.1666665}*cell[7] + V{0.333333}*cell[8] + x24 + x25) - x26*(V{0.333333}*cell[10] + V{0.333333}*cell[13] + V{0.333333}*cell[16] + V{0.1666665}*cell[17] + V{0.1666665}*cell[18] + V{0.1666665}*cell[8] + V{0.1666665}*cell[9] + x24 + x27) - x28*(V{0.333333}*cell[12] + V{0.1666665}*cell[13] + V{0.1666665}*cell[14] + V{0.333333}*cell[17] + V{0.1666665}*cell[4] + V{0.1666665}*cell[5] + V{0.333333}*cell[7] + x23 + x25 + x27 + V{0.1666665}) - x35*x45 - x37*x45 - x39*x45 + x43 + V{0.833332833333333});
auto x1 = -x19*(V{0.0185185185185185}*x128*x143 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x129 - x144 - x145 + x146*x37 + x163 + x169);
auto x2 = -x19*(V{0.0185185185185185}*x106*x143 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] + x107 + x146*x39 + x163 - x170 - x171 + x177);
auto x3 = -x19*(V{0.0185185185185185}*x113*x143 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] + x114 + x146*x35 + x155 + x158 + x160 + x169 + x177 - x178 - x179 + V{0.0555555555555555});
auto x4 = -x19*(V{0.00925925925925926}*x132*x143 + V{0.0277777777777778}) - x20*(x133 + x196);
auto x5 = -x19*(V{0.00925925925925926}*x123*x143 + V{0.0277777777777778}) - x20*(x197 + x198);
auto x6 = -x19*(V{0.00925925925925926}*x138*x143 + V{0.0277777777777778}) - x20*(x139 + x202);
auto x7 = -x19*(V{0.00925925925925926}*x126*x143 + V{0.0277777777777778}) - x20*(x204 + x205);
auto x8 = -x19*(V{0.00925925925925926}*x143*x64 + V{0.0277777777777778}) - x20*(x208 + x65);
auto x9 = -x19*(V{0.00925925925925926}*x143*x57 + V{0.0277777777777778}) - x20*(x210 + x211);
auto x10 = x19*(V{0.0185185185185185}*x216*x87 + V{-0.0555555555555556}) + x20*(V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] - x127 - x144 - x145 + V{0.0555555555555556}*x217*x37 - x220 - x221);
auto x11 = x19*(V{0.0185185185185185}*x216*x93 + V{-0.0555555555555556}) + x20*(V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] - x170 - x171 + V{0.0555555555555556}*x217*x39 - x220 - x222 - x94);
auto x12 = x19*(V{0.0185185185185185}*x111*x216 + V{-0.0555555555555556}) + x20*(V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] - x112 - x159 - x178 - x179 + V{0.0555555555555556}*x217*x35 - x218 - x219 - x221 - x222 + V{0.0555555555555555});
auto x13 = x19*(V{0.00925925925925926}*x216*x84 + V{-0.0277777777777778}) - x20*(x134 + x196);
auto x14 = -x19*(V{0.00925925925925926}*x117*x143 + V{0.0277777777777778}) - x20*(x135 + x198);
auto x15 = x19*(V{0.00925925925925926}*x104*x216 + V{-0.0277777777777778}) - x20*(x140 + x202);
auto x16 = -x19*(V{0.00925925925925926}*x143*x91 + V{0.0277777777777778}) - x20*(x141 + x205);
auto x17 = x19*(V{0.00925925925925926}*x216*x72 + V{-0.0277777777777778}) - x20*(x208 + x73);
auto x18 = -x19*(V{0.00925925925925926}*x143*x77 + V{0.0277777777777778}) - x20*(x211 + x78);
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
return { V{0.333333333333333}*x216, x35 + x37 + x39 };
}
};

}

}
