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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerCornerDensity3D<1, 1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerCornerStress3D<1, 1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{1});
auto x22 = V{0.1666665}*cell[0];
auto x23 = V{0.1666665}*cell[12] + V{0.333333}*cell[13] + V{0.1666665}*cell[3] + x22 + V{0.1666665};
auto x24 = V{0.1666665}*cell[11] + V{0.333333}*cell[15] + V{0.1666665}*cell[2];
auto x25 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{1});
auto x26 = V{0.1666665}*cell[10] + V{0.333333}*cell[17] + V{0.1666665}*cell[1];
auto x27 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{1});
auto x28 = cell[0] + cell[12] + V{2}*cell[13] + cell[3] + V{1};
auto x29 = cell[11] + V{2}*cell[15] + cell[2];
auto x30 = V{2}*cell[10] + V{2}*cell[14] + V{2}*cell[16] + cell[17] + cell[18] + cell[8] + cell[9] + x28 + x29;
auto x31 = cell[10] + V{2}*cell[17] + cell[1];
auto x32 = V{2}*cell[11] + cell[15] + cell[16] + V{2}*cell[18] + V{2}*cell[5] + cell[6] + cell[7] + x28 + x31;
auto x33 = cell[0] + V{2}*cell[12] + cell[13] + cell[14] + cell[4] + cell[5] + V{2}*cell[7] + V{2}*cell[9] + x29 + x31 + V{1};
auto x34 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x35 = V{1.5}*x34;
auto x36 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x37 = V{1.5}*x36;
auto x38 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x39 = V{1.5}*x38;
auto x40 = x37 + x39 + V{-1};
auto x41 = x35 + x40;
auto x42 = x41*(V{0.111111111111111}*x21*x30 + V{0.111111111111111}*x25*x32 + V{0.111111111111111}*x27*x33);
auto x43 = V{0.166666666666667}*x21*x30 + V{0.166666666666667}*x25*x32 + V{0.166666666666667}*x27*x33;
auto x44 = V{0.00925925925925926}*x21*x30 + V{0.00925925925925926}*x25*x32 + V{0.00925925925925926}*x27*x33;
auto x45 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x46 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x45;
auto x47 = -x46;
auto x48 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x49 = -x48;
auto x50 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x51 = x41 + x50;
auto x52 = x49 + x51;
auto x53 = x52 - V{4.5}*x47*x47;
auto x54 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x55 = V{4.5}*(x54*x54);
auto x56 = x48 + x51 - x55;
auto x57 = -x44*x56;
auto x58 = -x35;
auto x59 = V{1} - x37;
auto x60 = x58 + x59;
auto x61 = x50 + x60;
auto x62 = -x39;
auto x63 = x48 + x62;
auto x64 = x55 + x61 + x63;
auto x65 = x44*x64;
auto x66 = -x50;
auto x67 = -V{4.5}*x46*x46;
auto x68 = x41 + x48;
auto x69 = x66 + x67 + x68;
auto x70 = -x44*x69;
auto x71 = V{0.037037037037037}*x21*x30 + V{0.037037037037037}*x25*x32 + V{0.037037037037037}*x27*x33;
auto x72 = V{3}*x36;
auto x73 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x74 = x62 + x73;
auto x75 = x58 + x72 + x74 + V{1};
auto x76 = V{5.55111512312578e-17}*x21*x30 + V{5.55111512312578e-17}*x25*x32 + V{5.55111512312578e-17}*x27*x33;
auto x77 = V{0.0185185185185185}*x21*x30 + V{0.0185185185185185}*x25*x32 + V{0.0185185185185185}*x27*x33;
auto x78 = -x73;
auto x79 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x80 = -V{4.5}*x79*x79;
auto x81 = x51 + x78 + x80;
auto x82 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x45;
auto x83 = -V{4.5}*x82*x82;
auto x84 = x68 + x78 + x83;
auto x85 = V{7.40148683083438e-17}*x21*x30 + V{7.40148683083438e-17}*x25*x32 + V{7.40148683083438e-17}*x27*x33;
auto x86 = V{1.38777878078145e-17}*cell[0];
auto x87 = V{1.38777878078145e-17}*cell[12] + V{2.77555756156289e-17}*cell[13] + V{1.38777878078145e-17}*cell[3] + x86 + V{1.38777878078145e-17};
auto x88 = V{1.38777878078145e-17}*cell[11] + V{2.77555756156289e-17}*cell[15] + V{1.38777878078145e-17}*cell[2];
auto x89 = -x21*(V{2.77555756156289e-17}*cell[10] + V{2.77555756156289e-17}*cell[14] + V{2.77555756156289e-17}*cell[16] + V{1.38777878078145e-17}*cell[17] + V{1.38777878078145e-17}*cell[18] + V{1.38777878078145e-17}*cell[8] + V{1.38777878078145e-17}*cell[9] + x87 + x88);
auto x90 = V{1.38777878078145e-17}*cell[10] + V{2.77555756156289e-17}*cell[17] + V{1.38777878078145e-17}*cell[1];
auto x91 = -x25*(V{2.77555756156289e-17}*cell[11] + V{1.38777878078145e-17}*cell[15] + V{1.38777878078145e-17}*cell[16] + V{2.77555756156289e-17}*cell[18] + V{2.77555756156289e-17}*cell[5] + V{1.38777878078145e-17}*cell[6] + V{1.38777878078145e-17}*cell[7] + x87 + x90);
auto x92 = -x27*(V{2.77555756156289e-17}*cell[12] + V{1.38777878078145e-17}*cell[13] + V{1.38777878078145e-17}*cell[14] + V{1.38777878078145e-17}*cell[4] + V{1.38777878078145e-17}*cell[5] + V{2.77555756156289e-17}*cell[7] + V{2.77555756156289e-17}*cell[9] + x86 + x88 + x90 + V{1.38777878078145e-17});
auto x93 = V{3}*x38;
auto x94 = x61 + x93;
auto x95 = x77*x94;
auto x96 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x97 = V{4.5}*(x96*x96);
auto x98 = x60 + x63 + x73 + x97;
auto x99 = x35 + V{-1};
auto x100 = x37 + x50 - x93 + x99;
auto x101 = x100*x77;
auto x102 = -x42;
auto x103 = -x101 + x102 + x77*x98 + x89 + x91 + x92 + x95 + V{4.44089209850063e-16};
auto x104 = V{3}*x34;
auto x105 = x104 + x59 + x63;
auto x106 = x105*x77;
auto x107 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x108 = V{4.5}*(x107*x107);
auto x109 = x108 + x61 + x74;
auto x110 = -x104 + x40 + x48;
auto x111 = x110*x77;
auto x112 = V{3.70074341541719e-17}*x21*x30 + V{3.70074341541719e-17}*x25*x32 + V{3.70074341541719e-17}*x27*x33;
auto x113 = V{1.11022302462516e-16}*cell[12] + V{3.33066907387547e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{3.33066907387547e-16}*cell[5] + x106 + x109*x77 - x111 - x112*x34;
auto x114 = V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + x103 + x113 - x36*x85 - x38*x76 + x57 + x65 + x70 + x71*x75 - x77*x81 - x77*x84;
auto x115 = -x82;
auto x116 = x41 + x73;
auto x117 = x116 + x49;
auto x118 = x117 - V{4.5}*x115*x115;
auto x119 = -x79;
auto x120 = x116 + x66;
auto x121 = x120 - V{4.5}*x119*x119;
auto x122 = x68 + x73 - x97;
auto x123 = -x122*x44;
auto x124 = x44*x98;
auto x125 = -x44*x84;
auto x126 = x75*x77;
auto x127 = x39 - x72 + x73 + x99;
auto x128 = x127*x77;
auto x129 = V{1.66533453693773e-16}*cell[10] + V{1.66533453693773e-16}*cell[1] + x126 - x128 - x36*x76 + x64*x77;
auto x130 = V{2.22044604925031e-16}*cell[11] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{4.44089209850063e-16}*cell[17] + V{2.22044604925031e-16}*cell[18] + V{2.22044604925031e-16}*cell[2] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[8] + V{2.22044604925031e-16}*cell[9] + x102 + x113 + x123 + x124 + x125 + x129 - x38*x85 - x69*x77 + x71*x94 + x89 + x91 + x92 + V{4.44089209850063e-16};
auto x131 = -x108 + x51 + x73;
auto x132 = -x131*x44;
auto x133 = x109*x44;
auto x134 = -x44*x81;
auto x135 = V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.44089209850063e-16}*cell[13] + V{4.9960036108132e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{2.22044604925031e-16}*cell[17] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.44089209850063e-16}*cell[4] + V{4.9960036108132e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + V{2.22044604925031e-16}*cell[8] + x103 + x105*x71 - x112*x38 + x129 + x132 + x133 + x134 - x34*x85;
auto x136 = x21*(x114 - x44*x53) + x25*(-x118*x44 - x121*x77 + x130) + x27*(-x118*x77 - x121*x44 + x135 - x53*x77);
auto x137 = V{0.0555555555555556}*x21*x30 + V{0.0555555555555556}*x25*x32 + V{0.0555555555555556}*x27*x33;
auto x138 = V{3.46944695195361e-18}*cell[0];
auto x139 = V{3.46944695195361e-18}*cell[12] + V{6.93889390390723e-18}*cell[13] + V{3.46944695195361e-18}*cell[3] + x138 + V{3.46944695195361e-18};
auto x140 = V{3.46944695195361e-18}*cell[11] + V{6.93889390390723e-18}*cell[15] + V{3.46944695195361e-18}*cell[2];
auto x141 = x21*(V{6.93889390390723e-18}*cell[10] + V{6.93889390390723e-18}*cell[14] + V{6.93889390390723e-18}*cell[16] + V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x139 + x140);
auto x142 = V{3.46944695195361e-18}*cell[10] + V{6.93889390390723e-18}*cell[17] + V{3.46944695195361e-18}*cell[1];
auto x143 = x25*(V{6.93889390390723e-18}*cell[11] + V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[18] + V{6.93889390390723e-18}*cell[5] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + x139 + x142);
auto x144 = x27*(V{6.93889390390723e-18}*cell[12] + V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + V{6.93889390390723e-18}*cell[7] + V{6.93889390390723e-18}*cell[9] + x138 + x140 + x142 + V{3.46944695195361e-18});
auto x145 = V{0.0277777777777778}*x21*x30 + V{0.0277777777777778}*x25*x32 + V{0.0277777777777778}*x27*x33;
auto x146 = V{0.0833333333333333}*cell[12] - V{0.0833333333333334}*cell[13] - V{0.0833333333333334}*cell[14] + V{0.0833333333333333}*cell[3] - V{0.0833333333333334}*cell[4] - V{0.0833333333333334}*cell[5] + x141 + x143 + x144 - x145*x34 + V{0.0555555555555555};
auto x147 = V{0.0833333333333333}*cell[11] - V{0.0833333333333334}*cell[15] - V{0.0833333333333334}*cell[16] + V{0.0833333333333333}*cell[2] - V{0.0833333333333334}*cell[6] - V{0.0833333333333334}*cell[7] - x145*x38;
auto x148 = -V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] - V{0.166666666666667}*cell[1] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x137*x36 + x146 + x147;
auto x149 = V{0.0833333333333333}*cell[10] - V{0.0833333333333334}*cell[17] - V{0.0833333333333334}*cell[18] + V{0.0833333333333333}*cell[1] - V{0.0833333333333334}*cell[8] - V{0.0833333333333334}*cell[9] - x145*x36;
auto x150 = -V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] - V{0.166666666666667}*cell[2] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] + x137*x38 + x146 + x149;
auto x151 = -V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] - V{0.166666666666667}*cell[3] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] + x137*x34 + x141 + x143 + x144 + x147 + x149 + V{0.0555555555555555};
auto x152 = x21*x30 + x25*x32 + x27*x33;
auto x153 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x152;
auto x154 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x153;
auto x155 = V{0.0138888888888889}*x21*x30 + V{0.0138888888888889}*x25*x32 + V{0.0138888888888889}*x27*x33;
auto x156 = V{0.013888875}*cell[0];
auto x157 = V{0.013888875}*cell[12] + V{0.02777775}*cell[13] + V{0.013888875}*cell[3] + x156 + V{0.013888875};
auto x158 = V{0.013888875}*cell[11] + V{0.02777775}*cell[15] + V{0.013888875}*cell[2];
auto x159 = -x21*(V{0.02777775}*cell[10] + V{0.02777775}*cell[14] + V{0.02777775}*cell[16] + V{0.013888875}*cell[17] + V{0.013888875}*cell[18] + V{0.013888875}*cell[8] + V{0.013888875}*cell[9] + x157 + x158);
auto x160 = V{0.013888875}*cell[10] + V{0.02777775}*cell[17] + V{0.013888875}*cell[1];
auto x161 = -x25*(V{0.02777775}*cell[11] + V{0.013888875}*cell[15] + V{0.013888875}*cell[16] + V{0.02777775}*cell[18] + V{0.02777775}*cell[5] + V{0.013888875}*cell[6] + V{0.013888875}*cell[7] + x157 + x160);
auto x162 = -x27*(V{0.02777775}*cell[12] + V{0.013888875}*cell[13] + V{0.013888875}*cell[14] + V{0.013888875}*cell[4] + V{0.013888875}*cell[5] + V{0.02777775}*cell[7] + V{0.02777775}*cell[9] + x156 + x158 + x160 + V{0.013888875});
auto x163 = V{0.0277777777777778}*x21*x30 + V{0.0277777777777778}*x25*x32 + V{0.0277777777777778}*x27*x33;
auto x164 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x159 + x161 + x162 - x163*x36 + V{0.0138888472222222};
auto x165 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] - x163*x38;
auto x166 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x155*x34 + x164 + x165;
auto x167 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] - x154 + x166;
auto x168 = x120 + x80;
auto x169 = -x168*x44;
auto x170 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] + x154 + x166;
auto x171 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x153;
auto x172 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x163*x34;
auto x173 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x155*x38 + x164 + x172;
auto x174 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] - x171 + x173;
auto x175 = x117 + x83;
auto x176 = -x175*x44;
auto x177 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] + x171 + x173;
auto x178 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x152;
auto x179 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x155*x36 + x159 + x161 + x162 + x165 + x172 + V{0.0138888472222222};
auto x180 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x178 + x179;
auto x181 = x52 + x67;
auto x182 = -x181*x44;
auto x183 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x178 + x179;
auto x184 = x21*(x114 + x182) + x25*(x130 - x168*x77 + x176) + x27*(x135 + x169 - x175*x77 - x181*x77);
auto x0 = -x19*(V{0.111111111111111}*x136*x41 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] - x21*(V{0.333333}*cell[10] + V{0.333333}*cell[14] + V{0.333333}*cell[16] + V{0.1666665}*cell[17] + V{0.1666665}*cell[18] + V{0.1666665}*cell[8] + V{0.1666665}*cell[9] + x23 + x24) - x25*(V{0.333333}*cell[11] + V{0.1666665}*cell[15] + V{0.1666665}*cell[16] + V{0.333333}*cell[18] + V{0.333333}*cell[5] + V{0.1666665}*cell[6] + V{0.1666665}*cell[7] + x23 + x26) - x27*(V{0.333333}*cell[12] + V{0.1666665}*cell[13] + V{0.1666665}*cell[14] + V{0.1666665}*cell[4] + V{0.1666665}*cell[5] + V{0.333333}*cell[7] + V{0.333333}*cell[9] + x22 + x24 + x26 + V{0.1666665}) - x34*x43 - x36*x43 - x38*x43 + x42 + V{0.833332833333333});
auto x1 = -x19*(V{0.0185185185185185}*x127*x136 + V{0.0555555555555556}) + x20*(x128 + x148);
auto x2 = -x19*(V{0.0185185185185185}*x100*x136 + V{0.0555555555555556}) + x20*(x101 + x150);
auto x3 = -x19*(V{0.0185185185185185}*x110*x136 + V{0.0555555555555556}) + x20*(x111 + x151);
auto x4 = -x19*(V{0.00925925925925926}*x131*x136 + V{0.0277777777777778}) - x20*(x132 + x167);
auto x5 = -x19*(V{0.00925925925925926}*x121*x136 + V{0.0277777777777778}) - x20*(x169 + x170);
auto x6 = -x19*(V{0.00925925925925926}*x122*x136 + V{0.0277777777777778}) - x20*(x123 + x174);
auto x7 = -x19*(V{0.00925925925925926}*x118*x136 + V{0.0277777777777778}) - x20*(x176 + x177);
auto x8 = -x19*(V{0.00925925925925926}*x136*x56 + V{0.0277777777777778}) - x20*(x180 + x57);
auto x9 = -x19*(V{0.00925925925925926}*x136*x53 + V{0.0277777777777778}) - x20*(x182 + x183);
auto x10 = x19*(V{0.0185185185185185}*x184*x75 + V{-0.0555555555555556}) + x20*(-x126 + x148);
auto x11 = x19*(V{0.0185185185185185}*x184*x94 + V{-0.0555555555555556}) + x20*(x150 - x95);
auto x12 = x19*(V{0.0185185185185185}*x105*x184 + V{-0.0555555555555556}) + x20*(-x106 + x151);
auto x13 = x19*(V{0.00925925925925926}*x109*x184 + V{-0.0277777777777778}) - x20*(x133 + x167);
auto x14 = -x19*(V{0.00925925925925926}*x136*x81 + V{0.0277777777777778}) - x20*(x134 + x170);
auto x15 = x19*(V{0.00925925925925926}*x184*x98 + V{-0.0277777777777778}) - x20*(x124 + x174);
auto x16 = -x19*(V{0.00925925925925926}*x136*x84 + V{0.0277777777777778}) - x20*(x125 + x177);
auto x17 = x19*(V{0.00925925925925926}*x184*x64 + V{-0.0277777777777778}) - x20*(x180 + x65);
auto x18 = -x19*(V{0.00925925925925926}*x136*x69 + V{0.0277777777777778}) - x20*(x183 + x70);
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
return { V{0.333333333333333}*x184, x34 + x36 + x38 };
}
};

}

}
