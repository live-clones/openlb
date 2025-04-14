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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerCornerDensity3D<-1, -1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerCornerStress3D<-1, -1, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{-1};
auto x22 = V{1} / (x21);
auto x23 = V{0.1666665}*cell[0];
auto x24 = V{0.1666665}*cell[12] + V{0.1666665}*cell[3] + V{0.333333}*cell[4] + x23 + V{0.1666665};
auto x25 = V{0.1666665}*cell[11] + V{0.1666665}*cell[2] + V{0.333333}*cell[6];
auto x26 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{-1};
auto x27 = V{1} / (x26);
auto x28 = V{0.1666665}*cell[10] + V{0.1666665}*cell[1] + V{0.333333}*cell[8];
auto x29 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{-1};
auto x30 = V{1} / (x29);
auto x31 = cell[0] + cell[12] + cell[3] + V{2}*cell[4] + V{1};
auto x32 = cell[11] + cell[2] + V{2}*cell[6];
auto x33 = cell[17] + cell[18] + V{2}*cell[1] + V{2}*cell[5] + V{2}*cell[7] + cell[8] + cell[9] + x31 + x32;
auto x34 = cell[10] + cell[1] + V{2}*cell[8];
auto x35 = V{2}*cell[14] + cell[15] + cell[16] + V{2}*cell[2] + cell[6] + cell[7] + V{2}*cell[9] + x31 + x34;
auto x36 = cell[0] + cell[13] + cell[14] + V{2}*cell[16] + V{2}*cell[18] + V{2}*cell[3] + cell[4] + cell[5] + x32 + x34 + V{1};
auto x37 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x38 = V{1.5}*x37;
auto x39 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x40 = V{1.5}*x39;
auto x41 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x42 = V{1.5}*x41;
auto x43 = x40 + x42 + V{-1};
auto x44 = x38 + x43;
auto x45 = x44*(V{0.111111111111111}*x22*x33 + V{0.111111111111111}*x27*x35 + V{0.111111111111111}*x30*x36);
auto x46 = -1/x21;
auto x47 = -1/x26;
auto x48 = -1/x29;
auto x49 = V{0.166666666666667}*x33*x46 + V{0.166666666666667}*x35*x47 + V{0.166666666666667}*x36*x48;
auto x50 = V{0.00925925925925926}*x22*x33 + V{0.00925925925925926}*x27*x35 + V{0.00925925925925926}*x30*x36;
auto x51 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x52 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x53 = V{4.5}*(x52*x52);
auto x54 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x55 = x44 + x54;
auto x56 = x51 - x53 + x55;
auto x57 = x50*x56;
auto x58 = -x38;
auto x59 = V{1} - x40;
auto x60 = x58 + x59;
auto x61 = x54 + x60;
auto x62 = -x42;
auto x63 = x51 + x62;
auto x64 = x53 + x61 + x63;
auto x65 = -x50*x64;
auto x66 = -x54;
auto x67 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x68 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x67;
auto x69 = -V{4.5}*x68*x68;
auto x70 = x44 + x51;
auto x71 = x66 + x69 + x70;
auto x72 = x50*x71;
auto x73 = -x68;
auto x74 = -x51;
auto x75 = x55 + x74;
auto x76 = x75 - V{4.5}*x73*x73;
auto x77 = V{0.0185185185185185}*x22*x33 + V{0.0185185185185185}*x27*x35 + V{0.0185185185185185}*x30*x36;
auto x78 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x79 = -x78;
auto x80 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x81 = x44 + x80;
auto x82 = x66 + x81;
auto x83 = x82 - V{4.5}*x79*x79;
auto x84 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x67;
auto x85 = -x84;
auto x86 = x74 + x81;
auto x87 = x86 - V{4.5}*x85*x85;
auto x88 = V{0.037037037037037}*x22*x33 + V{0.037037037037037}*x27*x35 + V{0.037037037037037}*x30*x36;
auto x89 = V{3}*x39;
auto x90 = x38 + V{-1};
auto x91 = x42 + x80 - x89 + x90;
auto x92 = V{5.55111512312578e-17}*x33*x46 + V{5.55111512312578e-17}*x35*x47 + V{5.55111512312578e-17}*x36*x48;
auto x93 = V{3.70074341541719e-17}*x33*x46 + V{3.70074341541719e-17}*x35*x47 + V{3.70074341541719e-17}*x36*x48;
auto x94 = V{7.40148683083438e-17}*x33*x46 + V{7.40148683083438e-17}*x35*x47 + V{7.40148683083438e-17}*x36*x48;
auto x95 = V{3}*x41;
auto x96 = x40 + x54 + x90 - x95;
auto x97 = x77*x96;
auto x98 = x61 + x95;
auto x99 = -x77*x98;
auto x100 = V{1.38777878078145e-17}*cell[0];
auto x101 = V{1.38777878078145e-17}*cell[12] + V{1.38777878078145e-17}*cell[3] + V{2.77555756156289e-17}*cell[4] + x100 + V{1.38777878078145e-17};
auto x102 = V{1.38777878078145e-17}*cell[11] + V{1.38777878078145e-17}*cell[2] + V{2.77555756156289e-17}*cell[6];
auto x103 = x22*(V{1.38777878078145e-17}*cell[17] + V{1.38777878078145e-17}*cell[18] + V{2.77555756156289e-17}*cell[1] + V{2.77555756156289e-17}*cell[5] + V{2.77555756156289e-17}*cell[7] + V{1.38777878078145e-17}*cell[8] + V{1.38777878078145e-17}*cell[9] + x101 + x102);
auto x104 = V{1.38777878078145e-17}*cell[10] + V{1.38777878078145e-17}*cell[1] + V{2.77555756156289e-17}*cell[8];
auto x105 = x27*(V{2.77555756156289e-17}*cell[14] + V{1.38777878078145e-17}*cell[15] + V{1.38777878078145e-17}*cell[16] + V{2.77555756156289e-17}*cell[2] + V{1.38777878078145e-17}*cell[6] + V{1.38777878078145e-17}*cell[7] + V{2.77555756156289e-17}*cell[9] + x101 + x104);
auto x106 = x30*(V{1.38777878078145e-17}*cell[13] + V{1.38777878078145e-17}*cell[14] + V{2.77555756156289e-17}*cell[16] + V{2.77555756156289e-17}*cell[18] + V{2.77555756156289e-17}*cell[3] + V{1.38777878078145e-17}*cell[4] + V{1.38777878078145e-17}*cell[5] + x100 + x102 + x104 + V{1.38777878078145e-17});
auto x107 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x108 = V{4.5}*(x107*x107);
auto x109 = -x108 + x70 + x80;
auto x110 = x103 + x105 + x106 + x109*x77 + x45 + x97 + x99 + V{4.44089209850063e-16};
auto x111 = V{3}*x37;
auto x112 = -x111 + x43 + x51;
auto x113 = x112*x77;
auto x114 = x111 + x59 + x63;
auto x115 = -x114*x77;
auto x116 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x117 = V{4.5}*(x116*x116);
auto x118 = -x117 + x55 + x80;
auto x119 = V{3.33066907387547e-16}*cell[14] + V{3.33066907387547e-16}*cell[5] + x113 + x115 + x118*x77;
auto x120 = x109*x50;
auto x121 = x108 + x60 + x63 + x80;
auto x122 = -x121*x50;
auto x123 = -x80;
auto x124 = -V{4.5}*x84*x84;
auto x125 = x123 + x124 + x70;
auto x126 = x125*x50;
auto x127 = -V{4.5}*x78*x78;
auto x128 = x123 + x127 + x55;
auto x129 = x77*x91;
auto x130 = x62 + x80;
auto x131 = x130 + x58 + x89 + V{1};
auto x132 = -x131*x77;
auto x133 = V{1.66533453693773e-16}*cell[10] + V{2.22044604925031e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{2.22044604925031e-16}*cell[9] + x129 + x132 - x39*x92 + x56*x77;
auto x134 = x118*x50;
auto x135 = x117 + x130 + x61;
auto x136 = -x135*x50;
auto x137 = x128*x50;
auto x138 = x22*(V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{3.33066907387547e-16}*cell[13] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + x110 + x119 - x37*x93 - x39*x94 - x41*x92 + x50*x76 + x57 + x65 + x72 + x77*x83 + x77*x87 + x88*x91) + x27*(V{2.22044604925031e-16}*cell[11] + V{9.71445146547012e-17}*cell[12] + V{1.11022302462516e-16}*cell[13] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{2.22044604925031e-16}*cell[17] + V{2.22044604925031e-16}*cell[2] + V{9.71445146547012e-17}*cell[3] + V{1.11022302462516e-16}*cell[4] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{2.22044604925031e-16}*cell[8] + x103 + x105 + x106 + x119 + x120 + x122 + x126 + x128*x77 + x133 - x37*(V{3.23815048849004e-17}*x33*x46 + V{3.23815048849004e-17}*x35*x47 + V{3.23815048849004e-17}*x36*x48) - x41*x94 + x45 + x50*x87 + x76*x77 + x88*x96 + V{4.44089209850063e-16}) + x30*(V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.71844785465692e-16}*cell[13] + V{5.27355936696949e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.71844785465692e-16}*cell[4] + V{5.27355936696949e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + x110 + x112*x88 + x125*x77 + x133 + x134 + x136 + x137 - x37*x94 - x41*x93 + x50*x83 + x71*x77);
auto x139 = x33*x46 + x35*x47 + x36*x48;
auto x140 = V{3.46944695195361e-18}*cell[0];
auto x141 = V{3.46944695195361e-18}*cell[12] + V{3.46944695195361e-18}*cell[3] + V{6.93889390390723e-18}*cell[4] + x140 + V{3.46944695195361e-18};
auto x142 = V{3.46944695195361e-18}*cell[11] + V{3.46944695195361e-18}*cell[2] + V{6.93889390390723e-18}*cell[6];
auto x143 = x22*(V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{6.93889390390723e-18}*cell[1] + V{6.93889390390723e-18}*cell[5] + V{6.93889390390723e-18}*cell[7] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x141 + x142);
auto x144 = V{3.46944695195361e-18}*cell[10] + V{3.46944695195361e-18}*cell[1] + V{6.93889390390723e-18}*cell[8];
auto x145 = x27*(V{6.93889390390723e-18}*cell[14] + V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[2] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + V{6.93889390390723e-18}*cell[9] + x141 + x144);
auto x146 = x30*(V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{6.93889390390723e-18}*cell[16] + V{6.93889390390723e-18}*cell[18] + V{6.93889390390723e-18}*cell[3] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + x140 + x142 + x144 + V{3.46944695195361e-18});
auto x147 = V{0.0277777777777778}*x33*x46 + V{0.0277777777777778}*x35*x47 + V{0.0277777777777778}*x36*x48;
auto x148 = -V{0.0833333333333333}*cell[12] + V{0.0833333333333334}*cell[13] + V{0.0833333333333334}*cell[14] - V{0.0833333333333333}*cell[3] + V{0.0833333333333334}*cell[4] + V{0.0833333333333334}*cell[5] + x143 + x145 + x146 + x147*x37 + V{-0.0555555555555555};
auto x149 = -V{0.0833333333333333}*cell[11] + V{0.0833333333333334}*cell[15] + V{0.0833333333333334}*cell[16] - V{0.0833333333333333}*cell[2] + V{0.0833333333333334}*cell[6] + V{0.0833333333333334}*cell[7] + x147*x41;
auto x150 = V{0.166666666666667}*cell[10] - V{0.166666666666667}*cell[17] - V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[1] - V{0.166666666666667}*cell[8] - V{0.166666666666667}*cell[9] - V{0.0555555555555556}*x139*x39 + x148 + x149;
auto x151 = -V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] - V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] + x147*x39;
auto x152 = V{0.166666666666667}*cell[11] - V{0.166666666666667}*cell[15] - V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[2] - V{0.166666666666667}*cell[6] - V{0.166666666666667}*cell[7] - V{0.0555555555555556}*x139*x41 + x148 + x151;
auto x153 = V{0.166666666666667}*cell[12] - V{0.166666666666667}*cell[13] - V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[3] - V{0.166666666666667}*cell[4] - V{0.166666666666667}*cell[5] - V{0.0555555555555556}*x139*x37 + x143 + x145 + x146 + x149 + x151 + V{-0.0555555555555555};
auto x154 = x22*x33 + x27*x35 + x30*x36;
auto x155 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x154;
auto x156 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x155;
auto x157 = V{0.0138888888888889}*x22*x33 + V{0.0138888888888889}*x27*x35 + V{0.0138888888888889}*x30*x36;
auto x158 = V{0.013888875}*cell[0];
auto x159 = V{0.013888875}*cell[12] + V{0.013888875}*cell[3] + V{0.02777775}*cell[4] + x158 + V{0.013888875};
auto x160 = V{0.013888875}*cell[11] + V{0.013888875}*cell[2] + V{0.02777775}*cell[6];
auto x161 = x22*(V{0.013888875}*cell[17] + V{0.013888875}*cell[18] + V{0.02777775}*cell[1] + V{0.02777775}*cell[5] + V{0.02777775}*cell[7] + V{0.013888875}*cell[8] + V{0.013888875}*cell[9] + x159 + x160);
auto x162 = V{0.013888875}*cell[10] + V{0.013888875}*cell[1] + V{0.02777775}*cell[8];
auto x163 = x27*(V{0.02777775}*cell[14] + V{0.013888875}*cell[15] + V{0.013888875}*cell[16] + V{0.02777775}*cell[2] + V{0.013888875}*cell[6] + V{0.013888875}*cell[7] + V{0.02777775}*cell[9] + x159 + x162);
auto x164 = x30*(V{0.013888875}*cell[13] + V{0.013888875}*cell[14] + V{0.02777775}*cell[16] + V{0.02777775}*cell[18] + V{0.02777775}*cell[3] + V{0.013888875}*cell[4] + V{0.013888875}*cell[5] + x158 + x160 + x162 + V{0.013888875});
auto x165 = V{0.0277777777777778}*x22*x33 + V{0.0277777777777778}*x27*x35 + V{0.0277777777777778}*x30*x36;
auto x166 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x161 + x163 + x164 + x165*x39 + V{0.0138888472222222};
auto x167 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] + x165*x41;
auto x168 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] - x157*x37 + x166 + x167;
auto x169 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] + x156 + x168;
auto x170 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] - x156 + x168;
auto x171 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x155;
auto x172 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] + x165*x37;
auto x173 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] - x157*x41 + x166 + x172;
auto x174 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] + x171 + x173;
auto x175 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] - x171 + x173;
auto x176 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x154;
auto x177 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] - x157*x39 + x161 + x163 + x164 + x167 + x172 + V{0.0138888472222222};
auto x178 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] + x176 + x177;
auto x179 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] - x176 + x177;
auto x0 = -x19*(-V{0.111111111111111}*x138*x44 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x22*(V{0.1666665}*cell[17] + V{0.1666665}*cell[18] + V{0.333333}*cell[1] + V{0.333333}*cell[5] + V{0.333333}*cell[7] + V{0.1666665}*cell[8] + V{0.1666665}*cell[9] + x24 + x25) + x27*(V{0.333333}*cell[14] + V{0.1666665}*cell[15] + V{0.1666665}*cell[16] + V{0.333333}*cell[2] + V{0.1666665}*cell[6] + V{0.1666665}*cell[7] + V{0.333333}*cell[9] + x24 + x28) + x30*(V{0.1666665}*cell[13] + V{0.1666665}*cell[14] + V{0.333333}*cell[16] + V{0.333333}*cell[18] + V{0.333333}*cell[3] + V{0.1666665}*cell[4] + V{0.1666665}*cell[5] + x23 + x25 + x28 + V{0.1666665}) - x37*x49 - x39*x49 - x41*x49 - x45 + V{0.833332833333333});
auto x1 = -x19*(-V{0.0185185185185185}*x138*x91 + V{0.0555555555555556}) + x20*(-x129 - x150);
auto x2 = -x19*(-V{0.0185185185185185}*x138*x96 + V{0.0555555555555556}) + x20*(-x152 - x97);
auto x3 = -x19*(-V{0.0185185185185185}*x112*x138 + V{0.0555555555555556}) + x20*(-x113 - x153);
auto x4 = -x19*(-V{0.00925925925925926}*x118*x138 + V{0.0277777777777778}) - x20*(x134 + x169);
auto x5 = -x19*(-V{0.00925925925925926}*x138*x83 + V{0.0277777777777778}) - x20*(x170 + x50*(x127 + x82));
auto x6 = -x19*(-V{0.00925925925925926}*x109*x138 + V{0.0277777777777778}) - x20*(x120 + x174);
auto x7 = -x19*(-V{0.00925925925925926}*x138*x87 + V{0.0277777777777778}) - x20*(x175 + x50*(x124 + x86));
auto x8 = -x19*(-V{0.00925925925925926}*x138*x56 + V{0.0277777777777778}) - x20*(x178 + x57);
auto x9 = -x19*(-V{0.00925925925925926}*x138*x76 + V{0.0277777777777778}) - x20*(x179 + x50*(x69 + x75));
auto x10 = -x19*(V{0.0185185185185185}*x131*x138 + V{0.0555555555555556}) + x20*(-x132 - x150);
auto x11 = -x19*(V{0.0185185185185185}*x138*x98 + V{0.0555555555555556}) + x20*(-x152 - x99);
auto x12 = -x19*(V{0.0185185185185185}*x114*x138 + V{0.0555555555555556}) + x20*(-x115 - x153);
auto x13 = -x19*(V{0.00925925925925926}*x135*x138 + V{0.0277777777777778}) - x20*(x136 + x169);
auto x14 = -x19*(-V{0.00925925925925926}*x128*x138 + V{0.0277777777777778}) - x20*(x137 + x170);
auto x15 = -x19*(V{0.00925925925925926}*x121*x138 + V{0.0277777777777778}) - x20*(x122 + x174);
auto x16 = -x19*(-V{0.00925925925925926}*x125*x138 + V{0.0277777777777778}) - x20*(x126 + x175);
auto x17 = -x19*(V{0.00925925925925926}*x138*x64 + V{0.0277777777777778}) - x20*(x178 + x65);
auto x18 = -x19*(-V{0.00925925925925926}*x138*x71 + V{0.0277777777777778}) - x20*(x179 + x72);
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
return { -V{0.333333333333333}*x138, x37 + x39 + x41 };
}
};

}

}
