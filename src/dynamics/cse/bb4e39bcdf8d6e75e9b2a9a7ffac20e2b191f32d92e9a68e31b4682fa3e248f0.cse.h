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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerCornerDensity3D<1, -1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerCornerStress3D<1, -1, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{-1};
auto x22 = V{1} / (x21);
auto x23 = V{0.1666665}*cell[0];
auto x24 = V{0.1666665}*cell[10] + V{0.1666665}*cell[1] + V{0.333333}*cell[8] + x23 + V{0.1666665};
auto x25 = V{0.1666665}*cell[12] + V{0.333333}*cell[14] + V{0.1666665}*cell[3];
auto x26 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{-1};
auto x27 = V{1} / (x26);
auto x28 = V{0.1666665}*cell[11] + V{0.333333}*cell[16] + V{0.1666665}*cell[2];
auto x29 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{1});
auto x30 = cell[0] + cell[12] + V{2}*cell[14] + cell[3] + V{1};
auto x31 = cell[11] + V{2}*cell[16] + cell[2];
auto x32 = V{2}*cell[10] + V{2}*cell[13] + V{2}*cell[15] + cell[17] + cell[18] + cell[8] + cell[9] + x30 + x31;
auto x33 = cell[10] + cell[1] + V{2}*cell[8];
auto x34 = cell[15] + cell[16] + V{2}*cell[2] + V{2}*cell[4] + cell[6] + cell[7] + V{2}*cell[9] + x30 + x33;
auto x35 = cell[0] + cell[13] + cell[14] + V{2}*cell[18] + V{2}*cell[3] + cell[4] + cell[5] + V{2}*cell[6] + x31 + x33 + V{1};
auto x36 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x37 = V{1.5}*x36;
auto x38 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x39 = V{1.5}*x38;
auto x40 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x41 = V{1.5}*x40;
auto x42 = x39 + x41 + V{-1};
auto x43 = x37 + x42;
auto x44 = x43*(-V{0.111111111111111}*x22*x34 - V{0.111111111111111}*x27*x35 + V{0.111111111111111}*x29*x32);
auto x45 = -1/x21;
auto x46 = -1/x26;
auto x47 = V{0.166666666666667}*x29*x32 + V{0.166666666666667}*x34*x45 + V{0.166666666666667}*x35*x46;
auto x48 = V{1.38777878078145e-17}*cell[0];
auto x49 = V{1.38777878078145e-17}*cell[10] + V{1.38777878078145e-17}*cell[1] + V{2.77555756156289e-17}*cell[8] + x48 + V{1.38777878078145e-17};
auto x50 = V{1.38777878078145e-17}*cell[12] + V{2.77555756156289e-17}*cell[14] + V{1.38777878078145e-17}*cell[3];
auto x51 = x22*(V{1.38777878078145e-17}*cell[15] + V{1.38777878078145e-17}*cell[16] + V{2.77555756156289e-17}*cell[2] + V{2.77555756156289e-17}*cell[4] + V{1.38777878078145e-17}*cell[6] + V{1.38777878078145e-17}*cell[7] + V{2.77555756156289e-17}*cell[9] + x49 + x50);
auto x52 = V{1.38777878078145e-17}*cell[11] + V{2.77555756156289e-17}*cell[16] + V{1.38777878078145e-17}*cell[2];
auto x53 = x27*(V{1.38777878078145e-17}*cell[13] + V{1.38777878078145e-17}*cell[14] + V{2.77555756156289e-17}*cell[18] + V{2.77555756156289e-17}*cell[3] + V{1.38777878078145e-17}*cell[4] + V{1.38777878078145e-17}*cell[5] + V{2.77555756156289e-17}*cell[6] + x49 + x52);
auto x54 = -x29*(V{2.77555756156289e-17}*cell[10] + V{2.77555756156289e-17}*cell[13] + V{2.77555756156289e-17}*cell[15] + V{1.38777878078145e-17}*cell[17] + V{1.38777878078145e-17}*cell[18] + V{1.38777878078145e-17}*cell[8] + V{1.38777878078145e-17}*cell[9] + x48 + x50 + x52 + V{1.38777878078145e-17});
auto x55 = -x44;
auto x56 = -V{0.00925925925925926}*x22*x34 - V{0.00925925925925926}*x27*x35 + V{0.00925925925925926}*x29*x32;
auto x57 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x58 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x59 = V{4.5}*(x58*x58);
auto x60 = -x37;
auto x61 = V{1} - x39;
auto x62 = x60 + x61;
auto x63 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x64 = -x41;
auto x65 = x63 + x64;
auto x66 = x57 + x59 + x62 + x65;
auto x67 = x43 + x63;
auto x68 = x57 - x59 + x67;
auto x69 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x70 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x69;
auto x71 = -x70;
auto x72 = -x63;
auto x73 = x43 + x57;
auto x74 = x72 + x73;
auto x75 = x74 - V{4.5}*x71*x71;
auto x76 = -x57;
auto x77 = -V{4.5}*x70*x70;
auto x78 = x67 + x76 + x77;
auto x79 = -V{0.0185185185185185}*x22*x34 - V{0.0185185185185185}*x27*x35 + V{0.0185185185185185}*x29*x32;
auto x80 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x81 = V{4.5}*(x80*x80);
auto x82 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x83 = x43 + x82;
auto x84 = x57 - x81 + x83;
auto x85 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x69;
auto x86 = -x85;
auto x87 = x72 + x83;
auto x88 = x87 - V{4.5}*x86*x86;
auto x89 = -V{0.037037037037037}*x22*x34 - V{0.037037037037037}*x27*x35 + V{0.037037037037037}*x29*x32;
auto x90 = V{3}*x40;
auto x91 = x37 + V{-1};
auto x92 = x39 + x82 - x90 + x91;
auto x93 = V{7.40148683083438e-17}*x29*x32 + V{7.40148683083438e-17}*x34*x45 + V{7.40148683083438e-17}*x35*x46;
auto x94 = V{3}*x38;
auto x95 = x41 + x57 + x91 - x94;
auto x96 = -x79*x95;
auto x97 = x57 + x64;
auto x98 = x60 + x94 + x97 + V{1};
auto x99 = x79*x98;
auto x100 = V{5.55111512312578e-17}*x29*x32 + V{5.55111512312578e-17}*x34*x45 + V{5.55111512312578e-17}*x35*x46;
auto x101 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x102 = V{4.5}*(x101*x101);
auto x103 = -x102 + x63 + x83;
auto x104 = V{1.66533453693773e-16}*cell[10] + V{2.22044604925031e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{2.22044604925031e-16}*cell[9] - x100*x38 - x103*x79 + x96 + x99;
auto x105 = V{3}*x36;
auto x106 = -x105 + x42 + x63;
auto x107 = -x106*x79;
auto x108 = x105 + x61 + x65;
auto x109 = x108*x79;
auto x110 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x111 = -V{4.5}*x110*x110;
auto x112 = x111 + x76 + x83;
auto x113 = V{3.33066907387547e-16}*cell[14] + V{3.33066907387547e-16}*cell[5] + x107 + x109 - x112*x79;
auto x114 = x62 + x82;
auto x115 = x114 + x81 + x97;
auto x116 = -x110;
auto x117 = -x82;
auto x118 = x117 + x73;
auto x119 = x118 - V{4.5}*x116*x116;
auto x120 = -V{4.5}*x85*x85;
auto x121 = x117 + x120 + x67;
auto x122 = V{3.70074341541719e-17}*x29*x32 + V{3.70074341541719e-17}*x34*x45 + V{3.70074341541719e-17}*x35*x46;
auto x123 = -x79*x92;
auto x124 = x114 + x90;
auto x125 = x124*x79;
auto x126 = x123 + x125 + x51 + x53 + x54 + x55 - x78*x79 + V{4.44089209850063e-16};
auto x127 = x102 + x114 + x65;
auto x128 = x22*(V{2.22044604925031e-16}*cell[11] + V{9.71445146547012e-17}*cell[12] + V{1.11022302462516e-16}*cell[13] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{2.22044604925031e-16}*cell[17] + V{2.22044604925031e-16}*cell[2] + V{9.71445146547012e-17}*cell[3] + V{1.11022302462516e-16}*cell[4] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{2.22044604925031e-16}*cell[8] + x104 + x113 - x36*(V{3.23815048849004e-17}*x29*x32 + V{3.23815048849004e-17}*x34*x45 + V{3.23815048849004e-17}*x35*x46) - x40*x93 + x51 + x53 + x54 + x55 + x56*x66 - x56*x68 - x56*x75 - x56*x78 - x79*x84 - x79*x88 - x89*x92 + V{4.44089209850063e-16}) + x27*(V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.71844785465692e-16}*cell[13] + V{5.27355936696949e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.71844785465692e-16}*cell[4] + V{5.27355936696949e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + x104 - x106*x89 - x112*x56 + x115*x56 - x119*x56 - x121*x79 - x122*x40 + x126 - x36*x93 - x56*x84 - x68*x79) - x29*(V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{3.33066907387547e-16}*cell[13] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] - x100*x40 - x103*x56 + x113 + x115*x79 - x121*x56 - x122*x36 + x126 + x127*x56 - x38*x93 - x56*x88 + x66*x79 + x89*x98);
auto x129 = -x128;
auto x130 = x29*x32 + x34*x45 + x35*x46;
auto x131 = V{3.46944695195361e-18}*cell[0];
auto x132 = V{3.46944695195361e-18}*cell[12] + V{6.93889390390723e-18}*cell[14] + V{3.46944695195361e-18}*cell[3] + x131 + V{3.46944695195361e-18};
auto x133 = V{3.46944695195361e-18}*cell[10] + V{3.46944695195361e-18}*cell[1] + V{6.93889390390723e-18}*cell[8];
auto x134 = x22*(V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[2] + V{6.93889390390723e-18}*cell[4] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + V{6.93889390390723e-18}*cell[9] + x132 + x133);
auto x135 = V{3.46944695195361e-18}*cell[11] + V{6.93889390390723e-18}*cell[16] + V{3.46944695195361e-18}*cell[2];
auto x136 = x27*(V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{6.93889390390723e-18}*cell[18] + V{6.93889390390723e-18}*cell[3] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + V{6.93889390390723e-18}*cell[6] + x131 + x133 + x135 + V{3.46944695195361e-18});
auto x137 = -x29*(V{6.93889390390723e-18}*cell[10] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[15] + V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x132 + x135);
auto x138 = V{0.0277777777777778}*x29*x32 + V{0.0277777777777778}*x34*x45 + V{0.0277777777777778}*x35*x46;
auto x139 = -V{0.0833333333333333}*cell[12] + V{0.0833333333333334}*cell[13] + V{0.0833333333333334}*cell[14] - V{0.0833333333333333}*cell[3] + V{0.0833333333333334}*cell[4] + V{0.0833333333333334}*cell[5] + x134 + x136 + x137 + x138*x36 + V{-0.0555555555555555};
auto x140 = -V{0.0833333333333333}*cell[11] + V{0.0833333333333334}*cell[15] + V{0.0833333333333334}*cell[16] - V{0.0833333333333333}*cell[2] + V{0.0833333333333334}*cell[6] + V{0.0833333333333334}*cell[7] + x138*x40;
auto x141 = V{0.166666666666667}*cell[10] - V{0.166666666666667}*cell[17] - V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[1] - V{0.166666666666667}*cell[8] - V{0.166666666666667}*cell[9] - V{0.0555555555555556}*x130*x38 + x139 + x140;
auto x142 = -V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] - V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] + x138*x38;
auto x143 = V{0.166666666666667}*cell[11] - V{0.166666666666667}*cell[15] - V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[2] - V{0.166666666666667}*cell[6] - V{0.166666666666667}*cell[7] - V{0.0555555555555556}*x130*x40 + x139 + x142;
auto x144 = V{0.166666666666667}*cell[12] - V{0.166666666666667}*cell[13] - V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[3] - V{0.166666666666667}*cell[4] - V{0.166666666666667}*cell[5] - V{0.0555555555555556}*x130*x36 + x134 + x136 + x137 + x140 + x142 + V{-0.0555555555555555};
auto x145 = V{0.00925925925925926}*x22*x34 + V{0.00925925925925926}*x27*x35 - V{0.00925925925925926}*x29*x32;
auto x146 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x147 = x22*x34 + x27*x35 - x29*x32;
auto x148 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x147;
auto x149 = x146*x148;
auto x150 = -V{0.0138888888888889}*x22*x34 - V{0.0138888888888889}*x27*x35 + V{0.0138888888888889}*x29*x32;
auto x151 = V{0.013888875}*cell[0];
auto x152 = V{0.013888875}*cell[10] + V{0.013888875}*cell[1] + V{0.02777775}*cell[8] + x151 + V{0.013888875};
auto x153 = V{0.013888875}*cell[12] + V{0.02777775}*cell[14] + V{0.013888875}*cell[3];
auto x154 = x22*(V{0.013888875}*cell[15] + V{0.013888875}*cell[16] + V{0.02777775}*cell[2] + V{0.02777775}*cell[4] + V{0.013888875}*cell[6] + V{0.013888875}*cell[7] + V{0.02777775}*cell[9] + x152 + x153);
auto x155 = V{0.013888875}*cell[11] + V{0.02777775}*cell[16] + V{0.013888875}*cell[2];
auto x156 = x27*(V{0.013888875}*cell[13] + V{0.013888875}*cell[14] + V{0.02777775}*cell[18] + V{0.02777775}*cell[3] + V{0.013888875}*cell[4] + V{0.013888875}*cell[5] + V{0.02777775}*cell[6] + x152 + x155);
auto x157 = -x29*(V{0.02777775}*cell[10] + V{0.02777775}*cell[13] + V{0.02777775}*cell[15] + V{0.013888875}*cell[17] + V{0.013888875}*cell[18] + V{0.013888875}*cell[8] + V{0.013888875}*cell[9] + x151 + x153 + x155 + V{0.013888875});
auto x158 = -V{0.0277777777777778}*x22*x34 - V{0.0277777777777778}*x27*x35 + V{0.0277777777777778}*x29*x32;
auto x159 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x154 + x156 + x157 - x158*x38 + V{0.0138888472222222};
auto x160 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] - x158*x40;
auto x161 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x150*x36 + x159 + x160;
auto x162 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] + x149 + x161;
auto x163 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] - x149 + x161;
auto x164 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x148;
auto x165 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x158*x36;
auto x166 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x150*x40 + x159 + x165;
auto x167 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] + x164 + x166;
auto x168 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] - x164 + x166;
auto x169 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x146*x147;
auto x170 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x150*x38 + x154 + x156 + x157 + x160 + x165 + V{0.0138888472222222};
auto x171 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x169 + x170;
auto x172 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x169 + x170;
auto x0 = -x19*(V{0.111111111111111}*x129*x43 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x22*(V{0.1666665}*cell[15] + V{0.1666665}*cell[16] + V{0.333333}*cell[2] + V{0.333333}*cell[4] + V{0.1666665}*cell[6] + V{0.1666665}*cell[7] + V{0.333333}*cell[9] + x24 + x25) + x27*(V{0.1666665}*cell[13] + V{0.1666665}*cell[14] + V{0.333333}*cell[18] + V{0.333333}*cell[3] + V{0.1666665}*cell[4] + V{0.1666665}*cell[5] + V{0.333333}*cell[6] + x24 + x28) - x29*(V{0.333333}*cell[10] + V{0.333333}*cell[13] + V{0.333333}*cell[15] + V{0.1666665}*cell[17] + V{0.1666665}*cell[18] + V{0.1666665}*cell[8] + V{0.1666665}*cell[9] + x23 + x25 + x28 + V{0.1666665}) - x36*x47 - x38*x47 - x40*x47 + x44 + V{0.833332833333333});
auto x1 = -x19*(V{0.0185185185185185}*x129*x95 + V{0.0555555555555556}) + x20*(-x141 - x96);
auto x2 = -x19*(V{0.0185185185185185}*x129*x92 + V{0.0555555555555556}) + x20*(-x123 - x143);
auto x3 = -x19*(V{0.0185185185185185}*x106*x129 + V{0.0555555555555556}) + x20*(-x107 - x144);
auto x4 = -x19*(V{0.00925925925925926}*x129*x84 + V{0.0277777777777778}) - x20*(x145*x84 + x162);
auto x5 = -x19*(V{0.00925925925925926}*x119*x129 + V{0.0277777777777778}) - x20*(x145*(x111 + x118) + x163);
auto x6 = -x19*(V{0.00925925925925926}*x129*x68 + V{0.0277777777777778}) - x20*(x145*x68 + x167);
auto x7 = -x19*(V{0.00925925925925926}*x129*x75 + V{0.0277777777777778}) - x20*(x145*(x74 + x77) + x168);
auto x8 = -x19*(V{0.00925925925925926}*x103*x129 + V{0.0277777777777778}) - x20*(x103*x145 + x171);
auto x9 = -x19*(V{0.00925925925925926}*x129*x88 + V{0.0277777777777778}) - x20*(x145*(x120 + x87) + x172);
auto x10 = -x19*(V{0.0185185185185185}*x128*x98 + V{0.0555555555555556}) + x20*(-x141 - x99);
auto x11 = -x19*(V{0.0185185185185185}*x124*x128 + V{0.0555555555555556}) + x20*(-x125 - x143);
auto x12 = -x19*(V{0.0185185185185185}*x108*x128 + V{0.0555555555555556}) + x20*(-x109 - x144);
auto x13 = -x19*(V{0.00925925925925926}*x115*x128 + V{0.0277777777777778}) - x20*(-x115*x145 + x162);
auto x14 = -x19*(V{0.00925925925925926}*x112*x129 + V{0.0277777777777778}) - x20*(x112*x145 + x163);
auto x15 = -x19*(V{0.00925925925925926}*x128*x66 + V{0.0277777777777778}) - x20*(-x145*x66 + x167);
auto x16 = -x19*(V{0.00925925925925926}*x129*x78 + V{0.0277777777777778}) - x20*(x145*x78 + x168);
auto x17 = -x19*(V{0.00925925925925926}*x127*x128 + V{0.0277777777777778}) - x20*(-x127*x145 + x171);
auto x18 = -x19*(V{0.00925925925925926}*x121*x129 + V{0.0277777777777778}) - x20*(x121*x145 + x172);
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
return { V{0.333333333333333}*x129, x36 + x38 + x40 };
}
};

}

}
