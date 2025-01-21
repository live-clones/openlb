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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<0, 1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<0, 1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{2}*cell[17];
auto x22 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{1});
auto x23 = V{0.24999975}*cell[0] + V{0.24999975}*cell[10] + V{0.4999995}*cell[17] + V{0.24999975}*cell[1] + V{0.24999975};
auto x24 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{1});
auto x25 = cell[0] + cell[10] + cell[1] + x21 + V{1};
auto x26 = V{2}*cell[11] + cell[12] + V{2}*cell[13] + cell[15] + cell[16] + V{2}*cell[18] + cell[3] + V{2}*cell[5] + cell[6] + cell[7] + x25;
auto x27 = cell[11] + V{2}*cell[12] + cell[13] + cell[14] + V{2}*cell[15] + cell[2] + cell[4] + cell[5] + V{2}*cell[7] + V{2}*cell[9] + x25;
auto x28 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x29 = V{1.5}*x28;
auto x30 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x31 = V{1.5}*x30;
auto x32 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x33 = V{1.5}*x32;
auto x34 = x31 + x33 + V{-1};
auto x35 = x29 + x34;
auto x36 = x35*(V{0.166666666666667}*x22*x26 + V{0.166666666666667}*x24*x27);
auto x37 = V{0.25}*x22*x26 + V{0.25}*x24*x27;
auto x38 = V{0.0138888888888889}*x22*x26 + V{0.0138888888888889}*x24*x27;
auto x39 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x40 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x41 = V{4.5}*(x40*x40);
auto x42 = -x29;
auto x43 = V{1} - x33;
auto x44 = x42 + x43;
auto x45 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x46 = -x31;
auto x47 = x45 + x46;
auto x48 = x39 + x41 + x44 + x47;
auto x49 = x35 + x45;
auto x50 = x39 - x41 + x49;
auto x51 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x52 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x51;
auto x53 = -x52;
auto x54 = -x45;
auto x55 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x56 = x35 + x55;
auto x57 = x54 + x56;
auto x58 = x57 - V{4.5}*x53*x53;
auto x59 = V{0.0277777777777778}*x22*x26 + V{0.0277777777777778}*x24*x27;
auto x60 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x61 = -x60;
auto x62 = -x39;
auto x63 = x56 + x62;
auto x64 = x63 - V{4.5}*x61*x61;
auto x65 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x66 = V{4.5}*(x65*x65);
auto x67 = x49 + x55 - x66;
auto x68 = -x38*x67;
auto x69 = x44 + x55;
auto x70 = x47 + x66 + x69;
auto x71 = x38*x70;
auto x72 = -x55;
auto x73 = -V{4.5}*x52*x52;
auto x74 = x49 + x72 + x73;
auto x75 = -x38*x74;
auto x76 = V{3}*x28;
auto x77 = x43 + x47 + x76;
auto x78 = x59*x77;
auto x79 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x80 = V{4.5}*(x79*x79);
auto x81 = x39 + x46;
auto x82 = x69 + x80 + x81;
auto x83 = V{0.0555555555555556}*x22*x26 + V{0.0555555555555556}*x24*x27;
auto x84 = V{3}*x32;
auto x85 = x42 + x81 + x84 + V{1};
auto x86 = V{5.55111512312578e-17}*x22*x26 + V{5.55111512312578e-17}*x24*x27;
auto x87 = V{1.11022302462516e-16}*x22*x26 + V{1.11022302462516e-16}*x24*x27;
auto x88 = x34 + x45 - x76;
auto x89 = x59*x88;
auto x90 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x51;
auto x91 = -V{4.5}*x90*x90;
auto x92 = x49 + x62 + x91;
auto x93 = V{3}*x30;
auto x94 = x29 + V{-1};
auto x95 = x33 + x55 - x93 + x94;
auto x96 = -x59*x95;
auto x97 = x69 + x93;
auto x98 = x59*x97;
auto x99 = V{8.32667268468867e-17}*cell[0] + V{8.32667268468867e-17}*cell[10] + V{1.66533453693773e-16}*cell[17] + V{8.32667268468867e-17}*cell[1] + V{8.32667268468867e-17};
auto x100 = V{1.66533453693773e-16}*cell[10] + V{1.66533453693773e-16}*cell[1] - x22*(V{1.66533453693773e-16}*cell[11] + V{8.32667268468867e-17}*cell[12] + V{1.66533453693773e-16}*cell[13] + V{8.32667268468867e-17}*cell[15] + V{8.32667268468867e-17}*cell[16] + V{1.66533453693773e-16}*cell[18] + V{8.32667268468867e-17}*cell[3] + V{1.66533453693773e-16}*cell[5] + V{8.32667268468867e-17}*cell[6] + V{8.32667268468867e-17}*cell[7] + x99) - x24*(V{8.32667268468867e-17}*cell[11] + V{1.66533453693773e-16}*cell[12] + V{8.32667268468867e-17}*cell[13] + V{8.32667268468867e-17}*cell[14] + V{1.66533453693773e-16}*cell[15] + V{8.32667268468867e-17}*cell[2] + V{8.32667268468867e-17}*cell[4] + V{8.32667268468867e-17}*cell[5] + V{1.66533453693773e-16}*cell[7] + V{1.66533453693773e-16}*cell[9] + x99) - x30*(V{8.32667268468867e-17}*x22*x26 + V{8.32667268468867e-17}*x24*x27) - x36 + x96 + x98 + V{4.44089209850063e-16};
auto x101 = V{2.22044604925031e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{3.33066907387547e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{8.88178419700125e-16}*cell[17] + V{2.22044604925031e-16}*cell[18] + V{2.22044604925031e-16}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{3.33066907387547e-16}*cell[5] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{2.22044604925031e-16}*cell[9] + x100 - x28*x86 - x32*x87 + x48*(V{0.0277777777777778}*x22*x26 + V{0.0277777777777778}*x24*x27) - x50*(V{6.16790569236198e-18}*x22*x26 + V{6.16790569236198e-18}*x24*x27) + x59*x82 - x59*x92 + x68 + x71 + x75 + x78 + x83*x85 - x89;
auto x102 = -x90;
auto x103 = x35 + x39;
auto x104 = x103 + x54;
auto x105 = x104 - V{4.5}*x102*x102;
auto x106 = x39 + x56 - x80;
auto x107 = -x106*x38;
auto x108 = x38*x82;
auto x109 = -V{4.5}*x60*x60;
auto x110 = x103 + x109 + x72;
auto x111 = -x110*x38;
auto x112 = x59*x85;
auto x113 = x31 + x39 - x84 + x94;
auto x114 = x113*x59;
auto x115 = V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.44089209850063e-16}*cell[13] + V{4.9960036108132e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{4.44089209850063e-16}*cell[17] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.44089209850063e-16}*cell[4] + V{4.9960036108132e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + x100 + x107 + x108 + x111 + x112 - x114 - x28*x87 - x32*x86 + x48*(V{0.0277777777777778}*x22*x26 + V{0.0277777777777778}*x24*x27) - x50*(V{3.08395284618099e-18}*x22*x26 + V{3.08395284618099e-18}*x24*x27) + x59*x70 + x77*x83;
auto x116 = x22*(x101 - x38*x58 - x59*x64) + x24*(-x105*x59 + x115 - x38*x64 - x58*x59);
auto x117 = V{0.0833333333333334}*cell[13];
auto x118 = V{0.0833333333333334}*cell[14];
auto x119 = V{0.0833333333333334}*cell[15];
auto x120 = V{0.0833333333333334}*cell[16];
auto x121 = V{0.0833333333333334}*cell[4];
auto x122 = V{0.0833333333333334}*cell[5];
auto x123 = V{0.0833333333333334}*cell[6];
auto x124 = V{0.0833333333333334}*cell[7];
auto x125 = V{0.0833333333333333}*cell[11];
auto x126 = V{0.0833333333333333}*cell[12];
auto x127 = V{0.0833333333333333}*cell[2];
auto x128 = V{0.0833333333333333}*cell[3];
auto x129 = V{6.93889390390723e-18}*cell[0] + V{6.93889390390723e-18}*cell[10] + V{1.38777878078145e-17}*cell[17] + V{6.93889390390723e-18}*cell[1] + V{6.93889390390723e-18};
auto x130 = x22*(V{1.38777878078145e-17}*cell[11] + V{6.93889390390723e-18}*cell[12] + V{1.38777878078145e-17}*cell[13] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[16] + V{1.38777878078145e-17}*cell[18] + V{6.93889390390723e-18}*cell[3] + V{1.38777878078145e-17}*cell[5] + V{6.93889390390723e-18}*cell[6] + V{6.93889390390723e-18}*cell[7] + x129);
auto x131 = x24*(V{6.93889390390723e-18}*cell[11] + V{1.38777878078145e-17}*cell[12] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[14] + V{1.38777878078145e-17}*cell[15] + V{6.93889390390723e-18}*cell[2] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[5] + V{1.38777878078145e-17}*cell[7] + V{1.38777878078145e-17}*cell[9] + x129);
auto x132 = V{0.00231481481481481}*x22*x26 + V{0.00231481481481481}*x24*x27;
auto x133 = V{0.0416666666666667}*x22*x26 + V{0.0416666666666667}*x24*x27;
auto x134 = x133*x32;
auto x135 = x133*x28;
auto x136 = V{0.0833333333333333}*x22*x26 + V{0.0833333333333333}*x24*x27;
auto x137 = V{0.166666666666667}*cell[10] - V{0.333333333333333}*cell[17] - V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[1] - V{0.166666666666667}*cell[9] + x117 + x118 + x119 + x120 + x121 + x122 + x123 + x124 - x125 - x126 - x127 - x128 - x130 - x131 + x132*x48 + x132*x50 + x134 + x135 - x136*x30 + V{-0.0555555555555555};
auto x138 = V{0.00115740740740741}*x22*x26 + V{0.00115740740740741}*x24*x27;
auto x139 = V{0.0833333333333333}*cell[10] - V{0.166666666666667}*cell[17] - V{0.0833333333333334}*cell[18] + V{0.0833333333333333}*cell[1] - V{0.0833333333333334}*cell[9] + x130 + x131 - x133*x30 + x138*x48 + x138*x50 + V{0.0555555555555555};
auto x140 = -V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] - V{0.166666666666667}*cell[2] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] - x117 - x118 - x121 - x122 + x126 + x128 - x135 + x136*x32 + x139;
auto x141 = -V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] - V{0.166666666666667}*cell[3] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] - x119 - x120 - x123 - x124 + x125 + x127 - x134 + x136*x28 + x139;
auto x142 = x22*x26 + x24*x27;
auto x143 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x142;
auto x144 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x143;
auto x145 = V{0.0208333333333333}*x22*x26 + V{0.0208333333333333}*x24*x27;
auto x146 = V{0.0208333125}*cell[0] + V{0.0208333125}*cell[10] + V{0.041666625}*cell[17] + V{0.0208333125}*cell[1] + V{0.0208333125};
auto x147 = -x22*(V{0.041666625}*cell[11] + V{0.0208333125}*cell[12] + V{0.041666625}*cell[13] + V{0.0208333125}*cell[15] + V{0.0208333125}*cell[16] + V{0.041666625}*cell[18] + V{0.0208333125}*cell[3] + V{0.041666625}*cell[5] + V{0.0208333125}*cell[6] + V{0.0208333125}*cell[7] + x146);
auto x148 = -x24*(V{0.0208333125}*cell[11] + V{0.041666625}*cell[12] + V{0.0208333125}*cell[13] + V{0.0208333125}*cell[14] + V{0.041666625}*cell[15] + V{0.0208333125}*cell[2] + V{0.0208333125}*cell[4] + V{0.0208333125}*cell[5] + V{0.041666625}*cell[7] + V{0.041666625}*cell[9] + x146);
auto x149 = V{0.0416666666666667}*x22*x26 + V{0.0416666666666667}*x24*x27;
auto x150 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] + x147 + x148 - x149*x32 + V{0.0138888472222222};
auto x151 = V{0.000578703703703704}*x22*x26 + V{0.000578703703703704}*x24*x27;
auto x152 = V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[9] - x149*x30 - x151*x48 - x151*x50;
auto x153 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x145*x28 + x150 + x152;
auto x154 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] - x144 + x153;
auto x155 = x109 + x63;
auto x156 = -x155*x38;
auto x157 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] + x144 + x153;
auto x158 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x143;
auto x159 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x149*x28;
auto x160 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x145*x32 + x147 + x148 + x152 + x159 + V{0.0138888472222222};
auto x161 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] - x158 + x160;
auto x162 = x57 + x73;
auto x163 = -x162*x38;
auto x164 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] + x158 + x160;
auto x165 = V{0.00578703703703704}*x22*x26 + V{0.00578703703703704}*x24*x27;
auto x166 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x142;
auto x167 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x145*x30 + x150 + x159;
auto x168 = V{0.833333333333333}*cell[17] - V{0.0833333333333333}*cell[18] - V{0.0833333333333333}*cell[9] - x166 + x167;
auto x169 = x104 + x91;
auto x170 = V{0.00115740740740741}*x22*x26 + V{0.00115740740740741}*x24*x27;
auto x171 = -V{0.166666666666667}*cell[17] + V{0.416666666666667}*cell[18] + V{0.416666666666667}*cell[9] + x166 + x167 + x170*x48 + x170*x50;
auto x172 = x22*(x101 - x155*x59 + x163) + x24*(x115 + x156 - x162*x59 - x169*x59);
auto x0 = -x19*(V{0.166666666666667}*x116*x35 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[9] + x21 - x22*(V{0.4999995}*cell[11] + V{0.24999975}*cell[12] + V{0.4999995}*cell[13] + V{0.24999975}*cell[15] + V{0.24999975}*cell[16] + V{0.4999995}*cell[18] + V{0.24999975}*cell[3] + V{0.4999995}*cell[5] + V{0.24999975}*cell[6] + V{0.24999975}*cell[7] + x23) - x24*(V{0.24999975}*cell[11] + V{0.4999995}*cell[12] + V{0.24999975}*cell[13] + V{0.24999975}*cell[14] + V{0.4999995}*cell[15] + V{0.24999975}*cell[2] + V{0.24999975}*cell[4] + V{0.24999975}*cell[5] + V{0.4999995}*cell[7] + V{0.4999995}*cell[9] + x23) - x28*x37 - x30*x37 - x32*x37 + x36 - x38*x48 - x38*x50 + V{0.833332833333333});
auto x1 = -x19*(V{0.0277777777777778}*x116*x95 + V{0.0555555555555556}) - x20*(x137 + x96);
auto x2 = -x19*(V{0.0277777777777778}*x113*x116 + V{0.0555555555555556}) + x20*(x114 + x140);
auto x3 = -x19*(V{0.0277777777777778}*x116*x88 + V{0.0555555555555556}) + x20*(x141 + x89);
auto x4 = -x19*(V{0.0138888888888889}*x106*x116 + V{0.0277777777777778}) - x20*(x107 + x154);
auto x5 = -x19*(V{0.0138888888888889}*x116*x64 + V{0.0277777777777778}) - x20*(x156 + x157);
auto x6 = -x19*(V{0.0138888888888889}*x116*x67 + V{0.0277777777777778}) - x20*(x161 + x68);
auto x7 = -x19*(V{0.0138888888888889}*x116*x58 + V{0.0277777777777778}) - x20*(x163 + x164);
auto x8 = -x19*(V{0.0138888888888889}*x116*x50 + V{0.0277777777777778}) - x20*(-x165*x48 + x168 - x50*(V{0.0196759259259259}*x22*x26 + V{0.0196759259259259}*x24*x27));
auto x9 = -x19*(V{0.0138888888888889}*x105*x116 + V{0.0277777777777778}) - x20*(-x169*x38 + x171);
auto x10 = x19*(V{0.0277777777777778}*x172*x97 + V{-0.0555555555555556}) - x20*(x137 + x98);
auto x11 = x19*(V{0.0277777777777778}*x172*x85 + V{-0.0555555555555556}) + x20*(-x112 + x140);
auto x12 = x19*(V{0.0277777777777778}*x172*x77 + V{-0.0555555555555556}) + x20*(x141 - x78);
auto x13 = x19*(V{0.0138888888888889}*x172*x82 + V{-0.0277777777777778}) - x20*(x108 + x154);
auto x14 = -x19*(V{0.0138888888888889}*x110*x116 + V{0.0277777777777778}) - x20*(x111 + x157);
auto x15 = x19*(V{0.0138888888888889}*x172*x70 + V{-0.0277777777777778}) - x20*(x161 + x71);
auto x16 = -x19*(V{0.0138888888888889}*x116*x74 + V{0.0277777777777778}) - x20*(x164 + x75);
auto x17 = x19*(V{0.0138888888888889}*x172*x48 + V{-0.0277777777777778}) - x20*(-x165*x50 + x168 + x48*(V{0.00810185185185185}*x22*x26 + V{0.00810185185185185}*x24*x27));
auto x18 = -x19*(V{0.0138888888888889}*x116*x92 + V{0.0277777777777778}) - x20*(x171 - x38*x92);
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
return { V{0.5}*x172, x28 + x30 + x32 };
}
};

}

}
