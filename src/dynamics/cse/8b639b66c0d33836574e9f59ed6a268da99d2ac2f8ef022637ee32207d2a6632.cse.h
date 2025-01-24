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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<1, 1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<1, 1, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{2}*cell[7];
auto x22 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{-1};
auto x23 = V{1} / (x22);
auto x24 = V{0.24999975}*cell[0] + V{0.24999975}*cell[11] + V{0.24999975}*cell[2] + V{0.4999995}*cell[7] + V{0.24999975};
auto x25 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{1});
auto x26 = cell[0] + cell[11] + cell[2] + x21 + V{1};
auto x27 = cell[10] + V{2}*cell[12] + cell[13] + cell[14] + V{2}*cell[15] + V{2}*cell[17] + cell[1] + cell[4] + cell[5] + V{2}*cell[9] + x26;
auto x28 = cell[12] + cell[17] + cell[18] + V{2}*cell[1] + cell[3] + V{2}*cell[4] + V{2}*cell[5] + V{2}*cell[6] + cell[8] + cell[9] + x26;
auto x29 = -V{0.0138888888888889}*x23*x28 + V{0.0138888888888889}*x25*x27;
auto x30 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x31 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x30;
auto x32 = -x31;
auto x33 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x34 = -x33;
auto x35 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x36 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x37 = V{1.5}*x36;
auto x38 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x39 = V{1.5}*x38;
auto x40 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x41 = V{1.5}*x40;
auto x42 = x39 + x41 + V{-1};
auto x43 = x37 + x42;
auto x44 = x35 + x43;
auto x45 = x34 + x44;
auto x46 = x45 - V{4.5}*x32*x32;
auto x47 = x43*(-V{0.166666666666667}*x23*x28 + V{0.166666666666667}*x25*x27);
auto x48 = -1/x22;
auto x49 = V{0.25}*x25*x27 + V{0.25}*x28*x48;
auto x50 = -V{4.5}*x31*x31;
auto x51 = -x35 + x43;
auto x52 = x33 + x50 + x51;
auto x53 = V{5.55111512312578e-17}*x25*x27;
auto x54 = V{5.55111512312578e-17}*x28*x48 + x53;
auto x55 = V{1.11022302462516e-16}*x25*x27;
auto x56 = V{1.11022302462516e-16}*x28*x48 + x55;
auto x57 = V{8.32667268468867e-17}*x25*x27;
auto x58 = V{8.32667268468867e-17}*x28*x48 + x57;
auto x59 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x60 = -x59;
auto x61 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x62 = -x61;
auto x63 = x44 + x62;
auto x64 = x63 - V{4.5}*x60*x60;
auto x65 = -V{0.0277777777777778}*x23*x28 + V{0.0277777777777778}*x25*x27;
auto x66 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x30;
auto x67 = -x66;
auto x68 = x43 + x61;
auto x69 = x34 + x68;
auto x70 = x69 - V{4.5}*x67*x67;
auto x71 = -V{0.0277777777777778}*x23*x28 + V{0.0277777777777778}*x25*x27;
auto x72 = V{1.66533453693773e-16}*cell[1];
auto x73 = V{3}*x38;
auto x74 = -x37;
auto x75 = V{1} - x41;
auto x76 = x74 + x75;
auto x77 = x35 + x76;
auto x78 = x73 + x77;
auto x79 = x65*x78;
auto x80 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x81 = V{4.5}*(x80*x80);
auto x82 = x44 + x61 - x81;
auto x83 = -x29*x82;
auto x84 = -x39;
auto x85 = x61 + x84;
auto x86 = x77 + x81 + x85;
auto x87 = x29*x86;
auto x88 = -V{4.5}*x59*x59;
auto x89 = x51 + x61 + x88;
auto x90 = -x29*x89;
auto x91 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x92 = V{4.5}*(x91*x91);
auto x93 = x33 + x84;
auto x94 = x77 + x92 + x93;
auto x95 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x96 = V{4.5}*(x95*x95);
auto x97 = x33 + x76 + x85 + x96;
auto x98 = -V{0.0555555555555556}*x23*x28 + V{0.0555555555555556}*x25*x27;
auto x99 = V{3}*x40;
auto x100 = x74 + x93 + x99 + V{1};
auto x101 = x37 + V{-1};
auto x102 = x101 + x35 + x41 - x73;
auto x103 = x102*x65;
auto x104 = V{3}*x36;
auto x105 = -x104 + x42 + x61;
auto x106 = -x105*x65;
auto x107 = x104 + x75 + x85;
auto x108 = x107*x65;
auto x109 = V{8.32667268468867e-17}*cell[0] + V{8.32667268468867e-17}*cell[11] + V{8.32667268468867e-17}*cell[2] + V{1.66533453693773e-16}*cell[7] + V{8.32667268468867e-17};
auto x110 = x106 + x108 + x23*(V{8.32667268468867e-17}*cell[12] + V{8.32667268468867e-17}*cell[17] + V{8.32667268468867e-17}*cell[18] + V{8.32667268468867e-17}*cell[3] + V{1.66533453693773e-16}*cell[4] + V{1.66533453693773e-16}*cell[5] + V{1.66533453693773e-16}*cell[6] + V{8.32667268468867e-17}*cell[8] + V{8.32667268468867e-17}*cell[9] + x109 + x72) - x25*(V{8.32667268468867e-17}*cell[10] + V{1.66533453693773e-16}*cell[12] + V{8.32667268468867e-17}*cell[13] + V{8.32667268468867e-17}*cell[14] + V{1.66533453693773e-16}*cell[15] + V{1.66533453693773e-16}*cell[17] + V{8.32667268468867e-17}*cell[1] + V{8.32667268468867e-17}*cell[4] + V{8.32667268468867e-17}*cell[5] + V{1.66533453693773e-16}*cell[9] + x109) - x47 + V{4.44089209850063e-16};
auto x111 = V{1.66533453693773e-16}*cell[10] + V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.44089209850063e-16}*cell[13] + V{4.9960036108132e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{2.22044604925031e-16}*cell[17] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.44089209850063e-16}*cell[4] + V{4.9960036108132e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{6.66133814775094e-16}*cell[7] + V{2.22044604925031e-16}*cell[8] + x100*x98 - x103 + x110 - x52*(-V{4.62592926927149e-18}*x23*x28 + V{4.62592926927149e-18}*x25*x27) + x65*x94 + x65*x97 + x72 + x79 + x83 + x87 + x90;
auto x112 = -V{0.0277777777777778}*x23*x28 + V{0.0277777777777778}*x25*x27;
auto x113 = x100*x65;
auto x114 = x33 + x68 - x96;
auto x115 = -x114*x29;
auto x116 = x29*x97;
auto x117 = -V{4.5}*x66*x66;
auto x118 = x117 + x33 + x43 + x62;
auto x119 = -x118*x29;
auto x120 = x101 + x33 + x39 - x99;
auto x121 = x120*x65;
auto x122 = x33 + x44 - x92;
auto x123 = V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{3.33066907387547e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{4.44089209850063e-16}*cell[15] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{3.33066907387547e-16}*cell[5] + V{4.44089209850063e-16}*cell[6] + V{4.44089209850063e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] - x102*x98 + x110 + x113 + x115 + x116 + x119 - x121 - x122*x65 - x52*(-V{3.08395284618099e-18}*x23*x28 + V{3.08395284618099e-18}*x25*x27) - x65*x82;
auto x124 = -x23*(-x112*x46 + x123 - x29*x70 - x36*x58 - x38*x56 - x40*x54 - x64*x65) + x25*(x111 - x29*x64 - x36*x54 - x38*x58 - x40*x56 - x46*x71 - x65*x70);
auto x125 = V{0.166666666666667}*cell[10];
auto x126 = V{0.166666666666667}*cell[1];
auto x127 = V{0.0833333333333334}*cell[13];
auto x128 = V{0.0833333333333334}*cell[14];
auto x129 = V{0.0833333333333334}*cell[4];
auto x130 = V{0.0833333333333334}*cell[5];
auto x131 = V{0.0833333333333333}*x25*x27 + V{0.0833333333333333}*x28*x48;
auto x132 = V{0.0416666666666667}*x25*x27 + V{0.0416666666666667}*x28*x48;
auto x133 = x132*x40;
auto x134 = V{0.0833333333333334}*cell[15];
auto x135 = V{0.0833333333333334}*cell[6];
auto x136 = V{0.166666666666667}*cell[7];
auto x137 = V{6.93889390390723e-18}*cell[0] + V{6.93889390390723e-18}*cell[11] + V{6.93889390390723e-18}*cell[2] + V{1.38777878078145e-17}*cell[7] + V{6.93889390390723e-18};
auto x138 = V{6.93889390390723e-18}*cell[10] + V{1.38777878078145e-17}*cell[12] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[14] + V{1.38777878078145e-17}*cell[15] + V{1.38777878078145e-17}*cell[17] + V{6.93889390390723e-18}*cell[1] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[5] + V{1.38777878078145e-17}*cell[9] + x137;
auto x139 = x23*(V{6.93889390390723e-18}*cell[12] + V{6.93889390390723e-18}*cell[17] + V{6.93889390390723e-18}*cell[18] + V{1.38777878078145e-17}*cell[1] + V{6.93889390390723e-18}*cell[3] + V{1.38777878078145e-17}*cell[4] + V{1.38777878078145e-17}*cell[5] + V{1.38777878078145e-17}*cell[6] + V{6.93889390390723e-18}*cell[8] + V{6.93889390390723e-18}*cell[9] + x137);
auto x140 = -V{0.00115740740740741}*x23*x28 + V{0.00115740740740741}*x25*x27;
auto x141 = x140*x46;
auto x142 = x132*x36;
auto x143 = V{0.0833333333333333}*cell[11] + V{0.0833333333333333}*cell[2] - x134 - x135 - x136 + x138*x25 - x139 + x140*x52 - x141 - x142 + V{0.0555555555555555};
auto x144 = x25*x27;
auto x145 = x144 - x23*x28;
auto x146 = x144 + x28*x48;
auto x147 = -x138*x25;
auto x148 = -V{0.0833333333333333}*cell[12] - V{0.0833333333333333}*cell[3] + x127 + x128 + x129 + x130 + x133 + x139 + x147 + V{-0.0555555555555555};
auto x149 = V{0.0833333333333334}*cell[17];
auto x150 = V{0.0833333333333334}*cell[18];
auto x151 = V{0.0833333333333334}*cell[8];
auto x152 = V{0.0833333333333334}*cell[9];
auto x153 = x132*x38;
auto x154 = -V{0.0833333333333333}*cell[10] - V{0.0833333333333333}*cell[1] + x149 + x150 + x151 + x152 + x153;
auto x155 = V{0.166666666666667}*cell[11] - V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[2] - V{0.166666666666667}*cell[6] - V{0.333333333333333}*cell[7] - V{0.00231481481481481}*x145*x46 - V{0.0833333333333333}*x146*x36 + x148 + x154 + x52*(-V{0.00231481481481481}*x23*x28 + V{0.00231481481481481}*x25*x27);
auto x156 = V{0.166666666666667}*cell[12];
auto x157 = V{0.166666666666667}*cell[3];
auto x158 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x145;
auto x159 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x158;
auto x160 = -V{0.0208333333333333}*x23*x28 + V{0.0208333333333333}*x25*x27;
auto x161 = V{0.0208333125}*cell[0] + V{0.0208333125}*cell[11] + V{0.0208333125}*cell[2] + V{0.041666625}*cell[7] + V{0.0208333125};
auto x162 = x23*(V{0.0208333125}*cell[12] + V{0.0208333125}*cell[17] + V{0.0208333125}*cell[18] + V{0.041666625}*cell[1] + V{0.0208333125}*cell[3] + V{0.041666625}*cell[4] + V{0.041666625}*cell[5] + V{0.041666625}*cell[6] + V{0.0208333125}*cell[8] + V{0.0208333125}*cell[9] + x161);
auto x163 = -x25*(V{0.0208333125}*cell[10] + V{0.041666625}*cell[12] + V{0.0208333125}*cell[13] + V{0.0208333125}*cell[14] + V{0.041666625}*cell[15] + V{0.041666625}*cell[17] + V{0.0208333125}*cell[1] + V{0.0208333125}*cell[4] + V{0.0208333125}*cell[5] + V{0.041666625}*cell[9] + x161);
auto x164 = -V{0.0416666666666667}*x23*x28 + V{0.0416666666666667}*x25*x27;
auto x165 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x162 + x163 - x164*x38 + V{0.0138888472222222};
auto x166 = -V{0.000578703703703704}*x23*x28 + V{0.000578703703703704}*x25*x27;
auto x167 = x45 + x50;
auto x168 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0833333333333334}*cell[7] - x164*x36 + x166*x167 - x166*x52;
auto x169 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x160*x40 + x165 + x168;
auto x170 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] - x159 + x169;
auto x171 = x63 + x88;
auto x172 = -x171*x29;
auto x173 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] + x159 + x169;
auto x174 = -V{0.00115740740740741}*x23*x28 + V{0.00115740740740741}*x25*x27;
auto x175 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x158;
auto x176 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x164*x40;
auto x177 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x160*x36 + x165 + x176;
auto x178 = V{0.416666666666667}*cell[15] + V{0.416666666666667}*cell[6] - V{0.166666666666667}*cell[7] - x167*x174 + x174*x52 - x175 + x177;
auto x179 = -V{0.00578703703703704}*x23*x28 + V{0.00578703703703704}*x25*x27;
auto x180 = -V{0.0833333333333333}*cell[15] - V{0.0833333333333333}*cell[6] + V{0.833333333333333}*cell[7] + x175 + x177;
auto x181 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x145;
auto x182 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x160*x38 + x162 + x163 + x168 + x176 + V{0.0138888472222222};
auto x183 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x181 + x182;
auto x184 = x117 + x69;
auto x185 = -x184*x29;
auto x186 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x181 + x182;
auto x187 = -V{5.55111512312578e-17}*x23*x28 + x53;
auto x188 = -V{1.11022302462516e-16}*x23*x28 + x55;
auto x189 = -V{8.32667268468867e-17}*x23*x28 + x57;
auto x190 = -x23*(-x112*x167 + x123 - x171*x65 + x185 - x187*x40 - x188*x38 - x189*x36) + x25*(x111 - x167*x71 + x172 - x184*x65 - x187*x36 - x188*x40 - x189*x38);
auto x191 = -V{0.0833333333333333}*cell[11] - V{0.0833333333333333}*cell[2] + x134 + x135 + x136 + x141 + x142 - V{0.00115740740740741}*x145*x52;
auto x0 = -x19*(V{0.166666666666667}*x124*x43 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[8] + V{1}*cell[9] + x21 + x23*(V{0.24999975}*cell[12] + V{0.24999975}*cell[17] + V{0.24999975}*cell[18] + V{0.4999995}*cell[1] + V{0.24999975}*cell[3] + V{0.4999995}*cell[4] + V{0.4999995}*cell[5] + V{0.4999995}*cell[6] + V{0.24999975}*cell[8] + V{0.24999975}*cell[9] + x24) - x25*(V{0.24999975}*cell[10] + V{0.4999995}*cell[12] + V{0.24999975}*cell[13] + V{0.24999975}*cell[14] + V{0.4999995}*cell[15] + V{0.4999995}*cell[17] + V{0.24999975}*cell[1] + V{0.24999975}*cell[4] + V{0.24999975}*cell[5] + V{0.4999995}*cell[9] + x24) + x29*x46 - x29*x52 - x36*x49 - x38*x49 - x40*x49 + x47 + V{0.833332833333333});
auto x1 = -x19*(V{0.0277777777777778}*x102*x124 + V{0.0555555555555556}) + x20*(V{0.0833333333333333}*cell[12] + V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.0833333333333333}*cell[3] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x103 - x125 - x126 - x127 - x128 - x129 - x130 + x131*x38 - x133 + x143);
auto x2 = -x19*(V{0.0277777777777778}*x105*x124 + V{0.0555555555555556}) + x20*(-x106 - x155);
auto x3 = -x19*(V{0.0277777777777778}*x120*x124 + V{0.0555555555555556}) + x20*(V{0.0833333333333333}*cell[10] + V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] + V{0.0833333333333333}*cell[1] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] + x121 + x131*x40 + x143 - x149 - x150 - x151 - x152 - x153 - x156 - x157);
auto x4 = -x19*(V{0.0138888888888889}*x124*x82 + V{0.0277777777777778}) - x20*(x170 + x83);
auto x5 = -x19*(V{0.0138888888888889}*x124*x64 + V{0.0277777777777778}) - x20*(x172 + x173);
auto x6 = -x19*(V{0.0138888888888889}*x122*x124 + V{0.0277777777777778}) - x20*(-x122*x29 + x178);
auto x7 = -x19*(V{0.0138888888888889}*x124*x46 + V{0.0277777777777778}) - x20*(-x167*(-V{0.00810185185185185}*x23*x28 + V{0.00810185185185185}*x25*x27) - x179*x52 + x180);
auto x8 = -x19*(V{0.0138888888888889}*x114*x124 + V{0.0277777777777778}) - x20*(x115 + x183);
auto x9 = -x19*(V{0.0138888888888889}*x124*x70 + V{0.0277777777777778}) - x20*(x185 + x186);
auto x10 = x19*(V{0.0277777777777778}*x190*x78 + V{-0.0555555555555556}) + x20*(V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] - x125 - x126 + V{0.0833333333333333}*x146*x38 - x148 - x191 - x79);
auto x11 = x19*(V{0.0277777777777778}*x107*x190 + V{-0.0555555555555556}) + x20*(-x108 - x155);
auto x12 = x19*(V{0.0277777777777778}*x100*x190 + V{-0.0555555555555556}) + x20*(V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] - x113 - x139 + V{0.0833333333333333}*x146*x40 - x147 - x154 - x156 - x157 - x191 + V{0.0555555555555555});
auto x13 = x19*(V{0.0138888888888889}*x190*x86 + V{-0.0277777777777778}) - x20*(x170 + x87);
auto x14 = -x19*(V{0.0138888888888889}*x124*x89 + V{0.0277777777777778}) - x20*(x173 + x90);
auto x15 = x19*(V{0.0138888888888889}*x190*x94 + V{-0.0277777777777778}) - x20*(x178 + x29*x94);
auto x16 = -x19*(V{0.0138888888888889}*x124*x52 + V{0.0277777777777778}) - x20*(x167*x179 + x180 - x52*(-V{0.0196759259259259}*x23*x28 + V{0.0196759259259259}*x25*x27));
auto x17 = x19*(V{0.0138888888888889}*x190*x97 + V{-0.0277777777777778}) - x20*(x116 + x183);
auto x18 = -x19*(V{0.0138888888888889}*x118*x124 + V{0.0277777777777778}) - x20*(x119 + x186);
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
return { V{0.5}*x190, x36 + x38 + x40 };
}
};

}

}
