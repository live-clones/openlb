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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<1, -1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<1, -1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{2}*cell[16];
auto x22 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{-1};
auto x23 = V{1} / (x22);
auto x24 = V{0.24999975}*cell[0] + V{0.24999975}*cell[11] + V{0.4999995}*cell[16] + V{0.24999975}*cell[2] + V{0.24999975};
auto x25 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{1});
auto x26 = cell[0] + cell[11] + cell[2] + x21 + V{1};
auto x27 = V{2}*cell[10] + cell[12] + V{2}*cell[13] + V{2}*cell[14] + V{2}*cell[15] + cell[17] + cell[18] + cell[3] + cell[8] + cell[9] + x26;
auto x28 = cell[10] + cell[13] + cell[14] + V{2}*cell[18] + cell[1] + V{2}*cell[3] + cell[4] + cell[5] + V{2}*cell[6] + V{2}*cell[8] + x26;
auto x29 = -V{0.0138888888888889}*x23*x28 + V{0.0138888888888889}*x25*x27;
auto x30 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x31 = -x30;
auto x32 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x33 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x32;
auto x34 = -V{4.5}*x33*x33;
auto x35 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x36 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x37 = V{1.5}*x36;
auto x38 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x39 = V{1.5}*x38;
auto x40 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x41 = V{1.5}*x40;
auto x42 = x39 + x41 + V{-1};
auto x43 = x37 + x42;
auto x44 = x35 + x43;
auto x45 = x31 + x34 + x44;
auto x46 = x43*(-V{0.166666666666667}*x23*x28 + V{0.166666666666667}*x25*x27);
auto x47 = -1/x22;
auto x48 = V{0.25}*x25*x27 + V{0.25}*x28*x47;
auto x49 = -x33;
auto x50 = -x35 + x43;
auto x51 = x30 + x50;
auto x52 = x51 - V{4.5}*x49*x49;
auto x53 = -V{0.0277777777777778}*x23*x28 + V{0.0277777777777778}*x25*x27;
auto x54 = V{3}*x36;
auto x55 = V{1} - x39;
auto x56 = -x41;
auto x57 = x35 + x56;
auto x58 = x54 + x55 + x57;
auto x59 = x53*x58;
auto x60 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x61 = V{4.5}*(x60*x60);
auto x62 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x63 = -x37;
auto x64 = x55 + x63;
auto x65 = x62 + x64;
auto x66 = x57 + x61 + x65;
auto x67 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x68 = V{4.5}*(x67*x67);
auto x69 = x30 + x56;
auto x70 = x65 + x68 + x69;
auto x71 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x72 = V{4.5}*(x71*x71);
auto x73 = x30 + x57 + x64 + x72;
auto x74 = -V{0.0555555555555556}*x23*x28 + V{0.0555555555555556}*x25*x27;
auto x75 = V{3}*x38;
auto x76 = x63 + x69 + x75 + V{1};
auto x77 = V{5.55111512312578e-17}*x25*x27 + V{5.55111512312578e-17}*x28*x47;
auto x78 = V{1.11022302462516e-16}*x25*x27 + V{1.11022302462516e-16}*x28*x47;
auto x79 = V{8.32667268468867e-17}*x25*x27 + V{8.32667268468867e-17}*x28*x47;
auto x80 = x44 - x61 + x62;
auto x81 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x32;
auto x82 = -x81;
auto x83 = x50 + x62;
auto x84 = x83 - V{4.5}*x82*x82;
auto x85 = -x62;
auto x86 = -V{4.5}*x81*x81;
auto x87 = x44 + x85 + x86;
auto x88 = x35 + x42 - x54;
auto x89 = x53*x88;
auto x90 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x91 = -V{4.5}*x90*x90;
auto x92 = x43 + x62;
auto x93 = x31 + x91 + x92;
auto x94 = V{3}*x40;
auto x95 = x37 + V{-1};
auto x96 = x39 + x62 - x94 + x95;
auto x97 = -x53*x96;
auto x98 = x65 + x94;
auto x99 = x53*x98;
auto x100 = V{8.32667268468867e-17}*cell[0] + V{8.32667268468867e-17}*cell[11] + V{1.66533453693773e-16}*cell[16] + V{8.32667268468867e-17}*cell[2] + V{8.32667268468867e-17};
auto x101 = V{1.66533453693773e-16}*cell[10];
auto x102 = x23*(V{8.32667268468867e-17}*cell[10] + V{8.32667268468867e-17}*cell[13] + V{8.32667268468867e-17}*cell[14] + V{1.66533453693773e-16}*cell[18] + V{8.32667268468867e-17}*cell[1] + V{1.66533453693773e-16}*cell[3] + V{8.32667268468867e-17}*cell[4] + V{8.32667268468867e-17}*cell[5] + V{1.66533453693773e-16}*cell[6] + V{1.66533453693773e-16}*cell[8] + x100) - x25*(V{8.32667268468867e-17}*cell[12] + V{1.66533453693773e-16}*cell[13] + V{1.66533453693773e-16}*cell[14] + V{1.66533453693773e-16}*cell[15] + V{8.32667268468867e-17}*cell[17] + V{8.32667268468867e-17}*cell[18] + V{8.32667268468867e-17}*cell[3] + V{8.32667268468867e-17}*cell[8] + V{8.32667268468867e-17}*cell[9] + x100 + x101) - x46 + x97 + x99 + V{4.44089209850063e-16};
auto x103 = x25*(V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{3.33066907387547e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{4.44089209850063e-16}*cell[15] + V{4.44089209850063e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{3.33066907387547e-16}*cell[5] + V{4.44089209850063e-16}*cell[6] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + x102 + x29*x66 - x29*x80 - x29*x84 - x29*x87 - x36*x77 - x38*x78 - x40*x79 - x45*(-V{0.0277777777777778}*x23*x28 + V{0.0277777777777778}*x25*x27) - x52*(-V{3.08395284618099e-18}*x23*x28 + V{3.08395284618099e-18}*x25*x27) + x53*x70 + x53*x73 - x53*x93 + x59 + x74*x76 - x89);
auto x104 = x53*x76;
auto x105 = x30 - x68 + x92;
auto x106 = -x90;
auto x107 = x30 + x43 + x85;
auto x108 = x107 - V{4.5}*x106*x106;
auto x109 = x30 + x41 - x75 + x95;
auto x110 = x109*x53;
auto x111 = x30 + x44 - x72;
auto x112 = x23*(V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.71844785465692e-16}*cell[13] + V{5.27355936696949e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{6.66133814775094e-16}*cell[16] + V{2.22044604925031e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.71844785465692e-16}*cell[4] + V{5.27355936696949e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{2.22044604925031e-16}*cell[9] + x101 + x102 + x104 - x105*x29 - x108*x29 - x110 - x111*x53 + x29*x70 - x29*x93 - x36*x78 - x38*x79 - x40*x77 - x45*(-V{0.0277777777777778}*x23*x28 + V{0.0277777777777778}*x25*x27) - x52*(-V{4.62592926927149e-18}*x23*x28 + V{4.62592926927149e-18}*x25*x27) - x53*x80 - x53*x87 - x74*x88);
auto x113 = x103 - x112;
auto x114 = V{0.166666666666667}*cell[10];
auto x115 = V{0.166666666666667}*cell[1];
auto x116 = V{0.0833333333333334}*cell[13];
auto x117 = V{0.0833333333333334}*cell[14];
auto x118 = V{0.0833333333333334}*cell[4];
auto x119 = V{0.0833333333333334}*cell[5];
auto x120 = V{0.0833333333333333}*x25*x27 + V{0.0833333333333333}*x28*x47;
auto x121 = V{0.0416666666666667}*x25*x27 + V{0.0416666666666667}*x28*x47;
auto x122 = x121*x36;
auto x123 = V{0.0833333333333334}*cell[15];
auto x124 = V{0.0833333333333334}*cell[6];
auto x125 = V{0.166666666666667}*cell[16];
auto x126 = V{6.93889390390723e-18}*cell[0] + V{6.93889390390723e-18}*cell[11] + V{1.38777878078145e-17}*cell[16] + V{6.93889390390723e-18}*cell[2] + V{6.93889390390723e-18};
auto x127 = V{1.38777878078145e-17}*cell[10] + V{6.93889390390723e-18}*cell[12] + V{1.38777878078145e-17}*cell[13] + V{1.38777878078145e-17}*cell[14] + V{1.38777878078145e-17}*cell[15] + V{6.93889390390723e-18}*cell[17] + V{6.93889390390723e-18}*cell[18] + V{6.93889390390723e-18}*cell[3] + V{6.93889390390723e-18}*cell[8] + V{6.93889390390723e-18}*cell[9] + x126;
auto x128 = x23*(V{6.93889390390723e-18}*cell[10] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[14] + V{1.38777878078145e-17}*cell[18] + V{6.93889390390723e-18}*cell[1] + V{1.38777878078145e-17}*cell[3] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[5] + V{1.38777878078145e-17}*cell[6] + V{1.38777878078145e-17}*cell[8] + x126);
auto x129 = -V{0.00115740740740741}*x23*x28 + V{0.00115740740740741}*x25*x27;
auto x130 = x129*x45;
auto x131 = x121*x40;
auto x132 = V{0.0833333333333333}*cell[11] + V{0.0833333333333333}*cell[2] - x123 - x124 - x125 + x127*x25 - x128 + x129*x52 - x130 - x131 + V{0.0555555555555555};
auto x133 = x23*x28 - x25*x27;
auto x134 = -x133;
auto x135 = x25*x27 + x28*x47;
auto x136 = -x127*x25;
auto x137 = -V{0.0833333333333333}*cell[12] - V{0.0833333333333333}*cell[3] + x116 + x117 + x118 + x119 + x122 + x128 + x136 + V{-0.0555555555555555};
auto x138 = V{0.0833333333333334}*cell[17];
auto x139 = V{0.0833333333333334}*cell[18];
auto x140 = V{0.0833333333333334}*cell[8];
auto x141 = V{0.0833333333333334}*cell[9];
auto x142 = x121*x38;
auto x143 = -V{0.0833333333333333}*cell[10] - V{0.0833333333333333}*cell[1] + x138 + x139 + x140 + x141 + x142;
auto x144 = V{0.166666666666667}*cell[11] - V{0.166666666666667}*cell[15] - V{0.333333333333333}*cell[16] + V{0.166666666666667}*cell[2] - V{0.166666666666667}*cell[6] - V{0.00231481481481481}*x134*x45 - V{0.0833333333333333}*x135*x40 + x137 + x143 + x52*(-V{0.00231481481481481}*x23*x28 + V{0.00231481481481481}*x25*x27);
auto x145 = V{0.166666666666667}*cell[12];
auto x146 = V{0.166666666666667}*cell[3];
auto x147 = V{0.0138888888888889}*x23*x28 - V{0.0138888888888889}*x25*x27;
auto x148 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x149 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x133;
auto x150 = x148*x149;
auto x151 = -V{0.0208333333333333}*x23*x28 + V{0.0208333333333333}*x25*x27;
auto x152 = V{0.0208333125}*cell[0] + V{0.0208333125}*cell[11] + V{0.041666625}*cell[16] + V{0.0208333125}*cell[2] + V{0.0208333125};
auto x153 = x23*(V{0.0208333125}*cell[10] + V{0.0208333125}*cell[13] + V{0.0208333125}*cell[14] + V{0.041666625}*cell[18] + V{0.0208333125}*cell[1] + V{0.041666625}*cell[3] + V{0.0208333125}*cell[4] + V{0.0208333125}*cell[5] + V{0.041666625}*cell[6] + V{0.041666625}*cell[8] + x152);
auto x154 = -x25*(V{0.041666625}*cell[10] + V{0.0208333125}*cell[12] + V{0.041666625}*cell[13] + V{0.041666625}*cell[14] + V{0.041666625}*cell[15] + V{0.0208333125}*cell[17] + V{0.0208333125}*cell[18] + V{0.0208333125}*cell[3] + V{0.0208333125}*cell[8] + V{0.0208333125}*cell[9] + x152);
auto x155 = -V{0.0416666666666667}*x23*x28 + V{0.0416666666666667}*x25*x27;
auto x156 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x153 + x154 - x155*x38 + V{0.0138888472222222};
auto x157 = V{0.000578703703703704}*x23*x28 - V{0.000578703703703704}*x25*x27;
auto x158 = x34 + x51;
auto x159 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0833333333333334}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] - x155*x40 + x157*x158 - x157*x45;
auto x160 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x151*x36 + x156 + x159;
auto x161 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] + x150 + x160;
auto x162 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] - x150 + x160;
auto x163 = V{0.00115740740740741}*x23*x28 - V{0.00115740740740741}*x25*x27;
auto x164 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x149;
auto x165 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x155*x36;
auto x166 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x151*x40 + x156 + x165;
auto x167 = V{0.416666666666667}*cell[15] - V{0.166666666666667}*cell[16] + V{0.416666666666667}*cell[6] - x158*x163 + x163*x45 + x164 + x166;
auto x168 = V{0.00578703703703704}*x23*x28 - V{0.00578703703703704}*x25*x27;
auto x169 = -V{0.0833333333333333}*cell[15] + V{0.833333333333333}*cell[16] - V{0.0833333333333333}*cell[6] - x164 + x166;
auto x170 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x134*x148;
auto x171 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x151*x38 + x153 + x154 + x159 + x165 + V{0.0138888472222222};
auto x172 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x170 + x171;
auto x173 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x170 + x171;
auto x174 = -V{0.0833333333333333}*cell[11] - V{0.0833333333333333}*cell[2] + x123 + x124 + x125 + x130 + x131 - V{0.00115740740740741}*x134*x52;
auto x175 = -x103 + x112;
auto x0 = -x19*(V{0.166666666666667}*x113*x43 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[8] + V{1}*cell[9] + x21 + x23*(V{0.24999975}*cell[10] + V{0.24999975}*cell[13] + V{0.24999975}*cell[14] + V{0.4999995}*cell[18] + V{0.24999975}*cell[1] + V{0.4999995}*cell[3] + V{0.24999975}*cell[4] + V{0.24999975}*cell[5] + V{0.4999995}*cell[6] + V{0.4999995}*cell[8] + x24) - x25*(V{0.4999995}*cell[10] + V{0.24999975}*cell[12] + V{0.4999995}*cell[13] + V{0.4999995}*cell[14] + V{0.4999995}*cell[15] + V{0.24999975}*cell[17] + V{0.24999975}*cell[18] + V{0.24999975}*cell[3] + V{0.24999975}*cell[8] + V{0.24999975}*cell[9] + x24) + x29*x45 - x29*x52 - x36*x48 - x38*x48 - x40*x48 + x46 + V{0.833332833333333});
auto x1 = -x19*(V{0.0277777777777778}*x109*x113 + V{0.0555555555555556}) + x20*(V{0.0833333333333333}*cell[12] + V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.0833333333333333}*cell[3] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x110 - x114 - x115 - x116 - x117 - x118 - x119 + x120*x38 - x122 + x132);
auto x2 = -x19*(V{0.0277777777777778}*x113*x96 + V{0.0555555555555556}) + x20*(-x144 - x97);
auto x3 = -x19*(V{0.0277777777777778}*x113*x88 + V{0.0555555555555556}) + x20*(V{0.0833333333333333}*cell[10] + V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] + V{0.0833333333333333}*cell[1] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] + x120*x36 + x132 - x138 - x139 - x140 - x141 - x142 - x145 - x146 + x89);
auto x4 = -x19*(V{0.0138888888888889}*x105*x113 + V{0.0277777777777778}) - x20*(x105*x147 + x161);
auto x5 = -x19*(V{0.0138888888888889}*x108*x113 + V{0.0277777777777778}) - x20*(x147*(x107 + x91) + x162);
auto x6 = -x19*(V{0.0138888888888889}*x111*x113 + V{0.0277777777777778}) - x20*(x111*x147 + x167);
auto x7 = -x19*(V{0.0138888888888889}*x113*x52 + V{0.0277777777777778}) - x20*(x158*(V{0.0196759259259259}*x23*x28 - V{0.0196759259259259}*x25*x27) - x168*x45 + x169);
auto x8 = -x19*(V{0.0138888888888889}*x113*x80 + V{0.0277777777777778}) - x20*(x147*x80 + x172);
auto x9 = -x19*(V{0.0138888888888889}*x113*x84 + V{0.0277777777777778}) - x20*(x147*(x83 + x86) + x173);
auto x10 = -x19*(V{0.0277777777777778}*x175*x76 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] - x104 - x114 - x115 + V{0.0833333333333333}*x135*x38 - x137 - x174);
auto x11 = -x19*(V{0.0277777777777778}*x175*x98 + V{0.0555555555555556}) + x20*(-x144 - x99);
auto x12 = -x19*(V{0.0277777777777778}*x175*x58 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] - x128 + V{0.0833333333333333}*x135*x36 - x136 - x143 - x145 - x146 - x174 - x59 + V{0.0555555555555555});
auto x13 = -x19*(V{0.0138888888888889}*x175*x70 + V{0.0277777777777778}) - x20*(-x147*x70 + x161);
auto x14 = -x19*(V{0.0138888888888889}*x113*x93 + V{0.0277777777777778}) - x20*(x147*x93 + x162);
auto x15 = -x19*(V{0.0138888888888889}*x175*x73 + V{0.0277777777777778}) - x20*(-x147*x73 + x167);
auto x16 = -x19*(V{0.0138888888888889}*x113*x45 + V{0.0277777777777778}) - x20*(x158*x168 + x169 + x45*(V{0.00810185185185185}*x23*x28 - V{0.00810185185185185}*x25*x27));
auto x17 = -x19*(V{0.0138888888888889}*x175*x66 + V{0.0277777777777778}) - x20*(-x147*x66 + x172);
auto x18 = -x19*(V{0.0138888888888889}*x113*x87 + V{0.0277777777777778}) - x20*(x147*x87 + x173);
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
return { V{0.5}*x113, x36 + x38 + x40 };
}
};

}

}
