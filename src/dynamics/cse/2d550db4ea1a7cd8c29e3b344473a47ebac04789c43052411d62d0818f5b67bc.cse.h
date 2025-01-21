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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<0, -1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<0, -1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{2}*cell[9];
auto x22 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{-1};
auto x23 = V{1} / (x22);
auto x24 = V{0.24999975}*cell[0] + V{0.24999975}*cell[10] + V{0.24999975}*cell[1] + V{0.4999995}*cell[9] + V{0.24999975};
auto x25 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{1});
auto x26 = cell[0] + cell[10] + cell[1] + x21 + V{1};
auto x27 = cell[11] + V{2}*cell[12] + cell[13] + cell[14] + V{2}*cell[15] + V{2}*cell[17] + cell[2] + cell[4] + cell[5] + V{2}*cell[7] + x26;
auto x28 = cell[12] + V{2}*cell[14] + cell[15] + cell[16] + V{2}*cell[2] + cell[3] + V{2}*cell[4] + cell[6] + cell[7] + V{2}*cell[8] + x26;
auto x29 = -V{0.0138888888888889}*x23*x28 + V{0.0138888888888889}*x25*x27;
auto x30 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x31 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x30;
auto x32 = -x31;
auto x33 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x34 = -x33;
auto x35 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
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
auto x54 = V{1.11022302462516e-16}*x25*x27;
auto x55 = V{1.11022302462516e-16}*x28*x48 + x54;
auto x56 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x57 = -x56;
auto x58 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x59 = x51 + x58;
auto x60 = x59 - V{4.5}*x57*x57;
auto x61 = -V{0.0277777777777778}*x23*x28 + V{0.0277777777777778}*x25*x27;
auto x62 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x30;
auto x63 = -x62;
auto x64 = x43 + x58;
auto x65 = x34 + x64;
auto x66 = x65 - V{4.5}*x63*x63;
auto x67 = -V{0.0277777777777778}*x23*x28 + V{0.0277777777777778}*x25*x27;
auto x68 = V{8.32667268468867e-17}*x25*x27;
auto x69 = V{3}*x38;
auto x70 = x37 + V{-1};
auto x71 = x41 + x58 - x69 + x70;
auto x72 = -x61*x71;
auto x73 = -x37;
auto x74 = V{1} - x41;
auto x75 = x73 + x74;
auto x76 = x58 + x75;
auto x77 = x69 + x76;
auto x78 = x61*x77;
auto x79 = V{8.32667268468867e-17}*cell[0] + V{8.32667268468867e-17}*cell[10] + V{8.32667268468867e-17}*cell[1] + V{1.66533453693773e-16}*cell[9] + V{8.32667268468867e-17};
auto x80 = V{1.66533453693773e-16}*cell[10] + V{2.22044604925031e-16}*cell[17] + V{1.66533453693773e-16}*cell[1] + V{2.22044604925031e-16}*cell[8] + x23*(V{8.32667268468867e-17}*cell[12] + V{1.66533453693773e-16}*cell[14] + V{8.32667268468867e-17}*cell[15] + V{8.32667268468867e-17}*cell[16] + V{1.66533453693773e-16}*cell[2] + V{8.32667268468867e-17}*cell[3] + V{1.66533453693773e-16}*cell[4] + V{8.32667268468867e-17}*cell[6] + V{8.32667268468867e-17}*cell[7] + V{1.66533453693773e-16}*cell[8] + x79) - x25*(V{8.32667268468867e-17}*cell[11] + V{1.66533453693773e-16}*cell[12] + V{8.32667268468867e-17}*cell[13] + V{8.32667268468867e-17}*cell[14] + V{1.66533453693773e-16}*cell[15] + V{1.66533453693773e-16}*cell[17] + V{8.32667268468867e-17}*cell[2] + V{8.32667268468867e-17}*cell[4] + V{8.32667268468867e-17}*cell[5] + V{1.66533453693773e-16}*cell[7] + x79) - x47 + x72 + x78 + V{4.44089209850063e-16};
auto x81 = -x38*(V{8.32667268468867e-17}*x28*x48 + x68) + x80;
auto x82 = V{3}*x36;
auto x83 = -x39;
auto x84 = x35 + x83;
auto x85 = x74 + x82 + x84;
auto x86 = x61*x85;
auto x87 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x88 = V{4.5}*(x87*x87);
auto x89 = x44 + x58 - x88;
auto x90 = -x29*x89;
auto x91 = x76 + x84 + x88;
auto x92 = x29*x91;
auto x93 = -x58;
auto x94 = -V{4.5}*x56*x56;
auto x95 = x44 + x93 + x94;
auto x96 = -x29*x95;
auto x97 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x98 = V{4.5}*(x97*x97);
auto x99 = x33 + x83;
auto x100 = x76 + x98 + x99;
auto x101 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x102 = V{4.5}*(x101*x101);
auto x103 = x102 + x33 + x75 + x84;
auto x104 = -V{0.0555555555555556}*x23*x28 + V{0.0555555555555556}*x25*x27;
auto x105 = V{3}*x40;
auto x106 = x105 + x73 + x99 + V{1};
auto x107 = x35 + x42 - x82;
auto x108 = x107*x61;
auto x109 = V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.44089209850063e-16}*cell[13] + V{4.9960036108132e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.44089209850063e-16}*cell[4] + V{4.9960036108132e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + x100*x61 + x103*x61 + x104*x106 - x108 + x86 + x90 + x92 + x96;
auto x110 = V{4.85722573273506e-17}*x25*x27;
auto x111 = -V{0.0277777777777778}*x23*x28 + V{0.0277777777777778}*x25*x27;
auto x112 = x106*x61;
auto x113 = x33 + x64 - x98;
auto x114 = -x113*x29;
auto x115 = x100*x29;
auto x116 = -V{4.5}*x62*x62;
auto x117 = x116 + x33 + x43 + x93;
auto x118 = -x117*x29;
auto x119 = -x105 + x33 + x39 + x70;
auto x120 = x119*x61;
auto x121 = -x102 + x33 + x44;
auto x122 = V{2.22044604925031e-16}*cell[11] + V{9.71445146547012e-17}*cell[12] + V{1.11022302462516e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{2.22044604925031e-16}*cell[2] + V{9.71445146547012e-17}*cell[3] + V{1.11022302462516e-16}*cell[4] + V{3.33066907387547e-16}*cell[5] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[9] - x104*x107 + x112 + x114 + x115 + x118 - x120 - x121*x61 - x52*(-V{3.08395284618099e-18}*x23*x28 + V{3.08395284618099e-18}*x25*x27) - x61*x89 - x61*x95;
auto x123 = -x23*(-x111*x46 + x122 - x29*x66 - x36*x55 - x40*(x110 + V{4.85722573273506e-17}*x28*x48) + x81) + x25*(x109 - x29*x60 - x36*(V{5.55111512312578e-17}*x28*x48 + x53) - x40*x55 - x46*x67 - x61*x66 + x81);
auto x124 = -V{0.00231481481481481}*x23*x28 + V{0.00231481481481481}*x25*x27;
auto x125 = V{0.0416666666666667}*x25*x27;
auto x126 = x125 - V{0.0416666666666667}*x23*x28;
auto x127 = x45 + x50;
auto x128 = V{0.0833333333333333}*x25*x27;
auto x129 = V{0.0833333333333334}*cell[13];
auto x130 = V{0.0833333333333334}*cell[14];
auto x131 = V{0.0833333333333334}*cell[4];
auto x132 = V{0.0833333333333334}*cell[5];
auto x133 = V{0.0833333333333333}*cell[12];
auto x134 = V{0.0833333333333333}*cell[3];
auto x135 = V{6.93889390390723e-18}*cell[0] + V{6.93889390390723e-18}*cell[10] + V{6.93889390390723e-18}*cell[1] + V{1.38777878078145e-17}*cell[9] + V{6.93889390390723e-18};
auto x136 = x23*(V{6.93889390390723e-18}*cell[12] + V{1.38777878078145e-17}*cell[14] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[16] + V{1.38777878078145e-17}*cell[2] + V{6.93889390390723e-18}*cell[3] + V{1.38777878078145e-17}*cell[4] + V{6.93889390390723e-18}*cell[6] + V{6.93889390390723e-18}*cell[7] + V{1.38777878078145e-17}*cell[8] + x135);
auto x137 = x25*(V{6.93889390390723e-18}*cell[11] + V{1.38777878078145e-17}*cell[12] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[14] + V{1.38777878078145e-17}*cell[15] + V{1.38777878078145e-17}*cell[17] + V{6.93889390390723e-18}*cell[2] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[5] + V{1.38777878078145e-17}*cell[7] + x135);
auto x138 = -x137;
auto x139 = x129 + x130 + x131 + x132 - x133 - x134 + x136 + x138 + V{-0.0555555555555555};
auto x140 = V{0.0833333333333334}*cell[15];
auto x141 = V{0.0833333333333334}*cell[16];
auto x142 = V{0.0833333333333334}*cell[6];
auto x143 = V{0.0833333333333334}*cell[7];
auto x144 = V{0.0833333333333333}*cell[11];
auto x145 = V{0.0833333333333333}*cell[2];
auto x146 = x140 + x141 + x142 + x143 - x144 - x145;
auto x147 = V{0.166666666666667}*cell[10] - V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[1] - V{0.166666666666667}*cell[8] - V{0.333333333333333}*cell[9] - x124*x127 + x124*x52 + x126*x36 + x126*x40 + x139 + x146 - x38*(x128 - V{0.0833333333333333}*x23*x28);
auto x148 = V{0.166666666666667}*cell[11];
auto x149 = V{0.166666666666667}*cell[2];
auto x150 = x128 + V{0.0833333333333333}*x28*x48;
auto x151 = x125 + V{0.0416666666666667}*x28*x48;
auto x152 = x151*x40;
auto x153 = V{0.0833333333333334}*cell[17];
auto x154 = V{0.0833333333333334}*cell[8];
auto x155 = V{0.166666666666667}*cell[9];
auto x156 = -V{0.00115740740740741}*x23*x28 + V{0.00115740740740741}*x25*x27;
auto x157 = x156*x46;
auto x158 = x151*x38;
auto x159 = V{0.0833333333333333}*cell[10] + V{0.0833333333333333}*cell[1] - x136 + x137 - x153 - x154 - x155 + x156*x52 - x157 - x158 + V{0.0555555555555555};
auto x160 = V{0.166666666666667}*cell[12];
auto x161 = V{0.166666666666667}*cell[3];
auto x162 = x151*x36;
auto x163 = x25*x27;
auto x164 = x163 - x23*x28;
auto x165 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x164;
auto x166 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x165;
auto x167 = -V{0.0208333333333333}*x23*x28 + V{0.0208333333333333}*x25*x27;
auto x168 = V{0.0208333125}*cell[0] + V{0.0208333125}*cell[10] + V{0.0208333125}*cell[1] + V{0.041666625}*cell[9] + V{0.0208333125};
auto x169 = x23*(V{0.0208333125}*cell[12] + V{0.041666625}*cell[14] + V{0.0208333125}*cell[15] + V{0.0208333125}*cell[16] + V{0.041666625}*cell[2] + V{0.0208333125}*cell[3] + V{0.041666625}*cell[4] + V{0.0208333125}*cell[6] + V{0.0208333125}*cell[7] + V{0.041666625}*cell[8] + x168);
auto x170 = -x25*(V{0.0208333125}*cell[11] + V{0.041666625}*cell[12] + V{0.0208333125}*cell[13] + V{0.0208333125}*cell[14] + V{0.041666625}*cell[15] + V{0.041666625}*cell[17] + V{0.0208333125}*cell[2] + V{0.0208333125}*cell[4] + V{0.0208333125}*cell[5] + V{0.041666625}*cell[7] + x168);
auto x171 = -V{0.0416666666666667}*x23*x28 + V{0.0416666666666667}*x25*x27;
auto x172 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] + x169 + x170 - x171*x36 + V{0.0138888472222222};
auto x173 = -V{0.000578703703703704}*x23*x28 + V{0.000578703703703704}*x25*x27;
auto x174 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0833333333333334}*cell[9] + x127*x173 - x171*x38 - x173*x52;
auto x175 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x167*x40 + x172 + x174;
auto x176 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] - x166 + x175;
auto x177 = -x29*(x59 + x94);
auto x178 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] + x166 + x175;
auto x179 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x165;
auto x180 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x171*x40;
auto x181 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x167*x36 + x169 + x170 + x174 + x180 + V{0.0138888472222222};
auto x182 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] - x179 + x181;
auto x183 = x116 + x65;
auto x184 = -x183*x29;
auto x185 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] + x179 + x181;
auto x186 = -V{0.00115740740740741}*x23*x28 + V{0.00115740740740741}*x25*x27;
auto x187 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x164;
auto x188 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x167*x38 + x172 + x180;
auto x189 = V{0.416666666666667}*cell[17] + V{0.416666666666667}*cell[8] - V{0.166666666666667}*cell[9] - x127*x186 + x186*x52 - x187 + x188;
auto x190 = -V{0.00578703703703704}*x23*x28 + V{0.00578703703703704}*x25*x27;
auto x191 = -V{0.0833333333333333}*cell[17] - V{0.0833333333333333}*cell[8] + V{0.833333333333333}*cell[9] + x187 + x188;
auto x192 = -V{1.11022302462516e-16}*x23*x28 + x54;
auto x193 = -x38*(-V{8.32667268468867e-17}*x23*x28 + x68) + x80;
auto x194 = -x23*(-x111*x127 + x122 + x184 - x192*x36 + x193 - x40*(x110 - V{4.85722573273506e-17}*x23*x28)) + x25*(x109 - x127*x67 + x177 - x183*x61 - x192*x40 + x193 - x36*(-V{5.55111512312578e-17}*x23*x28 + x53));
auto x195 = x163 + x28*x48;
auto x196 = -V{0.0833333333333333}*cell[10] - V{0.0833333333333333}*cell[1] + x153 + x154 + x155 + x157 + x158 - V{0.00115740740740741}*x164*x52;
auto x0 = -x19*(V{0.166666666666667}*x123*x43 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + x21 + x23*(V{0.24999975}*cell[12] + V{0.4999995}*cell[14] + V{0.24999975}*cell[15] + V{0.24999975}*cell[16] + V{0.4999995}*cell[2] + V{0.24999975}*cell[3] + V{0.4999995}*cell[4] + V{0.24999975}*cell[6] + V{0.24999975}*cell[7] + V{0.4999995}*cell[8] + x24) - x25*(V{0.24999975}*cell[11] + V{0.4999995}*cell[12] + V{0.24999975}*cell[13] + V{0.24999975}*cell[14] + V{0.4999995}*cell[15] + V{0.4999995}*cell[17] + V{0.24999975}*cell[2] + V{0.24999975}*cell[4] + V{0.24999975}*cell[5] + V{0.4999995}*cell[7] + x24) + x29*x46 - x29*x52 - x36*x49 - x38*x49 - x40*x49 + x47 + V{0.833332833333333});
auto x1 = -x19*(V{0.0277777777777778}*x123*x71 + V{0.0555555555555556}) - x20*(x147 + x72);
auto x2 = -x19*(V{0.0277777777777778}*x107*x123 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] + x108 - x129 - x130 - x131 - x132 + x133 + x134 - x148 - x149 + x150*x36 - x152 + x159);
auto x3 = -x19*(V{0.0277777777777778}*x119*x123 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] + x120 - x140 - x141 - x142 - x143 + x144 + x145 + x150*x40 + x159 - x160 - x161 - x162);
auto x4 = -x19*(V{0.0138888888888889}*x123*x89 + V{0.0277777777777778}) - x20*(x176 + x90);
auto x5 = -x19*(V{0.0138888888888889}*x123*x60 + V{0.0277777777777778}) - x20*(x177 + x178);
auto x6 = -x19*(V{0.0138888888888889}*x113*x123 + V{0.0277777777777778}) - x20*(x114 + x182);
auto x7 = -x19*(V{0.0138888888888889}*x123*x66 + V{0.0277777777777778}) - x20*(x184 + x185);
auto x8 = -x19*(V{0.0138888888888889}*x121*x123 + V{0.0277777777777778}) - x20*(-x121*x29 + x189);
auto x9 = -x19*(V{0.0138888888888889}*x123*x46 + V{0.0277777777777778}) - x20*(-x127*(-V{0.00810185185185185}*x23*x28 + V{0.00810185185185185}*x25*x27) - x190*x52 + x191);
auto x10 = x19*(V{0.0277777777777778}*x194*x77 + V{-0.0555555555555556}) - x20*(x147 + x78);
auto x11 = x19*(V{0.0277777777777778}*x194*x85 + V{-0.0555555555555556}) + x20*(V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] - x139 - x148 - x149 - x152 + V{0.0833333333333333}*x195*x36 - x196 - x86);
auto x12 = x19*(V{0.0277777777777778}*x106*x194 + V{-0.0555555555555556}) + x20*(V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] - x112 - x136 - x138 - x146 - x160 - x161 - x162 + V{0.0833333333333333}*x195*x40 - x196 + V{0.0555555555555555});
auto x13 = x19*(V{0.0138888888888889}*x194*x91 + V{-0.0277777777777778}) - x20*(x176 + x92);
auto x14 = -x19*(V{0.0138888888888889}*x123*x95 + V{0.0277777777777778}) - x20*(x178 + x96);
auto x15 = x19*(V{0.0138888888888889}*x100*x194 + V{-0.0277777777777778}) - x20*(x115 + x182);
auto x16 = -x19*(V{0.0138888888888889}*x117*x123 + V{0.0277777777777778}) - x20*(x118 + x185);
auto x17 = x19*(V{0.0138888888888889}*x103*x194 + V{-0.0277777777777778}) - x20*(x103*x29 + x189);
auto x18 = -x19*(V{0.0138888888888889}*x123*x52 + V{0.0277777777777778}) - x20*(x127*x190 + x191 - x52*(-V{0.0196759259259259}*x23*x28 + V{0.0196759259259259}*x25*x27));
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
return { V{0.5}*x194, x36 + x38 + x40 };
}
};

}

}
