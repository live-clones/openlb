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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<0, -1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<0, -1, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{2}*cell[8];
auto x22 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{-1};
auto x23 = V{1} / (x22);
auto x24 = V{0.24999975}*cell[0] + V{0.24999975}*cell[10] + V{0.24999975}*cell[1] + V{0.4999995}*cell[8] + V{0.24999975};
auto x25 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{-1};
auto x26 = V{1} / (x25);
auto x27 = -1/x22;
auto x28 = cell[0] + cell[10] + cell[1] + x21 + V{1};
auto x29 = cell[12] + V{2}*cell[14] + cell[15] + cell[16] + V{2}*cell[2] + cell[3] + V{2}*cell[4] + cell[6] + cell[7] + V{2}*cell[9] + x28;
auto x30 = -1/x25;
auto x31 = cell[11] + cell[13] + cell[14] + V{2}*cell[16] + V{2}*cell[18] + cell[2] + V{2}*cell[3] + cell[4] + cell[5] + V{2}*cell[6] + x28;
auto x32 = V{0.25}*x27*x29 + V{0.25}*x30*x31;
auto x33 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x34 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x35 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x36 = V{0.0138888888888889}*x23*x29 + V{0.0138888888888889}*x26*x31;
auto x37 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x38 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x39 = V{4.5}*(x38*x38);
auto x40 = V{1.5}*x35;
auto x41 = -x40;
auto x42 = V{1.5}*x34;
auto x43 = V{1} - x42;
auto x44 = x41 + x43;
auto x45 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x46 = V{1.5}*x33;
auto x47 = -x46;
auto x48 = x45 + x47;
auto x49 = x37 + x39 + x44 + x48;
auto x50 = x42 + x46 + V{-1};
auto x51 = x40 + x50;
auto x52 = x45 + x51;
auto x53 = x37 - x39 + x52;
auto x54 = x51*(V{0.166666666666667}*x23*x29 + V{0.166666666666667}*x26*x31);
auto x55 = V{0.0277777777777778}*x23*x29 + V{0.0277777777777778}*x26*x31;
auto x56 = V{3}*x35;
auto x57 = x45 + x50 - x56;
auto x58 = x55*x57;
auto x59 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x60 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x61 = V{4.5}*(x60*x60);
auto x62 = x52 + x59 - x61;
auto x63 = x36*x62;
auto x64 = x44 + x59;
auto x65 = x48 + x61 + x64;
auto x66 = -x36*x65;
auto x67 = -x59;
auto x68 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x69 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x68;
auto x70 = -V{4.5}*x69*x69;
auto x71 = x52 + x67 + x70;
auto x72 = x36*x71;
auto x73 = -x69;
auto x74 = -x45;
auto x75 = x51 + x59;
auto x76 = x74 + x75;
auto x77 = x76 - V{4.5}*x73*x73;
auto x78 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x79 = V{4.5}*(x78*x78);
auto x80 = x37 + x75 - x79;
auto x81 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x82 = -V{4.5}*x81*x81;
auto x83 = x37 + x51;
auto x84 = x67 + x82 + x83;
auto x85 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x68;
auto x86 = -x85;
auto x87 = x74 + x83;
auto x88 = x87 - V{4.5}*x86*x86;
auto x89 = V{0.0555555555555556}*x23*x29 + V{0.0555555555555556}*x26*x31;
auto x90 = V{3}*x34;
auto x91 = x40 + V{-1};
auto x92 = x37 + x46 - x90 + x91;
auto x93 = V{1.11022302462516e-16}*x27*x29 + V{1.11022302462516e-16}*x30*x31;
auto x94 = x43 + x48 + x56;
auto x95 = x55*x94;
auto x96 = V{3}*x33;
auto x97 = x42 + x59 + x91 - x96;
auto x98 = x55*x97;
auto x99 = x64 + x96;
auto x100 = -x55*x99;
auto x101 = V{8.32667268468867e-17}*cell[0] + V{8.32667268468867e-17}*cell[10] + V{8.32667268468867e-17}*cell[1] + V{1.66533453693773e-16}*cell[8] + V{8.32667268468867e-17};
auto x102 = V{1.66533453693773e-16}*cell[10] + V{2.22044604925031e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{2.22044604925031e-16}*cell[9] + x100 + x23*(V{8.32667268468867e-17}*cell[12] + V{1.66533453693773e-16}*cell[14] + V{8.32667268468867e-17}*cell[15] + V{8.32667268468867e-17}*cell[16] + V{1.66533453693773e-16}*cell[2] + V{8.32667268468867e-17}*cell[3] + V{1.66533453693773e-16}*cell[4] + V{8.32667268468867e-17}*cell[6] + V{8.32667268468867e-17}*cell[7] + V{1.66533453693773e-16}*cell[9] + x101) + x26*(V{8.32667268468867e-17}*cell[11] + V{8.32667268468867e-17}*cell[13] + V{8.32667268468867e-17}*cell[14] + V{1.66533453693773e-16}*cell[16] + V{1.66533453693773e-16}*cell[18] + V{8.32667268468867e-17}*cell[2] + V{1.66533453693773e-16}*cell[3] + V{8.32667268468867e-17}*cell[4] + V{8.32667268468867e-17}*cell[5] + V{1.66533453693773e-16}*cell[6] + x101) - x33*(V{8.32667268468867e-17}*x27*x29 + V{8.32667268468867e-17}*x30*x31) + x54 + x98 + V{4.44089209850063e-16};
auto x103 = x55*x92;
auto x104 = x36*x80;
auto x105 = x37 + x47;
auto x106 = x105 + x64 + x79;
auto x107 = -x106*x36;
auto x108 = x36*x84;
auto x109 = -x81;
auto x110 = -x37;
auto x111 = x110 + x75;
auto x112 = x111 - V{4.5}*x109*x109;
auto x113 = -V{4.5}*x85*x85;
auto x114 = x110 + x113 + x52;
auto x115 = x105 + x41 + x90 + V{1};
auto x116 = x115*x55;
auto x117 = x23*(V{2.22044604925031e-16}*cell[11] + V{9.71445146547012e-17}*cell[12] + V{1.11022302462516e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{2.22044604925031e-16}*cell[2] + V{9.71445146547012e-17}*cell[3] + V{1.11022302462516e-16}*cell[4] + V{3.33066907387547e-16}*cell[5] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[8] + x102 - x34*x93 - x35*(V{4.85722573273506e-17}*x27*x29 + V{4.85722573273506e-17}*x30*x31) + x36*x77 - x49*(V{3.08395284618099e-18}*x23*x29 + V{3.08395284618099e-18}*x26*x31) + x53*(V{0.0277777777777778}*x23*x29 + V{0.0277777777777778}*x26*x31) + x55*x80 + x55*x84 + x55*x88 + x58 + x63 + x66 + x72 + x89*x92 - x95) + x26*(V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.71844785465692e-16}*cell[13] + V{5.27355936696949e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.71844785465692e-16}*cell[4] + V{5.27355936696949e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + x102 + x103 + x104 + x107 + x108 + x112*x36 + x114*x55 - x116 - x34*(V{5.55111512312578e-17}*x27*x29 + V{5.55111512312578e-17}*x30*x31) - x35*x93 + x53*(V{0.0277777777777778}*x23*x29 + V{0.0277777777777778}*x26*x31) + x55*x62 + x55*x71 + x57*x89);
auto x118 = V{0.00231481481481481}*x23*x29 + V{0.00231481481481481}*x26*x31;
auto x119 = V{0.0416666666666667}*x23*x29 + V{0.0416666666666667}*x26*x31;
auto x120 = V{0.0833333333333334}*cell[13];
auto x121 = V{0.0833333333333334}*cell[14];
auto x122 = V{0.0833333333333334}*cell[4];
auto x123 = V{0.0833333333333334}*cell[5];
auto x124 = V{0.0833333333333333}*cell[12];
auto x125 = V{0.0833333333333333}*cell[3];
auto x126 = V{6.93889390390723e-18}*cell[0] + V{6.93889390390723e-18}*cell[10] + V{6.93889390390723e-18}*cell[1] + V{1.38777878078145e-17}*cell[8] + V{6.93889390390723e-18};
auto x127 = x23*(V{6.93889390390723e-18}*cell[12] + V{1.38777878078145e-17}*cell[14] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[16] + V{1.38777878078145e-17}*cell[2] + V{6.93889390390723e-18}*cell[3] + V{1.38777878078145e-17}*cell[4] + V{6.93889390390723e-18}*cell[6] + V{6.93889390390723e-18}*cell[7] + V{1.38777878078145e-17}*cell[9] + x126);
auto x128 = x26*(V{6.93889390390723e-18}*cell[11] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[14] + V{1.38777878078145e-17}*cell[16] + V{1.38777878078145e-17}*cell[18] + V{6.93889390390723e-18}*cell[2] + V{1.38777878078145e-17}*cell[3] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[5] + V{1.38777878078145e-17}*cell[6] + x126);
auto x129 = x120 + x121 + x122 + x123 - x124 - x125 + x127 + x128 + V{-0.0555555555555555};
auto x130 = V{0.0833333333333334}*cell[15];
auto x131 = V{0.0833333333333334}*cell[16];
auto x132 = V{0.0833333333333334}*cell[6];
auto x133 = V{0.0833333333333334}*cell[7];
auto x134 = V{0.0833333333333333}*cell[11];
auto x135 = V{0.0833333333333333}*cell[2];
auto x136 = x130 + x131 + x132 + x133 - x134 - x135;
auto x137 = V{0.166666666666667}*cell[10] - V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[1] - V{0.333333333333333}*cell[8] - V{0.166666666666667}*cell[9] + x118*x49 + x118*x53 - x119*x34 - x119*x35 + x129 + x136 + x33*(V{0.0833333333333333}*x23*x29 + V{0.0833333333333333}*x26*x31);
auto x138 = V{0.166666666666667}*cell[11];
auto x139 = V{0.166666666666667}*cell[2];
auto x140 = V{0.0416666666666667}*x27*x29 + V{0.0416666666666667}*x30*x31;
auto x141 = x140*x35;
auto x142 = x27*x29 + x30*x31;
auto x143 = V{0.0833333333333334}*cell[18];
auto x144 = V{0.0833333333333334}*cell[9];
auto x145 = V{0.166666666666667}*cell[8];
auto x146 = x140*x33;
auto x147 = x23*x29 + x26*x31;
auto x148 = -V{0.0833333333333333}*cell[10] - V{0.0833333333333333}*cell[1] + x143 + x144 + x145 + x146 - V{0.00115740740740741}*x147*x49 - V{0.00115740740740741}*x147*x53;
auto x149 = V{0.166666666666667}*cell[12];
auto x150 = V{0.166666666666667}*cell[3];
auto x151 = x140*x34;
auto x152 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x147;
auto x153 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x152;
auto x154 = V{0.0208333333333333}*x23*x29 + V{0.0208333333333333}*x26*x31;
auto x155 = V{0.0208333125}*cell[0] + V{0.0208333125}*cell[10] + V{0.0208333125}*cell[1] + V{0.041666625}*cell[8] + V{0.0208333125};
auto x156 = x23*(V{0.0208333125}*cell[12] + V{0.041666625}*cell[14] + V{0.0208333125}*cell[15] + V{0.0208333125}*cell[16] + V{0.041666625}*cell[2] + V{0.0208333125}*cell[3] + V{0.041666625}*cell[4] + V{0.0208333125}*cell[6] + V{0.0208333125}*cell[7] + V{0.041666625}*cell[9] + x155);
auto x157 = x26*(V{0.0208333125}*cell[11] + V{0.0208333125}*cell[13] + V{0.0208333125}*cell[14] + V{0.041666625}*cell[16] + V{0.041666625}*cell[18] + V{0.0208333125}*cell[2] + V{0.041666625}*cell[3] + V{0.0208333125}*cell[4] + V{0.0208333125}*cell[5] + V{0.041666625}*cell[6] + x155);
auto x158 = V{0.0416666666666667}*x23*x29 + V{0.0416666666666667}*x26*x31;
auto x159 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] + x156 + x157 + x158*x34 + V{0.0138888472222222};
auto x160 = V{0.000578703703703704}*x23*x29 + V{0.000578703703703704}*x26*x31;
auto x161 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[8] + V{0.0416666666666667}*cell[9] + x158*x33 - x160*x49 - x160*x53;
auto x162 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] - x154*x35 + x159 + x161;
auto x163 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] + x153 + x162;
auto x164 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] - x153 + x162;
auto x165 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x152;
auto x166 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] + x158*x35;
auto x167 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] - x154*x34 + x156 + x157 + x161 + x166 + V{0.0138888472222222};
auto x168 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] + x165 + x167;
auto x169 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] - x165 + x167;
auto x170 = V{0.00578703703703704}*x23*x29 + V{0.00578703703703704}*x26*x31;
auto x171 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x147;
auto x172 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] - x154*x33 + x159 + x166;
auto x173 = -V{0.0833333333333333}*cell[18] + V{0.833333333333333}*cell[8] - V{0.0833333333333333}*cell[9] + x171 + x172;
auto x174 = V{0.00115740740740741}*x23*x29 + V{0.00115740740740741}*x26*x31;
auto x175 = V{0.416666666666667}*cell[18] - V{0.166666666666667}*cell[8] + V{0.416666666666667}*cell[9] - x171 + x172 + x174*x49 + x174*x53;
auto x176 = V{0.0833333333333333}*x27*x29 + V{0.0833333333333333}*x30*x31;
auto x177 = V{0.00115740740740741}*x23*x29 + V{0.00115740740740741}*x26*x31;
auto x178 = V{0.0833333333333333}*cell[10] + V{0.0833333333333333}*cell[1] - x127 - x128 - x143 - x144 - x145 - x146 + x177*x49 + x177*x53 + V{0.0555555555555555};
auto x0 = -x19*(-V{0.166666666666667}*x117*x51 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[9] + x21 + x23*(V{0.24999975}*cell[12] + V{0.4999995}*cell[14] + V{0.24999975}*cell[15] + V{0.24999975}*cell[16] + V{0.4999995}*cell[2] + V{0.24999975}*cell[3] + V{0.4999995}*cell[4] + V{0.24999975}*cell[6] + V{0.24999975}*cell[7] + V{0.4999995}*cell[9] + x24) + x26*(V{0.24999975}*cell[11] + V{0.24999975}*cell[13] + V{0.24999975}*cell[14] + V{0.4999995}*cell[16] + V{0.4999995}*cell[18] + V{0.24999975}*cell[2] + V{0.4999995}*cell[3] + V{0.24999975}*cell[4] + V{0.24999975}*cell[5] + V{0.4999995}*cell[6] + x24) - x32*x33 - x32*x34 - x32*x35 - x36*x49 - x36*x53 - x54 + V{0.833332833333333});
auto x1 = -x19*(-V{0.0277777777777778}*x117*x97 + V{0.0555555555555556}) - x20*(x137 + x98);
auto x2 = -x19*(-V{0.0277777777777778}*x117*x92 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] - x103 - x129 - x138 - x139 - x141 + V{0.0833333333333333}*x142*x34 - x148);
auto x3 = -x19*(-V{0.0277777777777778}*x117*x57 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] - x127 - x128 - x136 + V{0.0833333333333333}*x142*x35 - x148 - x149 - x150 - x151 - x58 + V{0.0555555555555555});
auto x4 = -x19*(-V{0.0138888888888889}*x117*x80 + V{0.0277777777777778}) - x20*(x104 + x163);
auto x5 = -x19*(-V{0.0138888888888889}*x112*x117 + V{0.0277777777777778}) - x20*(x164 + x36*(x111 + x82));
auto x6 = -x19*(-V{0.0138888888888889}*x117*x62 + V{0.0277777777777778}) - x20*(x168 + x63);
auto x7 = -x19*(-V{0.0138888888888889}*x117*x77 + V{0.0277777777777778}) - x20*(x169 + x36*(x70 + x76));
auto x8 = -x19*(-V{0.0138888888888889}*x117*x53 + V{0.0277777777777778}) - x20*(-x170*x49 + x173 + x53*(V{0.00810185185185185}*x23*x29 + V{0.00810185185185185}*x26*x31));
auto x9 = -x19*(-V{0.0138888888888889}*x117*x88 + V{0.0277777777777778}) - x20*(x175 + x36*(x113 + x87));
auto x10 = -x19*(V{0.0277777777777778}*x117*x99 + V{0.0555555555555556}) - x20*(x100 + x137);
auto x11 = -x19*(V{0.0277777777777778}*x115*x117 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] + x116 - x120 - x121 - x122 - x123 + x124 + x125 - x138 - x139 - x141 + x176*x34 + x178);
auto x12 = -x19*(V{0.0277777777777778}*x117*x94 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] - x130 - x131 - x132 - x133 + x134 + x135 - x149 - x150 - x151 + x176*x35 + x178 + x95);
auto x13 = -x19*(V{0.0138888888888889}*x106*x117 + V{0.0277777777777778}) - x20*(x107 + x163);
auto x14 = -x19*(-V{0.0138888888888889}*x117*x84 + V{0.0277777777777778}) - x20*(x108 + x164);
auto x15 = -x19*(V{0.0138888888888889}*x117*x65 + V{0.0277777777777778}) - x20*(x168 + x66);
auto x16 = -x19*(-V{0.0138888888888889}*x117*x71 + V{0.0277777777777778}) - x20*(x169 + x72);
auto x17 = -x19*(V{0.0138888888888889}*x117*x49 + V{0.0277777777777778}) - x20*(-x170*x53 + x173 - x49*(V{0.0196759259259259}*x23*x29 + V{0.0196759259259259}*x26*x31));
auto x18 = -x19*(-V{0.0138888888888889}*x114*x117 + V{0.0277777777777778}) - x20*(x114*x36 + x175);
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
return { -V{0.5}*x117, x33 + x34 + x35 };
}
};

}

}
