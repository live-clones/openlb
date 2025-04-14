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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<0, 1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<0, 1, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{2}*cell[18];
auto x22 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{-1};
auto x23 = V{1} / (x22);
auto x24 = V{0.24999975}*cell[0] + V{0.24999975}*cell[10] + V{0.4999995}*cell[18] + V{0.24999975}*cell[1] + V{0.24999975};
auto x25 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{1});
auto x26 = cell[0] + cell[10] + cell[1] + x21 + V{1};
auto x27 = V{2}*cell[11] + cell[12] + V{2}*cell[13] + cell[15] + cell[16] + V{2}*cell[17] + cell[3] + V{2}*cell[5] + cell[6] + cell[7] + x26;
auto x28 = cell[11] + cell[13] + cell[14] + V{2}*cell[16] + cell[2] + V{2}*cell[3] + cell[4] + cell[5] + V{2}*cell[6] + V{2}*cell[8] + x26;
auto x29 = -V{0.0138888888888889}*x23*x28 + V{0.0138888888888889}*x25*x27;
auto x30 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x31 = -x30;
auto x32 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x33 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x32;
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
auto x55 = V{1} - x41;
auto x56 = -x39;
auto x57 = x35 + x56;
auto x58 = x54 + x55 + x57;
auto x59 = x53*x58;
auto x60 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x61 = V{4.5}*(x60*x60);
auto x62 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x63 = -x37;
auto x64 = x55 + x63;
auto x65 = x62 + x64;
auto x66 = x57 + x61 + x65;
auto x67 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x68 = V{4.5}*(x67*x67);
auto x69 = x30 + x56;
auto x70 = x65 + x68 + x69;
auto x71 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x72 = V{4.5}*(x71*x71);
auto x73 = x30 + x57 + x64 + x72;
auto x74 = -V{0.0555555555555556}*x23*x28 + V{0.0555555555555556}*x25*x27;
auto x75 = V{3}*x40;
auto x76 = x63 + x69 + x75 + V{1};
auto x77 = V{5.55111512312578e-17}*x25*x27 + V{5.55111512312578e-17}*x28*x47;
auto x78 = V{1.11022302462516e-16}*x25*x27 + V{1.11022302462516e-16}*x28*x47;
auto x79 = x44 - x61 + x62;
auto x80 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x32;
auto x81 = -x80;
auto x82 = x50 + x62;
auto x83 = x82 - V{4.5}*x81*x81;
auto x84 = -x62;
auto x85 = -V{4.5}*x80*x80;
auto x86 = x44 + x84 + x85;
auto x87 = x35 + x42 - x54;
auto x88 = x53*x87;
auto x89 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x90 = -x89;
auto x91 = x43 + x62;
auto x92 = x31 + x91;
auto x93 = x92 - V{4.5}*x90*x90;
auto x94 = V{8.32667268468867e-17}*cell[0] + V{8.32667268468867e-17}*cell[10] + V{1.66533453693773e-16}*cell[18] + V{8.32667268468867e-17}*cell[1] + V{8.32667268468867e-17};
auto x95 = V{3}*x38;
auto x96 = x65 + x95;
auto x97 = x37 + V{-1};
auto x98 = x41 + x62 - x95 + x97;
auto x99 = V{1.66533453693773e-16}*cell[10] + V{4.44089209850063e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + x23*(V{8.32667268468867e-17}*cell[11] + V{8.32667268468867e-17}*cell[13] + V{8.32667268468867e-17}*cell[14] + V{1.66533453693773e-16}*cell[16] + V{8.32667268468867e-17}*cell[2] + V{1.66533453693773e-16}*cell[3] + V{8.32667268468867e-17}*cell[4] + V{8.32667268468867e-17}*cell[5] + V{1.66533453693773e-16}*cell[6] + V{1.66533453693773e-16}*cell[8] + x94) - x25*(V{1.66533453693773e-16}*cell[11] + V{8.32667268468867e-17}*cell[12] + V{1.66533453693773e-16}*cell[13] + V{8.32667268468867e-17}*cell[15] + V{8.32667268468867e-17}*cell[16] + V{1.66533453693773e-16}*cell[17] + V{8.32667268468867e-17}*cell[3] + V{1.66533453693773e-16}*cell[5] + V{8.32667268468867e-17}*cell[6] + V{8.32667268468867e-17}*cell[7] + x94) - x38*(V{8.32667268468867e-17}*x25*x27 + V{8.32667268468867e-17}*x28*x47) - x45*(-V{0.0277777777777778}*x23*x28 + V{0.0277777777777778}*x25*x27) - x46 - x52*(-V{3.08395284618099e-18}*x23*x28 + V{3.08395284618099e-18}*x25*x27) + x53*x96 - x53*x98 + V{4.44089209850063e-16};
auto x100 = x25*(V{2.22044604925031e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{3.33066907387547e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{4.44089209850063e-16}*cell[17] + V{2.22044604925031e-16}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{3.33066907387547e-16}*cell[5] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[8] + x29*x66 - x29*x79 - x29*x83 - x29*x86 - x36*x77 - x40*x78 + x53*x70 + x53*x73 - x53*x93 + x59 + x74*x76 - x88 + x99);
auto x101 = x53*x76;
auto x102 = x30 - x68 + x91;
auto x103 = -V{4.5}*x89*x89;
auto x104 = x103 + x30 + x43 + x84;
auto x105 = x30 + x39 - x75 + x97;
auto x106 = x105*x53;
auto x107 = x30 + x44 - x72;
auto x108 = x23*(V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.71844785465692e-16}*cell[13] + V{5.27355936696949e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.71844785465692e-16}*cell[4] + V{5.27355936696949e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + x101 - x102*x29 - x104*x29 - x106 - x107*x53 + x29*x70 - x29*x93 - x36*x78 - x40*x77 - x53*x79 - x53*x86 - x74*x87 + x99);
auto x109 = x100 - x108;
auto x110 = V{0.0277777777777778}*x23*x28 - V{0.0277777777777778}*x25*x27;
auto x111 = V{0.00231481481481481}*x23*x28 - V{0.00231481481481481}*x25*x27;
auto x112 = V{0.0416666666666667}*x25*x27;
auto x113 = x112 - V{0.0416666666666667}*x23*x28;
auto x114 = x34 + x51;
auto x115 = V{0.0833333333333333}*x25*x27;
auto x116 = V{0.0833333333333334}*cell[13];
auto x117 = V{0.0833333333333334}*cell[14];
auto x118 = V{0.0833333333333334}*cell[4];
auto x119 = V{0.0833333333333334}*cell[5];
auto x120 = V{0.0833333333333333}*cell[12];
auto x121 = V{0.0833333333333333}*cell[3];
auto x122 = V{6.93889390390723e-18}*cell[0] + V{6.93889390390723e-18}*cell[10] + V{1.38777878078145e-17}*cell[18] + V{6.93889390390723e-18}*cell[1] + V{6.93889390390723e-18};
auto x123 = x23*(V{6.93889390390723e-18}*cell[11] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[14] + V{1.38777878078145e-17}*cell[16] + V{6.93889390390723e-18}*cell[2] + V{1.38777878078145e-17}*cell[3] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[5] + V{1.38777878078145e-17}*cell[6] + V{1.38777878078145e-17}*cell[8] + x122);
auto x124 = x25*(V{1.38777878078145e-17}*cell[11] + V{6.93889390390723e-18}*cell[12] + V{1.38777878078145e-17}*cell[13] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[16] + V{1.38777878078145e-17}*cell[17] + V{6.93889390390723e-18}*cell[3] + V{1.38777878078145e-17}*cell[5] + V{6.93889390390723e-18}*cell[6] + V{6.93889390390723e-18}*cell[7] + x122);
auto x125 = -x124;
auto x126 = x116 + x117 + x118 + x119 - x120 - x121 + x123 + x125 + V{-0.0555555555555555};
auto x127 = V{0.0833333333333334}*cell[15];
auto x128 = V{0.0833333333333334}*cell[16];
auto x129 = V{0.0833333333333334}*cell[6];
auto x130 = V{0.0833333333333334}*cell[7];
auto x131 = V{0.0833333333333333}*cell[11];
auto x132 = V{0.0833333333333333}*cell[2];
auto x133 = x127 + x128 + x129 + x130 - x131 - x132;
auto x134 = V{0.166666666666667}*cell[10] - V{0.166666666666667}*cell[17] - V{0.333333333333333}*cell[18] + V{0.166666666666667}*cell[1] - V{0.166666666666667}*cell[8] - x111*x114 + x111*x45 + x113*x36 + x113*x40 + x126 + x133 - x38*(x115 - V{0.0833333333333333}*x23*x28);
auto x135 = V{0.166666666666667}*cell[11];
auto x136 = V{0.166666666666667}*cell[2];
auto x137 = x115 + V{0.0833333333333333}*x28*x47;
auto x138 = x112 + V{0.0416666666666667}*x28*x47;
auto x139 = x138*x36;
auto x140 = V{0.0833333333333334}*cell[17];
auto x141 = V{0.0833333333333334}*cell[8];
auto x142 = V{0.166666666666667}*cell[18];
auto x143 = -V{0.00115740740740741}*x23*x28 + V{0.00115740740740741}*x25*x27;
auto x144 = x143*x45;
auto x145 = x138*x38;
auto x146 = V{0.0833333333333333}*cell[10] + V{0.0833333333333333}*cell[1] - x123 + x124 - x140 - x141 - x142 + x143*x52 - x144 - x145 + V{0.0555555555555555};
auto x147 = V{0.166666666666667}*cell[12];
auto x148 = V{0.166666666666667}*cell[3];
auto x149 = x138*x40;
auto x150 = V{0.0138888888888889}*x23*x28 - V{0.0138888888888889}*x25*x27;
auto x151 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x152 = x23*x28 - x25*x27;
auto x153 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x152;
auto x154 = x151*x153;
auto x155 = -V{0.0208333333333333}*x23*x28 + V{0.0208333333333333}*x25*x27;
auto x156 = V{0.0208333125}*cell[0] + V{0.0208333125}*cell[10] + V{0.041666625}*cell[18] + V{0.0208333125}*cell[1] + V{0.0208333125};
auto x157 = x23*(V{0.0208333125}*cell[11] + V{0.0208333125}*cell[13] + V{0.0208333125}*cell[14] + V{0.041666625}*cell[16] + V{0.0208333125}*cell[2] + V{0.041666625}*cell[3] + V{0.0208333125}*cell[4] + V{0.0208333125}*cell[5] + V{0.041666625}*cell[6] + V{0.041666625}*cell[8] + x156);
auto x158 = -x25*(V{0.041666625}*cell[11] + V{0.0208333125}*cell[12] + V{0.041666625}*cell[13] + V{0.0208333125}*cell[15] + V{0.0208333125}*cell[16] + V{0.041666625}*cell[17] + V{0.0208333125}*cell[3] + V{0.041666625}*cell[5] + V{0.0208333125}*cell[6] + V{0.0208333125}*cell[7] + x156);
auto x159 = -V{0.0416666666666667}*x23*x28 + V{0.0416666666666667}*x25*x27;
auto x160 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] + x157 + x158 - x159*x40 + V{0.0138888472222222};
auto x161 = V{0.000578703703703704}*x23*x28 - V{0.000578703703703704}*x25*x27;
auto x162 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0833333333333334}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + x114*x161 - x159*x38 - x161*x45;
auto x163 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x155*x36 + x160 + x162;
auto x164 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] + x154 + x163;
auto x165 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] - x154 + x163;
auto x166 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x153;
auto x167 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x159*x36;
auto x168 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x155*x40 + x157 + x158 + x162 + x167 + V{0.0138888472222222};
auto x169 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] + x166 + x168;
auto x170 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] - x166 + x168;
auto x171 = V{0.00115740740740741}*x23*x28 - V{0.00115740740740741}*x25*x27;
auto x172 = -x152;
auto x173 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x151*x172;
auto x174 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x155*x38 + x160 + x167;
auto x175 = V{0.416666666666667}*cell[17] - V{0.166666666666667}*cell[18] + V{0.416666666666667}*cell[8] - x114*x171 + x171*x45 - x173 + x174;
auto x176 = V{0.00578703703703704}*x23*x28 - V{0.00578703703703704}*x25*x27;
auto x177 = -V{0.0833333333333333}*cell[17] + V{0.833333333333333}*cell[18] - V{0.0833333333333333}*cell[8] + x173 + x174;
auto x178 = -x100 + x108;
auto x179 = x25*x27 + x28*x47;
auto x180 = -V{0.0833333333333333}*cell[10] - V{0.0833333333333333}*cell[1] + x140 + x141 + x142 + x144 + x145 - V{0.00115740740740741}*x172*x52;
auto x0 = -x19*(V{0.166666666666667}*x109*x43 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + x21 + x23*(V{0.24999975}*cell[11] + V{0.24999975}*cell[13] + V{0.24999975}*cell[14] + V{0.4999995}*cell[16] + V{0.24999975}*cell[2] + V{0.4999995}*cell[3] + V{0.24999975}*cell[4] + V{0.24999975}*cell[5] + V{0.4999995}*cell[6] + V{0.4999995}*cell[8] + x24) - x25*(V{0.4999995}*cell[11] + V{0.24999975}*cell[12] + V{0.4999995}*cell[13] + V{0.24999975}*cell[15] + V{0.24999975}*cell[16] + V{0.4999995}*cell[17] + V{0.24999975}*cell[3] + V{0.4999995}*cell[5] + V{0.24999975}*cell[6] + V{0.24999975}*cell[7] + x24) + x29*x45 - x29*x52 - x36*x48 - x38*x48 - x40*x48 + x46 + V{0.833332833333333});
auto x1 = -x19*(V{0.0277777777777778}*x109*x98 + V{0.0555555555555556}) - x20*(x110*x98 + x134);
auto x2 = -x19*(V{0.0277777777777778}*x105*x109 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] + x106 - x116 - x117 - x118 - x119 + x120 + x121 - x135 - x136 + x137*x40 - x139 + x146);
auto x3 = -x19*(V{0.0277777777777778}*x109*x87 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] - x127 - x128 - x129 - x130 + x131 + x132 + x137*x36 + x146 - x147 - x148 - x149 + x88);
auto x4 = -x19*(V{0.0138888888888889}*x102*x109 + V{0.0277777777777778}) - x20*(x102*x150 + x164);
auto x5 = -x19*(V{0.0138888888888889}*x109*x93 + V{0.0277777777777778}) - x20*(x150*(x103 + x92) + x165);
auto x6 = -x19*(V{0.0138888888888889}*x109*x79 + V{0.0277777777777778}) - x20*(x150*x79 + x169);
auto x7 = -x19*(V{0.0138888888888889}*x109*x83 + V{0.0277777777777778}) - x20*(x150*(x82 + x85) + x170);
auto x8 = -x19*(V{0.0138888888888889}*x107*x109 + V{0.0277777777777778}) - x20*(x107*x150 + x175);
auto x9 = -x19*(V{0.0138888888888889}*x109*x52 + V{0.0277777777777778}) - x20*(x114*(V{0.0196759259259259}*x23*x28 - V{0.0196759259259259}*x25*x27) - x176*x45 + x177);
auto x10 = -x19*(V{0.0277777777777778}*x178*x96 + V{0.0555555555555556}) - x20*(-x110*x96 + x134);
auto x11 = -x19*(V{0.0277777777777778}*x178*x76 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] - x101 - x126 - x135 - x136 - x139 + V{0.0833333333333333}*x179*x40 - x180);
auto x12 = -x19*(V{0.0277777777777778}*x178*x58 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] - x123 - x125 - x133 - x147 - x148 - x149 + V{0.0833333333333333}*x179*x36 - x180 - x59 + V{0.0555555555555555});
auto x13 = -x19*(V{0.0138888888888889}*x178*x70 + V{0.0277777777777778}) - x20*(-x150*x70 + x164);
auto x14 = -x19*(V{0.0138888888888889}*x104*x109 + V{0.0277777777777778}) - x20*(x104*x150 + x165);
auto x15 = -x19*(V{0.0138888888888889}*x178*x66 + V{0.0277777777777778}) - x20*(-x150*x66 + x169);
auto x16 = -x19*(V{0.0138888888888889}*x109*x86 + V{0.0277777777777778}) - x20*(x150*x86 + x170);
auto x17 = -x19*(V{0.0138888888888889}*x178*x73 + V{0.0277777777777778}) - x20*(-x150*x73 + x175);
auto x18 = -x19*(V{0.0138888888888889}*x109*x45 + V{0.0277777777777778}) - x20*(x114*x176 + x177 + x45*(V{0.00810185185185185}*x23*x28 - V{0.00810185185185185}*x25*x27));
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
return { V{0.5}*x109, x36 + x38 + x40 };
}
};

}

}
