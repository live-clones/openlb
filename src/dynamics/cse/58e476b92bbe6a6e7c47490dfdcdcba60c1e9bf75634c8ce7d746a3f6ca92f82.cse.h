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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<2, -1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<2, -1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{2}*cell[5];
auto x22 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{-1};
auto x23 = V{1} / (x22);
auto x24 = V{0.24999975}*cell[0] + V{0.24999975}*cell[12] + V{0.24999975}*cell[3] + V{0.4999995}*cell[5] + V{0.24999975};
auto x25 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{1});
auto x26 = cell[0] + cell[12] + cell[3] + x21 + V{1};
auto x27 = cell[10] + V{2}*cell[11] + V{2}*cell[13] + cell[15] + cell[16] + V{2}*cell[17] + V{2}*cell[18] + cell[1] + cell[6] + cell[7] + x26;
auto x28 = cell[11] + cell[17] + cell[18] + V{2}*cell[1] + cell[2] + V{2}*cell[4] + V{2}*cell[6] + V{2}*cell[7] + cell[8] + cell[9] + x26;
auto x29 = -V{0.0138888888888889}*x23*x28 + V{0.0138888888888889}*x25*x27;
auto x30 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x31 = -x30;
auto x32 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x33 = -x32;
auto x34 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x35 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x36 = V{1.5}*x35;
auto x37 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x38 = V{1.5}*x37;
auto x39 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x40 = V{1.5}*x39;
auto x41 = x38 + x40 + V{-1};
auto x42 = x36 + x41;
auto x43 = x34 + x42;
auto x44 = x33 + x43;
auto x45 = x44 - V{4.5}*x31*x31;
auto x46 = x42*(-V{0.166666666666667}*x23*x28 + V{0.166666666666667}*x25*x27);
auto x47 = -1/x22;
auto x48 = V{0.25}*x25*x27 + V{0.25}*x28*x47;
auto x49 = -V{4.5}*x30*x30;
auto x50 = -x34 + x42;
auto x51 = x32 + x49 + x50;
auto x52 = V{1.11022302462516e-16}*x25*x27;
auto x53 = V{1.11022302462516e-16}*x28*x47 + x52;
auto x54 = V{8.32667268468867e-17}*x25*x27;
auto x55 = V{8.32667268468867e-17}*x28*x47 + x54;
auto x56 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x57 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x56;
auto x58 = -x57;
auto x59 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x60 = -x59;
auto x61 = x43 + x60;
auto x62 = x61 - V{4.5}*x58*x58;
auto x63 = V{5.55111512312578e-17}*x25*x27;
auto x64 = -V{0.0277777777777778}*x23*x28 + V{0.0277777777777778}*x25*x27;
auto x65 = -V{0.0277777777777778}*x23*x28 + V{0.0277777777777778}*x25*x27;
auto x66 = V{3}*x35;
auto x67 = x41 + x59 - x66;
auto x68 = -x65*x67;
auto x69 = V{1} - x40;
auto x70 = -x38;
auto x71 = x59 + x70;
auto x72 = x66 + x69 + x71;
auto x73 = x65*x72;
auto x74 = V{1.66533453693773e-16}*cell[1];
auto x75 = V{8.32667268468867e-17}*cell[0] + V{8.32667268468867e-17}*cell[12] + V{8.32667268468867e-17}*cell[3] + V{1.66533453693773e-16}*cell[5] + V{8.32667268468867e-17};
auto x76 = V{1.66533453693773e-16}*cell[11];
auto x77 = V{1.11022302462516e-16}*cell[12] + V{3.33066907387547e-16}*cell[13] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{6.66133814775094e-16}*cell[5] + x23*(V{8.32667268468867e-17}*cell[11] + V{8.32667268468867e-17}*cell[17] + V{8.32667268468867e-17}*cell[18] + V{8.32667268468867e-17}*cell[2] + V{1.66533453693773e-16}*cell[4] + V{1.66533453693773e-16}*cell[6] + V{1.66533453693773e-16}*cell[7] + V{8.32667268468867e-17}*cell[8] + V{8.32667268468867e-17}*cell[9] + x74 + x75) - x25*(V{8.32667268468867e-17}*cell[10] + V{1.66533453693773e-16}*cell[13] + V{8.32667268468867e-17}*cell[15] + V{8.32667268468867e-17}*cell[16] + V{1.66533453693773e-16}*cell[17] + V{1.66533453693773e-16}*cell[18] + V{8.32667268468867e-17}*cell[1] + V{8.32667268468867e-17}*cell[6] + V{8.32667268468867e-17}*cell[7] + x75 + x76) - x46 - x51*(-V{4.62592926927149e-18}*x23*x28 + V{4.62592926927149e-18}*x25*x27) + x68 + x73 + V{4.44089209850063e-16};
auto x78 = -x35*(V{5.55111512312578e-17}*x28*x47 + x63) - x45*x64 + x77;
auto x79 = V{3}*x37;
auto x80 = -x36;
auto x81 = x69 + x80;
auto x82 = x34 + x81;
auto x83 = x79 + x82;
auto x84 = x65*x83;
auto x85 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x86 = V{4.5}*(x85*x85);
auto x87 = x43 + x59 - x86;
auto x88 = -x29*x87;
auto x89 = x71 + x82 + x86;
auto x90 = x29*x89;
auto x91 = -V{4.5}*x57*x57;
auto x92 = x50 + x59 + x91;
auto x93 = -x29*x92;
auto x94 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x95 = V{4.5}*(x94*x94);
auto x96 = x32 + x70;
auto x97 = x82 + x95 + x96;
auto x98 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x99 = V{4.5}*(x98*x98);
auto x100 = x32 + x71 + x81 + x99;
auto x101 = -V{0.0555555555555556}*x23*x28 + V{0.0555555555555556}*x25*x27;
auto x102 = V{3}*x39;
auto x103 = x102 + x80 + x96 + V{1};
auto x104 = x36 + V{-1};
auto x105 = x104 + x34 + x40 - x79;
auto x106 = x105*x65;
auto x107 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x56;
auto x108 = -V{4.5}*x107*x107;
auto x109 = x42 + x59;
auto x110 = x108 + x109 + x33;
auto x111 = V{1.66533453693773e-16}*cell[10] + V{2.22044604925031e-16}*cell[11] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{4.44089209850063e-16}*cell[17] + V{2.22044604925031e-16}*cell[18] + V{2.22044604925031e-16}*cell[2] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[8] + V{2.22044604925031e-16}*cell[9] + x100*x65 + x101*x103 - x106 - x110*x65 + x65*x97 + x74 + x84 + x88 + x90 + x93;
auto x112 = -x107;
auto x113 = x32 + x42 + x60;
auto x114 = x113 - V{4.5}*x112*x112;
auto x115 = x103*x65;
auto x116 = x109 + x32 - x99;
auto x117 = -x116*x29;
auto x118 = x100*x29;
auto x119 = -x110*x29;
auto x120 = -x102 + x104 + x32 + x38;
auto x121 = x120*x65;
auto x122 = x32 + x43 - x95;
auto x123 = V{2.22044604925031e-16}*cell[10] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] - x101*x105 + x115 + x117 + x118 + x119 - x121 - x122*x65 - x65*x87 + x76;
auto x124 = -x23*(-x114*x29 + x123 - x37*x53 - x39*x55 - x62*x65 + x78) + x25*(x111 - x29*x62 - x37*x55 - x39*x53 + x78);
auto x125 = V{0.166666666666667}*cell[10];
auto x126 = V{0.166666666666667}*cell[1];
auto x127 = V{0.0833333333333334}*cell[15];
auto x128 = V{0.0833333333333334}*cell[16];
auto x129 = V{0.0833333333333334}*cell[6];
auto x130 = V{0.0833333333333334}*cell[7];
auto x131 = V{0.0833333333333333}*x25*x27 + V{0.0833333333333333}*x28*x47;
auto x132 = V{0.0416666666666667}*x25*x27 + V{0.0416666666666667}*x28*x47;
auto x133 = x132*x39;
auto x134 = V{0.0833333333333334}*cell[13];
auto x135 = V{0.0833333333333334}*cell[4];
auto x136 = V{0.166666666666667}*cell[5];
auto x137 = V{6.93889390390723e-18}*cell[0] + V{6.93889390390723e-18}*cell[12] + V{6.93889390390723e-18}*cell[3] + V{1.38777878078145e-17}*cell[5] + V{6.93889390390723e-18};
auto x138 = V{6.93889390390723e-18}*cell[10] + V{1.38777878078145e-17}*cell[11] + V{1.38777878078145e-17}*cell[13] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[16] + V{1.38777878078145e-17}*cell[17] + V{1.38777878078145e-17}*cell[18] + V{6.93889390390723e-18}*cell[1] + V{6.93889390390723e-18}*cell[6] + V{6.93889390390723e-18}*cell[7] + x137;
auto x139 = x23*(V{6.93889390390723e-18}*cell[11] + V{6.93889390390723e-18}*cell[17] + V{6.93889390390723e-18}*cell[18] + V{1.38777878078145e-17}*cell[1] + V{6.93889390390723e-18}*cell[2] + V{1.38777878078145e-17}*cell[4] + V{1.38777878078145e-17}*cell[6] + V{1.38777878078145e-17}*cell[7] + V{6.93889390390723e-18}*cell[8] + V{6.93889390390723e-18}*cell[9] + x137);
auto x140 = -V{0.00115740740740741}*x23*x28 + V{0.00115740740740741}*x25*x27;
auto x141 = x140*x45;
auto x142 = x132*x35;
auto x143 = V{0.0833333333333333}*cell[12] + V{0.0833333333333333}*cell[3] - x134 - x135 - x136 + x138*x25 - x139 + x140*x51 - x141 - x142 + V{0.0555555555555555};
auto x144 = V{0.166666666666667}*cell[11];
auto x145 = V{0.166666666666667}*cell[2];
auto x146 = V{0.0833333333333334}*cell[17];
auto x147 = V{0.0833333333333334}*cell[18];
auto x148 = V{0.0833333333333334}*cell[8];
auto x149 = V{0.0833333333333334}*cell[9];
auto x150 = x132*x37;
auto x151 = x25*x27;
auto x152 = x151 - x23*x28;
auto x153 = x151 + x28*x47;
auto x154 = -x138*x25;
auto x155 = -V{0.0833333333333333}*cell[11] - V{0.0833333333333333}*cell[2] + x127 + x128 + x129 + x130 + x133 + x139 + x154 + V{-0.0555555555555555};
auto x156 = -V{0.0833333333333333}*cell[10] - V{0.0833333333333333}*cell[1] + x146 + x147 + x148 + x149 + x150;
auto x157 = V{0.166666666666667}*cell[12] - V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[3] - V{0.166666666666667}*cell[4] - V{0.333333333333333}*cell[5] - V{0.00231481481481481}*x152*x45 - V{0.0833333333333333}*x153*x35 + x155 + x156 + x51*(-V{0.00231481481481481}*x23*x28 + V{0.00231481481481481}*x25*x27);
auto x158 = -V{0.00115740740740741}*x23*x28 + V{0.00115740740740741}*x25*x27;
auto x159 = x44 + x49;
auto x160 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x152;
auto x161 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x160;
auto x162 = -V{0.0208333333333333}*x23*x28 + V{0.0208333333333333}*x25*x27;
auto x163 = V{0.0208333125}*cell[0] + V{0.0208333125}*cell[12] + V{0.0208333125}*cell[3] + V{0.041666625}*cell[5] + V{0.0208333125};
auto x164 = x23*(V{0.0208333125}*cell[11] + V{0.0208333125}*cell[17] + V{0.0208333125}*cell[18] + V{0.041666625}*cell[1] + V{0.0208333125}*cell[2] + V{0.041666625}*cell[4] + V{0.041666625}*cell[6] + V{0.041666625}*cell[7] + V{0.0208333125}*cell[8] + V{0.0208333125}*cell[9] + x163);
auto x165 = -x25*(V{0.0208333125}*cell[10] + V{0.041666625}*cell[11] + V{0.041666625}*cell[13] + V{0.0208333125}*cell[15] + V{0.0208333125}*cell[16] + V{0.041666625}*cell[17] + V{0.041666625}*cell[18] + V{0.0208333125}*cell[1] + V{0.0208333125}*cell[6] + V{0.0208333125}*cell[7] + x163);
auto x166 = -V{0.0416666666666667}*x23*x28 + V{0.0416666666666667}*x25*x27;
auto x167 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x164 + x165 - x166*x37 + V{0.0138888472222222};
auto x168 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] - x166*x39;
auto x169 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x162*x35 + x167 + x168;
auto x170 = V{0.416666666666667}*cell[13] + V{0.416666666666667}*cell[4] - V{0.166666666666667}*cell[5] - x158*x159 + x158*x51 - x161 + x169;
auto x171 = -V{0.00578703703703704}*x23*x28 + V{0.00578703703703704}*x25*x27;
auto x172 = -V{0.0833333333333333}*cell[13] - V{0.0833333333333333}*cell[4] + V{0.833333333333333}*cell[5] + x161 + x169;
auto x173 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x160;
auto x174 = -V{0.000578703703703704}*x23*x28 + V{0.000578703703703704}*x25*x27;
auto x175 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0833333333333334}*cell[5] + x159*x174 - x166*x35 - x174*x51;
auto x176 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x162*x39 + x167 + x175;
auto x177 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] - x173 + x176;
auto x178 = x61 + x91;
auto x179 = -x178*x29;
auto x180 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] + x173 + x176;
auto x181 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x152;
auto x182 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x162*x37 + x164 + x165 + x168 + x175 + V{0.0138888472222222};
auto x183 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x181 + x182;
auto x184 = -x29*(x108 + x113);
auto x185 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x181 + x182;
auto x186 = -V{1.11022302462516e-16}*x23*x28 + x52;
auto x187 = -V{8.32667268468867e-17}*x23*x28 + x54;
auto x188 = -x159*x64 - x35*(-V{5.55111512312578e-17}*x23*x28 + x63) + x77;
auto x189 = -x23*(x123 - x178*x65 + x184 - x186*x37 - x187*x39 + x188) + x25*(x111 + x179 - x186*x39 - x187*x37 + x188);
auto x190 = -V{0.0833333333333333}*cell[12] - V{0.0833333333333333}*cell[3] + x134 + x135 + x136 + x141 + x142 - V{0.00115740740740741}*x152*x51;
auto x0 = -x19*(V{0.166666666666667}*x124*x42 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x21 + x23*(V{0.24999975}*cell[11] + V{0.24999975}*cell[17] + V{0.24999975}*cell[18] + V{0.4999995}*cell[1] + V{0.24999975}*cell[2] + V{0.4999995}*cell[4] + V{0.4999995}*cell[6] + V{0.4999995}*cell[7] + V{0.24999975}*cell[8] + V{0.24999975}*cell[9] + x24) - x25*(V{0.24999975}*cell[10] + V{0.4999995}*cell[11] + V{0.4999995}*cell[13] + V{0.24999975}*cell[15] + V{0.24999975}*cell[16] + V{0.4999995}*cell[17] + V{0.4999995}*cell[18] + V{0.24999975}*cell[1] + V{0.24999975}*cell[6] + V{0.24999975}*cell[7] + x24) + x29*x45 - x29*x51 - x35*x48 - x37*x48 - x39*x48 + x46 + V{0.833332833333333});
auto x1 = -x19*(V{0.0277777777777778}*x105*x124 + V{0.0555555555555556}) + x20*(V{0.0833333333333333}*cell[11] + V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.0833333333333333}*cell[2] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x106 - x125 - x126 - x127 - x128 - x129 - x130 + x131*x37 - x133 + x143);
auto x2 = -x19*(V{0.0277777777777778}*x120*x124 + V{0.0555555555555556}) + x20*(V{0.0833333333333333}*cell[10] + V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] + V{0.0833333333333333}*cell[1] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] + x121 + x131*x39 + x143 - x144 - x145 - x146 - x147 - x148 - x149 - x150);
auto x3 = -x19*(V{0.0277777777777778}*x124*x67 + V{0.0555555555555556}) + x20*(-x157 - x68);
auto x4 = -x19*(V{0.0138888888888889}*x122*x124 + V{0.0277777777777778}) - x20*(-x122*x29 + x170);
auto x5 = -x19*(V{0.0138888888888889}*x124*x45 + V{0.0277777777777778}) - x20*(-x159*(-V{0.00810185185185185}*x23*x28 + V{0.00810185185185185}*x25*x27) - x171*x51 + x172);
auto x6 = -x19*(V{0.0138888888888889}*x124*x87 + V{0.0277777777777778}) - x20*(x177 + x88);
auto x7 = -x19*(V{0.0138888888888889}*x124*x62 + V{0.0277777777777778}) - x20*(x179 + x180);
auto x8 = -x19*(V{0.0138888888888889}*x116*x124 + V{0.0277777777777778}) - x20*(x117 + x183);
auto x9 = -x19*(V{0.0138888888888889}*x114*x124 + V{0.0277777777777778}) - x20*(x184 + x185);
auto x10 = x19*(V{0.0277777777777778}*x189*x83 + V{-0.0555555555555556}) + x20*(V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] - x125 - x126 + V{0.0833333333333333}*x153*x37 - x155 - x190 - x84);
auto x11 = x19*(V{0.0277777777777778}*x103*x189 + V{-0.0555555555555556}) + x20*(V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] - x115 - x139 - x144 - x145 + V{0.0833333333333333}*x153*x39 - x154 - x156 - x190 + V{0.0555555555555555});
auto x12 = x19*(V{0.0277777777777778}*x189*x72 + V{-0.0555555555555556}) + x20*(-x157 - x73);
auto x13 = x19*(V{0.0138888888888889}*x189*x97 + V{-0.0277777777777778}) - x20*(x170 + x29*x97);
auto x14 = -x19*(V{0.0138888888888889}*x124*x51 + V{0.0277777777777778}) - x20*(x159*x171 + x172 - x51*(-V{0.0196759259259259}*x23*x28 + V{0.0196759259259259}*x25*x27));
auto x15 = x19*(V{0.0138888888888889}*x189*x89 + V{-0.0277777777777778}) - x20*(x177 + x90);
auto x16 = -x19*(V{0.0138888888888889}*x124*x92 + V{0.0277777777777778}) - x20*(x180 + x93);
auto x17 = x19*(V{0.0138888888888889}*x100*x189 + V{-0.0277777777777778}) - x20*(x118 + x183);
auto x18 = -x19*(V{0.0138888888888889}*x110*x124 + V{0.0277777777777778}) - x20*(x119 + x185);
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
return { V{0.5}*x189, x35 + x37 + x39 };
}
};

}

}
