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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<2, 1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<2, 1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{2}*cell[13];
auto x22 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{1});
auto x23 = V{0.24999975}*cell[0] + V{0.24999975}*cell[12] + V{0.4999995}*cell[13] + V{0.24999975}*cell[3] + V{0.24999975};
auto x24 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{1});
auto x25 = cell[0] + cell[12] + cell[3] + x21 + V{1};
auto x26 = V{2}*cell[10] + cell[11] + V{2}*cell[14] + V{2}*cell[15] + V{2}*cell[16] + cell[17] + cell[18] + cell[2] + cell[8] + cell[9] + x25;
auto x27 = cell[10] + V{2}*cell[11] + cell[15] + cell[16] + V{2}*cell[17] + V{2}*cell[18] + cell[1] + V{2}*cell[5] + cell[6] + cell[7] + x25;
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
auto x39 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x40 = V{4.5}*(x39*x39);
auto x41 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x42 = -x29;
auto x43 = V{1} - x31;
auto x44 = x42 + x43;
auto x45 = x41 + x44;
auto x46 = -x33;
auto x47 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x48 = x46 + x47;
auto x49 = x40 + x45 + x48;
auto x50 = x35 + x41;
auto x51 = -x40 + x47 + x50;
auto x52 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x53 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x52;
auto x54 = -x53;
auto x55 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x56 = -x55;
auto x57 = x50 + x56;
auto x58 = x57 - V{4.5}*x54*x54;
auto x59 = V{1.66533453693773e-16}*cell[11];
auto x60 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x61 = V{4.5}*(x60*x60);
auto x62 = x50 + x55 - x61;
auto x63 = -x38*x62;
auto x64 = x46 + x55;
auto x65 = x45 + x61 + x64;
auto x66 = x38*x65;
auto x67 = -x41;
auto x68 = -V{4.5}*x53*x53;
auto x69 = x35 + x55;
auto x70 = x67 + x68 + x69;
auto x71 = -x38*x70;
auto x72 = V{0.0277777777777778}*x22*x26 + V{0.0277777777777778}*x24*x27;
auto x73 = V{3}*x32;
auto x74 = x45 + x73;
auto x75 = x72*x74;
auto x76 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x77 = V{4.5}*(x76*x76);
auto x78 = x44 + x47 + x64 + x77;
auto x79 = V{0.0555555555555556}*x22*x26 + V{0.0555555555555556}*x24*x27;
auto x80 = V{3}*x30;
auto x81 = x42 + x48 + x80 + V{1};
auto x82 = V{1.11022302462516e-16}*x22*x26 + V{1.11022302462516e-16}*x24*x27;
auto x83 = V{8.32667268468867e-17}*x22*x26 + V{8.32667268468867e-17}*x24*x27;
auto x84 = x29 + V{-1};
auto x85 = x31 + x41 - x73 + x84;
auto x86 = x72*x85;
auto x87 = -x47;
auto x88 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x89 = -V{4.5}*x88*x88;
auto x90 = x50 + x87 + x89;
auto x91 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x52;
auto x92 = -V{4.5}*x91*x91;
auto x93 = x69 + x87 + x92;
auto x94 = V{3}*x28;
auto x95 = x34 + x55 - x94;
auto x96 = -x72*x95;
auto x97 = x43 + x64 + x94;
auto x98 = x72*x97;
auto x99 = V{1.66533453693773e-16}*cell[10];
auto x100 = V{8.32667268468867e-17}*cell[0] + V{8.32667268468867e-17}*cell[12] + V{1.66533453693773e-16}*cell[13] + V{8.32667268468867e-17}*cell[3] + V{8.32667268468867e-17};
auto x101 = V{1.11022302462516e-16}*cell[12] + V{6.66133814775094e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[5] - x22*(V{8.32667268468867e-17}*cell[11] + V{1.66533453693773e-16}*cell[14] + V{1.66533453693773e-16}*cell[15] + V{1.66533453693773e-16}*cell[16] + V{8.32667268468867e-17}*cell[17] + V{8.32667268468867e-17}*cell[18] + V{8.32667268468867e-17}*cell[2] + V{8.32667268468867e-17}*cell[8] + V{8.32667268468867e-17}*cell[9] + x100 + x99) - x24*(V{8.32667268468867e-17}*cell[10] + V{8.32667268468867e-17}*cell[15] + V{8.32667268468867e-17}*cell[16] + V{1.66533453693773e-16}*cell[17] + V{1.66533453693773e-16}*cell[18] + V{8.32667268468867e-17}*cell[1] + V{1.66533453693773e-16}*cell[5] + V{8.32667268468867e-17}*cell[6] + V{8.32667268468867e-17}*cell[7] + x100 + x59) - x28*(V{5.55111512312578e-17}*x22*x26 + V{5.55111512312578e-17}*x24*x27) - x36 + x49*(V{0.0277777777777778}*x22*x26 + V{0.0277777777777778}*x24*x27) - x51*(V{4.62592926927149e-18}*x22*x26 + V{4.62592926927149e-18}*x24*x27) + x96 + x98 + V{4.44089209850063e-16};
auto x102 = V{2.22044604925031e-16}*cell[10] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + x101 - x30*x82 - x32*x83 + x59 + x63 + x66 + x71 + x72*x78 - x72*x90 - x72*x93 + x75 + x79*x81 - x86;
auto x103 = -x91;
auto x104 = x35 + x47;
auto x105 = x104 + x56;
auto x106 = x105 - V{4.5}*x103*x103;
auto x107 = -x88;
auto x108 = x104 + x67;
auto x109 = x108 - V{4.5}*x107*x107;
auto x110 = x47 + x69 - x77;
auto x111 = -x110*x38;
auto x112 = x38*x78;
auto x113 = -x38*x93;
auto x114 = x72*x81;
auto x115 = x33 + x47 - x80 + x84;
auto x116 = x115*x72;
auto x117 = V{2.22044604925031e-16}*cell[11] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{4.44089209850063e-16}*cell[17] + V{2.22044604925031e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{2.22044604925031e-16}*cell[2] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[8] + V{2.22044604925031e-16}*cell[9] + x101 + x111 + x112 + x113 + x114 - x116 - x30*x83 - x32*x82 + x65*x72 - x70*x72 + x74*x79 + x99;
auto x118 = x22*(x102 - x38*x58) + x24*(-x106*x38 - x109*x72 + x117);
auto x119 = V{0.0833333333333334}*cell[15];
auto x120 = V{0.0833333333333334}*cell[16];
auto x121 = V{0.0833333333333334}*cell[6];
auto x122 = V{0.0833333333333334}*cell[7];
auto x123 = V{0.0833333333333333}*x22*x26 + V{0.0833333333333333}*x24*x27;
auto x124 = V{0.0416666666666667}*x22*x26 + V{0.0416666666666667}*x24*x27;
auto x125 = x124*x32;
auto x126 = V{6.93889390390723e-18}*cell[0] + V{6.93889390390723e-18}*cell[12] + V{1.38777878078145e-17}*cell[13] + V{6.93889390390723e-18}*cell[3] + V{6.93889390390723e-18};
auto x127 = V{1.38777878078145e-17}*cell[10] + V{6.93889390390723e-18}*cell[11] + V{1.38777878078145e-17}*cell[14] + V{1.38777878078145e-17}*cell[15] + V{1.38777878078145e-17}*cell[16] + V{6.93889390390723e-18}*cell[17] + V{6.93889390390723e-18}*cell[18] + V{6.93889390390723e-18}*cell[2] + V{6.93889390390723e-18}*cell[8] + V{6.93889390390723e-18}*cell[9] + x126;
auto x128 = V{6.93889390390723e-18}*cell[10] + V{1.38777878078145e-17}*cell[11] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[16] + V{1.38777878078145e-17}*cell[17] + V{1.38777878078145e-17}*cell[18] + V{6.93889390390723e-18}*cell[1] + V{1.38777878078145e-17}*cell[5] + V{6.93889390390723e-18}*cell[6] + V{6.93889390390723e-18}*cell[7] + x126;
auto x129 = V{0.00115740740740741}*x22*x26 + V{0.00115740740740741}*x24*x27;
auto x130 = V{0.0833333333333333}*cell[12] - V{0.166666666666667}*cell[13] - V{0.0833333333333334}*cell[14] + V{0.0833333333333333}*cell[3] - V{0.0833333333333334}*cell[5] - x124*x28 + x127*x22 + x128*x24 + x129*x49 + x129*x51 + V{0.0555555555555555};
auto x131 = -V{0.166666666666667}*cell[10] + V{0.0833333333333333}*cell[11] + V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] - V{0.166666666666667}*cell[1] + V{0.0833333333333333}*cell[2] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] - x119 - x120 - x121 - x122 + x123*x30 - x125 + x130;
auto x132 = V{0.0833333333333334}*cell[17];
auto x133 = V{0.0833333333333334}*cell[18];
auto x134 = V{0.0833333333333334}*cell[8];
auto x135 = V{0.0833333333333334}*cell[9];
auto x136 = x124*x30;
auto x137 = V{0.0833333333333333}*cell[10] - V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] + V{0.0833333333333333}*cell[1] - V{0.166666666666667}*cell[2] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] + x123*x32 + x130 - x132 - x133 - x134 - x135 - x136;
auto x138 = V{0.00231481481481481}*x22*x26 + V{0.00231481481481481}*x24*x27;
auto x139 = x22*x26 + x24*x27;
auto x140 = -V{0.0833333333333333}*cell[10] - V{0.0833333333333333}*cell[11] + V{0.166666666666667}*cell[12] - V{0.333333333333333}*cell[13] - V{0.166666666666667}*cell[14] - V{0.0833333333333333}*cell[1] - V{0.0833333333333333}*cell[2] + V{0.166666666666667}*cell[3] - V{0.166666666666667}*cell[5] + x119 + x120 + x121 + x122 + x125 - x127*x22 - x128*x24 + x132 + x133 + x134 + x135 + x136 + x138*x49 + x138*x51 - V{0.0833333333333333}*x139*x28 + V{-0.0555555555555555};
auto x141 = V{0.00578703703703704}*x22*x26 + V{0.00578703703703704}*x24*x27;
auto x142 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x139;
auto x143 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x142;
auto x144 = V{0.0208333333333333}*x22*x26 + V{0.0208333333333333}*x24*x27;
auto x145 = V{0.0208333125}*cell[0] + V{0.0208333125}*cell[12] + V{0.041666625}*cell[13] + V{0.0208333125}*cell[3] + V{0.0208333125};
auto x146 = -x22*(V{0.041666625}*cell[10] + V{0.0208333125}*cell[11] + V{0.041666625}*cell[14] + V{0.041666625}*cell[15] + V{0.041666625}*cell[16] + V{0.0208333125}*cell[17] + V{0.0208333125}*cell[18] + V{0.0208333125}*cell[2] + V{0.0208333125}*cell[8] + V{0.0208333125}*cell[9] + x145);
auto x147 = -x24*(V{0.0208333125}*cell[10] + V{0.041666625}*cell[11] + V{0.0208333125}*cell[15] + V{0.0208333125}*cell[16] + V{0.041666625}*cell[17] + V{0.041666625}*cell[18] + V{0.0208333125}*cell[1] + V{0.041666625}*cell[5] + V{0.0208333125}*cell[6] + V{0.0208333125}*cell[7] + x145);
auto x148 = V{0.0416666666666667}*x22*x26 + V{0.0416666666666667}*x24*x27;
auto x149 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x146 + x147 - x148*x30 + V{0.0138888472222222};
auto x150 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] - x148*x32;
auto x151 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x144*x28 + x149 + x150;
auto x152 = V{0.833333333333333}*cell[13] - V{0.0833333333333333}*cell[14] - V{0.0833333333333333}*cell[5] - x143 + x151;
auto x153 = x108 + x89;
auto x154 = V{0.00115740740740741}*x22*x26 + V{0.00115740740740741}*x24*x27;
auto x155 = -V{0.166666666666667}*cell[13] + V{0.416666666666667}*cell[14] + V{0.416666666666667}*cell[5] + x143 + x151 + x154*x49 + x154*x51;
auto x156 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x142;
auto x157 = V{0.000578703703703704}*x22*x26 + V{0.000578703703703704}*x24*x27;
auto x158 = V{0.0833333333333333}*cell[12] + V{0.0833333333333334}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[5] - x148*x28 - x157*x49 - x157*x51;
auto x159 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x144*x32 + x149 + x158;
auto x160 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] - x156 + x159;
auto x161 = -x38*(x105 + x92);
auto x162 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] + x156 + x159;
auto x163 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x139;
auto x164 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x144*x30 + x146 + x147 + x150 + x158 + V{0.0138888472222222};
auto x165 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x163 + x164;
auto x166 = -x38*(x57 + x68);
auto x167 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x163 + x164;
auto x168 = x22*(x102 + x166) + x24*(x117 - x153*x72 + x161);
auto x0 = -x19*(V{0.166666666666667}*x118*x35 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x21 - x22*(V{0.4999995}*cell[10] + V{0.24999975}*cell[11] + V{0.4999995}*cell[14] + V{0.4999995}*cell[15] + V{0.4999995}*cell[16] + V{0.24999975}*cell[17] + V{0.24999975}*cell[18] + V{0.24999975}*cell[2] + V{0.24999975}*cell[8] + V{0.24999975}*cell[9] + x23) - x24*(V{0.24999975}*cell[10] + V{0.4999995}*cell[11] + V{0.24999975}*cell[15] + V{0.24999975}*cell[16] + V{0.4999995}*cell[17] + V{0.4999995}*cell[18] + V{0.24999975}*cell[1] + V{0.4999995}*cell[5] + V{0.24999975}*cell[6] + V{0.24999975}*cell[7] + x23) - x28*x37 - x30*x37 - x32*x37 + x36 - x38*x49 - x38*x51 + V{0.833332833333333});
auto x1 = -x19*(V{0.0277777777777778}*x115*x118 + V{0.0555555555555556}) + x20*(x116 + x131);
auto x2 = -x19*(V{0.0277777777777778}*x118*x85 + V{0.0555555555555556}) + x20*(x137 + x86);
auto x3 = -x19*(V{0.0277777777777778}*x118*x95 + V{0.0555555555555556}) + x20*(-x140 - x96);
auto x4 = -x19*(V{0.0138888888888889}*x118*x51 + V{0.0277777777777778}) - x20*(-x141*x49 + x152 - x51*(V{0.0196759259259259}*x22*x26 + V{0.0196759259259259}*x24*x27));
auto x5 = -x19*(V{0.0138888888888889}*x109*x118 + V{0.0277777777777778}) - x20*(-x153*x38 + x155);
auto x6 = -x19*(V{0.0138888888888889}*x110*x118 + V{0.0277777777777778}) - x20*(x111 + x160);
auto x7 = -x19*(V{0.0138888888888889}*x106*x118 + V{0.0277777777777778}) - x20*(x161 + x162);
auto x8 = -x19*(V{0.0138888888888889}*x118*x62 + V{0.0277777777777778}) - x20*(x165 + x63);
auto x9 = -x19*(V{0.0138888888888889}*x118*x58 + V{0.0277777777777778}) - x20*(x166 + x167);
auto x10 = x19*(V{0.0277777777777778}*x168*x81 + V{-0.0555555555555556}) + x20*(-x114 + x131);
auto x11 = x19*(V{0.0277777777777778}*x168*x74 + V{-0.0555555555555556}) + x20*(x137 - x75);
auto x12 = x19*(V{0.0277777777777778}*x168*x97 + V{-0.0555555555555556}) + x20*(-x140 - x98);
auto x13 = x19*(V{0.0138888888888889}*x168*x49 + V{-0.0277777777777778}) - x20*(-x141*x51 + x152 + x49*(V{0.00810185185185185}*x22*x26 + V{0.00810185185185185}*x24*x27));
auto x14 = -x19*(V{0.0138888888888889}*x118*x90 + V{0.0277777777777778}) - x20*(x155 - x38*x90);
auto x15 = x19*(V{0.0138888888888889}*x168*x78 + V{-0.0277777777777778}) - x20*(x112 + x160);
auto x16 = -x19*(V{0.0138888888888889}*x118*x93 + V{0.0277777777777778}) - x20*(x113 + x162);
auto x17 = x19*(V{0.0138888888888889}*x168*x65 + V{-0.0277777777777778}) - x20*(x165 + x66);
auto x18 = -x19*(V{0.0138888888888889}*x118*x70 + V{0.0277777777777778}) - x20*(x167 + x71);
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
return { V{0.5}*x168, x28 + x30 + x32 };
}
};

}

}
