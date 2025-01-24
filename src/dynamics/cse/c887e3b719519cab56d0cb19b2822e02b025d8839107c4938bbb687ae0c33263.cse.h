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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell[10] + cell[14];
auto x22 = cell[12] + cell[7];
auto x23 = cell[0] + cell[11] + cell[13] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[8] + cell[9] + x21 + x22;
auto x24 = x23 + V{1};
auto x25 = V{1} / ((x24)*(x24));
auto x26 = x23 + V{1};
auto x27 = x25*x26;
auto x28 = V{0.5}*x27;
auto x29 = cell[13] - cell[4];
auto x30 = cell[15] - cell[6];
auto x31 = x29 + x30;
auto x32 = -cell[1];
auto x33 = cell[16] - cell[7];
auto x34 = x32 + x33;
auto x35 = -cell[5] + x21;
auto x36 = x31 + x34 + x35;
auto x37 = x36*x36;
auto x38 = cell[17] - cell[8];
auto x39 = x29 + x38;
auto x40 = cell[18] - cell[9];
auto x41 = -cell[2];
auto x42 = cell[11] - cell[14] + cell[5] + x41;
auto x43 = x39 + x40 + x42;
auto x44 = x43*x43;
auto x45 = x30 + x38;
auto x46 = -cell[3];
auto x47 = -cell[18] + cell[9];
auto x48 = x46 + x47;
auto x49 = -cell[16] + x22;
auto x50 = x45 + x48 + x49;
auto x51 = x50*x50;
auto x52 = V{1.5}*x25;
auto x53 = x37*x52;
auto x54 = x44*x52;
auto x55 = x51*x52;
auto x56 = x54 + x55 + V{-1};
auto x57 = x53 + x56;
auto x58 = -V{1.38777878078145e-17}*cell[0];
auto x59 = V{0.0833333333333333}*x27;
auto x60 = -V{0.0833333333333333}*cell[12] + V{0.0833333333333333}*cell[13] + V{0.0833333333333333}*cell[14] - V{0.0833333333333333}*cell[3] + V{0.0833333333333333}*cell[4] + V{0.0833333333333333}*cell[5] + x51*x59 + x58;
auto x61 = -V{0.0833333333333333}*cell[11] + V{0.0833333333333333}*cell[15] + V{0.0833333333333333}*cell[16] - V{0.0833333333333333}*cell[2] + V{0.0833333333333333}*cell[6] + V{0.0833333333333333}*cell[7] + x44*x59;
auto x62 = -V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] - V{0.166666666666667}*cell[1] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + V{0.166666666666667}*x25*x26*x37 - x60 - x61;
auto x63 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x64 = V{1} / (x24);
auto x65 = V{3}*cell[14];
auto x66 = V{3}*cell[16];
auto x67 = V{3}*cell[5];
auto x68 = V{3}*cell[7];
auto x69 = V{3}*cell[13] - V{3}*cell[4];
auto x70 = V{3}*cell[15] - V{3}*cell[6];
auto x71 = x64*(V{3}*cell[10] - V{3}*cell[1] + x65 + x66 - x67 - x68 + x69 + x70);
auto x72 = V{3}*x25;
auto x73 = x37*x72;
auto x74 = -V{0.0833333333333333}*cell[10] + V{0.0833333333333333}*cell[17] + V{0.0833333333333333}*cell[18] - V{0.0833333333333333}*cell[1] + V{0.0833333333333333}*cell[8] + V{0.0833333333333333}*cell[9] + x37*x59;
auto x75 = -V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] - V{0.166666666666667}*cell[2] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] + V{0.166666666666667}*x25*x26*x44 - x60 - x74;
auto x76 = V{3}*cell[18];
auto x77 = V{3}*cell[9];
auto x78 = V{3}*cell[17] - V{3}*cell[8];
auto x79 = x64*(V{3}*cell[11] - V{3}*cell[2] - x65 + x67 + x69 + x76 - x77 + x78);
auto x80 = x44*x72;
auto x81 = x53 + V{-1};
auto x82 = -V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] - V{0.166666666666667}*cell[3] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] + V{0.166666666666667}*x25*x26*x51 - x58 - x61 - x74;
auto x83 = x64*(V{3}*cell[12] - V{3}*cell[3] - x66 + x68 + x70 - x76 + x77 + x78);
auto x84 = x51*x72;
auto x85 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x86 = V{4.5}*x25;
auto x87 = cell[10] + x34;
auto x88 = cell[11] + V{2}*cell[13] - V{2}*cell[4] + x40 + x41 + x45 + x87;
auto x89 = x86*(x88*x88);
auto x90 = x57 + x71;
auto x91 = V{0.25}*x43;
auto x92 = x36*x64;
auto x93 = x91*x92;
auto x94 = V{0.0416666666666667}*x27;
auto x95 = -V{0.041666625}*cell[0];
auto x96 = V{0.0833333333333333}*x27;
auto x97 = V{0.0416667083333333}*cell[10] + V{0.0416667083333333}*cell[1] - x37*x96 + x95;
auto x98 = V{0.0416667083333333}*cell[11] + V{0.0416667083333333}*cell[2] - x44*x96;
auto x99 = -V{0.0833332916666667}*cell[12] + V{4.16666666690213e-08}*cell[15] + V{4.16666666690213e-08}*cell[16] + V{4.16666666690213e-08}*cell[17] + V{4.16666666690213e-08}*cell[18] - V{0.0833332916666667}*cell[3] + V{4.16666666690213e-08}*cell[6] + V{4.16666666690213e-08}*cell[7] + V{4.16666666690213e-08}*cell[8] + V{4.16666666690213e-08}*cell[9] + x51*x94 + x97 + x98;
auto x100 = x20*(V{0.375000041666667}*cell[13] - V{0.124999958333333}*cell[14] + V{0.375000041666667}*cell[4] - V{0.124999958333333}*cell[5] - x93 + x99) + V{0.0277777777777778};
auto x101 = -x79;
auto x102 = -cell[17] + cell[8];
auto x103 = -cell[11] + V{2}*cell[14] + cell[2] - V{2}*cell[5] + x102 + x30 + x47 + x87;
auto x104 = -x103;
auto x105 = x20*(-V{0.124999958333333}*cell[13] + V{0.375000041666667}*cell[14] - V{0.124999958333333}*cell[4] + V{0.375000041666667}*cell[5] + x93 + x99) + V{0.0277777777777778};
auto x106 = x32 + x35;
auto x107 = cell[12] + V{2}*cell[15] - V{2}*cell[6] + x106 + x39 + x48;
auto x108 = x86*(x107*x107);
auto x109 = V{0.25}*x50*x92;
auto x110 = V{0.0416667083333333}*cell[12] + V{4.16666666690213e-08}*cell[13] + V{4.16666666690213e-08}*cell[14] + V{0.0416667083333333}*cell[3] + V{4.16666666690213e-08}*cell[4] + V{4.16666666690213e-08}*cell[5] - x51*x96;
auto x111 = -V{0.0833332916666667}*cell[11] + V{4.16666666724907e-08}*cell[17] + V{4.16666666724907e-08}*cell[18] - V{0.0833332916666667}*cell[2] + V{4.16666666724907e-08}*cell[8] + V{4.16666666724907e-08}*cell[9] + x110 + x44*x94 + x97;
auto x112 = x20*(V{0.375000041666667}*cell[15] - V{0.124999958333333}*cell[16] + V{0.375000041666667}*cell[6] - V{0.124999958333333}*cell[7] - x109 + x111) + V{0.0277777777777778};
auto x113 = -x83;
auto x114 = -cell[12] + cell[3] + x29;
auto x115 = V{2}*cell[16] - V{2}*cell[7] + x102 + x106 + x114 + x40;
auto x116 = -x115;
auto x117 = x20*(-V{0.124999958333333}*cell[15] + V{0.375000041666667}*cell[16] - V{0.124999958333333}*cell[6] + V{0.375000041666667}*cell[7] + x109 + x111) + V{0.0277777777777778};
auto x118 = V{2}*cell[17] - V{2}*cell[8] + x31 + x42 + x46 + x49;
auto x119 = x86*(x118*x118);
auto x120 = x57 + x79;
auto x121 = x27*x50*x91;
auto x122 = -V{0.0833332916666667}*cell[10] + V{4.16666666724907e-08}*cell[15] + V{4.16666666724907e-08}*cell[16] - V{0.0833332916666667}*cell[1] + V{4.16666666724907e-08}*cell[6] + V{4.16666666724907e-08}*cell[7] + x110 + x37*x94 + x95 + x98;
auto x123 = x20*(V{0.375000041666667}*cell[17] - V{0.124999958333333}*cell[18] + V{0.375000041666667}*cell[8] - V{0.124999958333333}*cell[9] - x121 + x122) + V{0.0277777777777778};
auto x124 = -cell[15] + V{2}*cell[18] + cell[6] - V{2}*cell[9] + x114 + x33 + x42;
auto x125 = -x124;
auto x126 = x20*(-V{0.124999958333333}*cell[17] + V{0.375000041666667}*cell[18] - V{0.124999958333333}*cell[8] + V{0.375000041666667}*cell[9] + x121 + x122) + V{0.0277777777777778};
auto x127 = -x54;
auto x128 = V{1} - x55;
auto x129 = x127 + x128;
auto x130 = x129 + x71;
auto x131 = -x53;
auto x132 = x131 + x79;
auto x133 = x131 + x83;
auto x134 = -x71;
auto x135 = x57 + x83;
auto x0 = x20*(-V{0.4999995}*cell[0] + V{5.00000000028256e-07}*cell[10] + V{5.00000000028256e-07}*cell[11] + V{5.00000000028256e-07}*cell[12] + V{0.5000005}*cell[13] + V{0.5000005}*cell[14] + V{0.5000005}*cell[15] + V{0.5000005}*cell[16] + V{0.5000005}*cell[17] + V{0.5000005}*cell[18] + V{5.00000000028256e-07}*cell[1] + V{5.00000000028256e-07}*cell[2] + V{5.00000000028256e-07}*cell[3] + V{0.5000005}*cell[4] + V{0.5000005}*cell[5] + V{0.5000005}*cell[6] + V{0.5000005}*cell[7] + V{0.5000005}*cell[8] + V{0.5000005}*cell[9] - x28*x37 - x28*x44 - x28*x51) - x57*(V{0.333333333333333}*cell[0] + V{0.333333333333333}*cell[10] + V{0.333333333333333}*cell[11] + V{0.333333333333333}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] + V{0.333333333333333}) + V{-0.333333333333333};
auto x1 = x20*x62 - x63*(x56 + x71 - x73) + V{-0.0555555555555556};
auto x2 = x20*x75 - x63*(x55 + x79 - x80 + x81) + V{-0.0555555555555556};
auto x3 = x20*x82 - x63*(x54 + x81 + x83 - x84) + V{-0.0555555555555556};
auto x4 = -x100 - x85*(x79 - x89 + x90);
auto x5 = -(x105 + x85*(x101 - x86*x104*x104 + x90));
auto x6 = -x112 - x85*(-x108 + x83 + x90);
auto x7 = -(x117 + x85*(x113 - x86*x116*x116 + x90));
auto x8 = -x123 - x85*(-x119 + x120 + x83);
auto x9 = -(x126 + x85*(x113 + x120 - x86*x125*x125));
auto x10 = x20*x62 + x63*(x130 + x73) + V{-0.0555555555555556};
auto x11 = x20*x75 + x63*(x128 + x132 + x80) + V{-0.0555555555555556};
auto x12 = x20*x82 + x63*(x127 + x133 + x84 + V{1}) + V{-0.0555555555555556};
auto x13 = -x100 + V{0.0277777777777778}*x26*(x130 + x132 + x89);
auto x14 = -(x105 + x85*(x120 + x134 - x86*x103*x103));
auto x15 = -x112 + V{0.0277777777777778}*x26*(x108 + x130 + x133);
auto x16 = -(x117 + x85*(x134 + x135 - x86*x115*x115));
auto x17 = -x123 + V{0.0277777777777778}*x26*(x119 + x129 + x132 + x83);
auto x18 = -(x126 + x85*(x101 + x135 - x86*x124*x124));
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
return { x24, V{1}*x25*(x37 + x44 + x51) };
}
};

}

}
