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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::Wagner>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x22 = cell.template getFieldComponent<descriptors::OMEGA>(0);
auto x19 = cell.template getFieldComponent<descriptors::FORCE>(0);
auto x21 = cell.template getFieldComponent<descriptors::FORCE>(2);
auto x20 = cell.template getFieldComponent<descriptors::FORCE>(1);
auto x23 = cell.template getFieldComponent<descriptors::SCALAR>(0);
auto x24 = parameters.template get<descriptors::OMEGA>();
auto x25 = x24 + V{-1};
auto x26 = cell[10] + cell[14];
auto x27 = cell[12] + cell[7];
auto x28 = cell[0] + cell[11] + cell[13] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[8] + cell[9] + x26 + x27;
auto x29 = x28 + V{1};
auto x30 = V{0.4999995} - V{0.124999875}*x22;
auto x31 = x19*x19;
auto x32 = x20*x20;
auto x33 = x21*x21;
auto x34 = V{1} / (x29);
auto x35 = V{1}*x34;
auto x36 = cell[13] - cell[4];
auto x37 = cell[15] - cell[6];
auto x38 = x36 + x37;
auto x39 = -cell[1];
auto x40 = cell[16] - cell[7];
auto x41 = x39 + x40;
auto x42 = -cell[5] + x26;
auto x43 = x38 + x41 + x42;
auto x44 = x19*x43;
auto x45 = cell[17] - cell[8];
auto x46 = x36 + x45;
auto x47 = cell[18] - cell[9];
auto x48 = -cell[2];
auto x49 = cell[11] - cell[14] + cell[5] + x48;
auto x50 = x46 + x47 + x49;
auto x51 = x20*x50;
auto x52 = x37 + x45;
auto x53 = -cell[3];
auto x54 = -cell[18] + cell[9];
auto x55 = x53 + x54;
auto x56 = -cell[16] + x27;
auto x57 = x52 + x55 + x56;
auto x58 = x21*x57;
auto x59 = x22*x23*x34;
auto x60 = x28 + V{1};
auto x61 = V{1} / ((x29)*(x29));
auto x62 = V{1.5}*x61;
auto x63 = x43*x43;
auto x64 = x62*x63;
auto x65 = x50*x50;
auto x66 = x62*x65;
auto x67 = x57*x57;
auto x68 = x62*x67;
auto x69 = x66 + x68 + V{-1};
auto x70 = x64 + x69;
auto x71 = V{6}*cell[14];
auto x72 = V{6}*cell[16];
auto x73 = V{6}*cell[5];
auto x74 = V{6}*cell[7];
auto x75 = V{6}*cell[13] - V{6}*cell[4];
auto x76 = V{6}*cell[15] - V{6}*cell[6];
auto x77 = x34*(V{6}*cell[10] - V{6}*cell[1] + x71 + x72 - x73 - x74 + x75 + x76);
auto x78 = x77 + V{-3};
auto x79 = V{0.0555555555555556}*x19;
auto x80 = V{0.16666675} - V{0.0416666875}*x22;
auto x81 = V{0.08333325} - V{0.0208333125}*x22;
auto x82 = V{0.166666666666667}*x34;
auto x83 = x32*x81 + x51*x82;
auto x84 = x33*x81 + x58*x82 - V{2.08333333354744e-08}*x59;
auto x85 = -x31*x80 + x83 + x84;
auto x86 = V{0.0555555555555556}*x24;
auto x87 = V{3}*cell[14];
auto x88 = V{3}*cell[16];
auto x89 = V{3}*cell[5];
auto x90 = V{3}*cell[7];
auto x91 = V{3}*cell[13] - V{3}*cell[4];
auto x92 = V{3}*cell[15] - V{3}*cell[6];
auto x93 = x34*(V{3}*cell[10] - V{3}*cell[1] + x87 + x88 - x89 - x90 + x91 + x92);
auto x94 = V{3}*x61;
auto x95 = x63*x94;
auto x96 = V{6}*cell[18];
auto x97 = V{6}*cell[9];
auto x98 = V{6}*cell[17] - V{6}*cell[8];
auto x99 = x34*(V{6}*cell[11] - V{6}*cell[2] - x71 + x73 + x75 + x96 - x97 + x98);
auto x100 = x99 + V{-3};
auto x101 = V{0.0555555555555556}*x20;
auto x102 = x31*x81 + x44*x82;
auto x103 = x102 - x32*x80 + x84;
auto x104 = V{3}*cell[18];
auto x105 = V{3}*cell[9];
auto x106 = V{3}*cell[17] - V{3}*cell[8];
auto x107 = x34*(V{3}*cell[11] - V{3}*cell[2] + x104 - x105 + x106 - x87 + x89 + x91);
auto x108 = x65*x94;
auto x109 = x64 + V{-1};
auto x110 = x34*(V{6}*cell[12] - V{6}*cell[3] - x72 + x74 + x76 - x96 + x97 + x98);
auto x111 = x110 + V{-3};
auto x112 = V{0.0555555555555556}*x21;
auto x113 = x102 - x33*x80 - V{2.08333333347034e-08}*x59 + x83;
auto x114 = x34*(V{3}*cell[12] - V{3}*cell[3] - x104 + x105 + x106 - x88 + x90 + x92);
auto x115 = x67*x94;
auto x116 = V{9}*cell[18];
auto x117 = V{9}*cell[5];
auto x118 = V{9}*cell[14];
auto x119 = V{9}*cell[9];
auto x120 = V{9}*cell[13] - V{9}*cell[4];
auto x121 = V{9}*cell[17] - V{9}*cell[8];
auto x122 = x34*(V{9}*cell[11] - V{9}*cell[2] + x116 + x117 - x118 - x119 + x120 + x121);
auto x123 = V{0.0277777777777778}*x19;
auto x124 = V{9}*cell[16];
auto x125 = V{9}*cell[7];
auto x126 = V{9}*cell[15] - V{9}*cell[6];
auto x127 = x34*(V{9}*cell[10] - V{9}*cell[1] - x117 + x118 + x120 + x124 - x125 + x126);
auto x128 = V{0.0277777777777778}*x20;
auto x129 = V{0.25}*x22 + V{-1};
auto x130 = V{0.25}*x129*x19;
auto x131 = x130*x20;
auto x132 = V{0.02083334375}*x22 + V{-0.083333375};
auto x133 = x132*x32;
auto x134 = V{0.01041665625}*x22 + V{-0.041666625};
auto x135 = V{0.0833333333333333}*x34;
auto x136 = V{0.0104166770833333}*x59;
auto x137 = x132*x31 - x136;
auto x138 = x133 - x134*x33 + x135*x58 + x137;
auto x139 = x131 + x138;
auto x140 = V{0.0277777777777778}*x24;
auto x141 = V{4.5}*x61;
auto x142 = cell[10] + x41;
auto x143 = cell[11] + V{2}*cell[13] - V{2}*cell[4] + x142 + x47 + x48 + x52;
auto x144 = x141*(x143*x143);
auto x145 = x70 + x93;
auto x146 = V{3} - x77;
auto x147 = -x127;
auto x148 = x99 + V{3};
auto x149 = -x131 + x138;
auto x150 = -x107;
auto x151 = -cell[17] + cell[8];
auto x152 = -cell[11] + V{2}*cell[14] + cell[2] - V{2}*cell[5] + x142 + x151 + x37 + x54;
auto x153 = -x152;
auto x154 = x34*(V{9}*cell[12] - V{9}*cell[3] - x116 + x119 + x121 - x124 + x125 + x126);
auto x155 = V{0.0277777777777778}*x21;
auto x156 = x130*x21;
auto x157 = x132*x33;
auto x158 = -x134*x32 + x135*x51 + x137 + x157;
auto x159 = x156 + x158;
auto x160 = x39 + x42;
auto x161 = cell[12] + V{2}*cell[15] - V{2}*cell[6] + x160 + x46 + x55;
auto x162 = x141*(x161*x161);
auto x163 = x110 + V{3};
auto x164 = -x156 + x158;
auto x165 = -x114;
auto x166 = -cell[12] + cell[3] + x36;
auto x167 = V{2}*cell[16] - V{2}*cell[7] + x151 + x160 + x166 + x47;
auto x168 = -x167;
auto x169 = V{0.25}*x129*x20*x21;
auto x170 = -x133 + x134*x31 - x135*x44 + x136 - x157;
auto x171 = -x169 + x170;
auto x172 = V{2}*cell[17] - V{2}*cell[8] + x38 + x49 + x53 + x56;
auto x173 = x141*(x172*x172);
auto x174 = x107 + x70;
auto x175 = -x122;
auto x176 = V{3} - x99;
auto x177 = x169 + x170;
auto x178 = -cell[15] + V{2}*cell[18] + cell[6] - V{2}*cell[9] + x166 + x40 + x49;
auto x179 = -x178;
auto x180 = x77 + V{3};
auto x181 = -x66;
auto x182 = V{1} - x68;
auto x183 = x181 + x182;
auto x184 = x183 + x93;
auto x185 = -x64;
auto x186 = x107 + x185;
auto x187 = x114 + x185;
auto x188 = -x93;
auto x189 = V{3} - x110;
auto x190 = -x154;
auto x191 = x114 + x70;
auto x0 = -cell[0]*x25 - V{0.333333333333333}*x24*(x60*x70 + V{1}) - x29*(x30*x31 + x30*x32 + x30*x33 + x35*x44 + x35*x51 + x35*x58 + V{0.124999875}*x59);
auto x1 = -cell[1]*x25 - x29*(-x78*x79 + x85) - x86*(x60*(x69 + x93 - x95) + V{1});
auto x2 = -cell[2]*x25 - x29*(-x100*x101 + x103) - x86*(x60*(x107 - x108 + x109 + x68) + V{1});
auto x3 = -cell[3]*x25 - x29*(-x111*x112 + x113) - x86*(x60*(x109 + x114 - x115 + x66) + V{1});
auto x4 = -cell[4]*x25 - x140*(x60*(x107 - x144 + x145) + V{1}) - x29*(-x123*(x122 + x78) - x128*(x100 + x127) + x139);
auto x5 = -(cell[5]*x25 + x140*(x60*(-x141*x153*x153 + x145 + x150) + V{1}) + x29*(x123*(x122 + x146) - x128*(x147 + x148) + x149));
auto x6 = -cell[6]*x25 - x140*(x60*(x114 + x145 - x162) + V{1}) - x29*(-x123*(x154 + x78) - x155*(x111 + x127) + x159);
auto x7 = -(cell[7]*x25 + x140*(x60*(-x141*x168*x168 + x145 + x165) + V{1}) + x29*(x123*(x146 + x154) - x155*(x147 + x163) + x164));
auto x8 = -cell[8]*x25 - x140*(x60*(x114 - x173 + x174) + V{1}) + x29*(x128*(x100 + x154) + x155*(x111 + x122) + x171);
auto x9 = -(cell[9]*x25 + x140*(x60*(-x141*x179*x179 + x165 + x174) + V{1}) - x29*(-x128*(x154 + x176) + x155*(x163 + x175) + x177));
auto x10 = -cell[10]*x25 + V{0.0555555555555556}*x24*(x60*(x184 + x95) + V{-1}) - x29*(-x180*x79 + x85);
auto x11 = -cell[11]*x25 + V{0.0555555555555556}*x24*(x60*(x108 + x182 + x186) + V{-1}) - x29*(-x101*x148 + x103);
auto x12 = -cell[12]*x25 + V{0.0555555555555556}*x24*(x60*(x115 + x181 + x187 + V{1}) + V{-1}) - x29*(-x112*x163 + x113);
auto x13 = -cell[13]*x25 + V{0.0277777777777778}*x24*(x60*(x144 + x184 + x186) + V{-1}) - x29*(-x123*(x122 + x180) - x128*(x127 + x148) + x139);
auto x14 = -(cell[14]*x25 + x140*(x60*(-x141*x152*x152 + x174 + x188) + V{1}) + x29*(-x123*(x175 + x180) + x128*(x127 + x176) + x149));
auto x15 = -cell[15]*x25 + V{0.0277777777777778}*x24*(x60*(x162 + x184 + x187) + V{-1}) - x29*(-x123*(x154 + x180) - x155*(x127 + x163) + x159);
auto x16 = -(cell[16]*x25 + x140*(x60*(-x141*x167*x167 + x188 + x191) + V{1}) + x29*(-x123*(x180 + x190) + x155*(x127 + x189) + x164));
auto x17 = -cell[17]*x25 + x140*(x60*(x114 + x173 + x183 + x186) + V{-1}) + x29*(x128*(x148 + x154) + x155*(x122 + x163) + x171);
auto x18 = -(cell[18]*x25 + x140*(x60*(-x141*x178*x178 + x150 + x191) + V{1}) - x29*(x128*(x148 + x190) - x155*(x122 + x189) + x177));
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
return { x29, V{1}*x61*(x63 + x65 + x67) };
}
};

}

}
