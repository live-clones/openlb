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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::ShanChen>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x21 = cell.template getFieldComponent<descriptors::FORCE>(2);
auto x20 = cell.template getFieldComponent<descriptors::FORCE>(1);
auto x22 = parameters.template get<descriptors::OMEGA>();
auto x19 = cell.template getFieldComponent<descriptors::FORCE>(0);
auto x23 = x22 + V{-1};
auto x24 = cell[10] + cell[14] + cell[16];
auto x25 = cell[11] + cell[18] + cell[5];
auto x26 = cell[12] + cell[7] + cell[9];
auto x27 = cell[0] + cell[13] + cell[15] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + x24 + x25 + x26;
auto x28 = x27 + V{1};
auto x29 = V{1} / (x22);
auto x30 = x19*x29;
auto x31 = x27 + V{1};
auto x32 = V{1} / (x31);
auto x33 = V{1}*cell[14];
auto x34 = V{1}*cell[16];
auto x35 = V{1}*cell[5];
auto x36 = V{1}*cell[7];
auto x37 = V{1}*cell[13] - V{1}*cell[4];
auto x38 = V{1}*cell[15] - V{1}*cell[6];
auto x39 = x32*(V{1}*cell[10] - V{1}*cell[1] + x33 + x34 - x35 - x36 + x37 + x38);
auto x40 = x30 + x39;
auto x41 = V{1.5}*(x40*x40);
auto x42 = x20*x29;
auto x43 = V{1}*cell[18];
auto x44 = V{1}*cell[9];
auto x45 = V{1}*cell[17] - V{1}*cell[8];
auto x46 = x32*(V{1}*cell[11] - V{1}*cell[2] - x33 + x35 + x37 + x43 - x44 + x45);
auto x47 = x42 + x46;
auto x48 = V{1.5}*(x47*x47);
auto x49 = x21*x29;
auto x50 = x32*(V{1}*cell[12] - V{1}*cell[3] - x34 + x36 + x38 - x43 + x44 + x45);
auto x51 = x49 + x50;
auto x52 = V{1.5}*(x51*x51);
auto x53 = x41 + x48 + x52 + V{-1};
auto x54 = V{0.0555555555555556}*x22;
auto x55 = cell[13] - cell[4];
auto x56 = cell[15] - cell[6];
auto x57 = x30 + x32*(-cell[1] - cell[5] - cell[7] + x24 + x55 + x56);
auto x58 = V{4.5}*(x57*x57);
auto x59 = V{3}*cell[14];
auto x60 = V{3}*cell[16];
auto x61 = V{3}*cell[5];
auto x62 = V{3}*cell[7];
auto x63 = V{3}*cell[13] - V{3}*cell[4];
auto x64 = V{3}*cell[15] - V{3}*cell[6];
auto x65 = x32*(V{3}*cell[10] - V{3}*cell[1] + x59 + x60 - x61 - x62 + x63 + x64);
auto x66 = V{3}*x30;
auto x67 = x65 + x66;
auto x68 = x53 + x67;
auto x69 = cell[17] - cell[8];
auto x70 = x32*(-cell[14] - cell[2] - cell[9] + x25 + x55 + x69);
auto x71 = x42 + x70;
auto x72 = V{4.5}*(x71*x71);
auto x73 = V{3}*cell[18];
auto x74 = V{3}*cell[9];
auto x75 = V{3}*cell[17] - V{3}*cell[8];
auto x76 = x32*(V{3}*cell[11] - V{3}*cell[2] - x59 + x61 + x63 + x73 - x74 + x75);
auto x77 = V{3}*x42;
auto x78 = x76 + x77;
auto x79 = x53 + x78;
auto x80 = -cell[16] - cell[18] - cell[3] + x26 + x56 + x69;
auto x81 = x32*x80 + x49;
auto x82 = V{4.5}*(x81*x81);
auto x83 = x32*(V{3}*cell[12] - V{3}*cell[3] - x60 + x62 + x64 - x73 + x74 + x75);
auto x84 = V{3}*x49;
auto x85 = x83 + x84;
auto x86 = x53 + x85;
auto x87 = V{0.0277777777777778}*x22;
auto x88 = x57 + x71;
auto x89 = V{4.5}*(x88*x88);
auto x90 = -x42 + x57 - x70;
auto x91 = -x90;
auto x92 = -x76 - x77;
auto x93 = x57 + x81;
auto x94 = V{4.5}*(x93*x93);
auto x95 = -x21*x29 - x32*x80;
auto x96 = x57 + x95;
auto x97 = -x96;
auto x98 = -x83 - x84;
auto x99 = x71 + x81;
auto x100 = V{4.5}*(x99*x99);
auto x101 = x71 + x95;
auto x102 = -x101;
auto x103 = -x41 - x48 - x52 + V{1};
auto x104 = x103 + x67;
auto x105 = x103 + x78;
auto x106 = -x65 - x66;
auto x107 = V{0.5}*x19 + x39;
auto x108 = V{0.5}*x20 + x46;
auto x109 = V{0.5}*x21 + x50;
auto x0 = -cell[0]*x23 - V{0.333333333333333}*x22*(x28*x53 + V{1});
auto x1 = -cell[1]*x23 - x54*(x28*(-x58 + x68) + V{1});
auto x2 = -cell[2]*x23 - x54*(x28*(-x72 + x79) + V{1});
auto x3 = -cell[3]*x23 - x54*(x28*(-x82 + x86) + V{1});
auto x4 = -cell[4]*x23 - x87*(x28*(x68 + x78 - x89) + V{1});
auto x5 = -(cell[5]*x23 + x87*(x28*(x68 + x92 - V{4.5}*x91*x91) + V{1}));
auto x6 = -cell[6]*x23 - x87*(x28*(x68 + x85 - x94) + V{1});
auto x7 = -(cell[7]*x23 + x87*(x28*(x68 + x98 - V{4.5}*x97*x97) + V{1}));
auto x8 = -cell[8]*x23 - x87*(x28*(-x100 + x79 + x85) + V{1});
auto x9 = -(cell[9]*x23 + x87*(x28*(x79 + x98 - V{4.5}*x102*x102) + V{1}));
auto x10 = -cell[10]*x23 + x54*(x28*(x104 + x58) + V{-1});
auto x11 = -cell[11]*x23 + x54*(x28*(x105 + x72) + V{-1});
auto x12 = -cell[12]*x23 + x54*(x28*(x103 + x82 + x85) + V{-1});
auto x13 = -cell[13]*x23 + x87*(x28*(x104 + x78 + x89) + V{-1});
auto x14 = -(cell[14]*x23 + x87*(x28*(x106 + x79 - V{4.5}*x90*x90) + V{1}));
auto x15 = -cell[15]*x23 + x87*(x28*(x104 + x85 + x94) + V{-1});
auto x16 = -(cell[16]*x23 + x87*(x28*(x106 + x86 - V{4.5}*x96*x96) + V{1}));
auto x17 = -cell[17]*x23 + x87*(x28*(x100 + x105 + x85) + V{-1});
auto x18 = -(cell[18]*x23 + x87*(x28*(x86 + x92 - V{4.5}*x101*x101) + V{1}));
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
return { x31, x107*x107 + x108*x108 + x109*x109 };
}
};

}

}
