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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::SecondOrder, collision::BGK, forcing::MCGuo<momenta::Identity> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x22 = x21 + V{1};
auto x23 = cell.template getFieldComponent<descriptors::VELOCITY>(0)*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x24 = V{1.5}*x23;
auto x25 = cell.template getFieldComponent<descriptors::VELOCITY>(1)*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x26 = V{1.5}*x25;
auto x27 = cell.template getFieldComponent<descriptors::VELOCITY>(2)*cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x28 = V{1.5}*x27;
auto x29 = x26 + x28 + V{-1};
auto x30 = x24 + x29;
auto x31 = V{0.5}*x19 + V{-1};
auto x32 = cell.template getFieldComponent<descriptors::FORCE>(0)*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x33 = cell.template getFieldComponent<descriptors::FORCE>(1)*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x34 = cell.template getFieldComponent<descriptors::FORCE>(2)*cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x35 = x21 + V{1};
auto x36 = V{0.0555555555555556}*x19;
auto x37 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x38 = V{3}*x23;
auto x39 = V{6}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x40 = x39 + V{-3};
auto x41 = V{0.055556}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x42 = V{0.166667}*x33;
auto x43 = V{0.166667}*x34;
auto x44 = x42 + x43;
auto x45 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x46 = V{3}*x25;
auto x47 = x24 + V{-1};
auto x48 = V{6}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x49 = x48 + V{-3};
auto x50 = V{0.055556}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x51 = V{0.166667}*x32;
auto x52 = x43 + x51;
auto x53 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x54 = V{3}*x27;
auto x55 = V{6}*cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x56 = x55 + V{-3};
auto x57 = V{0.055556}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x58 = x42 + x51;
auto x59 = V{0.0277777777777778}*x19;
auto x60 = cell.template getFieldComponent<descriptors::VELOCITY>(0) + cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x61 = V{4.5}*(x60*x60);
auto x62 = x30 + x37;
auto x63 = V{9}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x64 = V{0.027778}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x65 = V{9}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x66 = V{0.027778}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x67 = V{0.083333}*x34;
auto x68 = -x67;
auto x69 = x31*x35;
auto x70 = V{1}*x69;
auto x71 = -x45;
auto x72 = cell.template getFieldComponent<descriptors::VELOCITY>(0) - cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x73 = -x72;
auto x74 = V{0.166668}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x75 = V{0.083334} - V{0.250002}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x76 = V{3} - x39;
auto x77 = cell.template getFieldComponent<descriptors::VELOCITY>(0) + cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x78 = V{4.5}*(x77*x77);
auto x79 = V{9}*cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x80 = V{0.027778}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x81 = V{0.083333}*x33;
auto x82 = -x81;
auto x83 = -x53;
auto x84 = -cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x85 = cell.template getFieldComponent<descriptors::VELOCITY>(0) + x84;
auto x86 = -x85;
auto x87 = V{0.166668}*cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x88 = cell.template getFieldComponent<descriptors::VELOCITY>(1) + cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x89 = V{4.5}*(x88*x88);
auto x90 = x30 + x45;
auto x91 = V{0.083333}*x32;
auto x92 = -x91;
auto x93 = cell.template getFieldComponent<descriptors::VELOCITY>(1) + x84;
auto x94 = -x93;
auto x95 = V{0.083334} - V{0.250002}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x96 = V{3} - x48;
auto x97 = -x26;
auto x98 = V{1} - x28;
auto x99 = x97 + x98;
auto x100 = x37 + x99;
auto x101 = x39 + V{3};
auto x102 = -x24;
auto x103 = x102 + x45;
auto x104 = x48 + V{3};
auto x105 = x102 + x53;
auto x106 = x55 + V{3};
auto x107 = -x37;
auto x108 = V{0.166668}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x109 = x30 + x53;
auto x110 = V{0.083334} - V{0.250002}*cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x111 = V{3} - x55;
auto x0 = -cell[0]*x20 - V{0.333333333333333}*x19*(x22*x30 + V{1}) + V{1}*x31*x35*(x32 + x33 + x34);
auto x1 = -cell[1]*x20 + x31*x35*(-x40*x41 + x44) - x36*(x22*(x29 + x37 - x38) + V{1});
auto x2 = -cell[2]*x20 + x31*x35*(-x49*x50 + x52) - x36*(x22*(x28 + x45 - x46 + x47) + V{1});
auto x3 = -cell[3]*x20 + x31*x35*(-x56*x57 + x58) - x36*(x22*(x26 + x47 + x53 - x54) + V{1});
auto x4 = -cell[4]*x20 - x59*(x22*(x45 - x61 + x62) + V{1}) - x70*(x64*(x40 + x63) + x66*(x49 + x65) + x68);
auto x5 = -(cell[5]*x20 - V{1}*x31*x35*(-cell.template getFieldComponent<descriptors::FORCE>(1)*(x74 + x75) + x64*(x63 + x76) + x67) + x59*(x22*(x62 + x71 - V{4.5}*x73*x73) + V{1}));
auto x6 = -cell[6]*x20 - x59*(x22*(x53 + x62 - x78) + V{1}) - x70*(x64*(x40 + x79) + x80*(x56 + x65) + x82);
auto x7 = -(cell[7]*x20 - V{1}*x31*x35*(-cell.template getFieldComponent<descriptors::FORCE>(2)*(x75 + x87) + x64*(x76 + x79) + x81) + x59*(x22*(x62 + x83 - V{4.5}*x86*x86) + V{1}));
auto x8 = -cell[8]*x20 - x59*(x22*(x53 - x89 + x90) + V{1}) - x70*(x66*(x49 + x79) + x80*(x56 + x63) + x92);
auto x9 = -(cell[9]*x20 - V{1}*x31*x35*(-cell.template getFieldComponent<descriptors::FORCE>(2)*(x87 + x95) + x66*(x79 + x96) + x91) + x59*(x22*(x83 + x90 - V{4.5}*x94*x94) + V{1}));
auto x10 = -cell[10]*x20 + x36*(x22*(x100 + x38) + V{-1}) + x69*(-x101*x41 + x44);
auto x11 = -cell[11]*x20 + x36*(x22*(x103 + x46 + x98) + V{-1}) + x69*(-x104*x50 + x52);
auto x12 = -cell[12]*x20 + x36*(x22*(x105 + x54 + x97 + V{1}) + V{-1}) + x69*(-x106*x57 + x58);
auto x13 = -cell[13]*x20 + V{0.0277777777777778}*x19*(x22*(x100 + x103 + x61) + V{-1}) - x70*(x64*(x101 + x63) + x66*(x104 + x65) + x68);
auto x14 = -(cell[14]*x20 - V{1}*x31*x35*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x108 + x95) + x66*(x65 + x96) + x67) + x59*(x22*(x107 + x90 - V{4.5}*x72*x72) + V{1}));
auto x15 = -cell[15]*x20 + V{0.0277777777777778}*x19*(x22*(x100 + x105 + x78) + V{-1}) - x70*(x64*(x101 + x79) + x80*(x106 + x65) + x82);
auto x16 = -(cell[16]*x20 - V{1}*x31*x35*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x108 + x110) + x80*(x111 + x65) + x81) + x59*(x22*(x107 + x109 - V{4.5}*x85*x85) + V{1}));
auto x17 = -cell[17]*x20 + V{0.0277777777777778}*x19*(x22*(x103 + x53 + x89 + x99) + V{-1}) - x70*(x66*(x104 + x79) + x80*(x106 + x63) + x92);
auto x18 = -(cell[18]*x20 - V{1}*x31*x35*(-cell.template getFieldComponent<descriptors::FORCE>(1)*(x110 + x74) + x80*(x111 + x63) + x91) + x59*(x22*(x109 + x71 - V{4.5}*x93*x93) + V{1}));
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
return { x35, x23 + x25 + x27 };
}
};

}

}
