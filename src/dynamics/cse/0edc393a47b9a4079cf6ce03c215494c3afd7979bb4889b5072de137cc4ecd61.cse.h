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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::SecondOrder, collision::BGK, forcing::Guo<momenta::Forced> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x22 = x21 + V{1};
auto x23 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x24 = cell.template getFieldComponent<descriptors::VELOCITY>(0) + x23;
auto x25 = x24*x24;
auto x26 = V{1.5}*x25;
auto x27 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x28 = cell.template getFieldComponent<descriptors::VELOCITY>(1) + x27;
auto x29 = x28*x28;
auto x30 = V{1.5}*x29;
auto x31 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x32 = cell.template getFieldComponent<descriptors::VELOCITY>(2) + x31;
auto x33 = x32*x32;
auto x34 = V{1.5}*x33;
auto x35 = x26 + x30 + x34 + V{-1};
auto x36 = V{0.5}*x19 + V{-1};
auto x37 = cell.template getFieldComponent<descriptors::FORCE>(0)*x24;
auto x38 = cell.template getFieldComponent<descriptors::FORCE>(1)*x28;
auto x39 = cell.template getFieldComponent<descriptors::FORCE>(2)*x32;
auto x40 = x21 + V{1};
auto x41 = V{0.0555555555555556}*x19;
auto x42 = V{1}*cell.template getFieldComponent<descriptors::VELOCITY>(0) + x23;
auto x43 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(0) + V{4.5}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x44 = x42*x43;
auto x45 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x46 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x47 = x45 + x46;
auto x48 = x35 + x47;
auto x49 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x50 = V{6}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x51 = x49 + x50;
auto x52 = x51 + V{-3};
auto x53 = V{0.055556}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x54 = V{0.166667}*x38;
auto x55 = V{0.166667}*x39;
auto x56 = x54 + x55;
auto x57 = V{1}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x58 = x27 + x57;
auto x59 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x60 = V{4.5}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x61 = x59 + x60;
auto x62 = x58*x61;
auto x63 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x64 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x65 = x63 + x64;
auto x66 = x35 + x65;
auto x67 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x68 = V{6}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x69 = x67 + x68;
auto x70 = x69 + V{-3};
auto x71 = V{0.055556}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x72 = V{0.166667}*x37;
auto x73 = x55 + x72;
auto x74 = V{1}*cell.template getFieldComponent<descriptors::VELOCITY>(2) + x31;
auto x75 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(2) + V{4.5}*cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x76 = x74*x75;
auto x77 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x78 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x79 = x77 + x78;
auto x80 = x35 + x79;
auto x81 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x82 = V{6}*cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x83 = x81 + x82;
auto x84 = x83 + V{-3};
auto x85 = V{0.055556}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x86 = x54 + x72;
auto x87 = V{0.0277777777777778}*x19;
auto x88 = (x42 + x58)*(x43 + x61);
auto x89 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(1) + V{9}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x90 = V{0.027778}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x91 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(0) + V{9}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x92 = V{0.027778}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x93 = V{0.083333}*x39;
auto x94 = -x93;
auto x95 = x36*x40;
auto x96 = V{1}*x95;
auto x97 = -x27 + x42 - x57;
auto x98 = x43 - x59 - x60;
auto x99 = -x63 - x64;
auto x100 = -V{0.125001}*cell.template getFieldComponent<descriptors::FORCE>(0) - V{0.250002}*cell.template getFieldComponent<descriptors::VELOCITY>(0) + V{0.083334};
auto x101 = V{0.083334}*cell.template getFieldComponent<descriptors::FORCE>(1) + V{0.166668}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x102 = x89 + V{3};
auto x103 = -x49 - x50;
auto x104 = (x42 + x74)*(x43 + x75);
auto x105 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(2) + V{9}*cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x106 = V{0.027778}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x107 = V{0.083333}*x38;
auto x108 = -x107;
auto x109 = -V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(2) - V{1}*cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x110 = x109 + x42;
auto x111 = -V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(2) - V{4.5}*cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x112 = x111 + x43;
auto x113 = -x77 - x78;
auto x114 = V{0.083334}*cell.template getFieldComponent<descriptors::FORCE>(2) + V{0.166668}*cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x115 = x105 + V{3};
auto x116 = (x58 + x74)*(x61 + x75);
auto x117 = V{0.083333}*x37;
auto x118 = -x117;
auto x119 = x109 + x58;
auto x120 = x111 + x61;
auto x121 = -V{0.125001}*cell.template getFieldComponent<descriptors::FORCE>(1) - V{0.250002}*cell.template getFieldComponent<descriptors::VELOCITY>(1) + V{0.083334};
auto x122 = -x67 - x68;
auto x123 = -x26 - x30 - x34 + V{1};
auto x124 = x123 + x47;
auto x125 = x51 + V{3};
auto x126 = x123 + x65;
auto x127 = x69 + V{3};
auto x128 = x83 + V{3};
auto x129 = -x45 - x46;
auto x130 = V{0.083334}*cell.template getFieldComponent<descriptors::FORCE>(0) + V{0.166668}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x131 = x91 + V{3};
auto x132 = -V{0.125001}*cell.template getFieldComponent<descriptors::FORCE>(2) - V{0.250002}*cell.template getFieldComponent<descriptors::VELOCITY>(2) + V{0.083334};
auto x133 = -x81 - x82;
auto x0 = -cell[0]*x20 - V{0.333333333333333}*x19*(x22*x35 + V{1}) + V{1}*x36*x40*(x37 + x38 + x39);
auto x1 = -cell[1]*x20 + x36*x40*(-x52*x53 + x56) - x41*(x22*(-x44 + x48) + V{1});
auto x2 = -cell[2]*x20 + x36*x40*(-x70*x71 + x73) - x41*(x22*(-x62 + x66) + V{1});
auto x3 = -cell[3]*x20 + x36*x40*(-x84*x85 + x86) - x41*(x22*(-x76 + x80) + V{1});
auto x4 = -cell[4]*x20 - x87*(x22*(x48 + x65 - x88) + V{1}) - x96*(x90*(x52 + x89) + x92*(x70 + x91) + x94);
auto x5 = -cell[5]*x20 + V{1}*x36*x40*(-cell.template getFieldComponent<descriptors::FORCE>(1)*(x100 + x101) + x90*(x102 + x103) + x93) - x87*(x22*(x48 - x97*x98 + x99) + V{1});
auto x6 = -cell[6]*x20 - x87*(x22*(-x104 + x48 + x79) + V{1}) - x96*(x106*(x84 + x91) + x108 + x90*(x105 + x52));
auto x7 = -cell[7]*x20 + V{1}*x36*x40*(-cell.template getFieldComponent<descriptors::FORCE>(2)*(x100 + x114) + x107 + x90*(x103 + x115)) - x87*(x22*(-x110*x112 + x113 + x48) + V{1});
auto x8 = -cell[8]*x20 - x87*(x22*(-x116 + x66 + x79) + V{1}) - x96*(x106*(x84 + x89) + x118 + x92*(x105 + x70));
auto x9 = -cell[9]*x20 + V{1}*x36*x40*(-cell.template getFieldComponent<descriptors::FORCE>(2)*(x114 + x121) + x117 + x92*(x115 + x122)) - x87*(x22*(x113 - x119*x120 + x66) + V{1});
auto x10 = -cell[10]*x20 + x41*(x22*(x124 + x44) + V{-1}) + x95*(-x125*x53 + x56);
auto x11 = -cell[11]*x20 + x41*(x22*(x126 + x62) + V{-1}) + x95*(-x127*x71 + x73);
auto x12 = -cell[12]*x20 + x41*(x22*(x123 + x76 + x79) + V{-1}) + x95*(-x128*x85 + x86);
auto x13 = -cell[13]*x20 + V{0.0277777777777778}*x19*(x22*(x124 + x65 + x88) + V{-1}) - x96*(x90*(x125 + x89) + x92*(x127 + x91) + x94);
auto x14 = -cell[14]*x20 + V{1}*x36*x40*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x121 + x130) + x92*(x122 + x131) + x93) - x87*(x22*(x129 + x66 - x97*x98) + V{1});
auto x15 = -cell[15]*x20 + V{0.0277777777777778}*x19*(x22*(x104 + x124 + x79) + V{-1}) - x96*(x106*(x128 + x91) + x108 + x90*(x105 + x125));
auto x16 = -cell[16]*x20 + V{1}*x36*x40*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x130 + x132) + x106*(x131 + x133) + x107) - x87*(x22*(-x110*x112 + x129 + x80) + V{1});
auto x17 = -cell[17]*x20 + V{0.0277777777777778}*x19*(x22*(x116 + x126 + x79) + V{-1}) - x96*(x106*(x128 + x89) + x118 + x92*(x105 + x127));
auto x18 = -cell[18]*x20 + V{1}*x36*x40*(-cell.template getFieldComponent<descriptors::FORCE>(1)*(x101 + x132) + x106*(x102 + x133) + x117) - x87*(x22*(-x119*x120 + x80 + x99) + V{1});
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
return { x40, x25 + x29 + x33 };
}
};

}

}
