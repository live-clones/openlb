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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::Guo<momenta::Forced> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = x9 + V{-1};
auto x11 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x12 = x11 + V{1};
auto x13 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x14 = x11 + V{1};
auto x15 = V{1} / (x14);
auto x16 = V{1}*cell[3];
auto x17 = V{1}*cell[1] - V{1}*cell[5];
auto x18 = V{1}*cell[2] - V{1}*cell[6] - V{1}*cell[7] + x16 + x17;
auto x19 = x13 - x15*x18;
auto x20 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x21 = x15*(-V{1}*cell[4] + V{1}*cell[7] + V{1}*cell[8] - x16 + x17);
auto x22 = x20 + x21;
auto x23 = x22*x22;
auto x24 = V{1.5}*x23;
auto x25 = x24 + V{-1};
auto x26 = x25 + V{1.5}*(x19*x19);
auto x27 = V{0.5}*x9 + V{-1};
auto x28 = x15*x18;
auto x29 = x13 - x28;
auto x30 = cell.template getFieldComponent<descriptors::FORCE>(0)*x29;
auto x31 = cell.template getFieldComponent<descriptors::FORCE>(1)*x22;
auto x32 = V{0.0277777777777778}*x9;
auto x33 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x34 = V{3}*cell[3];
auto x35 = V{3}*cell[1] - V{3}*cell[5];
auto x36 = x15*(-V{3}*cell[4] + V{3}*cell[7] + V{3}*cell[8] - x34 + x35);
auto x37 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x38 = V{4.5}*cell[3];
auto x39 = V{4.5}*cell[1] - V{4.5}*cell[5];
auto x40 = x15*(-V{4.5}*cell[4] + V{4.5}*cell[7] + V{4.5}*cell[8] - x38 + x39);
auto x41 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x42 = V{4.5}*cell[2] - V{4.5}*cell[6] - V{4.5}*cell[7] + x38 + x39;
auto x43 = -x15*x42 + x41;
auto x44 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x45 = V{3}*cell[2] - V{3}*cell[6] - V{3}*cell[7] + x34 + x35;
auto x46 = -x15*x45 + x26 + x44;
auto x47 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x48 = V{6}*cell[3];
auto x49 = V{6}*cell[7];
auto x50 = V{6}*cell[1] - V{6}*cell[5];
auto x51 = x15*(V{6}*cell[2] - V{6}*cell[6] + x48 - x49 + x50);
auto x52 = -x47 + x51;
auto x53 = x52 + V{3};
auto x54 = V{9}*cell[7];
auto x55 = V{9}*cell[3];
auto x56 = V{9}*cell[1] - V{9}*cell[5];
auto x57 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(1) + x15*(-V{9}*cell[4] + V{9}*cell[8] + x54 - x55 + x56);
auto x58 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(1) + x15*(-V{6}*cell[4] + V{6}*cell[8] - x48 + x49 + x50);
auto x59 = x58 + V{3};
auto x60 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x61 = x15*(V{9}*cell[2] - V{9}*cell[6] - x54 + x55 + x56);
auto x62 = -x60 + x61;
auto x63 = V{0.111111111111111}*x9;
auto x64 = V{0.111111}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x65 = -V{0.333333}*x31;
auto x66 = x14*x27;
auto x67 = V{1}*x66;
auto x68 = x37 + x40;
auto x69 = x33 + x36;
auto x70 = x47 - x51;
auto x71 = x57 + V{-3};
auto x72 = x58 + V{-3};
auto x73 = x60 - x61;
auto x74 = V{0.027778}*x66;
auto x75 = x22*x68;
auto x76 = V{0.333333}*x30;
auto x77 = V{0.111111}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x78 = x29*x29;
auto x79 = V{1.5}*x78;
auto x80 = x15*x45;
auto x81 = x15*x42;
auto x82 = x41 - x81;
auto x83 = -x24 - x79 + V{1};
auto x84 = x44 - x80 + x83;
auto x85 = x70 + V{3};
auto x0 = -cell[0]*x10 + V{1.333333}*x14*x27*(x30 + x31) - V{0.444444444444444}*x9*(x12*x26 + V{1});
auto x1 = -cell[1]*x10 + V{0.027778}*x14*x27*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x53 + x57) - cell.template getFieldComponent<descriptors::FORCE>(1)*(x59 + x62)) - x32*(x12*(-x33 - x36 + x46 - (-x19 + x20 + x21)*(x37 + x40 - x43)) + V{1});
auto x2 = -cell[2]*x10 - x63*(x12*(-x19*x43 + x46) + V{1}) - x67*(-x53*x64 + x65);
auto x3 = -cell[3]*x10 - x32*(x12*(x46 + x69 - (x19 + x22)*(x43 + x68)) + V{1}) - x74*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x70 + x71) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x72 + x73));
auto x4 = -cell[4]*x10 + V{1}*x14*x27*(-x72*x77 + x76) - x63*(x12*(x26 + x69 - x75) + V{1});
auto x5 = -cell[5]*x10 - x74*(cell.template getFieldComponent<descriptors::FORCE>(0)*(-x52 - x71) - cell.template getFieldComponent<descriptors::FORCE>(1)*(-x62 - x72)) + V{0.0277777777777778}*x9*(x12*(V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(0) - x25 - x69 - x79 - x80 + (x13 - x22 - x28)*(x41 - x68 - x81)) + V{-1});
auto x6 = -cell[6]*x10 - x67*(x64*x85 + x65) + V{0.111111111111111}*x9*(x12*(x29*x82 + x84) + V{-1});
auto x7 = -cell[7]*x10 - x74*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x57 + x85) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x59 + x73)) + V{0.0277777777777778}*x9*(x12*(x69 + x84 + (x22 + x29)*(x68 + x82)) + V{-1});
auto x8 = -cell[8]*x10 + x63*(x12*(x69 + x75 + x83) + V{-1}) + x67*(-x59*x77 + x76);
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x14, x23 + x78 };
}
};

}

}
