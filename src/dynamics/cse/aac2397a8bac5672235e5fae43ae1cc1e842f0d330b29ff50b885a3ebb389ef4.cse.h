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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::SecondOrder, collision::BGK, forcing::Guo<momenta::Forced> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = x9 + V{-1};
auto x11 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x12 = x11 + V{1};
auto x13 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x14 = cell.template getFieldComponent<descriptors::VELOCITY>(0) + x13;
auto x15 = x14*x14;
auto x16 = V{1.5}*x15;
auto x17 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x18 = cell.template getFieldComponent<descriptors::VELOCITY>(1) + x17;
auto x19 = x18*x18;
auto x20 = V{1.5}*x19;
auto x21 = x16 + x20 + V{-1};
auto x22 = V{0.5}*x9 + V{-1};
auto x23 = cell.template getFieldComponent<descriptors::FORCE>(0)*x14;
auto x24 = cell.template getFieldComponent<descriptors::FORCE>(1)*x18;
auto x25 = x11 + V{1};
auto x26 = V{0.0277777777777778}*x9;
auto x27 = V{1}*cell.template getFieldComponent<descriptors::VELOCITY>(0) + x13;
auto x28 = -V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(1) - V{1}*cell.template getFieldComponent<descriptors::VELOCITY>(1) + x27;
auto x29 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(0) + V{4.5}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x30 = -V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(1) - V{4.5}*cell.template getFieldComponent<descriptors::VELOCITY>(1) + x29;
auto x31 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(0) + V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x32 = x21 + x31;
auto x33 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x34 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x35 = -x33 - x34;
auto x36 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x37 = V{6}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x38 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x39 = V{9}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x40 = x38 + x39;
auto x41 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x42 = V{9}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x43 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x44 = V{6}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x45 = x43 + x44;
auto x46 = x45 + V{3};
auto x47 = V{0.111111111111111}*x9;
auto x48 = x27*x29;
auto x49 = x36 + x37;
auto x50 = x49 + V{-3};
auto x51 = V{0.111111}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x52 = -V{0.333333}*x24;
auto x53 = x22*x25;
auto x54 = V{1}*x53;
auto x55 = V{1}*cell.template getFieldComponent<descriptors::VELOCITY>(1) + x17;
auto x56 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(1) + V{4.5}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x57 = (x27 + x55)*(x29 + x56);
auto x58 = x33 + x34;
auto x59 = x45 + V{-3};
auto x60 = x41 + x42;
auto x61 = V{0.027778}*x53;
auto x62 = x55*x56;
auto x63 = V{0.333333}*x23;
auto x64 = V{0.111111}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x65 = -x16 - x20 + V{1};
auto x66 = x31 + x65;
auto x67 = x49 + V{3};
auto x0 = -cell[0]*x10 + V{1.333333}*x22*x25*(x23 + x24) - V{0.444444444444444}*x9*(x12*x21 + V{1});
auto x1 = -cell[1]*x10 + V{0.027778}*x22*x25*(cell.template getFieldComponent<descriptors::FORCE>(0)*(-x36 - x37 + x40 + V{3}) - cell.template getFieldComponent<descriptors::FORCE>(1)*(-x41 - x42 + x46)) - x26*(x12*(-x28*x30 + x32 + x35) + V{1});
auto x2 = -cell[2]*x10 - x47*(x12*(x32 - x48) + V{1}) - x54*(x50*x51 + x52);
auto x3 = -cell[3]*x10 - x26*(x12*(x32 - x57 + x58) + V{1}) - x61*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x40 + x50) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x59 + x60));
auto x4 = -cell[4]*x10 + V{1}*x22*x25*(-x59*x64 + x63) - x47*(x12*(x21 + x58 - x62) + V{1});
auto x5 = -cell[5]*x10 - x61*(cell.template getFieldComponent<descriptors::FORCE>(0)*(-x38 - x39 + x67) - cell.template getFieldComponent<descriptors::FORCE>(1)*(-x43 - x44 + x60 + V{3})) + V{0.0277777777777778}*x9*(x12*(x28*x30 + x35 + x66) + V{-1});
auto x6 = -cell[6]*x10 - x54*(x51*x67 + x52) + V{0.111111111111111}*x9*(x12*(x48 + x66) + V{-1});
auto x7 = -cell[7]*x10 - x61*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x40 + x67) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x46 + x60)) + V{0.0277777777777778}*x9*(x12*(x57 + x58 + x66) + V{-1});
auto x8 = -cell[8]*x10 + x47*(x12*(x58 + x62 + x65) + V{-1}) + x54*(-x46*x64 + x63);
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x25, x15 + x19 };
}
};

}

}
