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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, forcing::Guo<momenta::ForcedWithStress> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x11 = cell[2] + cell[3];
auto x12 = cell[7] + cell[8];
auto x13 = cell[0] + cell[1] + cell[4] + cell[5] + cell[6] + x11 + x12;
auto x14 = x13 + V{1};
auto x15 = V{1} / (x14);
auto x16 = x13 + V{1};
auto x17 = x16/((x14)*(x14));
auto x18 = cell[1] - cell[5];
auto x19 = -cell[6] - cell[7] + x11 + x18;
auto x20 = x15*x16;
auto x21 = -V{0.333333}*cell[0] + V{0.666667}*cell[1] + V{0.666667}*cell[3] + V{0.666667}*cell[5] + V{0.666667}*cell[7];
auto x22 = -cell.template getFieldComponent<descriptors::FORCE>(0)*x19*x20 + V{0.666667}*cell[2] - V{0.333333}*cell[4] + V{0.666667}*cell[6] - V{0.333333}*cell[8] - x17*x19*x19 + x21;
auto x23 = -cell[3] - cell[4] + x12 + x18;
auto x24 = cell.template getFieldComponent<descriptors::FORCE>(1)*x20*x23 - V{0.333333}*cell[2] + V{0.666667}*cell[4] - V{0.333333}*cell[6] + V{0.666667}*cell[8] - x17*x23*x23 + x21;
auto x25 = V{1}*cell[7];
auto x26 = V{1}*cell[1];
auto x27 = x20*(cell.template getFieldComponent<descriptors::FORCE>(0)*x23 - cell.template getFieldComponent<descriptors::FORCE>(1)*x19);
auto x28 = x17*x19*x23;
auto x29 = V{1}*cell[3];
auto x30 = V{1}*cell[5];
auto x31 = -x30;
auto x32 = x29 + x31;
auto x33 = V{1} / (V{2.52268963608289}*util::sqrt(x15*(x10*x10)*util::sqrt((x25 - x26 + V{0.5}*x27 + V{1}*x28 + x32)*(-V{2}*cell[1] + V{2}*cell[3] - V{2}*cell[5] + V{2}*cell[7] + V{1}*x27 + V{2}*x28) + V{1}*(x22*x22) + V{1}*(x24*x24)) + V{0.0392836979096202}/((x9)*(x9))) + V{0.5}/x9);
auto x34 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x35 = V{1}*cell[2] - V{1}*cell[6] - V{1}*cell[7] + x26 + x32;
auto x36 = -x15*x35 + x34;
auto x37 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x38 = V{1}*cell[8];
auto x39 = x15*(-V{1}*cell[4] + x25 + x26 - x29 + x31 + x38);
auto x40 = x37 + x39;
auto x41 = x40*x40;
auto x42 = V{1.5}*x41;
auto x43 = x42 + V{-1};
auto x44 = x43 + V{1.5}*(x36*x36);
auto x45 = V{1} - x33;
auto x46 = x15*x35;
auto x47 = x34 - x46;
auto x48 = cell.template getFieldComponent<descriptors::FORCE>(0)*x47;
auto x49 = cell.template getFieldComponent<descriptors::FORCE>(1)*x40;
auto x50 = V{1} - V{0.5}*x33;
auto x51 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778};
auto x52 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x53 = V{3}*cell[3];
auto x54 = V{3}*cell[1] - V{3}*cell[5];
auto x55 = x15*(-V{3}*cell[4] + V{3}*cell[7] + V{3}*cell[8] - x53 + x54);
auto x56 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x57 = V{4.5}*cell[3];
auto x58 = V{4.5}*cell[1] - V{4.5}*cell[5];
auto x59 = x15*(-V{4.5}*cell[4] + V{4.5}*cell[7] + V{4.5}*cell[8] - x57 + x58);
auto x60 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x61 = V{4.5}*cell[2] - V{4.5}*cell[6] - V{4.5}*cell[7] + x57 + x58;
auto x62 = -x15*x61 + x60;
auto x63 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x64 = V{3}*cell[2] - V{3}*cell[6] - V{3}*cell[7] + x53 + x54;
auto x65 = -x15*x64 + x44 + x63;
auto x66 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x67 = V{6}*cell[3];
auto x68 = V{6}*cell[7];
auto x69 = V{6}*cell[1] - V{6}*cell[5];
auto x70 = x15*(V{6}*cell[2] - V{6}*cell[6] + x67 - x68 + x69);
auto x71 = -x66 + x70;
auto x72 = x71 + V{3};
auto x73 = V{9}*cell[7];
auto x74 = V{9}*cell[3];
auto x75 = V{9}*cell[1] - V{9}*cell[5];
auto x76 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(1) + x15*(-V{9}*cell[4] + V{9}*cell[8] + x73 - x74 + x75);
auto x77 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(1) + x15*(-V{6}*cell[4] + V{6}*cell[8] - x67 + x68 + x69);
auto x78 = x77 + V{3};
auto x79 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x80 = x15*(V{9}*cell[2] - V{9}*cell[6] - x73 + x74 + x75);
auto x81 = -x79 + x80;
auto x82 = V{0.027778}*x14;
auto x83 = V{0.111111111111111}*cell[0] + V{0.111111111111111}*cell[1] + V{0.111111111111111}*cell[2] + V{0.111111111111111}*cell[3] + V{0.111111111111111}*cell[4] + V{0.111111111111111}*cell[5] + V{0.111111111111111}*cell[6] + V{0.111111111111111}*cell[7] + V{0.111111111111111}*cell[8] + V{0.111111111111111};
auto x84 = V{0.333333}*x49;
auto x85 = V{0.111111}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x86 = V{1}*x14;
auto x87 = x56 + x59;
auto x88 = x52 + x55;
auto x89 = x66 - x70;
auto x90 = x76 + V{-3};
auto x91 = x77 + V{-3};
auto x92 = x79 - x80;
auto x93 = x40*x87;
auto x94 = V{0.333333}*x48;
auto x95 = V{0.111111}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x96 = x47*x47;
auto x97 = V{1.5}*x96;
auto x98 = x15*x64;
auto x99 = x15*x61;
auto x100 = x60 - x99;
auto x101 = -x42 - x97 + V{1};
auto x102 = x101 + x63 - x98;
auto x103 = x89 + V{3};
auto x0 = V{1}*cell[0]*x45 - V{1.333333}*x14*x50*(x48 + x49) - x33*(x44*(V{0.444444444444444}*cell[0] + V{0.444444444444444}*cell[1] + V{0.444444444444444}*cell[2] + V{0.444444444444444}*cell[3] + V{0.444444444444444}*cell[4] + V{0.444444444444444}*cell[5] + V{0.444444444444444}*cell[6] + V{0.444444444444444}*cell[7] + V{0.444444444444444}*cell[8] + V{0.444444444444444}) + V{0.444444444444444});
auto x1 = V{1}*cell[1]*x45 - x33*(x51*(-x52 - x55 + x65 - (-x36 + x37 + x39)*(x56 + x59 - x62)) + V{0.0277777777777778}) - x50*x82*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x72 + x76) - cell.template getFieldComponent<descriptors::FORCE>(1)*(x78 + x81));
auto x2 = V{1}*cell[2]*x45 - x33*(x83*(-x36*x62 + x65) + V{0.111111111111111}) - x50*x86*(x72*x85 + x84);
auto x3 = x29*x45 - x33*(x51*(x65 + x88 - (x36 + x40)*(x62 + x87)) + V{0.0277777777777778}) + x50*x82*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x89 + x90) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x91 + x92));
auto x4 = V{1}*cell[4]*x45 - x33*(x83*(x44 + x88 - x93) + V{0.111111111111111}) - x50*x86*(-x91*x95 + x94);
auto x5 = x30*x45 - x33*(x51*(x43 - x63 + x88 + x97 + x98 - (x34 - x40 - x46)*(x60 - x87 - x99)) + V{0.0277777777777778}) + x50*x82*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x71 + x90) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x81 + x91));
auto x6 = V{1}*cell[6]*x45 + x33*(x83*(x100*x47 + x102) + V{-0.111111111111111}) + x50*x86*(x103*x85 - x84);
auto x7 = x25*x45 + x33*(x51*(x102 + x88 + (x100 + x87)*(x40 + x47)) + V{-0.0277777777777778}) + x50*x82*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x103 + x76) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x78 + x92));
auto x8 = x33*(x83*(x101 + x88 + x93) + V{-0.111111111111111}) + x38*x45 - x50*x86*(-x78*x95 + x94);
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x14, x41 + x96 };
}
};

}

}
