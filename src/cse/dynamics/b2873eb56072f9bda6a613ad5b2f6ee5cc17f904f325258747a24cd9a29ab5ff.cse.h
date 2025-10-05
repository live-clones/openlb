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
struct CSE<ZouHeDynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<1, 1>, momenta::BulkStress, momenta::DefineSeparately>, 1, 1>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(0);
auto x22 = V{1.5}*x21;
auto x23 = cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(2);
auto x24 = V{1.5}*x23;
auto x25 = x22 + x24;
auto x26 = V{1} / (cell.template getFieldComponent<momenta::FixedDensity::RHO>(0));
auto x27 = x26*(cell[0] + cell[10] + V{2}*cell[11] + cell[12] + V{2}*cell[13] + cell[15] + cell[16] + V{2}*cell[17] + V{2}*cell[18] + cell[1] + cell[3] + V{2}*cell[5] + cell[6] + cell[7] + V{1});
auto x28 = V{1} - x27;
auto x29 = -x28;
auto x30 = x29*x29;
auto x31 = V{1.5}*x30;
auto x32 = x31 + V{-1};
auto x33 = x25 + x32;
auto x34 = V{0.0555555555555556}*x19;
auto x35 = V{3}*cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(0);
auto x36 = V{3}*x21;
auto x37 = V{0.0555555555555556}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x38 = x28*x28;
auto x39 = x25 - V{3}*x38;
auto x40 = x26*(V{3}*cell[0] + V{3}*cell[10] + V{6}*cell[11] + V{3}*cell[12] + V{6}*cell[13] + V{3}*cell[15] + V{3}*cell[16] + V{6}*cell[17] + V{6}*cell[18] + V{3}*cell[1] + V{3}*cell[3] + V{6}*cell[5] + V{3}*cell[6] + V{3}*cell[7] + V{3});
auto x41 = V{2} - x40;
auto x42 = x40 + V{-4};
auto x43 = x25 + x42;
auto x44 = x31 + x43;
auto x45 = V{3}*cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(2);
auto x46 = V{3}*x23;
auto x47 = V{0.00694444444444444}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x48 = -x35;
auto x49 = cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(0) + x28;
auto x50 = -V{4.5}*x49*x49;
auto x51 = x48 + x50;
auto x52 = x44 + x51;
auto x53 = V{0.0208333333333333}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x54 = x27 + V{-1};
auto x55 = cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(0) + x54;
auto x56 = -x55;
auto x57 = x35 + x44 - V{4.5}*x56*x56;
auto x58 = -x49;
auto x59 = x25 + x41;
auto x60 = x31 + x59;
auto x61 = x35 + x60 - V{4.5}*x58*x58;
auto x62 = -V{4.5}*x55*x55;
auto x63 = x48 + x62;
auto x64 = x60 + x63;
auto x65 = V{0.25}*cell[7];
auto x66 = V{0.25}*cell[16];
auto x67 = V{0.25}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x68 = -V{0.25}*cell[15] + V{0.25}*cell[6];
auto x69 = cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(0)*x67 - V{0.25}*cell[10] + V{0.25}*cell[1] + x65 - x66 + x68;
auto x70 = V{0.0277777777777778}*x19;
auto x71 = cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(2);
auto x72 = V{4.5}*(x71*x71);
auto x73 = x33 + x35;
auto x74 = -x45;
auto x75 = cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(2);
auto x76 = -x75;
auto x77 = cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(2) + x28;
auto x78 = -V{4.5}*x77*x77;
auto x79 = x74 + x78;
auto x80 = x44 + x79;
auto x81 = cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(2) + x54;
auto x82 = -x81;
auto x83 = x44 + x45 - V{4.5}*x82*x82;
auto x84 = -x77;
auto x85 = x45 + x60 - V{4.5}*x84*x84;
auto x86 = -V{4.5}*x81*x81;
auto x87 = x74 + x86;
auto x88 = x60 + x87;
auto x89 = cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(2)*x67 - V{0.25}*cell[12] + V{0.25}*cell[3] - x65 + x66 + x68;
auto x90 = V{1.5}*x38;
auto x91 = x59 + x90;
auto x92 = x43 + x90;
auto x93 = V{1} - x90;
auto x94 = -x24 + x35 + x93;
auto x95 = -x22 + x45;
auto x0 = -cell[0]*x20 - V{0.333333333333333}*x19*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x33 + V{1});
auto x1 = -cell[1]*x20 - x34*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x24 + x32 + x35 - x36) + V{1});
auto x2 = -x20*(cell[11] + x37*(x39 + x41) - x37*(x39 + x42)) - x34*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(-V{4.5}*x38 + x44) + V{1});
auto x3 = -cell[3]*x20 - x34*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x22 + x32 + x45 - x46) + V{1});
auto x4 = x20*(-cell[13] + x47*x52 - x47*x61 + x53*x57 - x53*x64 + x69) - x70*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x57 + V{1});
auto x5 = -cell[5]*x20 - x70*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x61 + V{1});
auto x6 = -cell[6]*x20 - x70*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x45 - x72 + x73) + V{1});
auto x7 = -(cell[7]*x20 + x70*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x73 + x74 - V{4.5}*x76*x76) + V{1}));
auto x8 = x20*(-cell[17] + x47*x80 - x47*x85 + x53*x83 - x53*x88 + x89) - x70*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x83 + V{1});
auto x9 = -x20*(cell[18] + x47*(x87 + x91) - x47*(x45 + x86 + x92) - x53*(x79 + x92) + x53*(x45 + x78 + x91) + x89) - x70*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x80 + V{1});
auto x10 = -cell[10]*x20 + x34*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x36 + x94) + V{-1});
auto x11 = -cell[11]*x20 - x34*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(-V{3}*x30 + x59) + V{1});
auto x12 = -cell[12]*x20 + x34*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x46 + x93 + x95) + V{-1});
auto x13 = -cell[13]*x20 - x70*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x64 + V{1});
auto x14 = -x20*(cell[5] + x47*(x63 + x91) - x47*(x35 + x62 + x92) - x53*(x51 + x92) + x53*(x35 + x50 + x91) + x69) - x70*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x52 + V{1});
auto x15 = -cell[15]*x20 + x70*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x72 + x94 + x95) + V{-1});
auto x16 = -(cell[16]*x20 + x70*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x33 + x45 + x48 - V{4.5}*x75*x75) + V{1}));
auto x17 = -cell[17]*x20 - x70*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x88 + V{1});
auto x18 = -cell[18]*x20 - x70*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x85 + V{1});
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
return { cell.template getFieldComponent<momenta::FixedDensity::RHO>(0), x21 + x23 + V{1}*x38 };
}
};

}

}
