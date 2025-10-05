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
struct CSE<ZouHeDynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, 1, -1>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{-1};
auto x22 = cell[0] + cell[10] + cell[12] + V{2}*cell[14] + cell[15] + cell[16] + cell[1] + V{2}*cell[2] + cell[3] + V{2}*cell[4] + cell[6] + cell[7] + V{2}*cell[8] + V{2}*cell[9] + V{1};
auto x23 = -x22/x21;
auto x24 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x25 = V{1.5}*x24;
auto x26 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x27 = V{1.5}*x26;
auto x28 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x29 = V{1.5}*x28;
auto x30 = x27 + x29 + V{-1};
auto x31 = x25 + x30;
auto x32 = V{0.0555555555555556}*x19;
auto x33 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x34 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x35 = x25 + V{-1};
auto x36 = -V{3}*x26 + x29 + x34 + x35;
auto x37 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x38 = V{0.0277777777777778}*x19;
auto x39 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x40 = x39*x39;
auto x41 = x31 + x34;
auto x42 = x33 - V{4.5}*x40 + x41;
auto x43 = -x42;
auto x44 = V{0.25}*cell[10];
auto x45 = V{0.25}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x46 = V{0.00694444444444444}*x23;
auto x47 = V{0.0208333333333333}*x23;
auto x48 = -x33;
auto x49 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x50 = x41 + x48 - V{4.5}*x49*x49;
auto x51 = -x50;
auto x52 = -x27;
auto x53 = V{1} - x29;
auto x54 = x52 + x53;
auto x55 = x33 + x54;
auto x56 = -x25;
auto x57 = x34 + x56;
auto x58 = V{4.5}*x40 + x55 + x57;
auto x59 = -x34;
auto x60 = -x49;
auto x61 = x31 + x33;
auto x62 = x59 + x61 - V{4.5}*x60*x60;
auto x63 = -x62;
auto x64 = V{0.25}*cell[7];
auto x65 = V{0.25}*cell[16];
auto x66 = x64 - x65;
auto x67 = V{0.25}*cell[15];
auto x68 = V{0.25}*cell[6] - x67;
auto x69 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x70 = x69*x69;
auto x71 = -x37;
auto x72 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x73 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x72;
auto x74 = -x73;
auto x75 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x76 = x75*x75;
auto x77 = x37 + x41 - V{4.5}*x76;
auto x78 = -x77;
auto x79 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x72;
auto x80 = -x79;
auto x81 = x41 + x71 - V{4.5}*x80*x80;
auto x82 = -x81;
auto x83 = V{1} / (x21);
auto x84 = x22*x83;
auto x85 = V{0.0555555555555556}*x84;
auto x86 = V{3}*x26 + x53 + x57;
auto x87 = x37 + x56;
auto x88 = V{0.00694444444444444}*x84;
auto x89 = V{0.0208333333333333}*x84;
auto x90 = -V{0.25}*cell[6] + x67;
auto x91 = -x64 + x65;
auto x92 = x31 + x37;
auto x93 = V{0.25}*cell[12];
auto x94 = V{0.25}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x95 = x37 + x54 + x57 + V{4.5}*x76;
auto x96 = x59 + x92 - V{4.5}*x79*x79;
auto x97 = -x96;
auto x0 = -cell[0]*x20 - V{0.333333333333333}*x19*(x23*x31 + V{1});
auto x1 = -cell[1]*x20 - x32*(-x23*(V{3}*x24 - x30 - x33) + V{1});
auto x2 = -cell[2]*x20 - x32*(x23*x36 + V{1});
auto x3 = -cell[3]*x20 - x32*(-x23*(-x27 + V{3}*x28 - x35 - x37) + V{1});
auto x4 = -cell[4]*x20 - x38*(-x23*x43 + V{1});
auto x5 = x20*(-cell[14] + V{0.25}*cell[1] + x23*x45 + x43*x46 - x44 - x46*x58 + x47*x51 - x47*x63 + x66 + x68) - x38*(-x23*x63 + V{1});
auto x6 = -cell[6]*x20 - x38*(-x23*(-x37 - x61 + V{4.5}*x70) + V{1});
auto x7 = -(cell[7]*x20 + x38*(x23*(x61 + x71 - V{4.5}*x74*x74) + V{1}));
auto x8 = -cell[8]*x20 - x38*(-x23*x78 + V{1});
auto x9 = -cell[9]*x20 - x38*(-x23*x82 + V{1});
auto x10 = -cell[10]*x20 - x32*(-x23*(V{3}*x24 + x55) + V{1});
auto x11 = -x20*(cell[2] - x36*x85 - x85*x86) - x32*(-x23*x86 + V{1});
auto x12 = -cell[12]*x20 - x32*(-x23*(V{3}*x28 + x52 + x87 + V{1}) + V{1});
auto x13 = -x20*(V{0.25}*cell[1] + cell[4] + V{0.00694444444444444}*x22*x62*x83 - x42*x89 - x44 - x45*x84 - x50*x88 - x58*x89 - x90 - x91) - x38*(-x23*x58 + V{1});
auto x14 = -cell[14]*x20 - x38*(-x23*x51 + V{1});
auto x15 = -cell[15]*x20 - x38*(-x23*(x55 + V{4.5}*x70 + x87) + V{1});
auto x16 = -(cell[16]*x20 + x38*(x23*(x48 + x92 - V{4.5}*x73*x73) + V{1}));
auto x17 = -x20*(V{0.25}*cell[3] + cell[8] + V{0.00694444444444444}*x22*x83*x96 - x66 - x77*x89 - x81*x88 - x84*x94 - x89*x95 - x90 - x93) - x38*(-x23*x95 + V{1});
auto x18 = x20*(V{0.25}*cell[3] - cell[9] + x23*x94 + x46*x78 - x46*x95 + x47*x82 - x47*x97 + x68 + x91 - x93) - x38*(-x23*x97 + V{1});
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
return { -V{1}*x84, x24 + x26 + x28 };
}
};

}

}
