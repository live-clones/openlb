/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2008 Orestis Malaspinas
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

/* squareCavity2d.cpp
 * The reference is the paper in "Gaedtke, M., Wachter, S., Raedle, M., Nirschl, H., & Krause, M. J. (2018).
 * Application of a lattice Boltzmann method combined with a Smagorinsky turbulence model to spatially resolved heat flux inside a refrigerated vehicle.
 * Computers & Mathematics with Applications, 76(10), 2315-2329."
 */

// natural convection of air in a square cavity in 2D

#include <olb.h>
#include "referenceDataLaminar.cpp"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D2Q9<descriptors::FORCE>>,
  Temperature,  Lattice<double, descriptors::D2Q5<descriptors::VELOCITY>>
>;


namespace olb::parameters {

    struct T_HOT  : public descriptors::FIELD_BASE<1> { };
    struct T_COLD : public descriptors::FIELD_BASE<1> { };
    struct T_MEAN : public descriptors::FIELD_BASE<1> { };

    struct RESIDUUM : public descriptors::FIELD_BASE<1> { };
    struct CONVERGENCE_TEST : public descriptors::TYPED_FIELD_BASE<bool, 1> { };

}

/// @brief Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& params){
    using T = MyCase::value_t;
    Vector extent = params.get<parameters::DOMAIN_EXTENT>();
    std::vector<T> origin(2, T());
    IndicatorCuboid2D<T> cuboid(extent, origin);

    Mesh<T,MyCase::d> mesh(cuboid, params.get<parameters::PHYS_DELTA_X>(), singleton::mpi().getSize());
    mesh.setOverlap(params.get<parameters::OVERLAP>());

    return mesh;
}

int main(int argc, char* argv[]){
    initialize(&argc, &argv);

    MyCase::ParametersD myCaseParameters;
    {
        using namespace olb::parameters;
        //General Defaults
        myCaseParameters.set<MAX_PHYS_T>(1e4);
        myCaseParameters.set<RESOLUTION>(32.0);
        myCaseParameters.set<GRAVITATIONAL_CONST>(9.81);
        myCaseParameters.set<LATTICE_CHAR_VELOCITY>(1.0);

        myCaseParameters.set<RESIDUUM>(1e-4);
        myCaseParameters.set<CONVERGENCE_TEST>(false);

        //Fluid-related Defaults
        myCaseParameters.set<PHYS_DENSITY>(1.19);
        myCaseParameters.set<PHYS_KINEMATIC_VISCOSITY>(1.5126e-5);

        //Heat-related Defaults
        myCaseParameters.set<PRANDTL>(0.71);
        myCaseParameters.set<PRANDTL_TURB>(0.87);
        myCaseParameters.set<RAYLEIGH>(1e3);
        myCaseParameters.set<THERMAL_EXPANSION>(3.41e-3);
        myCaseParameters.set<THERMAL_DIFFUSIVITY>(25.684e-3);
        myCaseParameters.set<HEAT_CAPACITY>(1.01309e3);

        myCaseParameters.set<T_HOT>(285.15);
        myCaseParameters.set<T_COLD>(275.15);
        myCaseParameters.set<T_MEAN>(
            (myCaseParameters.get<T_HOT>() - myCaseParameters.get<T_COLD>())/2.0
        );

        //Computed Defaults
        //Size of Domain: ((Ra*nu^2)/(Pr*g*(T_hot - T_cold)*beta))^(1/3)
        float domain_l = util::pow((myCaseParameters.get<RAYLEIGH>() * util::pow(myCaseParameters.get<PHYS_KINEMATIC_VISCOSITY>(), 2) ) /
                                   (myCaseParameters.get<PRANDTL>() * myCaseParameters.get<GRAVITATIONAL_CONST>() * (myCaseParameters.get<T_HOT>() - myCaseParameters.get<T_COLD>()) * myCaseParameters.get<THERMAL_EXPANSION>()),
                                   (MyCase::value_t) 1/3
                                );
        myCaseParameters.set<DOMAIN_EXTENT>({domain_l, domain_l});
        myCaseParameters.set<PHYS_DELTA_X>(domain_l / myCaseParameters.get<RESOLUTION>());

        //Characteristic Velocity

        //LES-related Defaults
        myCaseParameters.set<SMAGORINTZKY_CONST>(0.1);
    }
    myCaseParameters.fromCLI(argc, argv);

    /// === Step 3: Create Mesh ===
    Mesh mesh = createMesh(myCaseParameters);

    /// === Step 4: Create Case ===
    MyCase myCase(myCaseParameters, mesh);
    
}