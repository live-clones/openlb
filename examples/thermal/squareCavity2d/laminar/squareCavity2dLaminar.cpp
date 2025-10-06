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

void prepareGeometry(MyCase& myCase){
    OstreamManager clout(std::cout, "prepraringGeometry");
    clout << "Prepare Geometry ..." << std::endl;

    using T = MyCase::value_t;
    auto& geometry = myCase.getGeometry();
    auto& params = myCase.getParameters();

    geometry.rename(0,4);

    Vector extent = params.get<parameters::DOMAIN_EXTENT>();
    std::vector<T> origin(2, T());
    IndicatorCuboid2D<T> cuboid2(extent, origin);

    geometry.rename(4, 1, cuboid2);

    std::vector<T> extendwallleft(2,T(0));
    extendwallleft[0] = params.get<parameters::PHYS_DELTA_X>();
    extendwallleft[1] = params.get<parameters::DOMAIN_EXTENT>()[0];
    std::vector<T> originwallleft(2,T(0));
    originwallleft[0] = 0.0;
    originwallleft[1] = 0.0;
    IndicatorCuboid2D<T> wallleft(extendwallleft, originwallleft);

      std::vector<T> extendwallright(2,T(0));
    extendwallright[0] = params.get<parameters::PHYS_DELTA_X>();
    extendwallright[1] = params.get<parameters::DOMAIN_EXTENT>()[0];
    std::vector<T> originwallright(2,T(0));
    originwallright[0] =  extendwallright[1] + 1.5*extendwallright[0];
    originwallright[1] = 0.0;
    IndicatorCuboid2D<T> wallright(extendwallright, originwallright);

    geometry.rename(4,2,1,wallleft);
    geometry.rename(4,3,1,wallright);

    /// Removes all not needed boundary voxels outside the surface
    geometry.clean();
    /// Removes all not needed boundary voxels inside the surface
    geometry.innerClean();
    geometry.checkForErrors();

    geometry.print();
    clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase){
    OstreamManager clout(std::cout,"prepareLattice");
    clout << "Prepare Lattice ..." << std::endl;

    using T = MyCase::value_t;
    auto& geometry = myCase.getGeometry();
    auto& params = myCase.getParameters();

    auto& NSElattice = myCase.getLattice(NavierStokes{});
    auto& ADlattice = myCase.getLattice(Temperature{});

    using NSEDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
    using TDESCRIPTOR = MyCase::descriptor_t_of<Temperature>;

    const T dx = params.get<parameters::PHYS_DELTA_X>();
    const T lx = params.get<parameters::DOMAIN_EXTENT>()[0];
    const T charU = params.get<parameters::LATTICE_CHAR_VELOCITY>();
    const T tau = params.get<parameters::LATTICE_RELAXATION_TIME>();

    const T nu = params.get<parameters::PHYS_KINEMATIC_VISCOSITY>();
    const T rho = params.get<parameters::PHYS_DENSITY>();

    const T alpha = params.get<parameters::THERMAL_DIFFUSIVITY>();
    const T beta = params.get<parameters::THERMAL_EXPANSION>();
    const T cp = params.get<parameters::HEAT_CAPACITY>();
    const T Pr = params.get<parameters::PRANDTL>();
    const T T_cold = params.get<parameters::T_COLD>();
    const T T_hot = params.get<parameters::T_HOT>();


    NSElattice.setUnitConverter<ThermalUnitConverter<T,NSEDESCRIPTOR,TDESCRIPTOR>>(
        (T) dx, //physDeltax
        (T) (tau - 0.5) / descriptors::invCs2<T,NSEDESCRIPTOR>() * util::pow(dx,2) / nu, //physDeltaT
        (T) lx, //charPhysLength
        (T) charU, //charPhysVelocity
        (T) nu, //physViscosity
        (T) rho, //physDensity
        (T) alpha, //physThermalConductivity
        (T) cp, //physSpecificHeatCapacity
        (T) beta, // physThermalExpansionCoefficient
        (T) T_cold, // charPhysLowTemperature
        (T) T_hot // charPhysHighTemperature
    );
    const auto& converter = NSElattice.getUnitConverter();
    converter.print();

    ADlattice.setUnitConverter(converter);

    NSElattice.defineDynamics<ForcedBGKdynamics>(geometry.getMaterialIndicator({1,2,3}));
    ADlattice.defineDynamics<AdvectionDiffusionBGKdynamics>(geometry.getMaterialIndicator({1,2,3}));

    clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase){
    OstreamManager clout(std::cout,"setInitialValues");
    clout << "Set initial values ..." << std::endl;

    using T = MyCase::value_t;
    auto& geometry = myCase.getGeometry();
    auto& params = myCase.getParameters();

    auto& NSElattice = myCase.getLattice(NavierStokes{});
    auto& ADlattice = myCase.getLattice(Temperature{});

    const auto& converter = NSElattice.getUnitConverter();
    T omega = converter.getLatticeRelaxationFrequency();
    T Tomega = converter.getLatticeThermalRelaxationFrequency();

    T Tcold = params.get<parameters::T_COLD>();
    T Thot = params.get<parameters::T_HOT>();
    T Tmean = params.get<parameters::T_MEAN>();

    boundary::set<boundary::BounceBack>(ADlattice, geometry, 4);
    boundary::set<boundary::BounceBack>(NSElattice, geometry, 4);

    /// sets boundary
    boundary::set<boundary::AdvectionDiffusionDirichlet>(ADlattice, geometry.getMaterialIndicator({2, 3}));
    boundary::set<boundary::LocalVelocity>(NSElattice, geometry.getMaterialIndicator({2, 3}));


    /// define initial conditions
    AnalyticalConst2D<T,T> rho(1.);
    AnalyticalConst2D<T,T> u0(0.0, 0.0);
    AnalyticalConst2D<T,T> T_cold(converter.getLatticeTemperature(Tcold));
    AnalyticalConst2D<T,T> T_hot(converter.getLatticeTemperature(Thot));
    AnalyticalConst2D<T,T> T_mean(converter.getLatticeTemperature(Tmean));


    /// for each material set Rho, U and the Equilibrium
    NSElattice.defineRhoU(geometry.getMaterialIndicator({1, 2, 3}), rho, u0);
    NSElattice.iniEquilibrium(geometry.getMaterialIndicator({1, 2, 3}), rho, u0);

    ADlattice.defineRho(geometry, 1, T_mean);
    ADlattice.iniEquilibrium(geometry, 1, T_mean, u0);
    ADlattice.defineRho(geometry, 2, T_hot);
    ADlattice.iniEquilibrium(geometry, 2, T_hot, u0);
    ADlattice.defineRho(geometry, 3, T_cold);
    ADlattice.iniEquilibrium(geometry, 3, T_cold, u0);

    NSElattice.setParameter<descriptors::OMEGA>(omega);
    ADlattice.setParameter<descriptors::OMEGA>(Tomega);

    /// Make the lattice ready for simulation
    NSElattice.initialize();
    ADlattice.initialize();

    clout << "Set initial values ... OK" << std::endl;

}




int main(int argc, char* argv[]){
    initialize(&argc, &argv);

    MyCase::ParametersD myCaseParameters;
    {
        using namespace olb::parameters;
        //General Defaults
        myCaseParameters.set<MAX_PHYS_T>(1e4);
        myCaseParameters.set<RESOLUTION>(64.0);
        myCaseParameters.set<GRAVITATIONAL_CONST>(9.81);
        myCaseParameters.set<LATTICE_CHAR_VELOCITY>(0.00797276);
        myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.65);

        //Fluid-related Defaults
        myCaseParameters.set<PHYS_DENSITY>(1.19);
        myCaseParameters.set<PHYS_KINEMATIC_VISCOSITY>(1.5126e-5);

        //Heat-related Defaults
        myCaseParameters.set<PRANDTL>(0.71);
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

    }
    myCaseParameters.fromCLI(argc, argv);

    /// === Step 3: Create Mesh ===
    Mesh mesh = createMesh(myCaseParameters);

    /// === Step 4: Create Case ===
    MyCase myCase(myCaseParameters, mesh);

    /// === Step 5: Prepare Geometry ===
    prepareGeometry(myCase);

    /// === Step 6: Prepare Lattice ===
    prepareLattice(myCase);

    /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
    setInitialValues(myCase);
    
}