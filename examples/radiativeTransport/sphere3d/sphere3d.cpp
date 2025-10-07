/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2014-2018 Albert Mink
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

/* sphere3d.cpp:
 * A 3D implementation of radiation being emitted from a sphere and
 * spreading in a bigger sphere through highly scattering media. The macroscopic
 * method is based on the diffusion approximation.
 * The theoretical background and validation of said method are detailed in
 * [A. Mink, G. Thäter, H. Nirschl, and M. J. Krause. “A 3D Lattice Boltzmann method
 * for light simulation in participating media”. In: Journal of Computational Science
 * 17.Part 2 (2016). Discrete Simulation of Fluid Dynamics 2015, pp. 431–437. DOI:
 * 10.1016/j.jocs.2016.03.014.]
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;

using MyCase = Case<
    Radiation, Lattice<double, descriptors::D3Q7<descriptors::VELOCITY>>
>;

namespace olb::parameters{
    struct DOMAIN_L : public descriptors::FIELD_BASE<1> { };

    struct ABSORPTION : public descriptors::FIELD_BASE<1> { };
    struct SCATTERING : public descriptors::FIELD_BASE<1> { };
    struct MU_EFF : public descriptors::FIELD_BASE<1> { };
}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
    using T = MyCase::value_t;

    const auto lx = params.get<parameters::DOMAIN_L>();
    const auto dx = params.get<parameters::PHYS_DELTA_X>();

    Vector<T,3> origin;
    IndicatorSphere3D<T> sphereAll(origin,  1 * lx);

    Mesh<T,MyCase::d> mesh(sphereAll, dx, 2*singleton::mpi().getSize());
    mesh.setOverlap(params.get<parameters::OVERLAP>());

    return mesh;
}

void prepareGeometry(MyCase& myCase){
    OstreamManager clout(std::cout, "prepraringGeometry");
    clout << "Prepare Geometry ..." << std::endl;

    using T = MyCase::value_t;
    auto& geometry = myCase.getGeometry();
    auto& params = myCase.getParameters();

    const T dx = params.get<parameters::PHYS_DELTA_X>();
    const T lx = params.get<parameters::DOMAIN_L>();

    Vector<T,3> origin;

    //Domain, Sphere of Size lx
    IndicatorSphere3D<T> sphereAll(origin,  1 * lx);
    geometry.rename(0, 2, sphereAll);
    geometry.rename(2, 1, {1,1,1});

    //Emitter Spehere of 10% domain:l
    IndicatorSphere3D<T> sphereEm(origin, 0.1 * lx);
    geometry.rename(1, 3, sphereEm);

    geometry.clean();
    geometry.innerClean();
    geometry.checkForErrors();

    geometry.print();
}

void prepareLattice(MyCase& myCase){
    OstreamManager clout(std::cout,"prepareLattice");
    clout << "Prepare Lattice ..." << std::endl;

    using T = MyCase::value_t;
    auto& geometry = myCase.getGeometry();
    auto& params = myCase.getParameters();

    auto& Rlattice = myCase.getLattice(Radiation{});

    using RDESCRIPTOR = MyCase::descriptor_t_of<Radiation>;

    Rlattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, RDESCRIPTOR>>(
        (int) params.get<parameters::RESOLUTION>(), // resolution: number of voxels per charPhysL
        (T)   params.get<parameters::LATTICE_RELAXATION_TIME>(), // latticeRelaxationTime: relaxation time, has to be greater than 0.5!
        (T)   1.0,  // charPhysLength: reference length of simulation geometry
        (T)   0.0,  // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__  || none for radiation
        (T)   1.0,  // physViscosity: physical kinematic viscosity in __m^2 / s__  || radiative diffusion coefficient D
        (T)   1.0   // physDensity: physical density in __kg / m^3__
    );

    const auto& converter = Rlattice.getUnitConverter();
    converter.print();

    const T omega = converter.getLatticeRelaxationFrequency();

    //define dynamics
    // Material=0 -->do nothing
    Rlattice.defineDynamics<NoDynamics<T, RDESCRIPTOR>>(geometry, 0);

    // Material=1 -->sourced advection diffusion dynamics OR poisson dynamics with sink term
    // when using poisson dynamics, setting the sink parameter is necessary
    // when using sourced advection diffusion dynamics, the sink term is handeled by the post processor
    auto bulkIndicator = geometry.getMaterialIndicator({1});
    // sLattice.defineDynamics<PoissonDynamics<T,DESCRIPTOR>>(bulkIndicator);
    Rlattice.defineDynamics<SourcedAdvectionDiffusionBGKdynamics<T, RDESCRIPTOR>>(bulkIndicator);

    // Material=2,3 -->first order equilibrium
    Rlattice.defineDynamics<EquilibriumBoundaryFirstOrder>(geometry, 2);
    Rlattice.defineDynamics<EquilibriumBoundaryFirstOrder>(geometry, 3);

    // set the parameters for the dynamics
    Rlattice.setParameter<descriptors::OMEGA>(omega);
    // Rlattice.setParameter<collision::Poisson::SINK>(latticeSink);

    clout << "Prepare Lattice ... OK" << std::endl;
    return;
}


int main(int argc, char* argv[]) {
    initialize(&argc, &argv);

    /// === Step 2: Set Parameters ===
    MyCase::ParametersD myCaseParameters;
    {
        using namespace olb::parameters;
        myCaseParameters.set<RESOLUTION>(40);
        myCaseParameters.set<MAX_PHYS_T>(100);
        myCaseParameters.set<DOMAIN_L>(1.0);
        myCaseParameters.set<ABSORPTION>(0.5);
        myCaseParameters.set<SCATTERING>(1.5);

        myCaseParameters.set<LATTICE_RELAXATION_TIME>(1.0);
        myCaseParameters.set<MU_EFF>([&] { return
            util::sqrt(3* myCaseParameters.get<ABSORPTION>() * (myCaseParameters.get<ABSORPTION>() + myCaseParameters.get<SCATTERING>()));
        });

        myCaseParameters.set<PHYS_DELTA_X>([&] { return
            myCaseParameters.get<DOMAIN_L>() / myCaseParameters.get<RESOLUTION>();
        });

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

}
