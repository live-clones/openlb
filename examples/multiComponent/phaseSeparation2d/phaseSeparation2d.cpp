/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Peter Weisbrod
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

/* phaseSeparation2d.cpp:
 * In this example the simulation is initialized with a given
 * density plus a small random number all over the domain. This
 * condition is unstable and leads to liquid-vapor phase separation.
 * Boundaries are assumed to be periodic. This example shows the
 * usage of multiphase flow.
 */


#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::graphics;

// === Step 1: Declarations ===
using MyCase = Case<
    NavierStokes, Lattice<double, descriptors::D2Q9<descriptors::VELOCITY, descriptors::EXTERNAL_FORCE, descriptors::STATISTIC>>
>;

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
    using T = MyCase::value_t;
    Vector extent = params.get<parameters::DOMAIN_EXTENT>();
    std::vector<T> origin(2,T());
    IndicatorCuboid2D<T> cuboid(extent, origin);

    T dx = 1;
    Mesh<T,MyCase::d> mesh(cuboid, dx, singleton::mpi().getSize());
    mesh.setOverlap(params.get<parameters::OVERLAP>());
    mesh.getCuboidDecomposition().setPeriodicity({true,true});

    return mesh;
}

void prepareGeometry(MyCase& myCase){
    auto& geometry = myCase.getGeometry();

    OstreamManager clout( std::cout,"prepareGeometry" );
    clout << "Prepare Geometry ..." << std::endl;

    // Sets material number for fluid
    geometry.rename( 0,1 );

    // Removes all not needed boundary voxels outside the surface
    geometry.clean();
    // Removes all not needed boundary voxels inside the surface
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
    auto& params   = myCase.getParameters();

    auto& NSElattice = myCase.getLattice(NavierStokes {});
    using NSEDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

    NSElattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,NSEDESCRIPTOR>>(
        (T) params.get<parameters::DOMAIN_EXTENT>()[0], //resolution
        (T) params.get<parameters::LATTICE_RELAXATION_TIME>(),
        (T) params.get<parameters::PHYS_CHAR_LENGTH>(),
        (T) params.get<parameters::PHYS_CHAR_VELOCITY>(),
        (T) params.get<parameters::PHYS_CHAR_VISCOSITY>(),
        (T) 1 //density
    );
    const auto& converter = NSElattice.getUnitConverter();
    //converter.print();

    using BulkDynamics = ForcedShanChenBGKdynamics<T, NSEDESCRIPTOR, momenta::ExternalVelocityTuple>;
    const T omega1 = 1.0;

    NSElattice.defineDynamics<BulkDynamics>(geometry, 1);

    // Initial conditions
    AnalyticalConst2D<T,T> noise( 2. );
    std::vector<T> v( 2,T() );
    AnalyticalConst2D<T,T> zeroVelocity( v );
    AnalyticalConst2D<T,T> oldRho( 199. );
    AnalyticalRandom2D<T,T> random;
    AnalyticalIdentity2D<T,T> newRho( random*noise+oldRho );

    // Initialize all values of distribution functions to their local equilibrium
    NSElattice.defineRhoU( geometry, 1, newRho, zeroVelocity );
    NSElattice.iniEquilibrium( geometry, 1, newRho, zeroVelocity );

    NSElattice.setParameter<descriptors::OMEGA>(omega1);

    NSElattice.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());
    {
        auto& communicator = NSElattice.getCommunicator(stage::Coupling());
        communicator.requestField<descriptors::STATISTIC>();
        communicator.requestOverlap(1);
        communicator.exchangeRequests();
    }

    using COUPLING = ShanChenForcedSingleComponentPostProcessor<T,NSEDESCRIPTOR,interaction::ShanChen94>;
    NSElattice.addPostProcessor<stage::Coupling>(meta::id<COUPLING>{});
    NSElattice.setParameter<COUPLING::G>(T(-120));
    NSElattice.setParameter<COUPLING::RHO0>(T(1));
    NSElattice.setParameter<COUPLING::OMEGA>(omega1);

    NSElattice.template addCustomTask<stage::PostStream>([&]() {
        NSElattice.executePostProcessors(stage::PreCoupling());
        NSElattice.executePostProcessors(stage::Coupling());
    });

    // Make the lattice ready for simulation
    NSElattice.initialize();
    clout << "Prepare Lattice ... OK" << std::endl;
}

int main(int argc, char *argv[]){
    initialize( &argc, &argv );

    /// === Step 2: Set Parameters ===
    MyCase::ParametersD myCaseParameters;
    {
        using namespace olb::parameters;
        myCaseParameters.set<MAX_LATTICE_T>(10000);
        myCaseParameters.set<LATTICE_STAT_ITER_T>(20);
        myCaseParameters.set<LATTICE_VTK_ITER_T>(20);

        myCaseParameters.set<LATTICE_RELAXATION_TIME>(1);
        myCaseParameters.set<PHYS_CHAR_LENGTH>(1e-5);
        myCaseParameters.set<PHYS_CHAR_VELOCITY>(1e-6);
        myCaseParameters.set<PHYS_CHAR_VISCOSITY>(0.1);

        myCaseParameters.set<DOMAIN_EXTENT>({200, 200});

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