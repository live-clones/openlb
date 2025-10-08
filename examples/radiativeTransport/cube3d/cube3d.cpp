/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2017 Albert Mink, Christopher McHardy
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

/* cube3d.cpp:
 * A 3D implementation of radiation being emitted from a square or one side of a cube and
 * spreading in a bigger cuboid through media. The model offers 35 different cases of
 * optical parameters to describe the participating medium and two methods to
 * choose from.
 * The theoretical background and validation of said methods are detailed in
 * [A. Mink, C. McHardy, L. Bressel, C. Rauh and M. J. Krause. “Radiative transfer lattice Boltzmann
 * methods: 3D models and their performance in different regimes of radiative transfer”. In: Journal of
 * Quantitative Spectroscopy & Radiative Transfer, Volume 243 (2020). DOI: 10.1016/j.jqsrt.2019.106810.]
*/

#include <olb.h>
#include "rtlbmDirectedBoundary.h"

#include <stdexcept>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;

using MyCase = Case<Radiation, Lattice<double, descriptors::D3Q27<tag::RTLBM>>>;

namespace olb::parameters {
    struct ANINOSOTROPY_FACTOR : public descriptors::FIELD_BASE<1> { };

    struct ABSORPTION : public descriptors::FIELD_BASE<1> { };
    struct SCATTERING : public descriptors::FIELD_BASE<1> { };
    struct MCVALUE : public descriptors::FIELD_BASE<1> { };
    struct TOTAL_ENERGY : public descriptors::FIELD_BASE<1> { };
    struct INTENSITY : public descriptors::FIELD_BASE<1> { };

    struct CASE_NUMBER : public descriptors::TYPED_FIELD_BASE<int, 1> { };
    struct DYNAMICS_NAME : public descriptors::TYPED_FIELD_BASE<std::string, 1> { };
    struct USE_MINK : public descriptors::TYPED_FIELD_BASE<bool, 1> { };
    struct USE_DIRECTED : public descriptors::FIELD_BASE<1> { };

    struct INLET_DIRICHLET : public descriptors::FIELD_BASE<1> { };

    struct BOUNDARY_SHIFT : public descriptors::FIELD_BASE<1> { };

} // namespace olb::parameters

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
    using T = MyCase::value_t;

    T boundaryShift = parameters.get<parameters::BOUNDARY_SHIFT>();
    T lx = parameters.get<parameters::DOMAIN_L>();
    T dx = parameters.get<parameters::PHYS_DELTA_X>();
    T originShift = lx/2;

    Vector<T,3> origin = {0.+boundaryShift, -originShift, -originShift};
    Vector<T,3> extent = {lx-boundaryShift, lx, lx};
    IndicatorCuboid3D<T> span(extent, origin);


    Mesh<T, MyCase::d> mesh(span, dx, 2*singleton::mpi().getSize());
    mesh.setOverlap(parameters.get<parameters::OVERLAP>());

    return mesh;
}

void prepareGeometry(MyCase& myCase){
    OstreamManager clout( std::cout, "prepareGeometryCube" );
    clout << "Prepare cubeGeometry ..." << std::endl;

    using T = MyCase::value_t;
    auto& geometry = myCase.getGeometry();
    auto& parameters = myCase.getParameters();

    T dx = parameters.get<parameters::PHYS_DELTA_X>();
    T boundaryShift = parameters.get<parameters::BOUNDARY_SHIFT>();
    T lx = parameters.get<parameters::DOMAIN_L>();
    T originShift = lx/2;

    Vector<T,3> origin = {0.+boundaryShift, -originShift, -originShift};
    Vector<T,3> extent = {lx-boundaryShift, lx, lx};
    IndicatorCuboid3D<T> span(extent, origin);

    origin = {0+boundaryShift,0,0};

    // select between a finite or infinite source of light
    // finite light beam source
    IndicatorCuboid3D<T> emittor(2*dx, 0.2, 0.2, origin);
    // infinite light beam source
    //IndicatorCuboid3D<T> emittor(2*conversionFactorLength, 1.0-2*conversionFactorLength, 1.0-2*conversionFactorLength, origin);

    // set material number, 0 outside, 1 domain, 2 outflow, 3 inflow
    geometry.rename(0, 2, span);
    geometry.rename(2, 1, {1,1,1});
    geometry.rename(2, 3, emittor);

    geometry.clean();
    geometry.innerClean();
    geometry.checkForErrors();

    geometry.print();

    clout << "Prepare cubeGeometry ... OK" << std::endl;
}

void attachUnitConverter(MyCase& myCase){
    OstreamManager clout( std::cout, "attachUnitConverter" );

    using T = MyCase::value_t;
    auto& params   = myCase.getParameters();

    auto& Rlattice = myCase.getLattice(Radiation {});
    using RDESCRIPTOR = MyCase::descriptor_t_of<Radiation>;

    Rlattice.setUnitConverter<RadiativeUnitConverter<T,RDESCRIPTOR>>(
        params.get<parameters::RESOLUTION>(),
        params.get<parameters::LATTICE_RELAXATION_TIME>(),
        params.get<parameters::ABSORPTION>(),
        params.get<parameters::SCATTERING>(),
        params.get<parameters::ANINOSOTROPY_FACTOR>()
    );

    const auto& converter = Rlattice.getUnitConverter();
    clout << "omega = " << converter.getLatticeRelaxationTime() << std::endl;
}

void prepareLatticeMink(MyCase& myCase){
    OstreamManager clout( std::cout, "prepareLattice" );
    clout << "Prepare Lattice for Mink ..." << std::endl;
    clout << "working with diffuse approximation" << std::endl;

    using T = MyCase::value_t;
    auto& geometry = myCase.getGeometry();
    auto& params   = myCase.getParameters();

    auto& Rlattice = myCase.getLattice(Radiation {});
    using RDESCRIPTOR = MyCase::descriptor_t_of<Radiation>;

    attachUnitConverter(myCase);
    const auto& converter = Rlattice.getUnitConverter();
    converter.print();
    converter.write();

    T latticeSink = converter.getLatticeAbsorption()/converter.getLatticeDiffusion()/8.;
    T inletDirichlet = params.get<parameters::INLET_DIRICHLET>();

    T latticeRelaxationFrequency = converter.getLatticeRelaxationFrequency();
    T latticeAbsorption = converter.getLatticeAbsorption();
    T latticeScattering = converter.getLatticeScattering();
    clout << "latticeSink= " << latticeSink << std::endl;

    Rlattice.defineDynamics<NoDynamics<T, RDESCRIPTOR>>(geometry, 0);
    Rlattice.defineDynamics<P1Dynamics<T, RDESCRIPTOR>>(geometry, 1);
    Rlattice.defineDynamics<EquilibriumBoundaryFirstOrder<T, RDESCRIPTOR>>(geometry, 2);
    Rlattice.defineDynamics<EquilibriumBoundaryFirstOrder<T, RDESCRIPTOR>>(geometry, 3);

    if(params.get<parameters::USE_DIRECTED>()){
        setRtlbmDirectedBoundary<T,RDESCRIPTOR>(Rlattice, geometry.getMaterialIndicator({3}), params.get<parameters::INTENSITY>());
        Rlattice.setParameter<>
    }


    Rlattice.setParameter<descriptors::OMEGA>(latticeRelaxationFrequency);
    Rlattice.setParameter<collision::P1::ABSORPTION>( latticeAbsorption );
    Rlattice.setParameter<collision::P1::SCATTERING>( latticeScattering );
    Rlattice.setParameter<collision::Poisson::SINK>( latticeSink );
    Rlattice.setParameter<parameters::INTENSITY>( inletDirichlet );

    clout << "Prepare Lattice ... OK" << std::endl;
    return;
}

void prepareLatticeMcHardy(MyCase& myCase){
    OstreamManager clout( std::cout, "prepareLattice" );
    clout << "Prepare Lattice for McHardy ..." << std::endl;
    clout << "working with direct discretization" << std::endl;

    using T = MyCase::value_t;
    auto& geometry = myCase.getGeometry();
    auto& params   = myCase.getParameters();

    auto& Rlattice = myCase.getLattice(Radiation {});
    using RDESCRIPTOR = MyCase::descriptor_t_of<Radiation>;

    attachUnitConverter(myCase);
    const auto& converter = Rlattice.getUnitConverter();
    converter.print();
    converter.write();

    T latticeRelaxationFrequency = converter.getLatticeRelaxationFrequency();
    T anisotropyFactor = params.get<parameters::ANINOSOTROPY_FACTOR>();

    int const q = descriptors::q<RDESCRIPTOR>() - 1;
    std::array<std::array<T,q>, q> phi;
    T solution [q*(q +1)/2];
    computeAnisotropyMatrix<RDESCRIPTOR>(1e-4, anisotropyFactor, solution, phi);

    T anisoMatrix[(q+1)*(q+1)] {};
    for ( int m = 0; m < q; m++ ) {
        for ( int n = 0; n < q; n++ ) {
        anisoMatrix[m+1+n+1] = phi[m][n];
        }
    }

    T latticeAbsorption = converter.getLatticeAbsorption();
    T latticeScattering = converter.getLatticeScattering();
    T inletDirichlet = params.get<parameters::INLET_DIRICHLET>();

    Rlattice.defineDynamics<NoDynamics<T, RDESCRIPTOR>>( geometry, 0 );
    // select between the mesoscopic method with or without a Runge Kutta scheme
    Rlattice.defineDynamics<RTLBMdynamicsMcHardyRK<T, RDESCRIPTOR>>( geometry, 1 );
    // Rlattice.defineDynamics<RTLBMdynamicsMcHardy<T, RDESCRIPTOR>>( geometry, 1 );
    Rlattice.defineDynamics<EquilibriumBoundaryFirstOrder>( geometry, 2 );
    Rlattice.defineDynamics<EquilibriumBoundaryFirstOrder>( geometry, 3 );


    Rlattice.setParameter<descriptors::OMEGA>( latticeRelaxationFrequency );
    Rlattice.setParameter<Light::ANISOMATRIX>( anisoMatrix );
    Rlattice.setParameter<Light::ABSORPTION>( latticeAbsorption );
    Rlattice.setParameter<Light::SCATTERING>( latticeScattering );
    Rlattice.setParameter<parameters::INTENSITY>( inletDirichlet );

    clout << "Prepare Lattice ... OK" << std::endl;
    return;
}

int main( int argc, char *argv[] ){
    OstreamManager clout(std::cout,"main");

    using T = MyCase::value_t;

    initialize(&argc, &argv);

    MyCase::ParametersD myCaseParameters;
    {
        using namespace olb::parameters;
        myCaseParameters.set<RESOLUTION>(50);
        myCaseParameters.set<DOMAIN_L>(1.);
        myCaseParameters.set<LATTICE_RELAXATION_TIME>(1.0);
        myCaseParameters.set<MAX_PHYS_T>(12);
        myCaseParameters.set<PHYS_SAVE_ITER>(2.0);
        myCaseParameters.set<INLET_DIRICHLET>(1.0);

        myCaseParameters.set<CASE_NUMBER>(1.);
        myCaseParameters.set<DYNAMICS_NAME>(std::string("mink"));
        myCaseParameters.set<USE_MINK>(true);
        myCaseParameters.set<USE_DIRECTED>(false);

        myCaseParameters.set<BOUNDARY_SHIFT>(0.);

        //Copmuted Values
        myCaseParameters.set<PHYS_DELTA_X>([&]{
            return myCaseParameters.get<DOMAIN_L>() / myCaseParameters.get<RESOLUTION>();
        });

    }
    myCaseParameters.fromCLI(argc, argv);

    //Check if given case number is valid
    int case_number = myCaseParameters.get<parameters::CASE_NUMBER>();
    if(case_number < 1 || case_number > 35){
        throw std::runtime_error("Please select a case number between 1 and 35");
    }else{
        //Load appropriate data from xml
        std::string fName("cube3d.xml");
        std::string case_name = "case" + std::to_string(case_number);
        clout << case_name << std::endl;
        XMLreader config(fName);

        T absorption, scattering, mcvalue, totalEnergy;
        std::stringstream xmlAbsorption(config["Application"][case_name].getAttribute("absorption"));
        xmlAbsorption >> absorption;
        myCaseParameters.set<parameters::ABSORPTION>(absorption);

        std::stringstream xmlScattering(config["Application"][case_name].getAttribute("scattering"));
        xmlScattering >> scattering;
        myCaseParameters.set<parameters::SCATTERING>(scattering);

        std::stringstream xmlMcValue(config["Application"][case_name].getAttribute("mcvalue"));
        xmlMcValue >> mcvalue;
        myCaseParameters.set<parameters::MCVALUE>(mcvalue);

        std::stringstream xmlTotalEnergy(config["Application"][case_name].getAttribute("totalEnergy"));
        xmlTotalEnergy >> totalEnergy;
        myCaseParameters.set<parameters::TOTAL_ENERGY>(totalEnergy);
    }

    //Check if given dynamics name is valid
    std::string s = myCaseParameters.get<parameters::DYNAMICS_NAME>();
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    clout << s << std::endl;
    if(!(std::string("MINK") == s || std::string("MCHARDY") == s)){
        throw std::runtime_error("Incompatible dynamics selected!\n\t\t Please choose between mink or mchardy");
    }
    if(std::string("MCHARDY") == s){myCaseParameters.set<parameters::USE_MINK>(false);}

    // Set output directory
    std::string caseName = "case" + std::to_string(case_number);
    std::string caseExtension = myCaseParameters.get<parameters::USE_MINK>() ? std::string("_mink/") : std::string("_mcHardy/");
    singleton::directories().setOutputDir( "./"+std::to_string(myCaseParameters.get<parameters::RESOLUTION>())+caseName+caseExtension );


    /// === Step 3: Create Mesh ===
    Mesh mesh = createMesh(myCaseParameters);

    /// === Step 4: Create Case ===
    MyCase myCase(myCaseParameters, mesh);

    /// === Step 5: Prepare Geometry ===
    prepareGeometry(myCase);

    /// === Step 6: Prepare Lattice ===
    myCaseParameters.get<parameters::USE_MINK>() ? prepareLatticeMink(myCase) : prepareLatticeMcHardy(myCase);

}
