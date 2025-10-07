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

    struct CASE_NUMBER : public descriptors::TYPED_FIELD_BASE<int, 1> { };
    struct DYNAMICS_NAME : public descriptors::TYPED_FIELD_BASE<std::string, 1> { };
    struct USE_MINK : public descriptors::TYPED_FIELD_BASE<bool, 1> { };

} // namespace olb::parameters



int main( int argc, char *argv[] ){
    OstreamManager clout(std::cout,"main");
    initialize(&argc, &argv);

    MyCase::ParametersD myCaseParameters;
    {
        using namespace olb::parameters;
        myCaseParameters.set<RESOLUTION>(50);
        myCaseParameters.set<LATTICE_RELAXATION_TIME>(1.0);
        myCaseParameters.set<MAX_PHYS_T>(12);
        myCaseParameters.set<PHYS_SAVE_ITER>(2.0);

        myCaseParameters.set<CASE_NUMBER>(1);
        myCaseParameters.set<DYNAMICS_NAME>(std::string("mink"));
        myCaseParameters.set<USE_MINK>(true);

    }
    myCaseParameters.fromCLI(argc, argv);

    //Check if given case number is valid
    int case_number = myCaseParameters.get<parameters::CASE_NUMBER>();
    if(case_number < 1 || case_number > 35){
        throw std::runtime_error("Please select a case number between 1 and 35");
    }

    //Check if given dynamics name is valid
    std::string s = myCaseParameters.get<parameters::DYNAMICS_NAME>();
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    clout << s << std::endl;
    if(!(std::string("MINK") == s || std::string("MCHARDY") == s)){
        throw std::runtime_error("Incompatible dynamics selected!\n\t\t Please choose between mink or mchardy");
    }
    if(std::string("MCHARDY") == s){myCaseParameters.set<parameters::USE_MINK>(false);}

    myCaseParameters.print();
}
