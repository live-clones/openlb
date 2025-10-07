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

}

int main(int argc, char* argv[]) {
    initialize(&argc, &argv);
}
