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

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;

using MyCase = Case<Radiation, Lattice<double, descriptors::D3Q27<tag::RTLBM>>>;

int main( int argc, char *argv[] ){
    initialize(&argc, &argv);
}
