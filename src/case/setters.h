/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Adrian Kummerlaender
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

#ifndef CASE_SETTERS_H
#define CASE_SETTERS_H

namespace olb {

namespace momenta {

template <typename T, typename DESCRIPTOR>
void setVelocity(SuperLattice<T,DESCRIPTOR>& sLattice,
                 FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& domainI,
                 AnalyticalF<DESCRIPTOR::d,T,T>& velocityF)
{
  const auto& converter = sLattice.getUnitConverter();
  AnalyticCalcMultiplication<DESCRIPTOR::d,T,T> scaledVelocityF(1/converter.getConversionFactorVelocity(),
                                                                velocityF);
  sLattice.defineU(std::move(domainI), scaledVelocityF);
}

}

}

#endif
