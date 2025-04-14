/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
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

#ifndef CUBOID_DECOMPOSITION_H
#define CUBOID_DECOMPOSITION_H

#include <optional>
#include <set>

#include "core/vector.h"
#include "utilities/aliases.h"

namespace olb {

/// Reduced interface for cuboid geometry (WIP, refinement)
template <typename T, unsigned D>
struct CuboidDecomposition {
  virtual int size() const = 0;

  virtual std::vector<Cuboid<T,D>>& cuboids() = 0;

  virtual void setPeriodicity(Vector<bool,D> periodicity) = 0;

  virtual const Cuboid<T,D>& get(int iC) const = 0;
  virtual const Cuboid<T,D>& getMotherCuboid() const = 0;

  virtual Cuboid<T,D>& get(int iC) = 0;

  virtual std::optional<int> getC(Vector<T,D> physR, int offset = 0) const = 0;

  virtual std::optional<LatticeR<D+1>> getLatticeR(Vector<T,D> physR) const = 0;
  virtual std::optional<LatticeR<D+1>> getFloorLatticeR(Vector<T,D> physR) const = 0;

  virtual Vector<T,D> getPhysR(LatticeR<D+1> latticeR) const = 0;

  virtual T getDeltaR() const = 0;

  /// Returns true iff physR is covered by cuboid decomposition
  virtual bool isInside(Vector<T,D> physR) const = 0;

  virtual std::set<int> getNeighborhood(int iCglob, int overlap = 0) const = 0;

  virtual void print() const = 0;

};

}

#endif
