/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007-2014 Mathias J. Krause
 *                2024 Adrian Kummerlaender
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

/** \file
 * The description of a vector of 3D cuboid -- header file.
 */

#ifndef CUBOID_GEOMETRY_H
#define CUBOID_GEOMETRY_H

#include <vector>
#include <fstream>

#include "cuboidDecomposition.h"

#include "core/singleton.h"
#include "core/vector.h"

#include "geometry/cuboid.h"

#include "io/ostreamManager.h"
#include "io/xmlReader.h"

#include "cuboid.h"
#include "cuboidGeometryMinimizer.h"

namespace olb {

template <typename T> class LoadBalancer;

/// A cuboid geometry represents a voxel mesh
/** A cuboid geometry is given by a number of cuboids. To represent
 * a connected domain it is required that the distance between two
 * neighbouring cuboids is less than the smallest delta of them.
 *
 * By the class, access is provided to the cuboids. Core methods of
 * the class are transforming lattice to physical positions in the
 * corresponding unit systems.
 */
template<typename T, unsigned D>
class CuboidGeometry : public CuboidDecomposition<T,D> {
private:
  /// Cuboid which contains all other cuboids
  Cuboid<T,D> _motherCuboid;
  /// Vector of the cuboids
  std::vector<Cuboid<T,D>> _cuboids;
  /// Periodicity flag
  Vector<bool,D> _periodicityOn;
  /// class specific ostream
  mutable OstreamManager clout;

  /// Returns the minimum of the ratio nX/NY in the structure
  T getMinRatio() const;
  /// Returns the maximum of the ratio nX/NY in the structure
  T getMaxRatio() const;
  /// Returns the minimum coordinate in the structure
  Vector<T,D> getMinPhysR() const;
  /// Returns the maximum coordinate in the structure
  Vector<T,D> getMaxPhysR() const;
  /// Returns the minimum volume in the structure
  T getMinPhysVolume() const;
  /// Returns the maximum volume in the structure
  T getMaxPhysVolume() const;
  /// Returns the minimum number of nodes in the structure
  std::size_t getMinLatticeVolume() const;
  /// Returns the maximum number of nodes in the structure
  std::size_t getMaxLatticeVolume() const;
  /// Returns the minimum number of nodes in the structure inside the indicator
  std::size_t getMinLatticeWeight() const;
  /// Returns the maximum number of nodes in the structure inside the indicator
  std::size_t getMaxLatticeWeight() const;
  /// Returns the minimum delta in the structure
  T getMinDeltaR() const;
  /// Returns the maximum delta in the structure
  T getMaxDeltaR() const;

public:
  //CuboidGeometry();
  /// Constructs cuboid decomposition of cuboid with origin and extent
  CuboidGeometry(Vector<T,D> origin, T deltaR, Vector<int,D> extent, int nC=1);
  /// Construction from an alrealy existing mother cuboid
  CuboidGeometry(const Cuboid<T,D>& motherCuboid, int nC);
  /// Constructs a cuboid structure with a uniform spacing of voxelSize which consists of  nC cuboids, the cuboids not needed are removed and too big ones are shrinked
  CuboidGeometry(IndicatorF<T,D>& indicatorF, T voxelSize, int nC=1);
  /// Constructs a cuboid structure with a uniform spacing of voxelSize which consists of nC cuboids, the cuboids not needed are removed and too big ones are shrinked. Uses an iterative method: Halves the largest child cuboid until the set number of cuboids is reached. Largest cuboid is determined by either volume or weight as choosen in minimizeBy.
  CuboidGeometry(IndicatorF<T,D>& indicatorF, T voxelSize, int nC, std::string minimizeBy);

  int size() const override {
    return static_cast<int>(_cuboids.size());
  }

  const Cuboid<T,D>& get(int iC) const override {
    return _cuboids.at(iC);
  }

  std::optional<int> getC(Vector<T,D> physR, int offset = 0) const override {
    for (int iC = 0; iC < size(); ++iC) {
      if (get(iC).isInside(physR, offset)) {
        return iC;
      }
    }
    return std::nullopt;
  }

  Vector<T,D> getPhysR(LatticeR<D+1> latticeR) const override;

  std::optional<LatticeR<D+1>> getLatticeR(Vector<T,D> physR) const override {
    if (auto iC = getC(physR, 0)) {
      const auto& c = get(*iC);
      return LatticeR<D>{util::floor((physR - c.getOrigin()) / c.getDeltaR() + 0.5)}.withPrefix(*iC);
    }
    return std::nullopt;
  }

  std::optional<LatticeR<D+1>> getFloorLatticeR(Vector<T,D> physR) const override {
    if (auto iC = getC(physR, 0)) {
      const auto& c = get(*iC);
      return LatticeR<D>{util::floor((physR - c.getOrigin()) / c.getDeltaR())}.withPrefix(*iC);
    }
    return std::nullopt;
  }

  T getDeltaR() const override {
    return getMinDeltaR();
  }

  bool isInside(Vector<T,D> physR) const override {
    for (unsigned i = 0; i < _cuboids.size(); i++) {
      if (_cuboids[i].isInside(physR, 1)) {
        return true;
      }
    }
    return false;
  }

  /// Read and write access to a single cuboid
  Cuboid<T,D>& get(int iC) override;
  /// Returns the smallest cuboid that includes all cuboids of the structure
  const Cuboid<T,D>& getMotherCuboid() const override;

  /// Returns the smallest cuboid that includes all cuboids of the structure
  Cuboid<T,D>& getMotherCuboid();

  /// Set flag to enable/disable periodicity depending of direction. Be aware that not all directions are true to ensure boundary conditions like for velocity are not disturbed.
  void setPeriodicity(Vector<bool,D> periodicity) override;

  std::vector<Cuboid<T,D>>& cuboids() override {
    return _cuboids;
  }

  /// Returns set of neighbors to cuboid iCglob within overlap
  std::set<int> getNeighborhood(int iCglob, int overlap = 0) const override;

  /// Compares two CuboidGeometries
  bool operator==(CuboidGeometry<T,D>& rhs);

  /// Sets the number of full cells of each cuboid
  void setWeights(IndicatorF<T,D>& indicatorF) requires (D == 3);

  /// Splits cuboid iC, removes it and adds p cuboids of same volume
  void split(int iC, int p);
  /// Splits cuboid iC, removes it, adds approx. width^3 sized new cuboids
  void splitRegular(int iC, int width);
  /// Splits cuboid iC, removes it and adds p cuboids of same weight
  void splitByWeight(int iC, int p, IndicatorF<T,D>& indicatorF) requires (D == 3);
  /// Splits cuboid iC along dimension iD into cuboids of fractions
  void splitFractional(int iC, int iD, std::vector<T> fractions);

  /// Removes the cuboid iC
  void remove(int iC);
  /// Removes all cuboids where indicatorF = 0
  void remove(IndicatorF<T,D>& indicatorF);
  /// Removes all cuboids where weight = 0
  void removeByWeight();

  /// Shrink cuboid iC so that no empty planes are left
  void shrink(int iC, IndicatorF<T,D>& indicatorF);
  /// Shrink all cuboids so that no empty planes are left
  void shrink(IndicatorF<T,D>& indicatorF);

  /// Refines mesh by splitting each cell into factor^3 cells
  void refine(int factor);
  /// Tries to refine mesh to given deltaR
  bool tryRefineTo(T deltaR);

  /// Prints cuboid geometry details
  void print() const override;
  /// Prints cuboid geometry details plus details of all cuboids
  void printExtended();

  /// Save CuboidGeometry into an existing XML File
  void writeToExistingFile(std::string completeFileName, LoadBalancer<T>& loadBalancer);
  /// Save CuboidGeometry into XML File
  void writeToFile(std::string fileName, LoadBalancer<T>& loadBalancer);

private:
  /// Helper Function to create cuboid parameters for XML tag
  std::string _cuboidParameters(Cuboid<T,D> const& cub);

};


/// Helper Function to retrieve nData-dimensional std::vector of type S from space separated tag
template<typename S>
std::vector<S> getDataFromTag(XMLreader const& reader, std::string attrName, int nData)
{
  std::vector<S> values(nData, S());
  std::stringstream extstr(reader.getAttribute(attrName));
  for (auto& valueI: values) {
    extstr >> valueI;
  }
  return values;
}


/// Load CuboidGeometry from XML File
template<typename T, unsigned D>
CuboidGeometry<T,D>* createCuboidGeometry(std::string fileName)
{
  OstreamManager clout("saveCuboidGeometry");
  XMLreader reader(fileName);

  std::vector<T> origin = getDataFromTag<T>(reader["CuboidGeometry"], "origin", 3);
  std::vector<int> extent = getDataFromTag<int>(reader["CuboidGeometry"], "extent", 3);
  T deltaR = getDataFromTag<T>(reader["CuboidGeometry"], "deltaR", 1)[0];
  std::size_t weight = getDataFromTag<size_t>(reader["CuboidGeometry"], "weight", 1)[0];

  CuboidGeometry<T,D>* cGeo = new CuboidGeometry<T,D> (origin, deltaR, extent);
  cGeo->getMotherCuboid().setWeight(weight);
  cGeo->cuboids().clear();

  for ( XMLreader* cub: reader["CuboidGeometry"] ) {
    origin = getDataFromTag<T>(*cub, "origin", 3);
    extent = getDataFromTag<int>(*cub, "extent", 3);
    deltaR = getDataFromTag<T>(*cub, "deltaR", 1)[0];
    weight = getDataFromTag<int>(*cub, "weight", 1)[0];

    cGeo->cuboids().emplace_back(Cuboid<T,D>(origin, deltaR, extent));
    cGeo->get(cGeo->size() - 1).setWeight(weight);
  }

  return cGeo;
}

template <typename T>
using CuboidGeometry2D = CuboidGeometry<T,2>;
template <typename T>
using CuboidGeometry3D = CuboidGeometry<T,3>;

}  // namespace olb

#endif
