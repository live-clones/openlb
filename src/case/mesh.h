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

#ifndef CASE_MESH_H
#define CASE_MESH_H

#include "io/stlReader.h"

namespace olb {

template <typename T, unsigned D>
class Mesh {
private:
  std::unique_ptr<CuboidDecomposition<T,D>> _decomposition;
  std::unique_ptr<LoadBalancer<T>> _balancer;
  std::optional<unsigned> _overlap;

  /// Arbitrary indicators related to the mesh
  std::unordered_map<std::string, std::shared_ptr<IndicatorF<T,D>>> _indicators;

public:
  Mesh(std::unique_ptr<CuboidDecomposition<T, D>> decomposition,
       std::unique_ptr<LoadBalancer<T>> balancer)
    : _decomposition(std::move(decomposition)),
      _balancer(std::move(balancer))
  {
    OLB_PRECONDITION(_decomposition);
    OLB_PRECONDITION(_balancer);
  }

  template <typename... ARGS>
  explicit Mesh(ARGS&&... args)
    : _decomposition(new CuboidDecomposition<T,D>(std::forward<ARGS>(args)...)),
      _balancer(new HeuristicLoadBalancer<T>(*_decomposition))
  { }

  CuboidDecomposition<T,D>& getCuboidDecomposition() {
    return *_decomposition;
  }

  LoadBalancer<T>& getLoadBalancer() {
    return *_balancer;
  }

  void setOverlap(unsigned overlap) {
    _overlap = overlap;
  }

  unsigned getOverlap() const {
    if (_overlap) {
      return *_overlap;
    } else {
      throw std::logic_error("Mesh overlap not specified");
    }
  }

  T getDeltaX() const {
    return _decomposition->getDeltaX();
  }

  /// Stores indicator under name
  void addIndicator(std::string name, std::shared_ptr<IndicatorF<T,D>> indicatorF) {
    _indicators[name] = indicatorF;
  }

  /// Reads STL with the given settings and stores under its path
  std::shared_ptr<IndicatorF<T,D>> readSTL(std::string path, T deltaX, T scaling) {
    std::shared_ptr<IndicatorF<T,D>> stlI(new STLreader<T>(path, deltaX, scaling));
    addIndicator(path, stlI);
    return stlI;
  }

  /// Return indicator by name
  std::shared_ptr<IndicatorF<T,D>> getIndicator(std::string name) {
    return _indicators.at(name);
  }

  /// Easy construction of mesh from STL
  static Mesh fromSTL(std::string path, T physDeltaX, T scaling) {
    std::shared_ptr<STLreader<T>> stlI(new STLreader<T>(path, physDeltaX, scaling));
    IndicatorLayer3D<T> extendedDomain(*stlI, physDeltaX);
    Mesh mesh(extendedDomain, physDeltaX, singleton::mpi().getSize());
    mesh.addIndicator(path, stlI);
    return mesh;
  }

};

}

#endif
