/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
 *
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

#ifndef REFINEMENT_ALGORITHM_LAGRAVA_H
#define REFINEMENT_ALGORITHM_LAGRAVA_H

/// Implementation of the refinement algorithm by Lagrava et al.
/**
 * DOI: 10.1016/j.jcp.2012.03.015
 **/
namespace olb::refinement::lagrava {

struct HalfTimeCoarseToFineO {
  static constexpr refinement::OperatorScope scope = refinement::OperatorScope::PerFineCell;

  using parameters = meta::list<descriptors::TAU>;

  using data = meta::list<fields::refinement::NORMAL,
                          fields::refinement::PREV_RHO,
                          fields::refinement::PREV_U,
                          fields::refinement::PREV_FNEQ,
                          fields::refinement::CONTEXT_NEIGHBORS>;

  template <typename COARSE_CELL, typename FINE_CELL, typename DATA, typename PARAMETERS>
  void apply(COARSE_CELL& cCellPtr, FINE_CELL& fCell, DATA& data, PARAMETERS& params) any_platform {
    using V = typename COARSE_CELL::value_t;
    using DESCRIPTOR = typename COARSE_CELL::descriptor_t;
    if (cCellPtr) {
      auto cCell = *cCellPtr;

      auto rhoPrev  = data->template getField<fields::refinement::PREV_RHO>();
      auto uPrev    = data->template getField<fields::refinement::PREV_U>();
      auto fNeqPrev = data->template getField<fields::refinement::PREV_FNEQ>();

      V rhoCurr{};
      Vector<V,DESCRIPTOR::d> uCurr{};
      Vector<V,DESCRIPTOR::q> fNeqCurr{};
      lbm<DESCRIPTOR>::computeRhoU(cCell, rhoCurr, uCurr);
      lbm<DESCRIPTOR>::computeFneq(cCell, fNeqCurr, rhoCurr, uCurr);

      auto rho = V{0.5}*(rhoPrev + rhoCurr);
      auto u = V{0.5}*(uPrev + uCurr);
      auto fNeq = V{0.5}*(fNeqPrev + fNeqCurr);
      V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);

      V coarseTau = params.template get<descriptors::TAU>();
      V scalingFactor = (coarseTau - V{0.25}) / coarseTau;

      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        fCell[iPop] = equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr) + scalingFactor*fNeq[iPop];
      }
    } else {
      auto normal = data->template getField<fields::refinement::NORMAL>();

      V rho{};
      Vector<V,DESCRIPTOR::d> u{};
      Vector<V,DESCRIPTOR::q> fNeq{};

      int normalneighbors=0;
      int coarseneighbors=0;
      for (unsigned iN=0; iN < fields::refinement::CONTEXT_NEIGHBORS::count<DESCRIPTOR>(); ++iN) {
        auto n = fields::refinement::CONTEXT_NEIGHBORS::c<DESCRIPTOR>(iN);  // direction (vector, length d)
        if (n*normal == 0) {  // normal to edge --> there should be neighbor
          normalneighbors +=1;
          if ( auto ncCellPtr = cCellPtr.neighbor(n) ) {  // this is a coarse neighbor (not only fine)
            coarseneighbors +=1;
            auto ncCell = *ncCellPtr;
            auto nData = data.neighbor(iN);  // new data object with data in neighbor direction

            auto rhoPrev  = nData->template getField<fields::refinement::PREV_RHO>();
            auto uPrev    = nData->template getField<fields::refinement::PREV_U>();
            auto fNeqPrev = nData->template getField<fields::refinement::PREV_FNEQ>();

            V rhoCurr{};
            Vector<V,DESCRIPTOR::d> uCurr{};
            Vector<V,DESCRIPTOR::q> fNeqCurr{};
            lbm<DESCRIPTOR>::computeRhoU(ncCell, rhoCurr, uCurr);
            lbm<DESCRIPTOR>::computeFneq(ncCell, fNeqCurr, rhoCurr, uCurr);

            // time-interpoltation
            bool cartesian = false;
            if ( DESCRIPTOR::d == 2 ) cartesian = true;
            if ( DESCRIPTOR::d == 3 ) if ( abs(n[0])+abs(n[1])+abs(n[2]) == V(1.) ) cartesian = true;
            if ( cartesian ) {  // neighbor lies in cartesian direction -> 2nd order interpolation in 1D
              rho += V(9./16.)*V{0.5}*(rhoPrev + rhoCurr);
              u += V(9./16.)*V{0.5}*(uPrev + uCurr);
              fNeq += V(9./16.)*V{0.5}*(fNeqPrev + fNeqCurr);

              ncCellPtr = cCellPtr.neighbor(3*n);
              ncCell = *ncCellPtr;
              nData = data.neighbor(iN).neighbor(iN).neighbor(iN);

              rhoPrev  = nData->template getField<fields::refinement::PREV_RHO>();
              if (rhoPrev==0) std::cout << "rhoPrev==0! iN=" << iN << "; n=" << n << std::endl;
              uPrev    = nData->template getField<fields::refinement::PREV_U>();
              fNeqPrev = nData->template getField<fields::refinement::PREV_FNEQ>();

              lbm<DESCRIPTOR>::computeRhoU(ncCell, rhoCurr, uCurr);
              lbm<DESCRIPTOR>::computeFneq(ncCell, fNeqCurr, rhoCurr, uCurr);

              rho -= V(1./16.)*V{0.5}*(rhoPrev + rhoCurr);
              u -= V(1./16.)*V{0.5}*(uPrev + uCurr);
              fNeq -= V(1./16.)*V{0.5}*(fNeqPrev + fNeqCurr);
            } else {  // neighbor in diagonal direction -> 2nd order interpolation in 2D
              // Lagrava weights: 81/256, -9/256, -9/256, +1/256 (using 16 points)
              // alternative: 5/16 direct neighbors, -1/32 further neighbors - applied directly to populations?? (using 12 points)
              rho += V(81./256.)*V{0.5}*(rhoPrev + rhoCurr);
              u += V(81./256.)*V{0.5}*(uPrev + uCurr);
              fNeq += V(81./256.)*V{0.5}*(fNeqPrev + fNeqCurr);

              for ( size_t dim = 0; dim < 3; dim++ ) {
                if ( n[dim] != V(0.) ) {
                  auto n2 = n;
                  n2[dim] *= 2;
                  unsigned iN2=0;
                  while ( iN2 < fields::refinement::CONTEXT_NEIGHBORS::count<DESCRIPTOR>() ) {
                    auto n3 = fields::refinement::CONTEXT_NEIGHBORS::c<DESCRIPTOR>(iN2);  // direction (vector, length d)
                    if ( n+n3 == n2 ) break;
                    ++iN2;
                  }
                  n2[dim] = n[dim]*3;

                  ncCellPtr = cCellPtr.neighbor(n2);
                  ncCell = *ncCellPtr;
                  nData = data.neighbor(iN).neighbor(iN2).neighbor(iN2);

                  rhoPrev  = nData->template getField<fields::refinement::PREV_RHO>();
                  uPrev    = nData->template getField<fields::refinement::PREV_U>();
                  fNeqPrev = nData->template getField<fields::refinement::PREV_FNEQ>();

                  lbm<DESCRIPTOR>::computeRhoU(ncCell, rhoCurr, uCurr);
                  lbm<DESCRIPTOR>::computeFneq(ncCell, fNeqCurr, rhoCurr, uCurr);
                  
                  rho -= V(9./256.)*V{0.5}*(rhoPrev + rhoCurr);
                  u -= V(9./256.)*V{0.5}*(uPrev + uCurr);
                  fNeq -= V(9./256.)*V{0.5}*(fNeqPrev + fNeqCurr);
                }
              }

              ncCellPtr = cCellPtr.neighbor(3*n);
              ncCell = *ncCellPtr;
              nData = data.neighbor(iN).neighbor(iN).neighbor(iN);  // new data object with data in 3* neighbor direction

              rhoPrev  = nData->template getField<fields::refinement::PREV_RHO>();
              uPrev    = nData->template getField<fields::refinement::PREV_U>();
              fNeqPrev = nData->template getField<fields::refinement::PREV_FNEQ>();

              lbm<DESCRIPTOR>::computeRhoU(ncCell, rhoCurr, uCurr);
              lbm<DESCRIPTOR>::computeFneq(ncCell, fNeqCurr, rhoCurr, uCurr);
              
              rho += V(1./256.)*V{0.5}*(rhoPrev + rhoCurr);
              u += V(1./256.)*V{0.5}*(uPrev + uCurr);
              fNeq += V(1./256.)*V{0.5}*(fNeqPrev + fNeqCurr);
            }
          }
        }
      }

      V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);

      V coarseTau = params.template get<descriptors::TAU>();
      V scalingFactor = (coarseTau - V{0.25}) / coarseTau;

      if ( rho > 1. ) {
        std::cout << "rho=" << rho << "; neighbors=" << fields::refinement::CONTEXT_NEIGHBORS::count<DESCRIPTOR>() 
        << "; normalneighbors=" << normalneighbors << "; DESCRIPTOR::d=" << DESCRIPTOR::d << ""
        << std::endl;
      }
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        fCell[iPop] = equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr) + scalingFactor*fNeq[iPop];
      }
      lbm<DESCRIPTOR>::computeRhoU(fCell, rho, u);
      if (rho>1.) std::cout << "rho="<<rho<<",3"<<std::endl;
    }
  }
};

struct FullTimeCoarseToFineO {
  static constexpr refinement::OperatorScope scope = refinement::OperatorScope::PerFineCell;

  using parameters = meta::list<descriptors::TAU>;

  using data = meta::list<fields::refinement::PREV_RHO,
                          fields::refinement::PREV_U,
                          fields::refinement::PREV_FNEQ,
                          fields::refinement::NORMAL>;

  template <typename COARSE_CELL, typename FINE_CELL, typename DATA, typename PARAMETERS>
  void apply(COARSE_CELL& cCellPtr, FINE_CELL& fCell, DATA& data, PARAMETERS& params) any_platform {
    using V = typename COARSE_CELL::value_t;
    using DESCRIPTOR = typename COARSE_CELL::descriptor_t;
    if (cCellPtr) {
      auto cCell = *cCellPtr;

      V rho{};
      Vector<V,DESCRIPTOR::d> u{};
      Vector<V,DESCRIPTOR::q> fNeq{};
      lbm<DESCRIPTOR>::computeRhoU(cCell, rho, u);
      lbm<DESCRIPTOR>::computeFneq(cCell, fNeq, rho, u);
      V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);

      data->template setField<fields::refinement::PREV_RHO>(rho);
      data->template setField<fields::refinement::PREV_U>(u);
      data->template setField<fields::refinement::PREV_FNEQ>(fNeq);

      V coarseTau = params.template get<descriptors::TAU>();
      V scalingFactor = (coarseTau - V{0.25}) / coarseTau;

      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        fCell[iPop] = equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr) + scalingFactor*fNeq[iPop];
      }
    } else {
      auto normal = data->template getField<fields::refinement::NORMAL>();

      V rho{};
      Vector<V,DESCRIPTOR::d> u{};
      Vector<V,DESCRIPTOR::q> fNeq{};

      for (unsigned iN=0; iN < fields::refinement::CONTEXT_NEIGHBORS::count<DESCRIPTOR>(); ++iN) {
        auto n = fields::refinement::CONTEXT_NEIGHBORS::c<DESCRIPTOR>(iN);  // direction (vector, length d)
        if (n*normal == 0) {  // normal to edge --> there should be neighbor
          if ( auto ncCellPtr = cCellPtr.neighbor(n) ) {  // this is a coarse neighbor (not only fine)
            auto ncCell = *ncCellPtr;

            V rhoCurr{};
            Vector<V,DESCRIPTOR::d> uCurr{};
            Vector<V,DESCRIPTOR::q> fNeqCurr{};
            lbm<DESCRIPTOR>::computeRhoU(ncCell, rhoCurr, uCurr);
            lbm<DESCRIPTOR>::computeFneq(ncCell, fNeqCurr, rhoCurr, uCurr);

            // time-interpoltation
            bool cartesian = false;
            if ( DESCRIPTOR::d == 2 ) cartesian = true;
            if ( DESCRIPTOR::d == 3 ) if ( abs(n[0])+abs(n[1])+abs(n[2]) == V(1.) ) cartesian = true;
            if ( cartesian ) {  // neighbor lies in cartesian direction -> 2nd order interpolation in 1D
              rho += V(9./16.)*rhoCurr;
              u += V(9./16.)*uCurr;
              fNeq += V(9./16.)*fNeqCurr;

              ncCellPtr = cCellPtr.neighbor(3*n);
              ncCell = *ncCellPtr;

              lbm<DESCRIPTOR>::computeRhoU(ncCell, rhoCurr, uCurr);
              lbm<DESCRIPTOR>::computeFneq(ncCell, fNeqCurr, rhoCurr, uCurr);

              rho -= V(1./16.)*rhoCurr;
              u -= V(1./16.)*uCurr;
              fNeq -= V(1./16.)*fNeqCurr;
            } else {  // neighbor in diagonal direction -> 2nd order interpolation in 2D
              // Lagrava weights: 81/256, -9/256, -9/256, +1/256 (using 16 points)
              // alternative: 5/16 direct neighbors, -1/32 further neighbors - applied directly to populations?? (using 12 points)
              rho += V(81./256.)*rhoCurr;
              u += V(81./256.)*uCurr;
              fNeq += V(81./256.)*fNeqCurr;

              for ( size_t dim = 0; dim < 3; dim++ ) {
                if ( n[dim] != V(0.) ) {
                  auto n2 = n;
                  n2[dim] *= 2;
                  unsigned iN2=0;
                  while ( iN2 < fields::refinement::CONTEXT_NEIGHBORS::count<DESCRIPTOR>() ) {
                    auto n3 = fields::refinement::CONTEXT_NEIGHBORS::c<DESCRIPTOR>(iN2);  // direction (vector, length d)
                    if ( n+n3 == n2 ) break;
                    ++iN2;
                  }
                  n2[dim] = n[dim]*3;

                  ncCellPtr = cCellPtr.neighbor(n2);
                  ncCell = *ncCellPtr;

                  lbm<DESCRIPTOR>::computeRhoU(ncCell, rhoCurr, uCurr);
                  lbm<DESCRIPTOR>::computeFneq(ncCell, fNeqCurr, rhoCurr, uCurr);
                  
                  rho -= V(9./256.)*rhoCurr;
                  u -= V(9./256.)*uCurr;
                  fNeq -= V(9./256.)*fNeqCurr;
                }
              }

              ncCellPtr = cCellPtr.neighbor(3*n);
              ncCell = *ncCellPtr;

              lbm<DESCRIPTOR>::computeRhoU(ncCell, rhoCurr, uCurr);
              lbm<DESCRIPTOR>::computeFneq(ncCell, fNeqCurr, rhoCurr, uCurr);
              
              rho += V(1./256.)*rhoCurr;
              u += V(1./256.)*uCurr;
              fNeq += V(1./256.)*fNeqCurr;
            }
          }
        }
      }

      V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);

      V coarseTau = params.template get<descriptors::TAU>();
      V scalingFactor = (coarseTau - V{0.25}) / coarseTau;

      if (rho>1.) std::cout << "rho="<<rho<<",1"<<std::endl;
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        fCell[iPop] = equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr) + scalingFactor*fNeq[iPop];
      }
      lbm<DESCRIPTOR>::computeRhoU(fCell, rho, u);
      if (rho>1.) std::cout << "rho="<<rho<<",2"<<std::endl;

    }
  }

  template <typename COARSE_CELL, typename FINE_CELL, typename DATA, typename PARAMETERS>
  void initialize_prev(COARSE_CELL& cCellPtr, FINE_CELL& fCell, DATA& data, PARAMETERS& params) any_platform {
    using V = typename COARSE_CELL::value_t;
    using DESCRIPTOR = typename COARSE_CELL::descriptor_t;
    if (cCellPtr) {
      auto cCell = *cCellPtr;

      V rho{};
      Vector<V,DESCRIPTOR::d> u{};
      Vector<V,DESCRIPTOR::q> fNeq{};
      lbm<DESCRIPTOR>::computeRhoU(cCell, rho, u);
      lbm<DESCRIPTOR>::computeFneq(cCell, fNeq, rho, u);

      data->template setField<fields::refinement::PREV_RHO>(rho);
      data->template setField<fields::refinement::PREV_U>(u);
      data->template setField<fields::refinement::PREV_FNEQ>(fNeq);

      auto normal = data->template getField<fields::refinement::NORMAL>();
      for (unsigned iN=0; iN < fields::refinement::CONTEXT_NEIGHBORS::count<DESCRIPTOR>(); ++iN) {
        auto n = fields::refinement::CONTEXT_NEIGHBORS::c<DESCRIPTOR>(iN);  // direction (vector, length d)
        if (n*normal == 0) {  // normal to edge --> there should be neighbor
          if ( auto ncCellPtr = cCellPtr.neighbor(3*n) ) {  // this is a coarse neighbor hopefulle
            auto ncCell = *ncCellPtr;
            auto nData = data.neighbor(iN).neighbor(iN);
            V rho{};
            Vector<V,DESCRIPTOR::d> u{};
            Vector<V,DESCRIPTOR::q> fNeq{};
            lbm<DESCRIPTOR>::computeRhoU(ncCell, rho, u);
            lbm<DESCRIPTOR>::computeFneq(ncCell, fNeq, rho, u);
            data->template setField<fields::refinement::PREV_RHO>(rho);
            data->template setField<fields::refinement::PREV_U>(u);
            data->template setField<fields::refinement::PREV_FNEQ>(fNeq);
          }
        }
      }
    }
  }
};

struct FineToCoarseO {
  static constexpr refinement::OperatorScope scope = refinement::OperatorScope::PerCoarseCell;

  using parameters = meta::list<descriptors::TAU>;

  using data = meta::list<>;

  template <typename COARSE_CELL, typename FINE_CELL, typename DATA, typename PARAMETERS>
  void apply(COARSE_CELL& cCell, FINE_CELL& fCell, DATA& data, PARAMETERS& params) any_platform {
    using V = typename COARSE_CELL::value_t;
    using DESCRIPTOR = typename COARSE_CELL::descriptor_t;

    V rho{};
    Vector<V,DESCRIPTOR::d> u{};
    Vector<V,DESCRIPTOR::q> fNeq{};
    lbm<DESCRIPTOR>::computeRhoU(fCell, rho, u);
    lbm<DESCRIPTOR>::computeFneq(fCell, fNeq, rho, u);
    V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);

    for (int jPop=1; jPop < DESCRIPTOR::q; ++jPop) {
      auto fNeighbor = fCell.neighbor(descriptors::c<DESCRIPTOR>(jPop));
      Vector<V,DESCRIPTOR::q> fNeq_{};
      lbm<DESCRIPTOR>::computeFneq(fNeighbor, fNeq_);
      fNeq += fNeq_;
    }
    fNeq /= descriptors::q<DESCRIPTOR>();

    V coarseTau = params.template get<descriptors::TAU>();
    V scalingFactor = coarseTau / (coarseTau - V{0.25});

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cCell[iPop] = equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr) + scalingFactor*fNeq[iPop];
    }
  }
};

template <typename T, typename DESCRIPTOR>
std::unique_ptr<SuperLatticeRefinement<T,DESCRIPTOR>> makeCoarseToFineCoupler(
  SuperLattice<T,DESCRIPTOR>& sLatticeCoarse,
  SuperGeometry<T,DESCRIPTOR::d>& sGeometryCoarse,
  SuperLattice<T,DESCRIPTOR>& sLatticeFine,
  SuperGeometry<T,DESCRIPTOR::d>& sGeometryFine)
{
  auto& loadBalancerFine = sLatticeFine.getLoadBalancer();
  auto& cGeometryFine = sLatticeFine.getCuboidGeometry();
  const auto& converterCoarse = sLatticeCoarse.getConverter();

  SuperIndicatorDomainFrontierDistanceF<T,DESCRIPTOR> c2fFrontierI(cGeometryFine,
                                                                   sGeometryFine,
                                                                   0);
  SuperIndicatorDomainFrontierDistanceF<T,DESCRIPTOR> c2fInsideI(cGeometryFine,
                                                                 sGeometryFine,
                                                                 1);
  SuperIndicatorDomainFrontierDistanceF<T,DESCRIPTOR> c2fOutsideI(cGeometryFine,
                                                                  sGeometryFine,
                                                                  -1);

  auto coarseToFine = std::make_unique<SuperLatticeRefinement<T,DESCRIPTOR>>(sLatticeCoarse,
                                                                             sLatticeFine,
                                                                             loadBalancerFine);
  for (int iC = 0; iC < loadBalancerFine.size(); ++iC) {
    auto& fBlock = sLatticeFine.getBlock(iC);

    fBlock.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> fineLatticeR) {
      if (c2fFrontierI.getBlockIndicatorF(iC)(fineLatticeR)) {
        coarseToFine->getBlock(iC).add(fineLatticeR);
      }
    });

    auto& insideI  = c2fInsideI.getBlockIndicatorF(iC);
    auto& outsideI = c2fOutsideI.getBlockIndicatorF(iC);

    fBlock.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> fineLatticeR) {
      if (auto index = coarseToFine->getBlock(iC).getDataIndex(fineLatticeR)) {
        auto [type, normal] = computeBoundaryTypeAndNormal(insideI, outsideI, fineLatticeR);
        coarseToFine->getBlock(iC).getData()
                     .template getField<fields::refinement::NORMAL>().set(*index, normal);
      }
    });
    fBlock.defineDynamics(outsideI, meta::id<NoDynamics<T,DESCRIPTOR>>{});

    coarseToFine->getBlock(iC).getData()
                 .template getData<OperatorParameters<refinement::lagrava::HalfTimeCoarseToFineO>>()
                 .template set<descriptors::TAU>(converterCoarse.getLatticeRelaxationTime());
    coarseToFine->getBlock(iC).getData()
                 .template getData<OperatorParameters<refinement::lagrava::FullTimeCoarseToFineO>>()
                 .template set<descriptors::TAU>(converterCoarse.getLatticeRelaxationTime());

  }

  coarseToFine->setProcessingContext(ProcessingContext::Simulation);

  return coarseToFine;
}

template <typename T, typename DESCRIPTOR>
std::unique_ptr<SuperLatticeRefinement<T,DESCRIPTOR>> makeFineToCoarseCoupler(
  SuperLattice<T,DESCRIPTOR>& sLatticeCoarse,
  SuperGeometry<T,DESCRIPTOR::d>& sGeometryCoarse,
  SuperLattice<T,DESCRIPTOR>& sLatticeFine,
  SuperGeometry<T,DESCRIPTOR::d>& sGeometryFine)
{
  auto& loadBalancerFine = dynamic_cast<RefinedLoadBalancer<T,DESCRIPTOR::d>&>(
    sLatticeFine.getLoadBalancer());
  auto& cGeometryFine = sLatticeFine.getCuboidGeometry();
  const auto& converterCoarse = sLatticeCoarse.getConverter();

  SuperIndicatorDomainFrontierDistanceF<T,DESCRIPTOR> f2cFrontierI(cGeometryFine,
                                                                   sGeometryFine,
                                                                   2);

  auto fineToCoarse = std::make_unique<SuperLatticeRefinement<T,DESCRIPTOR>>(
    sLatticeCoarse, sLatticeFine, loadBalancerFine);
  for (int iC = 0; iC < loadBalancerFine.size(); ++iC) {
    auto& cBlock = sLatticeCoarse.getBlock(loadBalancerFine.cloc(iC));

    cBlock.forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> coarseLatticeR) {
      if (f2cFrontierI.getBlockIndicatorF(iC)(2*coarseLatticeR)) {
        fineToCoarse->getBlock(iC).add(2*coarseLatticeR);
      }
    });
    fineToCoarse->getBlock(iC).getData()
                 .template getData<OperatorParameters<refinement::lagrava::FineToCoarseO>>()
                 .template set<descriptors::TAU>(converterCoarse.getLatticeRelaxationTime());
  }

  // for (int iC = 0; iC < sLatticeCoarse.getLoadBalancer().size(); ++iC) {
  //   auto& cBlock = sLatticeCoarse.getBlock(iC);
  //   cBlock.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> coarseLatticeR) {
  //     // if (f2cFrontierI.getBlockIndicatorF(iC)(2*coarseLatticeR)) {
  //       std::cout << "coarseLatticeR=" << coarseLatticeR << std::endl;
  //       auto cCell = cBlock.get(coarseLatticeR);
  //       T rho{};
  //       Vector<T,DESCRIPTOR::d> u{};
  //       Vector<T,DESCRIPTOR::q> fNeq{};
  //       lbm<DESCRIPTOR>::computeRhoU(cCell, rho, u);
  //       lbm<DESCRIPTOR>::computeFneq(cCell, fNeq, rho, u);
  //       cCell.template setField<fields::refinement::PREV_RHO>(rho);
  //       cCell.template setField<fields::refinement::PREV_U>(u);
  //       cCell.template setField<fields::refinement::PREV_FNEQ>(fNeq);
  //     // }
  //   });
  // }

  fineToCoarse->setProcessingContext(ProcessingContext::Simulation);

  return fineToCoarse;
}

}

#endif
