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

#ifndef FSI_SUPER_POROUS_ELEMENT_FORCE_REDUCTION_O_H
#define FSI_SUPER_POROUS_ELEMENT_FORCE_REDUCTION_O_H

#include "fields.h"
#include "operators.h"

#include "core/superD.h"

namespace olb {

namespace stage {

namespace fsi {

struct CollectFromFluid { };

}

}

template <typename T, typename DESCRIPTOR>
class SuperPorousElementForceReductionO final {
private:
  /// Fluid lattice
  SuperLattice<T,DESCRIPTOR>& _sLattice;
  /// Reduced element data obtained from the fluid lattice
  std::unique_ptr<SuperD<T,descriptors::fsi::REDUCED_ELEMENTS<DESCRIPTOR::d>>> _superReducedElementsD;

  /// Reduced per-element surface force
  std::map<unsigned, FieldD<T,DESCRIPTOR,fields::fsi::ELEMENT_FORCE>>  _forces;
  /// Reduced per-element torque
  std::map<unsigned, FieldD<T,DESCRIPTOR,fields::fsi::ELEMENT_TORQUE>> _torques;

  std::size_t _nElements;

  bool _rankDoesFSI;

  #ifdef PARALLEL_MODE_MPI
  MPI_Comm _mpiCommunicator;
  #endif

public:
  SuperPorousElementForceReductionO(SuperLattice<T,DESCRIPTOR>& sLattice,
                                    FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicatorF):
    _sLattice(sLattice),
    _superReducedElementsD(new SuperD<
      T,descriptors::fsi::REDUCED_ELEMENTS<DESCRIPTOR::d>>(sLattice.getLoadBalancer())),
    _rankDoesFSI{singleton::mpi().isMainProcessor()}
  {
    OstreamManager clout(std::cout, "SuperPorousElementForceReductionO");
    auto& load = _sLattice.getLoadBalancer();

    //clout << "Set up operators and communication" << std::endl;

    // Schedule boundary force collection for FSI domain + overlap for block migration
    for (int iC=0; iC < _sLattice.getLoadBalancer().size(); ++iC) {
      auto& block = _sLattice.getBlock(iC);
      auto& blockF = indicatorF->getBlockIndicatorF(iC);
      if (!blockF.isEmpty()) {
        block.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> latticeR) {
          if (blockF(latticeR)) {
                   if (block.isInsideCore(latticeR)) {
              block.addPostProcessor(typeid(stage::fsi::CollectFromFluid),
                                     latticeR,
                                     meta::id<CollectPorousBoundaryForceO>{});
            } else if (block.isPadding(latticeR, 2)) {
              block.addPostProcessor(typeid(stage::fsi::CollectFromFluid),
                                     latticeR,
                                     meta::id<GrowPaddingLayerO>{});
            }
          }
        });
      }
    }
    sLattice.template addPostProcessor<stage::fsi::CollectFromFluid>(
      std::forward<decltype(indicatorF)>(indicatorF),
      meta::id<IntegratePorousBoundaryForceO>{});

    {
      auto& c = sLattice.getCommunicator(stage::fsi::CollectFromFluid{});
      c.template requestField<descriptors::POPULATION>();
      c.template requestField<descriptors::POROSITY>();
      c.template requestField<fields::fsi::ELEMENT_TAG>();
      c.requestOverlap(1, std::forward<decltype(indicatorF)>(indicatorF));
      c.exchangeRequests();
    }

    #ifdef PARALLEL_MODE_MPI
    for (int iC=0; iC < load.size(); ++iC) {
      _rankDoesFSI |= sLattice.getBlock(iC).hasPostProcessor(typeid(stage::fsi::CollectFromFluid),
                                                             meta::id<IntegratePorousBoundaryForceO>{});
    }
    MPI_Comm_split(MPI_COMM_WORLD, _rankDoesFSI ? 0 : MPI_UNDEFINED, singleton::mpi().getRank(), &_mpiCommunicator);
    #endif

    sLattice.template addPostProcessor<stage::Evaluation>(
      std::forward<decltype(indicatorF)>(indicatorF),
      meta::id<CollectPorousBoundaryForceO>{});
    {
      auto& c = sLattice.getCommunicator(stage::Evaluation{});
      c.template requestField<descriptors::POROSITY>();
      c.requestOverlap(1, std::forward<decltype(indicatorF)>(indicatorF));
      c.exchangeRequests();
    }

    {
      auto& c = sLattice.getCommunicator(stage::Full{});
      c.template requestField<descriptors::POROSITY>();
      c.template requestField<fields::fsi::ELEMENT_TAG>();
      c.template requestField<fields::fsi::BOUNDARY_FORCE>();
      c.template requestField<fields::fsi::BOUNDARY_TORQUE>();
      c.exchangeRequests();
    }

    //clout << "Set operator parameters" << std::endl;

    for (int iC=0; iC < load.size(); ++iC) {
      auto& block = _sLattice.getBlock(iC);
      auto& elementsBlock = _superReducedElementsD->getBlock(iC);
      block.template setParameter<fields::array_of<fields::fsi::REDUCED_ELEMENT_TAG>>(
        elementsBlock.template getField<fields::fsi::ELEMENT_TAG>());
      block.template setParameter<fields::array_of<fields::fsi::ELEMENT_FORCE>>(
        elementsBlock.template getField<fields::fsi::ELEMENT_FORCE>());
      block.template setParameter<fields::array_of<fields::fsi::ELEMENT_TORQUE>>(
        elementsBlock.template getField<fields::fsi::ELEMENT_TORQUE>());
    }
  }

  /// Resize reduction target buffer must be >= total global number of elements
  void resize(std::size_t nElements) {
    _nElements = nElements;
    for (int iC=0; iC < _sLattice.getLoadBalancer().size(); ++iC) {
      auto& elementsBlock = _superReducedElementsD->getBlock(iC);
      elementsBlock.template getField<fields::fsi::ELEMENT_TAG>().resize(nElements);
      elementsBlock.template getField<fields::fsi::ELEMENT_FORCE>().resize(nElements);
      elementsBlock.template getField<fields::fsi::ELEMENT_TORQUE>().resize(nElements);

      auto& block = _sLattice.getBlock(iC);
      block.template setParameter<fields::array_of<fields::fsi::REDUCED_ELEMENT_TAG>>(
        elementsBlock.template getField<fields::fsi::ELEMENT_TAG>());
      block.template setParameter<fields::array_of<fields::fsi::ELEMENT_FORCE>>(
        elementsBlock.template getField<fields::fsi::ELEMENT_FORCE>());
      block.template setParameter<fields::array_of<fields::fsi::ELEMENT_TORQUE>>(
        elementsBlock.template getField<fields::fsi::ELEMENT_TORQUE>());
    }
    _superReducedElementsD->setProcessingContext(ProcessingContext::Simulation);
  }

  bool rankDoesFSI() const
  {
    return _rankDoesFSI;
  }

  /// Globally integrate all element surface forces
  void apply()
  {
    OstreamManager clout(std::cout, "SuperPorousBoundaryForceReductionO");

    _sLattice.executePostProcessors(stage::fsi::CollectFromFluid{});

    if (!_rankDoesFSI) {
      return;
    }

    _forces.clear();
    _torques.clear();

    /// Gather local element forces
    int maxElement{};
    for (int iC=0; iC < _sLattice.getLoadBalancer().size(); ++iC) {
      auto& block = _superReducedElementsD->getBlock(iC);
      auto& blockTags = block.template getField<fields::fsi::ELEMENT_TAG>();
      blockTags.setProcessingContext(ProcessingContext::Evaluation);
      auto& blockForces = block.template getField<fields::fsi::ELEMENT_FORCE>();
      blockForces.setProcessingContext(ProcessingContext::Evaluation);
      auto& blockTorques = block.template getField<fields::fsi::ELEMENT_TORQUE>();
      blockTorques.setProcessingContext(ProcessingContext::Evaluation);

      auto& blockLattice = _sLattice.getBlock(iC);
      const std::size_t nElements = blockLattice.template getData<OperatorParameters<IntegratePorousBoundaryForceO>>()
                                                .template get<fields::fsi::REDUCED_ELEMENTS_COUNT>();
      for (std::size_t i=0; i < nElements; ++i) {
        const int iElement = blockTags[0][i];
        if (iElement > maxElement) {
          maxElement = iElement;
        }
        if (iElement > 0) {
          FieldD<T,DESCRIPTOR,fields::fsi::ELEMENT_FORCE> f([&](std::size_t iD) {
            return blockForces[iD][i];
          });
          _forces[iElement] += f;
          FieldD<T,DESCRIPTOR,fields::fsi::ELEMENT_TORQUE> t([&](std::size_t iD) {
            return blockTorques[iD][i];
          });
          _torques[iElement] += t;
        }
      }
    }

    #ifdef PARALLEL_MODE_MPI
    int nTotalElements{};
    singleton::mpi().allreduce(&maxElement, &nTotalElements, 1, MPI_MAX, _mpiCommunicator);
    if (nTotalElements >= 0) {
      _nElements = nTotalElements;
    }

    /// Aggregate global element forces
    std::vector<T> localForces(DESCRIPTOR::d*(_nElements+1), T{});
    for (auto [tag, f] : _forces) {
      for (std::size_t iD=0; iD < DESCRIPTOR::template size<fields::fsi::ELEMENT_FORCE>(); ++iD) {
        localForces[DESCRIPTOR::d*tag+iD] = f[iD];
      }
    }
    std::vector<T> globalForces(DESCRIPTOR::d*(_nElements+1), T{});
    singleton::mpi().allreduce(localForces.data(), globalForces.data(), globalForces.size(), MPI_SUM, _mpiCommunicator);
    for (std::size_t iElement=0; iElement < _nElements; ++iElement) {
      _forces[iElement+1] = globalForces.data() + DESCRIPTOR::d*(iElement+1);
    }
    std::vector<T> localTorques(DESCRIPTOR::template size<fields::fsi::ELEMENT_TORQUE>()*(_nElements+1), T{});
    for (auto [tag, f] : _torques) {
      for (std::size_t iD=0; iD < DESCRIPTOR::template size<fields::fsi::ELEMENT_TORQUE>(); ++iD) {
        localTorques[DESCRIPTOR::template size<fields::fsi::ELEMENT_TORQUE>()*tag+iD] = f[iD];
      }
    }
    std::vector<T> globalTorques(DESCRIPTOR::template size<fields::fsi::ELEMENT_TORQUE>()*(_nElements+1), T{});
    singleton::mpi().allreduce(localTorques.data(), globalTorques.data(), globalTorques.size(), MPI_SUM, _mpiCommunicator);
    for (std::size_t iElement=0; iElement < _nElements; ++iElement) {
      _torques[iElement+1] = globalTorques.data() + DESCRIPTOR::template size<fields::fsi::ELEMENT_TORQUE>()*(iElement+1);
    }
    #endif
  }

  std::size_t getElementCount() const {
    return _nElements;
  }

  auto getForce(unsigned iElement) const {
    return _forces.at(iElement);
  }

  auto getTorque(unsigned iElement) const {
    return _torques.at(iElement);
  }

};

template <typename T, typename DESCRIPTOR>
class SuperPorousElementStressReductionO final {
private:
  /// Fluid lattice
  SuperLattice<T,DESCRIPTOR>& _sLattice;
  /// Reduced element data obtained from the fluid lattice
  std::unique_ptr<SuperD<T,descriptors::fsi::REDUCED_ELEMENTS<DESCRIPTOR::d>>> _superReducedElementsD;

  /// Reduced per-element stress
  std::map<unsigned, FieldD<T,DESCRIPTOR,fields::fsi::ELEMENT_STRESS>> _stresses;

  std::size_t _nElements;

  bool _rankDoesFSI;

  #ifdef PARALLEL_MODE_MPI
  MPI_Comm _mpiCommunicator;
  #endif

public:
  SuperPorousElementStressReductionO(SuperLattice<T,DESCRIPTOR>& sLattice,
                                    FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicatorF):
    _sLattice(sLattice),
    _superReducedElementsD(new SuperD<
      T,descriptors::fsi::REDUCED_ELEMENTS<DESCRIPTOR::d>>(sLattice.getLoadBalancer())),
    _rankDoesFSI{singleton::mpi().isMainProcessor()}
  {
    OstreamManager clout(std::cout, "SuperPorousElementStressReductionO");
    auto& load = _sLattice.getLoadBalancer();

    clout << "Set up operators and communication" << std::endl;

    // Schedule boundary force collection for FSI domain + overlap for block migration
    for (int iC=0; iC < _sLattice.getLoadBalancer().size(); ++iC) {
      auto& block = _sLattice.getBlock(iC);
      auto& blockF = indicatorF->getBlockIndicatorF(iC);
      if (!blockF.isEmpty()) {
        block.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> latticeR) {
          if (blockF(latticeR)) {
                   if (block.isInsideCore(latticeR)) {
              block.addPostProcessor(typeid(stage::fsi::CollectFromFluid),
                                     latticeR,
                                     meta::id<CollectPorousBoundaryStressO>{});
            } else if (block.isPadding(latticeR, 2)) {
              block.addPostProcessor(typeid(stage::fsi::CollectFromFluid),
                                     latticeR,
                                     meta::id<GrowPaddingLayerO>{});
            }
          }
        });
      }
    }
    sLattice.template addPostProcessor<stage::fsi::CollectFromFluid>(
      std::forward<decltype(indicatorF)>(indicatorF),
      meta::id<IntegratePorousBoundaryStressO>{});

    {
      auto& c = sLattice.getCommunicator(stage::fsi::CollectFromFluid{});
      c.template requestField<descriptors::POPULATION>();
      c.template requestField<descriptors::POROSITY>();
      c.template requestField<fields::fsi::ELEMENT_TAG>();
      c.requestOverlap(1, std::forward<decltype(indicatorF)>(indicatorF));
      c.exchangeRequests();
    }

    #ifdef PARALLEL_MODE_MPI
    for (int iC=0; iC < load.size(); ++iC) {
      _rankDoesFSI |= sLattice.getBlock(iC).hasPostProcessor(typeid(stage::fsi::CollectFromFluid),
                                                             meta::id<IntegratePorousBoundaryStressO>{});
    }
    MPI_Comm_split(MPI_COMM_WORLD, _rankDoesFSI ? 0 : MPI_UNDEFINED, singleton::mpi().getRank(), &_mpiCommunicator);
    #endif

    sLattice.template addPostProcessor<stage::Evaluation>(
      std::forward<decltype(indicatorF)>(indicatorF),
      meta::id<CollectPorousBoundaryStressO>{});
    {
      auto& c = sLattice.getCommunicator(stage::Evaluation{});
      c.template requestField<descriptors::POROSITY>();
      c.requestOverlap(1, std::forward<decltype(indicatorF)>(indicatorF));
      c.exchangeRequests();
    }

    {
      auto& c = sLattice.getCommunicator(stage::Full{});
      c.template requestField<descriptors::POROSITY>();
      c.template requestField<fields::fsi::ELEMENT_TAG>();
      c.template requestField<fields::fsi::BOUNDARY_STRESS>();
      c.exchangeRequests();
    }

    //clout << "Set operator parameters" << std::endl;

    for (int iC=0; iC < load.size(); ++iC) {
      auto& block = _sLattice.getBlock(iC);
      auto& elementsBlock = _superReducedElementsD->getBlock(iC);
      block.template setParameter<fields::array_of<fields::fsi::REDUCED_ELEMENT_TAG>>(
        elementsBlock.template getField<fields::fsi::ELEMENT_TAG>());
      block.template setParameter<fields::array_of<fields::fsi::ELEMENT_STRESS>>(
        elementsBlock.template getField<fields::fsi::ELEMENT_STRESS>());
    }
  }

  /// Resize reduction target buffer must be >= total global number of elements
  void resize(std::size_t nElements) {
    _nElements = nElements;
    for (int iC=0; iC < _sLattice.getLoadBalancer().size(); ++iC) {
      auto& elementsBlock = _superReducedElementsD->getBlock(iC);
      elementsBlock.template getField<fields::fsi::ELEMENT_TAG>().resize(nElements);
      elementsBlock.template getField<fields::fsi::ELEMENT_STRESS>().resize(nElements);

      auto& block = _sLattice.getBlock(iC);
      block.template setParameter<fields::array_of<fields::fsi::REDUCED_ELEMENT_TAG>>(
        elementsBlock.template getField<fields::fsi::ELEMENT_TAG>());
      block.template setParameter<fields::array_of<fields::fsi::ELEMENT_STRESS>>(
        elementsBlock.template getField<fields::fsi::ELEMENT_STRESS>());
    }
    _superReducedElementsD->setProcessingContext(ProcessingContext::Simulation);
  }

  bool rankDoesFSI() const
  {
    return _rankDoesFSI;
  }

  /// Globally integrate all element surface stresses
  void apply()
  {
    OstreamManager clout(std::cout, "SuperPorousBoundaryStressReductionO");

    _sLattice.executePostProcessors(stage::fsi::CollectFromFluid{});

    if (!_rankDoesFSI) {
      return;
    }

    _stresses.clear();

    /// Gather local element stresses
    int maxElement{};
    for (int iC=0; iC < _sLattice.getLoadBalancer().size(); ++iC) {
      auto& block = _superReducedElementsD->getBlock(iC);
      auto& blockTags = block.template getField<fields::fsi::ELEMENT_TAG>();
      blockTags.setProcessingContext(ProcessingContext::Evaluation);
      auto& blockStress = block.template getField<fields::fsi::ELEMENT_STRESS>();
      blockStress.setProcessingContext(ProcessingContext::Evaluation);

      auto& blockLattice = _sLattice.getBlock(iC);
      const std::size_t nElements = blockLattice.template getData<OperatorParameters<IntegratePorousBoundaryStressO>>()
                                                .template get<fields::fsi::REDUCED_ELEMENTS_COUNT>();
      for (std::size_t i=0; i < nElements; ++i) {
        const int iElement = blockTags[0][i];
        if (iElement > maxElement) {
          maxElement = iElement;
        }
        if (iElement > 0) {
          FieldD<T,DESCRIPTOR,fields::fsi::ELEMENT_STRESS> f([&](std::size_t iD) {
            return blockStress[iD][i];
          });
          _stresses[iElement] += f;
        }
      }
    }

    #ifdef PARALLEL_MODE_MPI
    int nTotalElements{};
    singleton::mpi().allreduce(&maxElement, &nTotalElements, 1, MPI_MAX, _mpiCommunicator);
    if (nTotalElements >= 0) {
      _nElements = nTotalElements;
    }

    /// Aggregate global element forces
    std::vector<T> localStresses(DESCRIPTOR::template size<fields::fsi::ELEMENT_STRESS>()*(_nElements+1), T{});
    for (auto [tag, f] : _stresses) {
      for (std::size_t iD=0; iD < DESCRIPTOR::template size<fields::fsi::ELEMENT_STRESS>(); ++iD) {
        localStresses[DESCRIPTOR::template size<fields::fsi::ELEMENT_STRESS>()*tag+iD] = f[iD];
      }
    }
    std::vector<T> globalStresses(DESCRIPTOR::template size<fields::fsi::ELEMENT_STRESS>()*(_nElements+1), T{});
    singleton::mpi().allreduce(localStresses.data(), globalStresses.data(), globalStresses.size(), MPI_SUM, _mpiCommunicator);
    for (std::size_t iElement=0; iElement < _nElements; ++iElement) {
      _stresses[iElement+1] = globalStresses.data() + DESCRIPTOR::template size<fields::fsi::ELEMENT_STRESS>()*(iElement+1);
    }
    #endif
  }

  std::size_t getElementCount() const {
    return _nElements;
  }

  auto getStress(unsigned iElement) const {
    return _stresses.at(iElement);
  }

};

}

#endif
