/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Fedor Bukreev
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

#ifndef SET_TURBULENT_WALL_MODEL_3D_H
#define SET_TURBULENT_WALL_MODEL_3D_H

#include "postprocessor/turbulentWallModelPostProcessor3D.h"

namespace olb {

//======================================================================
// ======== Dynamics for WMLES ======//
//======================================================================

template <typename T, typename DESCRIPTOR>
using ForcedWMHRRdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
    momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>,
    momenta::BulkStress,
    momenta::DefineToNEq
  >,
  equilibria::ThirdOrder,
  collision::ParameterFromCell<collision::HYBRID,collision::LocalSmagorinskyEffectiveOmega<collision::HRR>>,
  forcing::WFGuoThirdOrder<momenta::ForcedWithStress>
>;

template <typename T, typename DESCRIPTOR>
using WMHRRdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
    momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>,
    momenta::BulkStress,
    momenta::DefineToNEq
  >,
  equilibria::ThirdOrder,
  collision::ParameterFromCell<collision::HYBRID,collision::LocalSmagorinskyEffectiveOmega<collision::HRR>>
>;

template <typename T, typename DESCRIPTOR>
using ForcedVanDriestWMHRRdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
    momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>,
    momenta::BulkStress,
    momenta::DefineToNEq
  >,
  equilibria::ThirdOrder,
  collision::ParameterFromCell<collision::HYBRID,collision::LocalVanDriestSmagorinskyEffectiveOmega<collision::HRR>>,
  forcing::WFGuoThirdOrder<momenta::ForcedWithStress>
>;

template <typename T, typename DESCRIPTOR>
using VanDriestWMHRRdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
    momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>,
    momenta::BulkStress,
    momenta::DefineToNEq
  >,
  equilibria::ThirdOrder,
  collision::ParameterFromCell<collision::HYBRID,collision::LocalVanDriestSmagorinskyEffectiveOmega<collision::HRR>>
>;

template <typename T, typename DESCRIPTOR>
using ForcedVanDriestExternalRhoWMHRRdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
    momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>,
    momenta::BulkStress,
    momenta::DefineToNEq
  >,
  equilibria::ThirdOrder,
  collision::ParameterFromCell<collision::HYBRID,
  collision::ParameterFromCell<collision::HYBRID_RHO,
  collision::LocalVanDriestSmagorinskyEffectiveOmega<collision::ExternalRhoHRR>>>,
  forcing::WFGuoThirdOrder<momenta::ForcedWithStress>
>;

template <typename T, typename DESCRIPTOR>
using VanDriestExternalRhoWMHRRdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
    momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>,
    momenta::BulkStress,
    momenta::DefineToNEq
  >,
  equilibria::ThirdOrder,
  collision::ParameterFromCell<collision::HYBRID,
  collision::ParameterFromCell<collision::HYBRID_RHO,
  collision::LocalVanDriestSmagorinskyEffectiveOmega<collision::ExternalRhoHRR>>>
>;

template <typename T, typename DESCRIPTOR>
using ForcedExternalRhoWMHRRdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
    momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>,
    momenta::BulkStress,
    momenta::DefineToNEq
  >,
  equilibria::ThirdOrder,
  collision::ParameterFromCell<collision::HYBRID,
  collision::ParameterFromCell<collision::HYBRID_RHO,
  collision::LocalSmagorinskyEffectiveOmega<collision::ExternalRhoHRR>>>,
  forcing::WFGuoThirdOrder<momenta::ForcedWithStress>
>;

template <typename T, typename DESCRIPTOR>
using ExternalRhoWMHRRdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
    momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>,
    momenta::BulkStress,
    momenta::DefineToNEq
  >,
  equilibria::ThirdOrder,
  collision::ParameterFromCell<collision::HYBRID,
  collision::ParameterFromCell<collision::HYBRID_RHO,
  collision::LocalSmagorinskyEffectiveOmega<collision::ExternalRhoHRR>>>
>;

//======================================================================
// ======== Wall model parameters ======//
//======================================================================

template <typename T>
struct WallModelParameters{
  /*  Used method for density reconstruction
   *  0: use local density
   *  1: extrapolation (Guo)
   *  2: constant (rho = 1.)
   */
  int rhoMethod = 0;

  /*  Used method for non-equilibrium population reconstruction
   *  0: extrapolation NEQ (Guo Zhaoli)
   *  1: first or second order finite differnce (Malaspinas)
   *  2: equilibrium scheme (no fNeq)
   */
  int fNeqMethod = 1;

  /*  Used wall profile
   *  0: power law profile
   *  1: Spalding profile
   */
  int wallFunctionProfile = 1;

  /// check if descriptor with body force is used
  bool bodyForce = true;

  /// interpolate sampling velocity along given normal between lattice voxels
  bool interpolateSampleVelocity = false;

  /// use van Driest damping function for turbulent viscosity in boundary cell
  bool useVanDriest = false;

  ///  distance from cell to real wall in lattice units if no geometry indicator is given as input
  T latticeWallDistance = T(0.5);

  ///  distance from cell to velocity sampling point in lattice units
  T samplingCellDistance = T(3.);
};

//======================================================================
// Set wall model dynamics depending on chosen paramertes.
// Use before the boundary setter.
//======================================================================

template<typename T, typename DESCRIPTOR>
void setTurbulentWallModelDynamics(SuperLattice<T, DESCRIPTOR>& sLattice,
                                   FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& bulkIndicator,
                                   WallModelParameters<T>& wallModelParameters)
{
  if( wallModelParameters.bodyForce ) {
    if( wallModelParameters.useVanDriest ) {
      if( wallModelParameters.rhoMethod != 0 ) {
        sLattice.template defineDynamics<ForcedVanDriestExternalRhoWMHRRdynamics>(std::move(bulkIndicator));
      } else {
        sLattice.template defineDynamics<ForcedVanDriestWMHRRdynamics>(std::move(bulkIndicator));
      }
    }
    else {
      if( wallModelParameters.rhoMethod != 0 ) {
        sLattice.template defineDynamics<ForcedExternalRhoWMHRRdynamics>(std::move(bulkIndicator));
      } else {
        sLattice.template defineDynamics<ForcedWMHRRdynamics>(std::move(bulkIndicator));
      }
    }
  }
  else {
    if( wallModelParameters.useVanDriest ) {
      if( wallModelParameters.rhoMethod != 0 ) {
        sLattice.template defineDynamics<VanDriestExternalRhoWMHRRdynamics>(std::move(bulkIndicator));
      } else {
        sLattice.template defineDynamics<VanDriestWMHRRdynamics>(std::move(bulkIndicator));
      }
    }
    else {
      if( wallModelParameters.rhoMethod != 0 ) {
        sLattice.template defineDynamics<ExternalRhoWMHRRdynamics>(std::move(bulkIndicator));
      } else {
        sLattice.template defineDynamics<WMHRRdynamics>(std::move(bulkIndicator));
      }
    }
  }

  AnalyticalConst3D<T, T> one(T(1.));
  AnalyticalConst3D<T, T> zero(T(0.));
  AnalyticalConst3D<T, T> zeroVector(T(0.), T(0.), T(0.));
  std::vector<T> zeroStrainVec = {T(0), T(0), T(0), T(0), T(0), T(0)};
  AnalyticalConst3D<T, T> zeroStrain(zeroStrainVec);
  sLattice.template defineField<descriptors::STRAINRATE>(std::move(bulkIndicator), zeroStrain);
  sLattice.template defineField<descriptors::OMEGA>(std::move(bulkIndicator), one);
  sLattice.template defineField<descriptors::POROSITY>(std::move(bulkIndicator), one);
  sLattice.template defineField<descriptors::U_TAU>(std::move(bulkIndicator), zero);
  sLattice.template defineField<collision::HYBRID>(std::move(bulkIndicator), one);
  if( wallModelParameters.rhoMethod != 0) {
    sLattice.template defineField<collision::HYBRID_RHO>(std::move(bulkIndicator), one);
    sLattice.template defineField<descriptors::DENSITY>(std::move(bulkIndicator), one);
  }
  sLattice.template defineField<descriptors::Y1>(std::move(bulkIndicator), zeroVector);
  sLattice.template defineField<descriptors::VELOCITY>(std::move(bulkIndicator), zeroVector);
  if( wallModelParameters.useVanDriest ) {
    sLattice.template defineField<descriptors::VISCOSITY>(std::move(bulkIndicator), zero);
  }

  /*{
    auto& communicator = sLattice.getCommunicator(stage::PostStream());
    if( wallModelParameters.fNeqMethod == 1 ) {
      communicator.template requestField<descriptors::VELOCITY>();
    }
    //communicator.template requestField<descriptors::OMEGA>();
    //communicator.template requestField<descriptors::Y1>();
    //communicator.template requestField<descriptors::STRAINRATE>();
    communicator.template requestField<descriptors::POPULATION>();
    //communicator.template requestField<collision::HYBRID>();
    //communicator.template requestField<descriptors::POROSITY>();
    if( wallModelParameters.useVanDriest ) {
      communicator.template requestField<descriptors::VISCOSITY>();
    }
    if( wallModelParameters.rhoMethod != 0 ) {
      communicator.template requestField<descriptors::DENSITY>();
    }
    communicator.requestOverlap(sLattice.getOverlap());
    communicator.exchangeRequests();
  }*/
}

template<typename T, typename DESCRIPTOR>
void setTurbulentWallModelDynamics(SuperLattice<T, DESCRIPTOR>& sLattice,
                                   SuperGeometry<T,DESCRIPTOR::d>& superGeometry,
                                   int materialOfBulk,
                                   WallModelParameters<T>& wallModelParameters)
{
  // Getting the indicators by material numbers and calling the superLattice method via the indicators:
  setTurbulentWallModelDynamics<T,DESCRIPTOR>(sLattice,
                        FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>(superGeometry.getMaterialIndicator(materialOfBulk)),
                        wallModelParameters);
}

//======================================================================
// Setter of the wall model boundary.
// Can be used with or without geometry indicator input.
//======================================================================

template<typename T, typename DESCRIPTOR>
void setTurbulentWallModel(SuperLattice<T, DESCRIPTOR>& sLattice,
                           FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& boundaryIndicator,
                           FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& bulkIndicator,
                           WallModelParameters<T>& wallModelParameters,
                           IndicatorF<T,DESCRIPTOR::d>* indicatorAnalyticalBoundary = nullptr)
{
  int _overlap = 3;
  OstreamManager clout(std::cout, "TurbulentWallModelSetter");

  AnalyticalConst3D<T, T> one(T(1.));
  AnalyticalConst3D<T, T> zero(T(0.));
  AnalyticalConst3D<T, T> zeroVector(T(0.), T(0.), T(0.));
  sLattice.template defineField<descriptors::OMEGA>(std::move(boundaryIndicator), one);
  sLattice.template defineField<descriptors::POROSITY>(std::move(boundaryIndicator), zero);
  sLattice.template defineField<descriptors::U_TAU>(std::move(boundaryIndicator), zero);
  sLattice.template defineField<collision::HYBRID>(std::move(boundaryIndicator), one);
  if( wallModelParameters.rhoMethod != 0) {
    sLattice.template defineField<collision::HYBRID_RHO>(std::move(boundaryIndicator), one);
    sLattice.template defineField<descriptors::DENSITY>(std::move(boundaryIndicator), one);
  }
  sLattice.template defineField<descriptors::Y1>(std::move(boundaryIndicator), zeroVector);
  sLattice.template defineField<descriptors::VELOCITY>(std::move(boundaryIndicator), zeroVector);
  if( wallModelParameters.useVanDriest ) {
    sLattice.template defineField<descriptors::VISCOSITY>(std::move(bulkIndicator), zero);
  }

  auto& load = sLattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    setTurbulentWallModel<T,DESCRIPTOR>(sLattice.getBlock(iC),
                          (bulkIndicator->getBlockIndicatorF(iC)).getBlockGeometry(),
                          boundaryIndicator->getBlockIndicatorF(iC),
                          bulkIndicator->getBlockIndicatorF(iC),
                          wallModelParameters,
                          indicatorAnalyticalBoundary);
  }
  addPoints2CommBC<T,DESCRIPTOR>(sLattice, std::forward<decltype(boundaryIndicator)>(boundaryIndicator), _overlap);
  clout.setMultiOutput(false);
  sLattice.template setParameter<descriptors::SAMPLING_DISTANCE>( wallModelParameters.samplingCellDistance );
}

/// Set parameters of the turbulent wall model after Bouzidi bounce-back on material cells of sLattice
template<typename T, typename DESCRIPTOR>
void setTurbulentWallModel(SuperLattice<T, DESCRIPTOR>& sLattice,
                           SuperGeometry<T,DESCRIPTOR::d>& superGeometry,
                           int materialOfSolidObstacle,
                           WallModelParameters<T>& wallModelParameters,
                           IndicatorF<T,DESCRIPTOR::d>* indicatorAnalyticalBoundary = nullptr,
                           std::vector<int> bulkMaterials = std::vector<int>(1,1))
{
  // Getting the indicators by material numbers and calling the superLattice method via the indicators:
  setTurbulentWallModel<T,DESCRIPTOR>(sLattice,
                      FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>(superGeometry.getMaterialIndicator(materialOfSolidObstacle)),
                      FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>(superGeometry.getMaterialIndicator(std::move(bulkMaterials))),
                      wallModelParameters,
                      indicatorAnalyticalBoundary);
}

/// Set parameters of the turbulent wall model after Bouzidi bounce-back on indicated cells of block lattice
template<typename T, typename DESCRIPTOR>
void setTurbulentWallModel(BlockLattice<T,DESCRIPTOR>& block,
                           BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
                           BlockIndicatorF<T,DESCRIPTOR::d>& boundaryIndicator,
                           BlockIndicatorF<T,DESCRIPTOR::d>& bulkIndicator,
                           WallModelParameters<T>& wallModelParameters,
                           IndicatorF<T,DESCRIPTOR::d>* indicatorAnalyticalBoundary = nullptr)
{
  OstreamManager clout(std::cout, "TurbulentWallModelSetter");
  clout.setMultiOutput(true);

  // Defining boundary distance y1 and distance to exchange location y2 in wall normal direction for every boundary cell
  const T deltaR = blockGeometry.getDeltaR();
  // for each solid cell:
  int kMax = 2;
  if( wallModelParameters.fNeqMethod == 1 ) {
    kMax = 3;
  }
  for(int k=1; k < kMax; k++){
    block.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> solidLatticeR) {
      // Check if cell is solid cell
      if (boundaryIndicator(solidLatticeR)) {
        for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
          Vector<T,DESCRIPTOR::d> boundaryLatticeR(solidLatticeR + k*descriptors::c<DESCRIPTOR>(iPop));
          if (blockGeometry.getNeighborhoodRadius(boundaryLatticeR) >= 1) {
            if (blockGeometry.isInside(boundaryLatticeR)) {
              Vector<T,DESCRIPTOR::d> boundaryPhysR{ };
              blockGeometry.getPhysR(boundaryPhysR,boundaryLatticeR);
              // check if neighbor is fluid cell
              if (bulkIndicator(boundaryLatticeR)) {
                if (indicatorAnalyticalBoundary) {
                  // Calculate surface normal
                  Vector<T,3> normal = indicatorAnalyticalBoundary->surfaceNormal(boundaryPhysR, deltaR);
                  T y1 = 0.;
                  // Calculate boundary distance y1 in surface normal direction
                  indicatorAnalyticalBoundary->distance(y1, boundaryPhysR, normal, blockGeometry.getIcGlob());
                  y1 /= deltaR; // y1 in lattice units
                  if ( util::abs(y1) > T(k) ) {
                    y1 = util::sign(y1) * T(k);
                  }
                  block.get(boundaryLatticeR).template setField<descriptors::Y1>(-y1*normal);
                } else {
                  Vector<T,3> normal(T(0),T(0),T(0));
                  for(int jPop = 0; jPop < DESCRIPTOR::q; jPop++){
                    Vector<T,DESCRIPTOR::d> boundaryLatticeR2(boundaryLatticeR + k*descriptors::c<DESCRIPTOR>(jPop));
                    if(boundaryIndicator(boundaryLatticeR2)){
                      normal[0] += k*descriptors::c<DESCRIPTOR>(jPop,0);
                      normal[1] += k*descriptors::c<DESCRIPTOR>(jPop,1);
                      normal[2] += k*descriptors::c<DESCRIPTOR>(jPop,2);
                    }
                  }
                  T normalNorm = util::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
                  if(normalNorm != T(0)) {
                    normal[0] /= normalNorm;
                    normal[1] /= normalNorm;
                    normal[2] /= normalNorm;
                  }
                  for(int i=0; i<3; i++) {
                    if(util::abs(normal[i]) < T(0.4)) {
                      normal[i] = T(0);
                    }
                  }
                  normalNorm = util::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
                  normal[0] /= normalNorm;
                  normal[1] /= normalNorm;
                  normal[2] /= normalNorm;
                  if(normal[0] != T(0) && normal[0] > T(0)) {
                    normal[0] /= normal[0];
                  }
                  else if(normal[0] != T(0) && normal[0] < T(0)) {
                    normal[0] /= normal[0];
                    normal[0] *= T(-1.);
                  }
                  if(normal[1] != T(0) && normal[1] > T(0)) {
                    normal[1] /= normal[1];
                  }
                  else if(normal[1] != T(0) && normal[1] < T(0)) {
                    normal[1] /= normal[1];
                    normal[1] *= T(-1.);
                  }
                  if(normal[2] != T(0) && normal[2] > T(0)) {
                    normal[2] /= normal[2];
                  }
                  else if(normal[2] != T(0) && normal[2] < T(0)) {
                    normal[2] /= normal[2];
                    normal[2] *= T(-1.);
                  }
                  if(util::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]) > T(1)) {
                    normal[0] = T(0);
                    normal[1] = T(0);
                    normal[2] = T(0);
                  }
                  if(normalNorm != T(0)) {
                    auto field = block.get(boundaryLatticeR).template getField<descriptors::Y1>();
                    if( util::norm<DESCRIPTOR::d>(field) == T(0)) {
                      T y1 = (k-1) + wallModelParameters.latticeWallDistance;
                      block.get(boundaryLatticeR).template setField<descriptors::Y1>(-y1*normal);
                    }
                  }
                }
                T pi[util::TensorVal<DESCRIPTOR>::n] {T(0)};
                block.get(boundaryLatticeR).template setField<descriptors::STRAINRATE>(pi);
                // Setting turbulent wall model post processor
                if( wallModelParameters.bodyForce ) {
                  if( wallModelParameters.interpolateSampleVelocity ) {
                    if( wallModelParameters.useVanDriest && k == 1 ) {
                      if( wallModelParameters.wallFunctionProfile == 0 ) {
                        block.addPostProcessor(typeid(stage::PostStream),
                                              boundaryLatticeR,
                                              meta::id<TurbulentWallModelPostProcessor<true,true,true,0>>{});
                      } else {
                        block.addPostProcessor(typeid(stage::PostStream),
                                              boundaryLatticeR,
                                              meta::id<TurbulentWallModelPostProcessor<true,true,true,1>>{});
                      }
                    } else {
                      if( wallModelParameters.wallFunctionProfile == 0 ) {
                        block.addPostProcessor(typeid(stage::PostStream),
                                              boundaryLatticeR,
                                              meta::id<TurbulentWallModelPostProcessor<true,true,false,0>>{});
                      } else {
                        block.addPostProcessor(typeid(stage::PostStream),
                                              boundaryLatticeR,
                                              meta::id<TurbulentWallModelPostProcessor<true,true,false,1>>{});
                      }
                    }
                  } else {
                    if( wallModelParameters.useVanDriest && k == 1 ) {
                      if( wallModelParameters.wallFunctionProfile == 0 ) {
                        block.addPostProcessor(typeid(stage::PostStream),
                                              boundaryLatticeR,
                                              meta::id<TurbulentWallModelPostProcessor<true,false,true,0>>{});
                      } else {
                        block.addPostProcessor(typeid(stage::PostStream),
                                              boundaryLatticeR,
                                              meta::id<TurbulentWallModelPostProcessor<true,false,true,1>>{});
                      }
                    } else {
                      if( wallModelParameters.wallFunctionProfile == 0 ) {
                        block.addPostProcessor(typeid(stage::PostStream),
                                              boundaryLatticeR,
                                              meta::id<TurbulentWallModelPostProcessor<true,false,false,0>>{});
                      } else {
                        block.addPostProcessor(typeid(stage::PostStream),
                                              boundaryLatticeR,
                                              meta::id<TurbulentWallModelPostProcessor<true,false,false,1>>{});
                      }
                    }
                  }
                } else {
                  if( wallModelParameters.interpolateSampleVelocity ) {
                    if( wallModelParameters.useVanDriest && k == 1 ) {
                      if( wallModelParameters.wallFunctionProfile == 0 ) {
                        block.addPostProcessor(typeid(stage::PostStream),
                                              boundaryLatticeR,
                                              meta::id<TurbulentWallModelPostProcessor<false,true,true,0>>{});
                      } else {
                        block.addPostProcessor(typeid(stage::PostStream),
                                              boundaryLatticeR,
                                              meta::id<TurbulentWallModelPostProcessor<false,true,true,1>>{});
                      }
                    } else {
                      if( wallModelParameters.wallFunctionProfile == 0 ) {
                        block.addPostProcessor(typeid(stage::PostStream),
                                              boundaryLatticeR,
                                              meta::id<TurbulentWallModelPostProcessor<false,true,false,0>>{});
                      } else {
                        block.addPostProcessor(typeid(stage::PostStream),
                                              boundaryLatticeR,
                                              meta::id<TurbulentWallModelPostProcessor<false,true,false,1>>{});
                      }
                    }
                  } else {
                    if( wallModelParameters.useVanDriest && k == 1 ) {
                      if( wallModelParameters.wallFunctionProfile == 0 ) {
                        block.addPostProcessor(typeid(stage::PostStream),
                                              boundaryLatticeR,
                                              meta::id<TurbulentWallModelPostProcessor<false,false,true,0>>{});
                      } else {
                        block.addPostProcessor(typeid(stage::PostStream),
                                              boundaryLatticeR,
                                              meta::id<TurbulentWallModelPostProcessor<false,false,true,1>>{});
                      }
                    } else {
                      if( wallModelParameters.wallFunctionProfile == 0 ) {
                        /*block.addPostProcessor(typeid(stage::PostStream),
                                              boundaryLatticeR,
                                              meta::id<TurbulentWallModelPostProcessor<false,false,false,0>>{});*/
                        block.addPostProcessor(typeid(stage::PostStream),
                                              boundaryLatticeR,
                                              meta::id<TurbulentSimpleWallModelPostProcessor>{});
                      } else {
                        block.addPostProcessor(typeid(stage::PostStream),
                                              boundaryLatticeR,
                                              meta::id<TurbulentWallModelPostProcessor<false,false,false,1>>{});
                      }
                    }
                  }
                }
                /*if( k == 1 ) {
                  if( wallModelParameters.fNeqMethod == 0 ) {
                    if( wallModelParameters.rhoMethod == 0 ) {
                      block.addPostProcessor(typeid(stage::PostStream),
                                            boundaryLatticeR,
                                            meta::id<TurbulentWallModelFneqGuoPostProcessor<0>>{});
                    }
                    else if( wallModelParameters.rhoMethod == 1 ) {
                      block.addPostProcessor(typeid(stage::PostStream),
                                            boundaryLatticeR,
                                            meta::id<TurbulentWallModelFneqGuoPostProcessor<1>>{});
                    }
                    else {
                      block.addPostProcessor(typeid(stage::PostStream),
                                            boundaryLatticeR,
                                            meta::id<TurbulentWallModelFneqGuoPostProcessor<2>>{});
                    }
                  }
                  else if( wallModelParameters.fNeqMethod == 1 ) {
                    if( wallModelParameters.rhoMethod == 0 ) {
                      block.addPostProcessor(typeid(stage::PostStream),
                                            boundaryLatticeR,
                                            meta::id<TurbulentWallModelFneqFDMPostProcessor<0>>{});
                    }
                    else if( wallModelParameters.rhoMethod == 1 ) {
                      block.addPostProcessor(typeid(stage::PostStream),
                                            boundaryLatticeR,
                                            meta::id<TurbulentWallModelFneqFDMPostProcessor<1>>{});
                    }
                    else {
                      block.addPostProcessor(typeid(stage::PostStream),
                                            boundaryLatticeR,
                                            meta::id<TurbulentWallModelFneqFDMPostProcessor<2>>{});
                    }
                  }
                  else {
                    if( wallModelParameters.rhoMethod == 0 ) {
                      block.addPostProcessor(typeid(stage::PostStream),
                                            boundaryLatticeR,
                                            meta::id<TurbulentWallModelFneqZeroPostProcessor<0>>{});
                    }
                    else if( wallModelParameters.rhoMethod == 1 ) {
                      block.addPostProcessor(typeid(stage::PostStream),
                                            boundaryLatticeR,
                                            meta::id<TurbulentWallModelFneqZeroPostProcessor<1>>{});
                    }
                    else {
                      block.addPostProcessor(typeid(stage::PostStream),
                                            boundaryLatticeR,
                                            meta::id<TurbulentWallModelFneqZeroPostProcessor<2>>{});
                    }
                  }
                }*/
              }
            }
          }
        }
      }
    });
  }
}
}
#endif
