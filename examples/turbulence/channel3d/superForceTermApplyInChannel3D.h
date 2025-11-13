/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Fedor Bukreev, Yuji (Sam) Shimojima
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
#ifndef SUPER_FORCE_TERM_APPLY_IN_CHANNEL3D_H
#define SUPER_FORCE_TERM_APPLY_IN_CHANNEL3D_H

namespace olb {
namespace descriptors {
namespace channel {
struct LATTICE_VOXCELSIZE : public descriptors::FIELD_BASE<1> {};
struct LATTICE_SURFACESIZE : public descriptors::FIELD_BASE<1> {};
struct CELL_LATTICEU_dVolume : public descriptors::FIELD_BASE<0, 1, 0> {};
struct CELL_UTAU_dSurface : public descriptors::FIELD_BASE<0, 1, 0> {};
struct TAGS_BULK : public TYPED_FIELD_BASE<bool, 1> {};
struct TAGS_BCBULK : public TYPED_FIELD_BASE<bool, 1> {};
struct LATTICE_TAU_W : public descriptors::FIELD_BASE<0, 1, 0> {};
struct LATTICE_U_TAU : public FIELD_BASE<0, 1, 0> {};

} // namespace channel
} // namespace descriptors

namespace stage::integral {
struct IntegralFrictionVelocity {};
struct IntegralVelocity {};
} // namespace stage::integral

struct CellSpatialCalculationVelocityField {
  static constexpr OperatorScope scope = OperatorScope::PerCell;
  int                            getPriority() const { return 2; }

  template <typename CELL>
  void apply(CELL& cell) any_platform
  {
    using DESCRIPTOR                   = typename CELL::descriptor_t;
    using V                            = typename CELL::value_t;
    V                        dV        = cell.template getField<descriptors::channel::LATTICE_VOXCELSIZE>();
    Vector<V, DESCRIPTOR::d> lattice_u = {};
    cell.computeU(lattice_u.data());
    /*
      UdV
    */
    Vector<V, DESCRIPTOR::d> lattice_u_dV = lattice_u * dV;
    cell.template setField<descriptors::channel::CELL_LATTICEU_dVolume>(lattice_u_dV);
  }
};

struct CellSpatialCalculationFrictionVelocityField {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters                     = meta::list<descriptors::OMEGA>;
  int getPriority() const { return 3; }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& params) any_platform
  {
    using DESCRIPTOR = typename CELL::descriptor_t;
    using V          = typename CELL::value_t;

    auto y1     = cell.template getField<descriptors::Y1>();
    V    y1Norm = util::norm<DESCRIPTOR::d>(y1);
    if (y1Norm != (V)0) {
      Vector<V, DESCRIPTOR::d> lattice_u = {};

      cell.computeU(lattice_u.data());
      Vector<V, 3> normal(y1[0] / y1Norm, y1[1] / y1Norm, y1[2] / y1Norm);
      auto u_t = lattice_u - (lattice_u[0] * normal[0] + lattice_u[1] * normal[1] + lattice_u[2] * normal[2]) * normal;
      auto tangent                   = u_t / util::norm<DESCRIPTOR::d>(u_t);
      const V                  omega = params.template get<descriptors::OMEGA>();
      const V                  nu    = ((V)1.0 / omega - (V)0.5) / descriptors::invCs2<V, DESCRIPTOR>();
      V                        u_1t = lattice_u[0] * tangent[0] + lattice_u[1] * tangent[1] + lattice_u[2] * tangent[2];
      Vector<V, DESCRIPTOR::d> u_tau = {};
      V                        rho   = cell.computeRho();
      u_tau[0]                       = util::sqrt((nu * (u_1t - (V)0.0) / y1Norm) / rho);
      V dS                           = cell.template getField<descriptors::channel::LATTICE_SURFACESIZE>();
      /*
         UTAUdS
      */
      Vector<V, DESCRIPTOR::d> u_tau_dS = u_tau * dS;
      cell.template setField<descriptors::channel::CELL_UTAU_dSurface>(u_tau_dS); // for reduce to calaculate forcing
      cell.template setField<descriptors::channel::LATTICE_U_TAU>(u_tau);
    }
  }
};

template <typename T, typename DESCRIPTOR>
class SuperForceTermApplyInChannel3d {
private:
  using _COUPLEE =
      meta::plain_map<meta::list<names::NavierStokes>, meta::list<descriptors::VALUED_DESCRIPTOR<T, DESCRIPTOR>>>;

  SuperLattice<T, DESCRIPTOR>&                             _sLattice;
  SuperLatticeCoupling<TurbulentChannelForce<T>, _COUPLEE> _coupling;

  SuperLatticeFieldReductionO<T, DESCRIPTOR, descriptors::channel::CELL_LATTICEU_dVolume, reduction::SumO,
                              reduction::checkTag<descriptors::channel::TAGS_BULK>>
      _sumUdV;

  SuperLatticeFieldReductionO<T, DESCRIPTOR, descriptors::channel::CELL_UTAU_dSurface, reduction::SumO,
                              reduction::checkTag<descriptors::channel::TAGS_BCBULK>>
      _sumUTAUdS;

  const T _ReTau;
  const T _charLength;
  const T _lengthX;
  const T _lengthY;
  const T _lengthZ;

  const int _lattice_channel_surfacearea;

public:
  FieldD<T, DESCRIPTOR, descriptors::channel::CELL_LATTICEU_dVolume> getVelocityIntegral()
  {
    /* ∫UdV */
    _sLattice.executePostProcessors(stage::integral::IntegralVelocity {});
    return _sumUdV.compute();
  }

  FieldD<T, DESCRIPTOR, descriptors::channel::CELL_LATTICEU_dVolume> getFrictionVelocityIntegral()
  {
    /* ∫UTAUdS */
    _sLattice.executePostProcessors(stage::integral::IntegralFrictionVelocity {});
    return _sumUTAUdS.compute();
  }

  void applyForing(FieldD<T, DESCRIPTOR, descriptors::channel::CELL_UTAU_dSurface> utauDs)
  {
    /*
      utauDs = ∫UTAUdS
      utau_ave = 1/S∫UTAUdS
    */
    const auto utau_ave = utauDs / _lattice_channel_surfacearea;
    _coupling.template setParameter<typename TurbulentChannelForce<T>::LATTICE_UTAU>(utau_ave[0]);
    _coupling.apply();
    return;
  }

  void applyForing()
  {
    /*
     In the coupling, 1/V ∫UdV will be calculated
   */
    _coupling.template setParameter<typename TurbulentChannelForce<T>::LATTICE_U_SUM>(getVelocityIntegral());
    _coupling.apply();
    return;
  }

  ~SuperForceTermApplyInChannel3d() {};

  /**
 * @brief Constructor for setting up the velocity integral and forcing coupling in a turbulent channel.
 * @details
 * Please look at `examples/turbulent/channel3d` for an example of how to use this class.
 * @param sLattice             SuperLattice
 * @param sGeometry            SuperGeometry.
 * @param converter            UnitConverter.
 * @param ReTau                Friction Reynolds number.
 * @param charLength           Calaculated characteristic physical.
 * @param lengthX              Domain physical length in the X (streamwise) direction.
 * @param lengthY              Domain physical length in the Y (spanwise) direction.
 * @param lengthZ              Domain physical length in the Z (wall-normal) direction.
 * @param latticeWallDistance  Distance from the wall in lattice units (default: 0.5).
 * @param boudaryMaterials     Material numbers identifying boundary regions (default: {2}).
 * @param bulkMaterials        Material numbers identifying bulk (fluid) regions (default: {1}).
 */
  SuperForceTermApplyInChannel3d(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T, 3>& sGeometry,
                                 UnitConverter<T, DESCRIPTOR> const& converter, const T ReTau, const T charLength,
                                 const T lengthX, const T lengthY, const T lengthZ, const T latticeWallDistance = 0.5,
                                 std::vector<int> boudaryMaterials = std::vector<int>(1, 2),
                                 std::vector<int> bulkMaterials    = std::vector<int>(1, 1))
      : _sLattice(sLattice)
      , _coupling(TurbulentChannelForce<T> {}, names::NavierStokes {}, sLattice)
      , _sumUdV(sLattice)
      , _sumUTAUdS(sLattice)
      , _ReTau(ReTau)
      , _charLength(charLength)
      , _lengthX(lengthX)
      , _lengthY(lengthY)
      , _lengthZ(lengthZ)
      , _lattice_channel_surfacearea(converter.getLatticeLength(lengthX) * converter.getLatticeLength(lengthY))
  {
    OstreamManager clout(std::cout, "SuperForceTermApplyInChannel3d");
    auto           bulkCopy = bulkMaterials;
    {
      clout << "set up the paramerters in the coupling between foucing calculation and LBM" << std::endl;

      _coupling.template setParameter<typename TurbulentChannelForce<T>::CHAR_LATTICE_U>(
          converter.getCharLatticeVelocity());
      _coupling.template setParameter<typename TurbulentChannelForce<T>::CHAR_LATTICE_L>(
          converter.getLatticeLength(_charLength));
      _coupling.template setParameter<typename TurbulentChannelForce<T>::LATTICE_UTAU>(
          converter.getLatticeVelocity(_ReTau * converter.getPhysViscosity() / ((T)0.5 * _lengthZ)));
      _coupling.template setParameter<typename TurbulentChannelForce<T>::LATTICE_CHANNEL_VOLUME>(
          converter.getLatticeLength(_lengthX) * converter.getLatticeLength(_lengthY) *
          converter.getLatticeLength(_lengthZ));
      _coupling.template setParameter<typename TurbulentChannelForce<T>::LATTICE_U_SUM>(
          Vector<T, DESCRIPTOR::d> {converter.getCharLatticeVelocity(), (T)0, (T)0});
      _coupling.restrictTo(sGeometry.getMaterialIndicator(std::move(bulkMaterials)));
    }
    {
      clout << "Set up the integral field and post processor" << std::endl;

      T                                               dx = converter.getPhysDeltaX(); // pyhsical length
      FunctorPtr<SuperIndicatorF<T, DESCRIPTOR::d>>&& superboundaryIndicator =
          FunctorPtr<SuperIndicatorF<T, DESCRIPTOR::d>>(sGeometry.getMaterialIndicator(std::move(boudaryMaterials)));
      FunctorPtr<SuperIndicatorF<T, DESCRIPTOR::d>>&& superbulkIndicator =
          FunctorPtr<SuperIndicatorF<T, DESCRIPTOR::d>>(sGeometry.getMaterialIndicator(std::move(bulkCopy)));

      auto& load = sLattice.getLoadBalancer();

      FunctorPtr<SuperIndicatorF<T, DESCRIPTOR::d>>&& outsideI =
          FunctorPtr<SuperIndicatorF<T, DESCRIPTOR::d>>(sGeometry.getMaterialIndicator(0));
      for (int iC = 0; iC < load.size(); ++iC) {
        auto&                              block             = sLattice.getBlock(iC);
        BlockIndicatorF<T, DESCRIPTOR::d>& boundaryIndicator = superboundaryIndicator->getBlockIndicatorF(iC);
        BlockIndicatorF<T, DESCRIPTOR::d>& bulkIndicator     = superbulkIndicator->getBlockIndicatorF(iC);
        BlockIndicatorF<T, DESCRIPTOR::d>& outsideIndicator  = outsideI->getBlockIndicatorF(iC);
        const auto& blockGeometry = superbulkIndicator->getBlockIndicatorF(iC).getBlockGeometry();
        const T     deltaR        = blockGeometry.getDeltaR(); // pyhsical voxels length
        block.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
          T    voxcelSize = 0.0; // lattice length
          auto cell       = sLattice.getBlock(iC).get(loc);
          if (bulkIndicator(loc)) {
            voxcelSize = (deltaR * deltaR * deltaR) / (dx * dx * dx);
            block.addPostProcessor(typeid(stage::integral::IntegralVelocity), loc,
                                   meta::id<CellSpatialCalculationVelocityField> {});
            cell.template setField<descriptors::channel::TAGS_BULK>(true);
          }
          cell.template setField<descriptors::channel::LATTICE_VOXCELSIZE>(voxcelSize);
        });
        //set voxel size for boundary cells
        block.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
          if (boundaryIndicator(loc)) {
            T voxcelSize = 0.0; // lattice length
            for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
              LatticeR<DESCRIPTOR::d> boundaryBulkLatticeR(loc + descriptors::c<DESCRIPTOR>(iPop));
              if (blockGeometry.isInside(boundaryBulkLatticeR)) {
                Vector<T, DESCRIPTOR::d> boundaryBulkPhysR = {};
                blockGeometry.getPhysR(boundaryBulkPhysR, boundaryBulkLatticeR);
                // check if neighbor is fluid cell
                if (bulkIndicator(boundaryBulkLatticeR)) {
                  auto cellBC = sLattice.getBlock(iC).get(boundaryBulkLatticeR);

                  const auto [normalType, normal] = computeBoundaryTypeAndNormal(bulkIndicator, outsideIndicator, loc);
                  T normalNorm                    = util::norm<DESCRIPTOR::d>(normal);
                  if (normalNorm != T(0)) {
                    T y1          = latticeWallDistance * deltaR;
                    voxcelSize    = (deltaR * deltaR * util::abs(y1)) / (dx * dx * dx);
                    T sarfaceSize = (deltaR * deltaR) / (dx * dx);
                    cellBC.template setField<descriptors::channel::LATTICE_VOXCELSIZE>(voxcelSize);
                    cellBC.template setField<descriptors::channel::LATTICE_SURFACESIZE>(sarfaceSize);
                    cellBC.template setField<descriptors::channel::TAGS_BCBULK>(true);
                    block.addPostProcessor(typeid(stage::integral::IntegralFrictionVelocity), boundaryBulkLatticeR,
                                           meta::id<CellSpatialCalculationFrictionVelocityField> {});
                  }
                }
              }
            }
          }
        });
      }
    }
    sLattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

    sLattice.setProcessingContext(ProcessingContext::Simulation);
  }
};
} // namespace olb
#endif // SUPER_FORCE_TERM_APPLY_IN_CHANNEL3D_H
