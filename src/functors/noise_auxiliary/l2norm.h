
namespace olb {

template<unsigned int ndim, typename T, typename DESCRIPTOR>

T L2Norm(SuperLattice<T,DESCRIPTOR>& sLattice,
         SuperGeometry<T,ndim>& sGeometry,
         UnitConverter<T,DESCRIPTOR> const& converter,
         int fluidMaterial
         ) {
  OstreamManager clout(std::cout, "L2Norm");
  SuperLatticePhysPressure3D<T,DESCRIPTOR> pressurenD( sLattice, converter );
  AnalyticalConst<ndim,T,T> rho0nD( 0. );
  SuperAbsoluteErrorL2Norm3D<T> absPressureErrorL2Norm( pressurenD, rho0nD, sGeometry.getMaterialIndicator( fluidMaterial ) );
  T result[ndim];
  int tmp[] = {int()};
  absPressureErrorL2Norm( result, tmp );
  return result[0];
};

}