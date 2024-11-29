namespace olb {

template <unsigned ndim, typename T, typename DESCRIPTOR>
class DampingTerm : public AnalyticalF<ndim,T,T> {
protected:
  Vector<T,ndim> x_0, Lx;
  int A=4;
  int n=4;
  T _boundary_depth_pu;
  T _damping_strength;
  UnitConverter<T, DESCRIPTOR> _converter;
public:
  DampingTerm(UnitConverter<T, DESCRIPTOR> converter, int boundary_depth_lu, Vector<T,ndim> domain_lengths, T damping_strength = 1. )
      : AnalyticalF<ndim,T,T>(1), _converter(converter), _damping_strength(damping_strength)
      {
        _boundary_depth_pu = converter.getPhysLength(boundary_depth_lu);
        for (size_t d=0; d<ndim; d++) {
          Lx[d] = domain_lengths[d]/2;  // Lx is usually half the domain
          x_0[d] = Lx[d] - _boundary_depth_pu;  // subtract boundary depth from domain length
          // x_0[d] /= Lx[d];  // normalize x_0 (Lx is still needed in operator() to normalize x; Lx will be replaced by 1)
        }
      };

  bool operator()(T output[], const T input[]) override
  {
    // normalize x to [0,1]
    Vector<T,ndim> x, distance_from_border;
    bool is_boundary = false;
    for (size_t d=0; d<ndim; d++) {
      x[d] = ( std::abs(input[d]) - x_0[d] ) / _boundary_depth_pu;  // x_0 becomes 0
      x[d] = std::max(x[d], T(0));  // if x[d]<0, ignore for X; if x[d]>1, set to one to avoid negative values of sigma
      if ( x[d] > 0 ) {
        is_boundary = true;
      }
      distance_from_border[d] = Lx[d] - std::abs(input[d]);
    }
    if ( !is_boundary ) { output[0] = 0; return true; }

    T X=0.;
    for ( size_t d=0; d<ndim; d++ ) {
      if ( distance_from_border[d] == *std::min_element(std::begin(distance_from_border), std::end(distance_from_border)) ) {
        X = x[d];
      }
    }
    T sigma;
    sigma = A * ( ( (std::pow(X, n))*(1-X)*(n+1)*(n+2) ) / std::pow(_boundary_depth_pu,n+2) );  // x_0 left out
    sigma /= 9.8304 / std::pow(_boundary_depth_pu, 6);  // calculated as max of function (at X=0.8) for A=n=4
    sigma = std::max(sigma, T(0));
    sigma = std::min(sigma, T(1));  // hard set 0<sigma<1
    sigma *= _damping_strength;
    output[0] = sigma;
    return true;
  };
};
 
}