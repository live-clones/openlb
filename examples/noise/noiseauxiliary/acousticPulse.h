namespace olb {

template <unsigned ndim, typename T>
class AcousticPulse : public AnalyticalF<ndim,T,T> {
protected:
  T rho0;
  T amplitude;
  T alpha;
  Vector<T,ndim> x0;
public:
  AcousticPulse(T rho0, T amplitude, T alpha, Vector<T,ndim> x0=Vector<T,ndim>(0.) ) : AnalyticalF<ndim,T,T>(1),
    rho0(rho0), amplitude(amplitude), alpha(alpha), x0(x0) {};

  bool operator()(T output[], const T input[]) override
  {
    T distance=0;
    for ( unsigned d=0; d<ndim; d++ ) {
      T distance_i = input[d] - x0[d];
      distance += distance_i*distance_i;  // actually, distance would be square root, but would be squared in exponent
    }
    output[0] = rho0+amplitude*util::exp(-alpha*distance);
    // if ( output[0] > rho0 ) std::cout << "input=[" << input[0] << "," << input[1] << "," << input[2] << "], distance=" << distance << ", output=" << output[0] << std::endl;
    return true;
  };
};

}