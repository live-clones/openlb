namespace olb {

template <unsigned ndim, typename T>
class AcousticPulse : public AnalyticalF<ndim,T,T> {
protected:
  T rho0;
  T amplitude;
  T alpha;
  T scaling=2500;  // shock is too large otherwise -> scaling it to 1/50 extend
public:
  AcousticPulse(T rho0, T amplitude, T alpha ) : AnalyticalF<ndim,T,T>(1),
    rho0(rho0), amplitude(amplitude), alpha(alpha) {};

  bool operator()(T output[], const T input[]) override
  {
    T distance=0;
    for ( unsigned d=0; d<ndim; d++ ) {
      distance += input[d]*input[d]*scaling;
    }
    output[0] = rho0+amplitude*util::exp(-alpha*distance);
    return true;
  };
};

}