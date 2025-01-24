
#include "olb3D.h"
#include "olb3D.hh" // use only generic version!
#include <cmath>
#include <vector>

namespace olb {

template <typename T, typename NSDESCRIPTOR, typename TDESCRIPTOR>
class ThermalCreep3d : public AnalyticalF3D<T, T> {
protected:
  T                                                         _Kn;
  T                                                         _sigma;
  T                                                         _Dp;
  T                                                         _Rs;
  T                                                         _DT;
  T                                                         _L;
  const ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR>& _converter;
  std::vector<T>                                            _axisPoint;
  std::vector<T>                                            _axisDirection;

public:
  /**
     * \param Kn
     *        Knudsen number (default is 0)
     * \param sigma
     *        tangential momentum accommodation coefficient default is 1
     * \param Dp
     *        pressure difference
     *        ly of the pipe, equal to h/2
     * \param Rs
     *        Specific gas constant
     * \param DT
     *        Temperature difference
     * \param L
     *        Length of the pipe
     * \param axisPoint
     *        The origin point for the coordinate system
     * \param axisDirection
     *        The orientation of the axes
     * \param isforced
     *
    **/
  ThermalCreep3d(
      std::vector<T> axisPoint, std::vector<T> axisDirection, T Dp, T Rs, T DT,
      T Kn, T L,
      const ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR>& converter,
      T                                                         sigma = T(1))
      : AnalyticalF3D<T, T>(3)
      , _converter(converter)
  {
    this->getName() = "velocity_profile";
    _Rs             = Rs;
    _Dp             = Dp;
    _DT             = DT;
    _L              = L;
    _Kn             = Kn;
    _sigma          = sigma;
    _axisPoint      = axisPoint;
    _axisDirection  = util::normalize(axisDirection);
  }
  ThermalCreep3d(
      T Dp, T Rs, T DT, T Kn, T L,
      const ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR>& converter,
      T                                                         sigma = T(1))
      : ThermalCreep3d({0, 0, 0}, {1, 0, 0}, Dp, Rs, DT, Kn, L, converter,
                       sigma)
  {}

  bool operator()(T output[], const T x[]) override
  {
    T Ks = (2 - _sigma) / _sigma;
    T alpha =
        0.; // 1 for second degree in Knudsen number, 0 for only first degree
    T radius      = _converter.getCharPhysLength() / 2;
    T mu          = _converter.getPhysViscosity() * _converter.getPhysDensity();
    T maxVelocity = (_Dp / _L) * (radius * radius / (4 * mu));
    T r2          = (x[1] - _axisPoint[1]) * (x[1] - _axisPoint[1]) +
           (x[2] - _axisPoint[2]) * (x[2] - _axisPoint[2]);
    T creep =
        (3. / 4.) * (_DT / _L) * mu * (_Rs / _converter.getCharPhysPressure());

    output[0] =
        ((-1) * maxVelocity * _axisDirection[0] *
         (1. - r2 / (radius * radius) + 4 * (_Kn - alpha * _Kn * _Kn) * Ks)) +
        creep;
    output[1] = (-1) * maxVelocity * _axisDirection[1] *
                    (1. -
                     ((x[0] - _axisPoint[0]) * (x[0] - _axisPoint[0]) +
                      (x[2] - _axisPoint[2]) * (x[2] - _axisPoint[2])) /
                         (radius * radius) +
                     4 * (_Kn - alpha * _Kn * _Kn) * Ks) +
                creep;

    output[2] = (-1) * maxVelocity * _axisDirection[2] *
                    (1. -
                     ((x[0] - _axisPoint[0]) * (x[0] - _axisPoint[0]) +
                      (x[1] - _axisPoint[1]) * (x[1] - _axisPoint[1])) /
                         (radius * radius) +
                     4 * (_Kn - alpha * _Kn * _Kn) * Ks) +
                creep;

    if (_Kn >= 10) {
      output[0] = maxVelocity * _axisDirection[0];
      output[1] = maxVelocity * _axisDirection[1];
      output[2] = maxVelocity * _axisDirection[2];
    }
    return true;
  }
};

} // namespace olb