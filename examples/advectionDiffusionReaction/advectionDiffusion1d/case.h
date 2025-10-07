#include <olb.h>

using namespace olb;
//using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::names;

using MyCase = Case<
  AdvectionDiffusion, Lattice<double, descriptors::D2Q5<descriptors::VELOCITY>>
>;

namespace olb::parameters {

  struct RUNS : public descriptors::FIELD_BASE<1> {};
  struct N0 : public descriptors::FIELD_BASE<1> {};

}

template <typename T, typename TDESCRIPTOR>
class AdePhysTemp1D : public AnalyticalF2D<T,T> {

protected:
  MyCase myCase;
  T t;
public:
  AdePhysTemp1D(T time, MyCase& _myCase) : AnalyticalF2D<T,T>(1), myCase(_myCase),
    t(time) {};

  bool operator()(T output[], const T input[]) override
  {
    using T = MyCase::value_t;
    auto& lattice = myCase.getLattice(AdvectionDiffusion{});
    const auto& converter = lattice.getUnitConverter();
    const T mue = converter.getPhysDiffusivity();
    const T uMag = converter.getCharPhysVelocity();
    const T res = converter.getResolution();

    T x = input[0];
    T gf = res/(res+1.);

    // initial condition (pseudo 1D)
    output[0] = util::sin(M_PI*(x-uMag*t)*gf) * util::exp(-mue*util::pow(M_PI,2)*t*util::pow(gf,2));

    return true;
  };
};

void prepareLattice(MyCase& myCase)
{
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<AdvectionDiffusion>;
  using namespace olb::parameters;

  auto& ADlattice = myCase.getLattice(AdvectionDiffusion{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const T physDeltaX = parameters.get<PHYS_CHAR_LENGTH>() / parameters.get<RESOLUTION>();
  const T physDeltaT = physDeltaX * physDeltaX;
  const T physLength = parameters.get<PHYS_CHAR_LENGTH>();
  const T mue = parameters.get<PHYS_DIFFUSIVITY>();
  const T physVelocity = parameters.get<PECLET>() * mue / physLength;
  const T physDensity = parameters.get<PHYS_CHAR_DENSITY>();

  ADlattice.setUnitConverter<AdeUnitConverter<T,DESCRIPTOR>>(
    physDeltaX,    // physDeltaX
    physDeltaT,    // physDeltaT (diffusive scaling)
    physLength,    // charPhysLength
    physVelocity,  // charPhysVelocity from Peclet
    mue,           // physDiffusivity
    physDensity    // physDensity,
  );

  const auto& converter = ADlattice.getUnitConverter();
  converter.print();

  ADlattice.defineDynamics<AdvectionDiffusionBGKdynamics>(geometry.getMaterialIndicator({1}));
}

void setInitialValues(MyCase& myCase) {
  using T = myCase::value_t;
  auto& ADlattice = myCase.getLattice();
  auto& geometry = myCase.getGeometry();
  const auto& converter = ADlattice.getUnitConverter();

  //TODO make 1D?
  AnalyticalConst1D<T,T> u0( converter.getCharLatticeVelocity());
  AdePhysTemp1D<T> Tinit( 0.0, converter );

  auto bulkIndicator = geometry.getMaterialIndicator({0,1});

  ADlattice.defineField<descriptors::VELOCITY>( bulkIndicator, u0 );
  ADlattice.defineRho( bulkIndicator, Tinit );
  ADlattice.iniEquilibrium( bulkIndicator, Tinit, u0 );

  ADlattice.setParameter<descriptors::OMEGA>( converter.getLatticeAdeRelaxationFrequency() );

  /// Make the lattice ready for simulation
  ADlattice.initialize();
}
