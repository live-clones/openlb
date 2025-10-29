/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Max Gaedtke, Albert Mink
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

/** \file
 * Unit conversion handling -- header file.
 */

#ifndef UNITCONVERTER_H
#define UNITCONVERTER_H


#include "utilities/omath.h"
#include "io/ostreamManager.h"
#include "core/util.h"
#include "io/xmlReader.h"

#include "descriptor/fields.h"

// known design issues
//    1. How can we prevent abuse of constructur by mixing up parameters?
//    2. physical problems may have different names for viscosity, e.g. diffusity,  temperature conductivity
//    4. Feedback about stability or comment the chosen discretization
//    5. Explain why Desctiptor as template
//    6. Is it worth to introduce invConversionDensity to avoid division


namespace olb {

namespace fields {

namespace converter {

struct PHYS_VELOCITY : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value >= 0;
  }
};
struct PHYS_FORCE : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value >= 0;
  }
};
struct PHYS_LENGTH : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value > 0;
  }
 };
struct PHYS_DELTA_X : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value > 0;
  }
};
struct PHYS_DELTA_T : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value > 0;
  }
};
struct PHYS_PRESSURE : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value >= 0;
  }
};

struct LATTICE_TIME : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value >= 0;
  }
};

struct LATTICE_VISCOSITY : public descriptors::FIELD_BASE<1> { };

}

}


struct UnitConverterBase {
  virtual ~UnitConverterBase() = default;

  virtual void print() const = 0;
  virtual void print(std::ostream& fout) const = 0;
  virtual void write(std::string const& fileName = "unitConverter") const = 0;
};


/** Conversion between physical and lattice units, as well as discretization.
* Be aware of the nomenclature:
* We distingish between physical (dimensioned) and lattice (dimensionless) values.
* A specific conversion factor maps the two different scopes,
* e.g. __physLength = conversionLength * latticeLength__
*
* For pressure and temperature we first shift the physical values by a characteristic value to asure a lattice pressure and lattice temperature between 0 and 1, e.g. __physPressure - charPhysPressure = conversionPressure * latticePressure__
*
*  \param latticeRelaxationTime   relaxation time, have to be greater than 0.5!
*  - - -
*  \param physViscosity         physical kinematic viscosity in __m^2 / s__
*  \param physDensity           physical density in __kg / m^3__
*  - - -
*  \param conversionLength      conversion factor for length __m__
*  \param conversionTime        conversion factor for time __s__
*  \param conversionMass        conversion factor for mass __kg__
*  - - -
*  \param conversionVelocity    conversion velocity __m / s__
*  \param conversionViscosity   conversion kinematic viscosity __m^2 / s__
*  \param conversionDensity     conversion density __kg / m^3__
*  \param conversionForce       conversion force __kg m / s^2__
*  \param conversionPressure    conversion pressure __kg / m s^2__
*  - - -
*  \param resolution            number of grid points per charPhysLength
*  - - -
*  \param charLatticeVelocity
*/
template <typename T, typename DESCRIPTOR>
class UnitConverter : public UnitConverterBase {
public:
  /** Documentation of constructor:
    *  \param physDeltaX              spacing between two lattice cells in __m__
    *  \param physDeltaT              time step in __s__
    *  \param charPhysLength          reference/characteristic length of simulation geometry in __m__
    *  \param charPhysVelocity        maximal or highest expected velocity during simulation in __m / s__
    *  \param physViscosity           physical kinematic viscosity in __m^2 / s__
    *  \param physDensity             physical density in __kg / m^3__
    *  \param charPhysPressure        reference/characteristic physical pressure in Pa = kg / m s^2
    */
  UnitConverter( T physDeltaX, T physDeltaT, T charPhysLength, T charPhysVelocity,
                           T physViscosity, T physDensity, T charPhysPressure = 0 )
    : _conversionLength(physDeltaX),
      _conversionTime(physDeltaT),
      _conversionVelocity(_conversionLength / _conversionTime),
      _conversionDensity(physDensity),
      _conversionMass( _conversionDensity * util::pow(_conversionLength, 3) ),
      _conversionViscosity(_conversionLength * _conversionLength / _conversionTime),
      _conversionForce( _conversionMass * _conversionLength / (_conversionTime * _conversionTime) ),
      _conversionTorque( _conversionMass * _conversionLength * _conversionLength / (_conversionTime * _conversionTime) ),
      _conversionPressure( _conversionForce / util::pow(_conversionLength, 2) ),
      _charPhysLength(charPhysLength),
      _charPhysVelocity(charPhysVelocity),
      _physViscosity(physViscosity),
      _physDensity(physDensity),
      _charPhysPressure(charPhysPressure),
      _resolution((size_t)(_charPhysLength / _conversionLength + 0.5)),
      _latticeRelaxationTime( (_physViscosity / _conversionViscosity * descriptors::invCs2<T,DESCRIPTOR>()) + 0.5 ),
      _charLatticeVelocity( _charPhysVelocity / _conversionVelocity ),
      clout(std::cout,"UnitConverter")
  {
  }

  virtual ~UnitConverter() = default;

  /// return resolution
  int getResolution(  ) const
  {
    return _resolution;
  }
  /// return relaxation time in lattice units
  T getLatticeRelaxationTime(  ) const
  {
    return _latticeRelaxationTime;
  }
  /// return relaxation frequency in lattice units
  T getLatticeRelaxationFrequency(  ) const
  {
    return 1./_latticeRelaxationTime;
  }
  /// return relaxation frequency in lattice units computed from given physical diffusivity in __m^2 / s__
  template <typename DESCRIPTOR_>
  T getLatticeRelaxationFrequencyFromDiffusivity( const T physDiffusivity ) const
  {
    return 1.0 / ( physDiffusivity / _conversionViscosity * descriptors::invCs2<T,DESCRIPTOR_>() + 0.5 );
  }
  /// return characteristic length in physical units
  T getCharPhysLength(  ) const
  {
    return _charPhysLength;
  }
  /// return characteristic velocity in physical units
  T getCharPhysVelocity(  ) const
  {
    return _charPhysVelocity;
  }
  /// return characteristic velocity in lattice units
  T getCharLatticeVelocity(  ) const
  {
    return _charLatticeVelocity;
  }
  /// return characteristic CFL number
  T getCharCFLnumber(  ) const
  {
    return _charLatticeVelocity;
  }
  /// return viscosity in physical units
  T getPhysViscosity(  ) const
  {
    return _physViscosity;
  }
  /// return density in physical units
  T getPhysDensity(  ) const
  {
    return _physDensity;
  }
  /// return characteristic pressure in physical units
  T getCharPhysPressure(  ) const
  {
    return _charPhysPressure;
  }
  /// return Reynolds number
  T getReynoldsNumber(  ) const
  {
    return _charPhysVelocity * _charPhysLength / _physViscosity;
  }
  /// return Mach number
  T getMachNumber(  ) const
  {
    return getCharLatticeVelocity() * util::sqrt(descriptors::invCs2<T,DESCRIPTOR>());
  }
  /// return Knudsen number
  virtual T getKnudsenNumber(  ) const
  {
    // This calculates the lattice Knudsen number.
    // See e.g. (7.22) in "The Lattice Boltzmann Method: Principles and Practice" [kruger2017lattice].
    return getMachNumber() / getReynoldsNumber();
  }
  /// conversion from lattice to  physical length
  T getPhysLength( int latticeLength ) const
  {
    return _conversionLength * latticeLength;
  }
  /// conversion from physical to lattice length, returns number of voxels for given physical length
  int getLatticeLength( T physLength ) const
  {
    return int( physLength / _conversionLength + 0.5 );
  }
  /// access (read-only) to private member variable
  T getConversionFactorLength() const
  {
    return _conversionLength;
  }
  /// returns grid spacing (voxel length) in __m__
  T getPhysDeltaX() const
  {
    return _conversionLength;
  }

  /// conversion from lattice to  physical time
  T getPhysTime( size_t latticeTime ) const
  {
    return _conversionTime * latticeTime;
  }
  /// conversion from physical to lattice time
  size_t getLatticeTime( T physTime ) const
  {
    return size_t(physTime / _conversionTime + 0.5);
  }
  /// access (read-only) to private member variable
  T getConversionFactorTime() const
  {
    return _conversionTime;
  }
  /// returns time spacing (timestep length) in __s__
  T getPhysDeltaT() const
  {
    return _conversionTime;
  }

  /// conversion from lattice to  physical velocity
  T getPhysVelocity( T latticeVelocity ) const
  {
    return _conversionVelocity * latticeVelocity;
  }
  /// conversion from physical to lattice velocity
  T getLatticeVelocity( T physVelocity ) const
  {
    return physVelocity / _conversionVelocity;
  }
  /// conversion from physical to lattice velocity
  template <unsigned D>
  Vector<T,D> getLatticeVelocity(Vector<T,D> physU) const
  {
    return Vector<T,D>([&](std::size_t iD) -> T {
      return this->getLatticeVelocity(physU[iD]);
    });
  }
  /// access (read-only) to private member variable
  T getConversionFactorVelocity() const
  {
    return _conversionVelocity;
  }

  /// conversion from lattice to  physical density
  T getPhysDensity( T latticeDensity ) const
  {
    return _conversionDensity * latticeDensity;
  }
  /// conversion from physical to lattice density
  T getLatticeDensity( T physDensity ) const
  {
    return physDensity / _conversionDensity;
  }
  T getLatticeDensityFromPhysPressure( T physPressure ) const
  {
    return util::densityFromPressure<T,DESCRIPTOR>( getLatticePressure( physPressure ) );
  }
  /// access (read-only) to private member variable
  T getConversionFactorDensity() const
  {
    return _conversionDensity;
  }

  /// conversion from lattice to  physical mass
  T getPhysMass( T latticeMass ) const
  {
    return _conversionMass * latticeMass;
  }
  /// conversion from physical to lattice mass
  T getLatticeMass( T physMass ) const
  {
    return physMass / _conversionMass;
  }
  /// access (read-only) to private member variable
  T getConversionFactorMass() const
  {
    return _conversionMass;
  }

  /// conversion from lattice to  physical viscosity
  T getPhysViscosity( T latticeViscosity ) const
  {
    return _conversionViscosity * latticeViscosity;
  }
  /// conversion from physical to lattice viscosity
  T getLatticeViscosity(  ) const
  {
    return _physViscosity / _conversionViscosity;
  }
  /// access (read-only) to private member variable
  T getConversionFactorViscosity() const
  {
    return _conversionViscosity;
  }

  /// conversion from lattice to  physical force
  T getPhysForce( T latticeForce ) const
  {
    return _conversionForce * latticeForce;
  }
  /// conversion from lattice to  physical force vector
  template <unsigned D>
  Vector<T,D> getPhysForce(Vector<T,D> latticeForce) const
  {
    return Vector<T,D>([&](std::size_t iD) -> T {
      return this->getPhysForce(latticeForce[iD]);
    });
  }
  /// conversion from physical to lattice force
  T getLatticeForce( T physForce ) const
  {
    return physForce / _conversionForce;
  }
  /// access (read-only) to private member variable
  T getConversionFactorForce() const
  {
    return _conversionForce;
  }

  /// conversion from lattice to  physical torque
  T getPhysTorque ( T latticeTorque ) const
  {
    return _conversionTorque * latticeTorque;
  }
  /// conversion from lattice to  physical force vector
  template <unsigned D>
  Vector<T,D> getPhysTorque(Vector<T,D> latticeTorque) const
  {
    return Vector<T,D>([&](std::size_t iD) -> T {
      return this->getPhysTorque(latticeTorque[iD]);
    });
  }
  /// conversion from physical to lattice torque
  T getLatticeTorque( T physTorque ) const
  {
    return physTorque / _conversionTorque;
  }
  /// access (read-only) to private member variable
  T getConversionFactorTorque() const
  {
    return _conversionTorque;
  }

  /// conversion from lattice to  physical pressure
  T getPhysPressure( T latticePressure ) const
  {
    return _conversionPressure * latticePressure + _charPhysPressure;
  }
  /// conversion from physical to lattice pressure
  T getLatticePressure( T physPressure ) const
  {
    return ( physPressure - _charPhysPressure ) / _conversionPressure;
  }
  /// access (read-only) to private member variable
  T getConversionFactorPressure() const
  {
    return _conversionPressure;
  }
  /// nice terminal output for conversion factors, characteristical and physical data
  virtual void print() const;
  void print(std::ostream& fout) const;

  virtual void write(std::string const& fileName = "unitConverter") const;


  // from thermalUnitConverter
  /// return thermal relaxation time in lattice units
  virtual T getLatticeThermalRelaxationTime(  ) const {
    throw std::logic_error("Undefined");
  };
  /// return thermal relaxation frequency in lattice units
  virtual T getLatticeThermalRelaxationFrequency() const {
    throw std::logic_error("Undefined");
  };
  /// return characteristic low temperature in physical units
  virtual T getCharPhysLowTemperature(  ) const {
    throw std::logic_error("Undefined");
  };
  /// return characteristic high temperature in physical units
  virtual T getCharPhysHighTemperature(  ) const {
    throw std::logic_error("Undefined");
  };
  /// return characteristic temperature difference in physical units
  virtual T getCharPhysTemperatureDifference() const {
    throw std::logic_error("Undefined");
  };
  /// return thermal expansion coefficient in physical units
  virtual T getPhysThermalExpansionCoefficient() const {
    throw std::logic_error("Undefined");
  };
  /// return thermal diffusivity in physical units
  virtual T getPhysThermalDiffusivity(  ) const {
    throw std::logic_error("Undefined");
  };
  /// return specific heat capacity in physical units
  virtual T getPhysSpecificHeatCapacity(  ) const {
    throw std::logic_error("Undefined");
  };
  /// return thermal conductivity in physical units
  virtual T getThermalConductivity(  ) const {
    throw std::logic_error("Undefined");
  };
  /// conversion from lattice to physical temperature
  virtual T getPhysTemperature( T latticeTemperature ) const {
    throw std::logic_error("Undefined");
  };
  /// conversion from physical to lattice temperature
  virtual T getLatticeTemperature( T physTemperature ) const {
    throw std::logic_error("Undefined");
  };
  /// conversion from lattice to physical thermal diffusivity
  virtual T getPhysThermalDiffusivity( T latticeThermalDiffusivity ) const {
    throw std::logic_error("Undefined");
  };
  /// conversion from physical to lattice thermal diffusivity
  virtual T getLatticeThermalDiffusivity( T physThermalDiffusivity ) const {
    throw std::logic_error("Undefined");
  };
  /// access (read-only) to private member variable
  virtual T getConversionFactorThermalDiffusivity() const {
    throw std::logic_error("Undefined");
  };
  /// conversion from lattice to physical specific heat capacity
  virtual T getPhysSpecificHeatCapacity( T latticeSpecificHeatCapacity ) const {
    throw std::logic_error("Undefined");
  };
  /// conversion from physical to lattice specific heat capacity
  virtual T getLatticeSpecificHeatCapacity( T physSpecificHeatCapacity ) const {
    throw std::logic_error("Undefined");
  };
  /// access (read-only) to private member variable
  virtual T getConversionFactorSpecificHeatCapacity() const {
    throw std::logic_error("Undefined");
  };
  /// conversion from lattice to physical thermal  conductivity
  virtual T getPhysThermalConductivity( T latticeThermalConductivity ) const {
    throw std::logic_error("Undefined");
  };
  /// conversion from physical to lattice thermal  conductivity
  virtual T getLatticeThermalConductivity( T physThermalConductivity ) const {
    throw std::logic_error("Undefined");
  };
  /// access (read-only) to private member variable
  virtual T getConversionFactorThermalConductivity() const {
    throw std::logic_error("Undefined");
  };
  /// conversion from lattice to physical heat flux
  virtual T getPhysHeatFlux( T latticeHeatFlux ) const {
    throw std::logic_error("Undefined");
  };
  /// conversion from physical to lattice heat flux
  virtual T getLatticeHeatFlux( T physHeatFlux ) const {
    throw std::logic_error("Undefined");
  };
  /// access (read-only) to private member variable
  virtual T getConversionFactorHeatFlux() const {
    throw std::logic_error("Undefined");
  };
  virtual T getPrandtlNumber() const {
    throw std::logic_error("Undefined");
  };
  virtual T getRayleighNumber() const {
    throw std::logic_error("Undefined");
  };

  // from multiPhaseUnitConverter
  /// return characteristic temperature in physical units
  virtual T getCharPhysTemperature() const {
    throw std::logic_error("Undefined");
  };
  /// return equation of state parameter a in physical units
  virtual T getPhysEoSa() const {
    throw std::logic_error("Undefined");
  };
  /// return characteristic temperature in physical units
  virtual T getPhysTemperature(  ) const  {
    throw std::logic_error("Undefined");
  };
  /// return equation of state parameter b in physical units
  virtual T getPhysEoSb() const {
    throw std::logic_error("Undefined");
  };
  /// return molar mass in physical units
  virtual T getPhysMolarMass() const {
    throw std::logic_error("Undefined");
  };
  /// return surface tension in physical units
  virtual T getPhysSurfaceTension() const {
    throw std::logic_error("Undefined");
  };
  /// access (read-only) to private member variable
  virtual T getConversionFactorEoSa() const {
    throw std::logic_error("Undefined");
  };
  /// access (read-only) to private member variable
  virtual T getConversionFactorEoSb() const {
    throw std::logic_error("Undefined");
  };
  /// access (read-only) to private member variable
  virtual T getConversionFactorMolarMass() const  {
    throw std::logic_error("Undefined");
  };
  /// access (read-only) to private member variable
  virtual T getConversionFactorGasConstant() const  {
    throw std::logic_error("Undefined");
  };
  /// access (read-only) to private member variable
  virtual T getConversionFactorTemperature() const  {
    throw std::logic_error("Undefined");
  };
  /// access (read-only) to private member variable
  virtual T getLatticeSurfaceTension() const  {
    throw std::logic_error("Undefined");
  };

  virtual T getConversionFactorSurfaceTension() const {
    throw std::logic_error("Undefined");
  };
  virtual T getConversionFactorChemicalPotential() const {
    throw std::logic_error("Undefined");
  };
  virtual T computeRelaxationTimefromPhysViscosity( T userViscosity ) const {
    throw std::logic_error("Undefined");
  };
  virtual T computeLatticeSurfaceTension( T userSurfaceTension ) const {
    throw std::logic_error("Undefined");
  };
  virtual T computeReynolds(  T userVelocity, T userLength, T userViscosity  ) const {
    throw std::logic_error("Undefined");
  };
  virtual T computeWeber( T userVelocity, T userLength, T userSurfaceTension ) const {
    throw std::logic_error("Undefined");
  };

  //from powerLawUnitConverter
  /// return consistency coefficient in physical units
  virtual T getPhysConsistencyCoeff(  ) const {
    throw std::logic_error("Undefined");
  };
  /// conversion from lattice to  physical consistency coefficient
  virtual T getPhysConsistencyCoeff( T latticeConsistencyCoeff ) const {
    throw std::logic_error("Undefined");
  };
  /// conversion from physical to lattice consistency coefficient
  virtual T getLatticeConsistencyCoeff(  ) const {
    throw std::logic_error("Undefined");
  };
  /// access (read-only) to private member variable
  virtual T getConversionFactorConsistencyCoeff(  ) const {
    throw std::logic_error("Undefined");
  };
  /// access (read-only) to private member variable
  virtual T getPowerLawIndex(  ) const {
    throw std::logic_error("Undefined");
  };



  // from adeUnitConverter
  virtual T getLatticeAdeRelaxationTime() const {
    throw std::logic_error("Undefined");
  };
  virtual T getLatticeAdeRelaxationFrequency() const {
    throw std::logic_error("Undefined");
  };
  virtual T getPhysDiffusivity() const {
    throw std::logic_error("Undefined");
  };
  virtual T getLatticeDiffusivity() const {
    throw std::logic_error("Undefined");
  };
  virtual T getConversionFactorDiffusivity() const {
    throw std::logic_error("Undefined");
  };
  virtual T getPecletNumber() const {
    throw std::logic_error("Undefined");
  };

  // from adsorptionUnitConverter
  virtual T getConversionFactorParticleDensity() const {
    throw std::logic_error("Undefined");
  };
  virtual T getPhysParticleConcentration(T c) const {
    throw std::logic_error("Undefined");
  };
  virtual T getPhysConcentration(T c) const {
    throw std::logic_error("Undefined");
  };
  virtual T getPhysLoading(T Cq) const {
    throw std::logic_error("Undefined");
  };
  virtual T getSchmidtNumber() const {
    throw std::logic_error("Undefined");
  };
  virtual T getFourierNumber() const {
    throw std::logic_error("Undefined");
  };

  // from radiativeUnitConverter
  virtual T getPhysAbsorption() const {
    throw std::logic_error("Undefined");
  };
  virtual T getPhysScattering() const {
    throw std::logic_error("Undefined");
  };
  virtual T getAnisotropyFactor() const {
    throw std::logic_error("Undefined");
  };
  virtual T getExtinction() const {
    throw std::logic_error("Undefined");
  };
  virtual T getScatteringAlbedo() const {
    throw std::logic_error("Undefined");
  };
  virtual T getPhysDiffusion() const {
    throw std::logic_error("Undefined");
  };
  virtual T getEffectiveAttenuation() const {
    throw std::logic_error("Undefined");
  };
  virtual T getLatticeAbsorption() const {
    throw std::logic_error("Undefined");
  };
  virtual T getLatticeScattering() const {
    throw std::logic_error("Undefined");
  };
  virtual T getLatticeDiffusion() const {
    throw std::logic_error("Undefined");
  };
  virtual T getRefractiveRelative() const {
    throw std::logic_error("Undefined");
  };

  // from LinElaUnitConverter
  virtual T getEpsilon() const {
    throw std::logic_error("Undefined");
  };
  virtual T getCharPhysDisplacement() const {
    throw std::logic_error("Undefined");
  };
  virtual T getLatticeShearModulus() const {
    throw std::logic_error("Undefined");
  };
  virtual T getPhysShearModulus() const {
    throw std::logic_error("Undefined");
  };
  virtual T getLatticeBulkModulus() const {
    throw std::logic_error("Undefined");
  };
  virtual T getPhysBulkModulus() const {
    throw std::logic_error("Undefined");
  };
  virtual T getLatticeLambda() const {
    throw std::logic_error("Undefined");
  };
  virtual T getPhysLambda() const {
    throw std::logic_error("Undefined");
  };
  virtual T getLatticeYoungsModulus() const {
    throw std::logic_error("Undefined");
  };
  virtual T getPhysYoungsModulus() const {
    throw std::logic_error("Undefined");
  };
  virtual T getDampingFactor() const {
    throw std::logic_error("Undefined");
  };
  virtual T getPoissonRatio() const {
    throw std::logic_error("Undefined");
  };

protected:
  // conversion factors
  const T _conversionLength;      // m
  const T _conversionTime;        // s
  const T _conversionVelocity;    // m / s
  const T _conversionDensity;     // kg / m^3
  const T _conversionMass;        // kg
  const T _conversionViscosity;   // m^2 / s
  const T _conversionForce;       // kg m / s^2
  const T _conversionTorque;      // kg m^2 / s^2
  const T _conversionPressure;    // kg / m s^2

  // physical units, e.g characteristic or reference values
  const T _charPhysLength;        // m
  const T _charPhysVelocity;      // m / s
  const T _physViscosity;         // m^2 / s
  const T _physDensity;           // kg / m^3
  const T _charPhysPressure;      // kg / m s^2

  // lattice units, discretization parameters
  const size_t _resolution;
  const T _latticeRelaxationTime;
  const T _charLatticeVelocity;   //

private:
  mutable OstreamManager clout;
};

template <typename T, typename DESCRIPTOR>
UnitConverter<T,DESCRIPTOR> convectivelyRefineUnitConverter(
  const UnitConverter<T,DESCRIPTOR>& converter,
  unsigned scale = 2)
{
  const T refinementFactor = T{1} / scale;
  return UnitConverter<T,DESCRIPTOR>(
    converter.getPhysDeltaX() * refinementFactor,
    converter.getPhysDeltaT() * refinementFactor,
    converter.getCharPhysLength(),
    converter.getCharPhysVelocity(),
    converter.getPhysViscosity(),
    converter.getPhysDensity(),
    converter.getCharPhysPressure()
  );
}

/// creator function with data given by a XML file
template <typename T, typename DESCRIPTOR>
UnitConverter<T, DESCRIPTOR>* createUnitConverter(XMLreader const& params);

template <typename T, typename DESCRIPTOR>
class UnitConverterFromResolutionAndRelaxationTime : public UnitConverter<T, DESCRIPTOR> {
public:
  UnitConverterFromResolutionAndRelaxationTime(
    size_t resolution,
    T latticeRelaxationTime,
    T charPhysLength,
    T charPhysVelocity,
    T physViscosity,
    T physDensity,
    T charPhysPressure = 0) : UnitConverter<T, DESCRIPTOR>(
        (charPhysLength/resolution),
        (latticeRelaxationTime - 0.5) / descriptors::invCs2<T,DESCRIPTOR>() * util::pow((charPhysLength/resolution),2) / physViscosity,
        charPhysLength,
        charPhysVelocity,
        physViscosity,
        physDensity,
        charPhysPressure)
  { }
};

template <typename T, typename DESCRIPTOR>
class UnitConverterFromResolutionAndLatticeVelocity : public UnitConverter<T, DESCRIPTOR> {
public:
  UnitConverterFromResolutionAndLatticeVelocity(
    size_t resolution,
    T charLatticeVelocity,
    T charPhysLength,
    T charPhysVelocity,
    T physViscosity,
    T physDensity,
    T charPhysPressure = 0) : UnitConverter<T, DESCRIPTOR>(
        (charPhysLength/resolution),
        (charLatticeVelocity / charPhysVelocity * charPhysLength / resolution),
        charPhysLength,
        charPhysVelocity,
        physViscosity,
        physDensity,
        charPhysPressure)
  { }
};

template <typename T, typename DESCRIPTOR>
class UnitConverterFromRelaxationTimeAndLatticeVelocity : public UnitConverter<T, DESCRIPTOR> {
public:
  UnitConverterFromRelaxationTimeAndLatticeVelocity(
    T latticeRelaxationTime,
    T charLatticeVelocity,
    T charPhysLength,
    T charPhysVelocity,
    T physViscosity,
    T physDensity,
    T charPhysPressure = 0) : UnitConverter<T, DESCRIPTOR>(
        (physViscosity * charLatticeVelocity / charPhysVelocity * descriptors::invCs2<T,DESCRIPTOR>() / (latticeRelaxationTime - 0.5)),
        (physViscosity * charLatticeVelocity * charLatticeVelocity / charPhysVelocity / charPhysVelocity * descriptors::invCs2<T,DESCRIPTOR>() / (latticeRelaxationTime - 0.5)),
        charPhysLength,
        charPhysVelocity,
        physViscosity,
        physDensity,
        charPhysPressure)
  {
  }
};

template<typename T, typename DESCRIPTOR, typename MEMBER, typename... ARGS>
std::unique_ptr<UnitConverter<T,DESCRIPTOR>> createUnitConverter(ARGS&&... args) {
  return std::make_unique<MEMBER>(std::forward<ARGS>(args)...);
}

}  // namespace olb

#include "thermalUnitConverter.h"
#include "multiPhaseUnitConverter.h"
#endif
