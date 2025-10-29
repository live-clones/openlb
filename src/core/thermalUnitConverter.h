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

#ifndef THERMALUNITCONVERTER_H
#define THERMALUNITCONVERTER_H


#include <math.h>
#include "io/ostreamManager.h"
#include "core/util.h"
#include "io/xmlReader.h"
#include "core/unitConverter.h"

// All OpenLB code is contained in this namespace.
namespace olb {



/** Conversion between physical and lattice units, as well as discretization specialized for thermal applications with boussinesq approximation.
* Be aware of the nomenclature:
* We distingish between physical (dimensioned) and lattice (dimensionless) values.
* A specific conversion factor maps the two different scopes,
* e.g. __physLength = conversionLength * latticeLength__
*
* For pressure and temperature we first shift the physical values by a characteristic value to asure a lattice pressure and between 0 and 1, e.g. __physPressure - charPhysPressure = conversionPressure * latticePressure__. For the temperature we set lattice values between 0.5 and 1.5 by __latticeTemperature = (physTemperature - charPhysLowTemperature) / conversionTemperature + 0.5 with conversionTemperature = charPhysHighTemperature - charPhysLowTemperature = charPhysTemperatureDifference
*
* TODO: Extend documentation for ThermalUnitConverter
*/
template <typename T, typename DESCRIPTOR, typename ThermalLattice>
class ThermalUnitConverter : public UnitConverter<T, DESCRIPTOR> {
public:
  /** Documentation of constructor:
    * TODO: Extend constructur documentation
    */
  ThermalUnitConverter(
    T physDeltaX,
    T physDeltaT,
    T charPhysLength,
    T charPhysVelocity,
    T physViscosity,
    T physDensity,
    T physThermalConductivity,
    T physSpecificHeatCapacity,
    T physThermalExpansionCoefficient,
    T charPhysLowTemperature,
    T charPhysHighTemperature,
    T charPhysPressure = 0 )
    : UnitConverter<T, DESCRIPTOR>(
        physDeltaX, physDeltaT, charPhysLength, charPhysVelocity,
        physViscosity, physDensity, charPhysPressure),
      _conversionTemperature(charPhysHighTemperature - charPhysLowTemperature),
      _conversionThermalDiffusivity(this->_conversionViscosity),
      _conversionSpecificHeatCapacity(this->_conversionVelocity * this->_conversionVelocity / _conversionTemperature),
      _conversionThermalConductivity(this->_conversionForce / this->_conversionTime / _conversionTemperature),
      _conversionHeatFlux(this->_conversionMass / util::pow(this->_conversionTime, 3)),
      _charPhysLowTemperature(charPhysLowTemperature),
      _charPhysHighTemperature(charPhysHighTemperature),
      _charPhysTemperatureDifference(charPhysHighTemperature - charPhysLowTemperature),
      _physThermalExpansionCoefficient(physThermalExpansionCoefficient),
      _physThermalDiffusivity(physThermalConductivity / (physDensity * physSpecificHeatCapacity)),
      _physSpecificHeatCapacity(physSpecificHeatCapacity),
      _physThermalConductivity(physThermalConductivity),
      _latticeThermalRelaxationTime( (_physThermalDiffusivity / _conversionThermalDiffusivity * descriptors::invCs2<T,ThermalLattice>()) + 0.5 ),
      clout(std::cout,"ThermalUnitConv")
  {
  };

  /// return thermal relaxation time in lattice units
  T getLatticeThermalRelaxationTime(  ) const override
  {
    return _latticeThermalRelaxationTime;
  };
  /// return thermal relaxation frequency in lattice units
  T getLatticeThermalRelaxationFrequency(  ) const override
  {
    return 1.0 / _latticeThermalRelaxationTime;
  };

  /// return characteristic low temperature in physical units
  T getCharPhysLowTemperature(  ) const override
  {
    return _charPhysLowTemperature;
  };
  /// return characteristic high temperature in physical units
  T getCharPhysHighTemperature(  ) const override
  {
    return _charPhysHighTemperature;
  };
  /// return characteristic temperature difference in physical units
  T getCharPhysTemperatureDifference(  ) const override
  {
    return _charPhysTemperatureDifference;
  };
  /// return thermal expansion coefficient in physical units
  T getPhysThermalExpansionCoefficient(  ) const override
  {
    return _physThermalExpansionCoefficient;
  };
  /// return thermal diffusivity in physical units
  T getPhysThermalDiffusivity(  ) const override
  {
    return _physThermalDiffusivity;
  };
  /// return specific heat capacity in physical units
  T getPhysSpecificHeatCapacity(  ) const override
  {
    return _physSpecificHeatCapacity;
  };
  /// return thermal conductivity in physical units
  T getThermalConductivity(  ) const override
  {
    return _physThermalConductivity;
  };

  /// conversion from lattice to physical temperature
  T getPhysTemperature( T latticeTemperature ) const override
  {
    return _conversionTemperature * (latticeTemperature - 0.5) + _charPhysLowTemperature;
  };
  /// conversion from physical to lattice temperature
  T getLatticeTemperature( T physTemperature ) const override
  {
    return (physTemperature - _charPhysLowTemperature) / _conversionTemperature + 0.5;
  };
  /// access (read-only) to private member variable
  T getConversionFactorTemperature() const override
  {
    return _conversionTemperature;
  };

  /// conversion from lattice to physical thermal diffusivity
  T getPhysThermalDiffusivity( T latticeThermalDiffusivity ) const override
  {
    return _conversionThermalDiffusivity * latticeThermalDiffusivity;
  };
  /// conversion from physical to lattice thermal diffusivity
  T getLatticeThermalDiffusivity( T physThermalDiffusivity ) const override
  {
    return physThermalDiffusivity / _conversionThermalDiffusivity;
  };
  /// access (read-only) to private member variable
  T getConversionFactorThermalDiffusivity() const override
  {
    return _conversionThermalDiffusivity;
  };


  /// conversion from lattice to physical specific heat capacity
  T getPhysSpecificHeatCapacity( T latticeSpecificHeatCapacity ) const override
  {
    return _conversionSpecificHeatCapacity * latticeSpecificHeatCapacity;
  };
  /// conversion from physical to lattice specific heat capacity
  T getLatticeSpecificHeatCapacity( T physSpecificHeatCapacity ) const override
  {
    return physSpecificHeatCapacity / _conversionSpecificHeatCapacity;
  };
  /// access (read-only) to private member variable
  T getConversionFactorSpecificHeatCapacity() const override
  {
    return _conversionSpecificHeatCapacity;
  };

  /// conversion from lattice to physical thermal  conductivity
  T getPhysThermalConductivity( T latticeThermalConductivity ) const override
  {
    return _conversionThermalConductivity * latticeThermalConductivity;
  };
  /// conversion from physical to lattice thermal  conductivity
  T getLatticeThermalConductivity( T physThermalConductivity ) const override
  {
    return physThermalConductivity / _conversionThermalConductivity;
  };
  /// access (read-only) to private member variable
  T getConversionFactorThermalConductivity() const override
  {
    return _conversionThermalConductivity;
  };

  /// conversion from lattice to physical heat flux
  T getPhysHeatFlux( T latticeHeatFlux ) const override
  {
    return _conversionHeatFlux * latticeHeatFlux;
  };
  /// conversion from physical to lattice heat flux
  T getLatticeHeatFlux( T physHeatFlux ) const override
  {
    return physHeatFlux / _conversionHeatFlux;
  };
  /// access (read-only) to private member variable
  T getConversionFactorHeatFlux() const override
  {
    return _conversionHeatFlux;
  };
  T getPrandtlNumber() const override
  {
    return this->_physViscosity/_physThermalDiffusivity;
  };
  T getRayleighNumber() const override
  {
    return 9.81 * _physThermalExpansionCoefficient/this->_physViscosity/_physThermalDiffusivity * (_charPhysHighTemperature - _charPhysLowTemperature) * util::pow(this->_charPhysLength,3);
  };
/// nice terminal output for conversion factors, characteristical and physical data
  void print() const override;



protected:
  // conversion factors
  const T _conversionTemperature; // K
  const T _conversionThermalDiffusivity; // m^2 / s
  const T _conversionSpecificHeatCapacity; // J / kg K = m^2 / s^2 K
  const T _conversionThermalConductivity; // W / m K = kg m / s^3 K
  const T _conversionHeatFlux; // W / m^2 = kg / s^3

  // physical units, e.g characteristic or reference values
  const T _charPhysLowTemperature; // K
  const T _charPhysHighTemperature; // K
  const T _charPhysTemperatureDifference; // K
  const T _physThermalExpansionCoefficient; // 1 / K
  const T _physThermalDiffusivity; // m^2 / s
  const T _physSpecificHeatCapacity; // J / kg K = m^2 / s^2 K
  const T _physThermalConductivity; // W / m K = kg m / s^3 K

  // lattice units, discretization parameters
  const T _latticeThermalRelaxationTime; // -

private:
  mutable OstreamManager clout;
};

}  // namespace olb

#endif
