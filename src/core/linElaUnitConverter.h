/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Louis Kronberg, Stephan Simonis
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
 * Unit conversion handling for Advection-Diffusion Problmes -- header file.
 */
#ifndef LinEla_UNITCONVERTER_H
#define LinEla_UNITCONVERTER_H

#include <math.h>
#include <fstream>
#include <iostream>
#include <unistd.h>

#include "io/ostreamManager.h"
#include "io/fileName.h"

#include "core/util.h"
#include "core/unitConverter.h"
#include "core/singleton.h"


namespace olb
{

  template <typename T, typename DESCRIPTOR>
  class LinElaUnitConverter : public UnitConverter<T, DESCRIPTOR>
  {
  public:

    constexpr LinElaUnitConverter(
      T physDeltaX,
      T physDeltaT,
      T charPhysLength,
      T charPhysDisplacement,
      T youngsModulus,
      T poissonRatio,
      T dampingFactor
    )

      :UnitConverter<T, DESCRIPTOR>(
        physDeltaX,
        physDeltaT,
        charPhysLength,
        charPhysDisplacement,
        youngsModulus,
        poissonRatio,
        dampingFactor
      ),

      _physDeltaX           ( physDeltaX ),
      _physDeltaT           ( physDeltaT ),
      _charPhysLength       ( charPhysLength ),
      _charPhysDisplacement     ( charPhysDisplacement ),
      _epsilon              ( physDeltaX / charPhysLength ),
      _charPhysTime         ( physDeltaT / (physDeltaX * physDeltaX / (charPhysLength * charPhysLength)) ),
      _youngsModulus        ( youngsModulus  * physDeltaT / (physDeltaX * physDeltaX * dampingFactor) ),
      _bulkModulus          ( (youngsModulus / (2. * (1. - poissonRatio)))  * physDeltaT / (physDeltaX * physDeltaX * dampingFactor) ), // K
      _shearModulus         ( (youngsModulus / (2. * (1. + poissonRatio))) * physDeltaT / (physDeltaX * physDeltaX * dampingFactor) ), // dynamic viscosity in fluid dynamics
      _lambda               ( (youngsModulus * poissonRatio / (1. - poissonRatio * poissonRatio)) * physDeltaT / (physDeltaX * physDeltaX * dampingFactor) ),
      _poissonRatio         ( poissonRatio ),
      _dampingFactor        ( dampingFactor ), // Kappa
      clout(std::cout, "LinElaUnitConv")
    {
    };

    constexpr int getResolution( ) const
    {
      return 1 / _physDeltaX;
    };

    constexpr T getEpsilon( ) const
    {
      return _physDeltaX / _charPhysLength;
    };

    constexpr T getConversionFactorLength( ) const
    {
      return _physDeltaX;
    }

    constexpr T getPhysDeltaX( ) const
    {
      return _physDeltaX;
    }

    constexpr T getConversionFactorTime( ) const
    {
      return _physDeltaT;
    }

    constexpr T getPhysDeltaT( ) const
    {
      return _physDeltaT;
    }

    constexpr T getCharPhysDisplacement( ) const
    {
      return _charPhysDisplacement;
    }

    constexpr T getCharPhysTime( ) const
    {
      return _charPhysTime;
    }

    constexpr T getCharPhysLength( ) const
    {
      return _charPhysLength;
    }

      constexpr T getLatticeShearModulus( ) const
    {
      return _shearModulus;
    };

    constexpr T getPhysShearModulus( ) const
    {
      return _shearModulus * (_physDeltaX * _physDeltaX * _dampingFactor) / _physDeltaT;
    };

    constexpr T getLatticeBulkModulus( ) const
    {
      return _bulkModulus;
    };

    constexpr T getPhysBulkModulus( ) const
    {
      return _bulkModulus * (_physDeltaX * _physDeltaX * _dampingFactor) / _physDeltaT;
    };

    constexpr T getLatticeLambda( ) const
    {
      return _lambda;
    };

    constexpr T getPhysLambda( ) const
    {
      return _lambda * (_physDeltaX * _physDeltaX * _dampingFactor) / _physDeltaT;
    };

    constexpr T getLatticeYoungsModulus() const
    {
      return _youngsModulus;
    };

    constexpr T getPhysYoungsModulus() const
    {
      return _youngsModulus * (_physDeltaX * _physDeltaX * _dampingFactor) / _physDeltaT;
    };

    constexpr T getDampingFactor() const
    {
      return _dampingFactor;
    };

    constexpr T getPoissonRatio() const
    {
      return _poissonRatio;
    };

    void print() const override;

    void write(std::string const& fileName = "LinElaUnitConverter") const override;

  protected:
    // lattice units, discretization parameters
    const T _charPhysLength;
    const T _charPhysTime;
    const T _charPhysDisplacement;
    const T _physDeltaX;
    const T _physDeltaT;
    const T _epsilon;
    const T _youngsModulus;
    const T _bulkModulus;
    const T _shearModulus;
    const T _lambda;
    const T _poissonRatio;
    const T _dampingFactor;
    // conversionfactors

  private:
    mutable OstreamManager clout;
  };

  template <typename T, class DESCRIPTOR>
  void LinElaUnitConverter<T, DESCRIPTOR>::print() const

  {
    clout << "----------------- UnitConverter information -----------------" << std::endl;
    clout << "-- Parameters:" << std::endl;
    clout << "Resolution:                           N=              " << this->getResolution() << std::endl;
    clout << "Characteristical length(m):           charL=          " << this->getCharPhysLength() << std::endl;
    clout << "Characteristical time(s):             charT=          " << this->getCharPhysTime() << std::endl;
    clout << "Characteristical displacement (m/s):  charU=          " << this->getCharPhysDisplacement() << std::endl;
    clout << "Smallness Parameter:              epsilon=        " << this->getEpsilon() << std::endl;
    clout << std::endl;

    clout << std::endl;
    clout << "-- Conversion factors:" << std::endl;
    clout << "Voxel length(m):                  physDeltaX=     " << this->getConversionFactorLength() << std::endl;
    clout << "Time step(s):                     physDeltaT=     " << this->getConversionFactorTime() << std::endl;
    clout << "Shear modulus mue (Phys): "                         << this->getPhysShearModulus() << std::endl;
    clout << "Shear modulus mue (Lattice): "                      << this->getLatticeShearModulus() << std::endl;
    clout << "Bulk modulus K (Lattice): "                         << this->getLatticeBulkModulus() << std::endl;
    clout << "Bulk modulus K (Phys): "                            << this->getPhysBulkModulus() << std::endl;
    clout << "Youngs modulus E (Phys): "                          << this->getPhysYoungsModulus() << std::endl;
    clout << "Youngs modulus E (Lattice): "                       << this->getLatticeYoungsModulus() << std::endl;
    clout << "Lambda (Phys): "                                    << this->getPhysLambda() << std::endl;
    clout << "Lambda (Lattice): "                                 << this->getLatticeLambda() << std::endl;
    clout << "Damping factor Kappa (Lattice): "                   << this->getDampingFactor() << std::endl;
    clout << "-------------------------------------------------------------" << std::endl;

  }

  template <typename T, typename DESCRIPTOR>
  void LinElaUnitConverter<T, DESCRIPTOR>::write(std::string const& fileName) const
  {
    std::string dataFile = singleton::directories().getLogOutDir() + fileName + ".dat";

    if (singleton::mpi().isMainProcessor()) {
      std::ofstream fout;
      fout.open(dataFile.c_str(), std::ios::trunc);

      fout << "UnitConverter information\n\n";
      fout << "----------------- UnitConverter information -----------------\n";
      fout << "-- Parameters:" << std::endl;
      fout << "Resolution:                       N=              " << this->getResolution()                     << "\n";
      fout << "Characteristical length(m):       charL=          " << this->getCharPhysLength()                 << "\n";
      fout << "Characteristical Displacement(m/s):      charU=          " << this->getCharPhysDisplacement()           << "\n";
      fout << "\n";
      fout << "-- Conversion factors:"                                                                          << "\n";
      fout << "Voxel length(m):                  physDeltaX=     " << this->getConversionFactorLength()         << "\n";
      fout << "Time step(s):                     physDeltaT=     " << this->getConversionFactorTime()           << "\n";

      fout << "-------------------------------------------------------------" << "\n";

      fout.close();
    }
  }

  template <typename T, typename DESCRIPTOR>
  class NewLinElaUnitConverter : public UnitConverter<T, DESCRIPTOR> {
  public:

    constexpr NewLinElaUnitConverter(
      T epsilon,
      T charPhysLength,
      T charPhysTime,
      T charPhysDisplacement,
      T youngsModulus,
      T poissonRatio,
      T dampingFactor
    )

      :UnitConverter<T, DESCRIPTOR>(
        epsilon,
        charPhysLength,
        charPhysTime,
        charPhysDisplacement,
        youngsModulus,
        poissonRatio,
        dampingFactor
        ),

        _epsilon              (  epsilon                                                                                                                                               ),
        _charPhysLength       (  charPhysLength                                                                                                                                        ),
        _charPhysTime         (  charPhysTime                                                                                                                                          ),
        _charPhysDisplacement     (  charPhysDisplacement                                                                                                                                      ),
        _physDeltaX           (  charPhysLength * epsilon                                                                                                                              ),
        _physDeltaT           (  charPhysTime   * epsilon * epsilon                                                                                                                    ),
        _youngsModulus        (  youngsModulus                                                                      * charPhysTime / (charPhysLength * charPhysLength * dampingFactor) ),
        _bulkModulus          ( (youngsModulus  / (2. * (1. - poissonRatio)))                                       * charPhysTime / (charPhysLength * charPhysLength * dampingFactor) ),
        _shearModulus         ( (youngsModulus  / (2. * (1. + poissonRatio)))                                       * charPhysTime / (charPhysLength * charPhysLength * dampingFactor) ),
        _lambda               ( (youngsModulus  * poissonRatio                / (1. - poissonRatio * poissonRatio)) * charPhysTime / (charPhysLength * charPhysLength * dampingFactor) ),
        _poissonRatio         (  poissonRatio                                                                                                                                          ),
        _dampingFactor        (  dampingFactor                                                                                                                                         ),
        clout(std::cout, "NewLinElaUnitConv")
      {
      };

    constexpr T getEpsilon( ) const
    {
      return _epsilon;
    };

    constexpr T getConversionFactorLength( ) const
    {
      return _physDeltaX;
    }

    constexpr T getConversionFactorTime( ) const
    {
      return _physDeltaT;
    }

    constexpr T getCharPhysDisplacement( ) const
    {
      return _charPhysDisplacement;
    }

    constexpr T getCharPhysTime( ) const
    {
      return _charPhysTime;
    }

    constexpr T getCharPhysLength( ) const
    {
      return _charPhysLength;
    }

      constexpr T getLatticeShearModulus( ) const
    {
      return _shearModulus;
    };

    constexpr T getPhysShearModulus( ) const
    {
      return _shearModulus * (_charPhysLength * _charPhysLength * _dampingFactor) / _charPhysTime;
    };

    constexpr T getLatticeBulkModulus( ) const
    {
      return _bulkModulus;
    };

      constexpr T getPhysBulkModulus( ) const
    {
      return _bulkModulus * (_charPhysLength * _charPhysLength * _dampingFactor) / _charPhysTime;
    };

    constexpr T getLatticeLambda( ) const
    {
      return _lambda;
    };

    constexpr T getPhysLambda( ) const
    {
      return _lambda * (_charPhysLength * _charPhysLength * _dampingFactor) / _charPhysTime;
    };

    constexpr T getLatticeYoungsModulus() const
    {
      return _youngsModulus;
    };

    constexpr T getPhysYoungsModulus() const
    {
      return _youngsModulus * (_charPhysLength * _charPhysLength * _dampingFactor) / _charPhysTime;
    };

    constexpr T getPoissonRatio() const
    {
      return _poissonRatio;
    }

    constexpr T getDampingFactor() const
    {
      return _dampingFactor;
    };

    void print() const override;

    void write(std::string const& fileName = "NewLinElaUnitConverter") const override;

  protected:
    // lattice units, discretization parameters
    const T _epsilon;
    const T _charPhysLength;
    const T _charPhysTime;
    const T _charPhysDisplacement;
    const T _physDeltaX;
    const T _physDeltaT;
    const T _youngsModulus;
    const T _bulkModulus;
    const T _shearModulus;
    const T _lambda;
    const T _poissonRatio;
    const T _dampingFactor;
    // conversionfactors

  private:
    mutable OstreamManager clout;
  };

  template <typename T, class DESCRIPTOR>
  void NewLinElaUnitConverter<T, DESCRIPTOR>::print() const

  {
    clout << "----------------- UnitConverter information -----------------" << std::endl;
    clout << "-- Parameters:" << std::endl;
    clout << "Epsilon:                                          " << this->getEpsilon() << std::endl;
    clout << "Characteristical time(s):         charT=          " << this->getCharPhysTime() << std::endl;
    clout << "Characteristical length(m):       charL=          " << this->getCharPhysLength() << std::endl;
    clout << "Characteristical Displacement(m/s):      charU=          " << this->getCharPhysDisplacement() << std::endl;
    clout << std::endl;

    clout << std::endl;
    clout << "-------------------- Conversion factors ----------------------" << std::endl;
    clout << "Voxel length(m):                  physDeltaX=       " << this->getConversionFactorLength() << std::endl;
    clout << "Time step(s):                     physDeltaT=       " << this->getConversionFactorTime() << std::endl;
    clout << "Displacement factor(m/s):         physDisplacement= " << this->getConversionFactorDisplacement() << std::endl;
    std::cout << std::endl;
    clout << "Youngs modulus (Phys):            physE=            " << this->getPhysYoungsModulus() << std::endl;
    clout << "Youngs modulus (Lattice):         latticeE=         " << this->getLatticeYoungsModulus() << std::endl;
    clout << "Shear modulus (Phys):             physMue=          " << this->getPhysShearModulus() << std::endl;
    clout << "Shear modulus (Lattice):          latticeMue=       " << this->getLatticeShearModulus() << std::endl;
    clout << "Bulk modulus (Phys):              physK=            " << this->getPhysBulkModulus() << std::endl;
    clout << "Bulk modulus (Lattice):           latticeK=         " << this->getLatticeBulkModulus() << std::endl;
    clout << "Lambda (Phys):                    physLambda=       " << this->getPhysLambda() << std::endl;
    clout << "Lambda (Lattice):                 latticeLambda=    " << this->getLatticeLambda() << std::endl;
    clout << "Damping factor:                   kappa=            " << this->getDampingFactor() << std::endl;
    clout << "Poisson ratio:                    nu=               " << this->getPoissonRatio() << std::endl;
    clout << "-------------------------------------------------------------" << std::endl;

  }

  template <typename T, typename DESCRIPTOR>
  void NewLinElaUnitConverter<T, DESCRIPTOR>::write(std::string const& fileName) const
  {
    std::string dataFile = singleton::directories().getLogOutDir() + fileName + ".dat";

    if (singleton::mpi().isMainProcessor())
    {
      std::ofstream fout;
      fout.open(dataFile.c_str(), std::ios::trunc);

      fout << "----------------- UnitConverter information -----------------" << std::endl;
      fout << "-- Parameters:" << std::endl;
      fout << "Epsilon:                                   " << this->getEpsilon() << std::endl;
      fout << "Characteristical time(s):           charT= " << this->getCharPhysTime() << std::endl;
      fout << "Characteristical length(m):         charL= " << this->getCharPhysLength() << std::endl;
      fout << "Characteristical Displacement(m/s): charU= " << this->getCharPhysDisplacement() << std::endl;
      fout << std::endl;

      fout << "-------------------- Conversion factors ----------------------" << std::endl;
      fout << "Voxel length(m):                  physDeltaX=        " << this->getConversionFactorLength() << std::endl;
      fout << "Time step(s):                     physDeltaT=        " << this->getConversionFactorTime() << std::endl;
      fout << "Displacement factor(m/s):         physDisplacement=  " << this->getConversionFactorDisplacement() << std::endl;
      fout << std::endl;

      fout << "-------------------- Material Properties ----------------------" << std::endl;
      fout << "Youngs modulus (Phys):            physE=          " << this->getPhysYoungsModulus() << std::endl;
      fout << "Youngs modulus (Lattice):         latticeE=       " << this->getLatticeYoungsModulus() << std::endl;
      fout << "Shear modulus (Phys):             physMue=        " << this->getPhysShearModulus() << std::endl;
      fout << "Shear modulus (Lattice):          latticeMue=     " << this->getLatticeShearModulus() << std::endl;
      fout << "Bulk modulus (Phys):              physK=          " << this->getPhysBulkModulus() << std::endl;
      fout << "Bulk modulus (Lattice):           latticeK=       " << this->getLatticeBulkModulus() << std::endl;
      fout << "Lambda (Phys):                    physLambda=     " << this->getPhysLambda() << std::endl;
      fout << "Lambda (Lattice):                 latticeLambda=  " << this->getLatticeLambda() << std::endl;
      fout << "Damping factor:                   kappa=          " << this->getDampingFactor() << std::endl;
      fout << "Poisson ratio:                    nu=             " << this->getPoissonRatio() << std::endl;
      fout << "-------------------------------------------------------------" << "\n";

      fout.close();
    }
  }

  template <typename T, typename DESCRIPTOR>
  class LinElaUnitConverterNoEpsilon : public NewLinElaUnitConverter<T, DESCRIPTOR>
  {
    public:
      constexpr LinElaUnitConverterNoEpsilon(
        T physDeltaX,
        T physDeltaT,
        T charPhysLength,
        T charPhysDisplacement,
        T youngsModulus,
        T poissonRatio,
        T dampingFactor) : NewLinElaUnitConverter<T, DESCRIPTOR>(
          (physDeltaX /charPhysLength),
          charPhysLength,
          (physDeltaT / (physDeltaX * physDeltaX / (charPhysLength * charPhysLength))),
          charPhysDisplacement,
          youngsModulus,
          poissonRatio,
          dampingFactor)
        {
        }
  };
}
#endif
