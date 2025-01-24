/*  This file is part of the OpenLB library
*
*  Copyright (C) 2012-2024 Mathias J. Krause, Benjamin Förster, Julius Jeßberger
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


#ifndef OPTI_CASE_DUAL_HH
#define OPTI_CASE_DUAL_HH

#include "optiCaseDual.h"
#include "serialization.h"

namespace olb {

namespace opti {

template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
void OptiCaseDual<S,SOLVER,C>::readFromXML(XMLreader const& xml)
{
  std::string type ("");
  xml.readOrWarn<std::string>("Optimization", "ControlType", "", type);
  _controlType = ((type == "Force") || (type == "force")) ? ForceControl : PorosityControl;
  xml.readOrWarn<std::string>("Optimization", "Projection", "", _projectionName);

  type = "";
  xml.readOrWarn<std::string>("Optimization", "StartValueType", "", type);
  if (type == "Porosity") {
    _startValueType = Porosity;
  } else if (type == "Permeability") {
    _startValueType = Permeability;
  } else if (type == "Control") {
    _startValueType = Control;
  } else {
    _startValueType = ProjectedControl;
  }
   xml.readOrWarn<bool>("Optimization", "ReferenceSolution", "", _computeReference);
}

template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
void OptiCaseDual<S,SOLVER,C>::initialize(XMLreader const& xml)
{
  _converter = createUnitConverter<S,descriptor>(xml);

  if (_controlType == ForceControl) {
    _fieldDim = dim;
  } else {
    _fieldDim = 1;
  }

  _primalSolver = createLbSolver <SOLVER<S,SolverMode::Primal>> (xml);
  _dualSolver = createLbSolver <SOLVER<S,SolverMode::Dual>> (xml);
  this->_postEvaluation = std::bind(&SOLVER<S,SolverMode::Primal>::postProcess, _primalSolver);

  _primalSolver->initialize();
  _controlIndicator = _primalSolver->parameters(names::Opti()).designDomain;
  _primalGeometry = _primalSolver->parameters(names::Results()).geometry;
  _refLattice = std::make_shared<SuperLattice<S,descriptor>>(*_primalGeometry);

  _serializer = std::make_shared<SimpleGeometrySerializer<S,dim>>(*_primalGeometry);
  _dimCtrl = _serializer->getNoCells() * _fieldDim;

  _controller = new Controller<S>(_dimCtrl);
}

template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
void OptiCaseDual<S,SOLVER,C>::initializeFields()
{
  _projection = projection::construct<S,descriptor>(*_converter, _projectionName);

  _projectedControl = std::make_shared<SuperLatticeSerialDataF<S,descriptor>>(
    *_refLattice,
    *_controller,
    _fieldDim,
    _serializer,
    [this](S x) { return _projection->project(x); });
  _dProjectionDcontrol = std::make_shared<SuperLatticeSerialDataF<S,descriptor>>(
    *_refLattice,
    *_controller,
    _fieldDim,
    _serializer,
    [this](S x) { return _projection->derivative(x); });
  _primalSolver->parameters(names::Opti()).controlledField = _projectedControl;
  _dualSolver->parameters(names::Opti()).controlledField = _projectedControl;
}


template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
S OptiCaseDual<S,SOLVER,C>::evaluateObjective(
  const C& control, unsigned optiStep)
{
  _controller->setControl(control, _dimCtrl);
  _primalSolver->parameters(names::OutputOpti()).counterOptiStep = optiStep;
  _primalSolver->solve();

  _objective->setPrimalSolver(_primalSolver);
  return _objective->evaluate();
}


template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
void OptiCaseDual<S,SOLVER,C>::computeDerivatives(
  const C& control, C& derivatives, unsigned optiStep)
{
  _controller->setControl(control, _dimCtrl);
  _dualSolver->parameters(names::OutputOpti()).counterOptiStep = optiStep;

  const auto& primalResults = _primalSolver->parameters(names::Results());
  auto& dualParams = _dualSolver->parameters(names::Opti());

  dualParams.fpop = std::make_shared<SuperLatticeFpop3D<S,descriptor>>(*(primalResults.lattice));
  dualParams.dObjectiveDf = _objective->derivativeByPopulations();

  _dualSolver->solve();

  derivativesFromDualSolution(derivatives);
}

template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
void OptiCaseDual<S,SOLVER,C>::derivativesFromDualSolution(
  C& derivatives)
{
  const auto& primalGeometry = _primalSolver->parameters(names::Results()).geometry;
  const S omega = _converter->getLatticeRelaxationFrequency();
  const int nC = primalGeometry->getCuboidGeometry().size();

  // add derivative of regularizing term
  for (int iC = 0; iC < nC; iC++) {

    const Vector<int,3> extend = primalGeometry->getCuboidGeometry().get(iC).getExtent();

    for (int iX = 0; iX < extend[0]; iX++) {
      for (int iY = 0; iY < extend[1]; iY++) {
        for (int iZ = 0; iZ < extend[2]; iZ++) {
          const LatticeR<dim+1> latticeR(iC, iX, iY, iZ);

          if (evaluateSuperIndicatorFglobally<dim,S>(*_controlIndicator, latticeR.data())) {
            C derivativesHelp(_fieldDim, 0);

            if (primalGeometry->getLoadBalancer().rank(iC)  == singleton::mpi().getRank()) {
              const S rho_f = _primalSolver->parameters(names::Results()).lattice->get(latticeR).computeRho();
              S dProjectionDcontrol[_fieldDim];
              (*_dProjectionDcontrol)(dProjectionDcontrol, latticeR.data());
              S dObjectiveDcontrol[_fieldDim];
              (*_objective->derivativeByControl())(dObjectiveDcontrol, latticeR.data());

              // compute derivative coming from dependence of objective on state -> use adjoint solution
              if ( _controlType == ForceControl ) {
                for (int iDim=0; iDim<_fieldDim; iDim++) {
                  for (int jPop=0; jPop < descriptor::q; ++jPop) {
                    const S phi_j = _dualSolver->parameters(names::Results()).lattice->get(latticeR)[jPop];
                    derivativesHelp[iDim] -= rho_f * descriptors::t<S,descriptor>(jPop)
                       * descriptors::invCs2<S,descriptor>() * descriptors::c<descriptor>(jPop,iDim) * phi_j;
                  }
                }
              }
              else if ( _controlType == PorosityControl ) {
                S u_f[3];
                _primalSolver->parameters(names::Results()).lattice->get(latticeR).computeU(u_f);
                const S d = _dualSolver->parameters(names::Results()).lattice->get(latticeR).template getField<descriptors::POROSITY>();

                for (int jPop = 0; jPop < descriptor::q; ++jPop) {
                  const S phi_j = _dualSolver->parameters(names::Results()).lattice->get(latticeR)[jPop];
                  const S feq_j = equilibrium<descriptor>::secondOrder(jPop, rho_f, u_f) + descriptors::t<S,descriptor>(jPop);
                  for (int iDim = 0; iDim < descriptor::d; iDim++) {
                    derivativesHelp[0] +=  phi_j*feq_j*( descriptors::c<descriptor>(jPop,iDim) - d*u_f[iDim] )*u_f[iDim]*dProjectionDcontrol[0];
                  }
                }
                derivativesHelp[0] *= -omega*descriptors::invCs2<S,descriptor>();
              }
              // add (partial) derivative coming from direct dependence of objective on control
              for (int iDim=0; iDim<_fieldDim; iDim++) {
                derivativesHelp[iDim] += dObjectiveDcontrol[iDim];  // * dProjectionDcontrol[iDim];  // todo uncomment if objective depends on projection
              }
            }

#ifdef PARALLEL_MODE_MPI
            singleton::mpi().bCast(&derivativesHelp[0], _fieldDim, primalGeometry->getLoadBalancer().rank(iC));
#endif
            for (int iDim=0; iDim<_fieldDim; ++iDim) {
              const int index = _serializer->getSerializedComponentIndex(latticeR, iDim, _fieldDim);
              derivatives[index] = derivativesHelp[iDim];
            }
          }
        }
      }
    }
  }
}


} // namespace opti

} // namespace olb

#endif
