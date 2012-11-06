/**
 * @file   material_standard_linear_solid_deviatoric.cc
 *
 * @author Vladislav Yastrebov <vladislav.yastrebov@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date   Wed Feb 08 16:53:31 2012
 *
 * @brief  Material Visco-elastic
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_standard_linear_solid_deviatoric.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialStandardLinearSolidDeviatoric<spatial_dimension>::MaterialStandardLinearSolidDeviatoric(SolidMechanicsModel & model, const ID & id)  :
  Material(model, id), MaterialElastic<spatial_dimension>(model, id),
  stress_dev("stress_dev", id),
  history_integral("history_integral", id) {
  AKANTU_DEBUG_IN();

  this->registerParam("Eta",  eta,   1., ParamAccessType(_pat_parsable | _pat_modifiable), "Viscosity");
  this->registerParam("Ev",   Ev,    1., ParamAccessType(_pat_parsable | _pat_modifiable), "Stiffness of the viscous element");
  this->registerParam("Einf", E_inf, 1., ParamAccessType(_pat_readable), "Stiffness of the elastic element");

  UInt stress_size = spatial_dimension * spatial_dimension;

  this->initInternalVector(this->stress_dev, stress_size);
  this->initInternalVector(this->history_integral, stress_size);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialStandardLinearSolidDeviatoric<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  updateInternalParameters();

  MaterialElastic<spatial_dimension>::initMaterial();

  this->resizeInternalVector(this->stress_dev);
  this->resizeInternalVector(this->history_integral);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialStandardLinearSolidDeviatoric<spatial_dimension>::updateInternalParameters() {
  MaterialElastic<spatial_dimension>::updateInternalParameters();
  E_inf = this->E - this->Ev;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialStandardLinearSolidDeviatoric<spatial_dimension>::setToSteadyState(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Vector<Real> & stress_dev_vect  = stress_dev(el_type, ghost_type);
  Vector<Real> & history_int_vect = history_integral(el_type, ghost_type);

  Vector<Real>::iterator<types::RMatrix> stress_d = stress_dev_vect.begin(spatial_dimension, spatial_dimension);
  Vector<Real>::iterator<types::RMatrix> history_int = history_int_vect.begin(spatial_dimension, spatial_dimension);

  /// Loop on all quadrature points
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  types::RMatrix & dev_s = *stress_d;
  types::RMatrix & h = *history_int;

  /// Compute the first invariant of strain
  Real Theta = grad_u.trace();

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j) {
      dev_s(i, j) = 2 * this->mu * (.5 * (grad_u(i,j) + grad_u(j,i)) - 1./3. * Theta *(i == j));
      h(i, j) = 0.;
    }

  /// Save the deviator of stress
  ++stress_d;
  ++history_int;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialStandardLinearSolidDeviatoric<spatial_dimension>::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real tau = 0.;
  // if(std::abs(Ev) > std::numeric_limits<Real>::epsilon())
  tau = eta / Ev;

  Vector<Real> & stress_dev_vect  = stress_dev(el_type, ghost_type);
  Vector<Real> & history_int_vect = history_integral(el_type, ghost_type);

  Vector<Real>::iterator<types::RMatrix> stress_d = stress_dev_vect.begin(spatial_dimension, spatial_dimension);
  Vector<Real>::iterator<types::RMatrix> history_int = history_int_vect.begin(spatial_dimension, spatial_dimension);

  types::RMatrix s(spatial_dimension, spatial_dimension);

  Real dt = this->model->getTimeStep();
  Real exp_dt_tau = exp( -dt/tau );
  Real exp_dt_tau_2 = exp( -.5*dt/tau );

  /// Loop on all quadrature points
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  types::RMatrix & dev_s = *stress_d;
  types::RMatrix & h = *history_int;

  s.clear();
  sigma.clear();

  /// Compute the first invariant of strain
  Real Theta = grad_u.trace();

  Real gamma_inf = E_inf / this->E;
  Real gamma_v   = Ev    / this->E;

  types::RMatrix U_rond_prim(spatial_dimension, spatial_dimension);
  U_rond_prim.eye(gamma_inf * this->kpa * Theta);

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j) {
      s(i, j) = 2 * this->mu * (.5 * (grad_u(i,j) + grad_u(j,i)) - 1./3. * Theta *(i == j));
      h(i, j)  = exp_dt_tau * h(i, j) + exp_dt_tau_2 * (s(i, j) - dev_s(i, j));
      dev_s(i, j) = s(i, j);
      sigma(i, j) =  U_rond_prim(i,j) + gamma_inf * s(i, j) + gamma_v * h(i, j);
    }

  /// Save the deviator of stress
  ++stress_d;
  ++history_int;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */


INSTANSIATE_MATERIAL(MaterialStandardLinearSolidDeviatoric);

__END_AKANTU__
