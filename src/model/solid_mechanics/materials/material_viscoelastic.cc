/**
 * @file   material_viscoelastic.hh
 * @author Vlad Yastrebov <vladislav.yastrebov@epfl.ch> 
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Thu Feb 7 2012
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
#include "material_viscoelastic.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialViscoElastic<spatial_dimension>::MaterialViscoElastic(SolidMechanicsModel & model, const ID & id)  :
  Material(model, id), MaterialElastic<spatial_dimension>(model, id),
  stress_dev("stress_dev", id),
  history_integral("history_integral", id) {
  AKANTU_DEBUG_IN();

  eta = 1.;
  E_inf = 1.;
  Ev  = 1.;

  UInt stress_size = spatial_dimension * spatial_dimension;

  this->initInternalVector(this->stress_dev, stress_size);
  this->initInternalVector(this->history_integral, stress_size);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialViscoElastic<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  E_inf = this->E - this->Ev;

  MaterialElastic<spatial_dimension>::initMaterial();

  this->resizeInternalVector(this->stress_dev);
  this->resizeInternalVector(this->history_integral);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialViscoElastic<spatial_dimension>::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real tau = 0.;
  // if(std::abs(Ev) > std::numeric_limits<Real>::epsilon())
  tau = eta / Ev;

  Vector<Real> & stress_dev_vect  = stress_dev(el_type, ghost_type);
  Vector<Real> & history_int_vect = history_integral(el_type, ghost_type);

  Vector<Real>::iterator<types::Matrix> stress_d = stress_dev_vect.begin(spatial_dimension, spatial_dimension);
  Vector<Real>::iterator<types::Matrix> history_int = history_int_vect.begin(spatial_dimension, spatial_dimension);

  types::Matrix e(spatial_dimension, spatial_dimension);
  types::Matrix s(spatial_dimension, spatial_dimension);

  types::Matrix theta_sp(spatial_dimension, spatial_dimension);

  Real dt = this->model->getTimeStep();
  Real exp_dt_tau = exp( -dt/tau );
  Real exp_dt_tau_2 = exp( -.5*dt/tau );

  Real plane_stress_coeff = this->nu / (this->nu - 1);

  /// Loop on all quadrature points
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;

  types::Matrix & dev_s = *stress_d;
  types::Matrix & h = *history_int;

  e.clear();
  s.clear();

  sigma.clear();

  /// Compute the first invariant of strain
  Real Theta = grad_u.trace();

  ///\todo correct the trace in case of plane stress
  if(spatial_dimension == 2 && this->plane_stress == true)
    Theta += plane_stress_coeff * (grad_u(0,0) + grad_u(1,1));

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j) {
      s(i, j) = 2 * this->mu * (.5 * (grad_u(i,j) + grad_u(j,i)) - 1./3. * Theta *(i == j));
    }
  
  /// Update the internal variable
  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j) {
      h(i, j)  = exp_dt_tau * h(i, j) + exp_dt_tau_2 * (s(i, j) - dev_s(i, j));
      dev_s(i, j) = s(i, j);
    }

  types::Matrix U_rond_prim(spatial_dimension, spatial_dimension);
  U_rond_prim.eye(this->kpa * Theta);

  Real gamma_inf = E_inf / this->E;
  Real gamma_v   = Ev    / this->E;

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      sigma(i, j) =  U_rond_prim(i,j) + gamma_inf * s(i, j) + gamma_v * h(i, j);

  /// Save the deviator of stress
  ++stress_d;
  ++history_int;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
bool MaterialViscoElastic<spatial_dimension>::setParam(const std::string & key, const std::string & value,
				      const ID & id) {
  std::stringstream sstr(value);
  if(key == "eta") { sstr >> eta; }
  else if(key == "Ev") { sstr >> Ev; }
  else { return MaterialElastic<spatial_dimension>::setParam(key, value, id); }
  return true;
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialViscoElastic<spatial_dimension>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "MaterialViscoElastic [" << std::endl;
  MaterialElastic<spatial_dimension>::printself(stream, indent + 1);
  stream << space << " + Eta : " << eta << std::endl;
  stream << space << " + Ev  : " << Ev << std::endl;
  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialViscoElastic);

__END_AKANTU__
