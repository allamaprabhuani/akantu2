/**
 * @file   material_viscoelastic.hh
 * @author Vlad Yastrebov <vladislav.yastrebov@epfl.ch>
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
MaterialViscoElastic::MaterialViscoElastic(Model & model, const ID & id)  :
  Material(model, id), MaterialElastic(model, id),
  stress_dev("stress_dev", id),
  history_integral("history_integral", id) {
  AKANTU_DEBUG_IN();

  eta = 1.;
  Ev  = 1.;

  UInt stress_size = spatial_dimension * spatial_dimension;

  initInternalVector(this->stress_dev, stress_size);
  initInternalVector(this->history_integral, stress_size);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialViscoElastic::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialElastic::initMaterial();

  resizeInternalVector(this->stress_dev);
  resizeInternalVector(this->history_integral);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialViscoElastic::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real tau = eta / Ev;
  Real K = E / ( 3 * (1 - 2*nu) ); // kpa ?

  Vector<Real> & stress_dev_vect  = stress_dev(el_type, ghost_type);
  Vector<Real> & history_int_vect = history_integral(el_type, ghost_type);

  Vector<Real>::iterator<types::Matrix> stress_d = stress_dev_vect.begin(spatial_dimension, spatial_dimension);
  Vector<Real>::iterator<types::Matrix> history_int = history_int_vect.begin(spatial_dimension, spatial_dimension);


  types::Matrix e(spatial_dimension, spatial_dimension);
  types::Matrix s(spatial_dimension, spatial_dimension);
  types::Matrix theta_sp(spatial_dimension, spatial_dimension);

  Real dt = model->getTimeStep();
  Real exp_dt_thau = exp( -dt/tau );
  Real exp_dt_thau_2 = exp( -.5*dt/tau );

  /// Loop on all quadrature points
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;

  types::Matrix sigma(stress_val, spatial_dimension, spatial_dimension);
  types::Matrix F(strain_val, spatial_dimension, spatial_dimension);

  e.clear();
  s.clear();
  sigma.clear();

  /// Compute the first invariant of strain
  Real Theta = F.trace();

  theta_sp.eye(1/double(spatial_dimension) * Theta);

  /// Compute the deviator of strain @f$ @f$
  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j) {
      e(i, j) = F(i, j) - theta_sp(i, j);
      s(i, j) = E / (1 + nu) * e(i, j);
    }


  /// Update the internal variable
  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      (*history_int)(i, j)  = exp_dt_thau * (*history_int)(i, j) + exp_dt_thau_2 * (s(i, j) - (*stress_d)(i, j));

  Real alpha = 2./3. * K * Theta;
  Real beta = 1. / ( Ev + E );
  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      sigma(i, j) = (i == j) * alpha + beta * (E * s(i, j) + Ev * (*history_int)(i, j));


  /// Save the deviator of stress
  stress_d->copy(s);


  ++stress_d;
  ++history_int;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialViscoElastic::computePotentialEnergy(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();
/*
  if(ghost_type != _not_ghost) return;

  Vector<Real> & stress_el = stress_elastic(el_type, ghost_type);
  Real * stress_el_val = stress_el.storage();

  Real * epot = potential_energy(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;

  MaterialElastic::computePotentialEnergy(strain_val, stress_el_val, epot);
  epot++;
  stress_el_val += spatial_dimension*spatial_dimension;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
*/
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool MaterialViscoElastic::setParam(const std::string & key, const std::string & value,
				      const ID & id) {
  std::stringstream sstr(value);
  if(key == "eta") { sstr >> eta; }
  else if(key == "Ev") { sstr >> Ev; }
  else { return MaterialElastic::setParam(key, value, id); }
  return true;
}


/* -------------------------------------------------------------------------- */
void MaterialViscoElastic::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "MaterialViscoElastic [" << std::endl;
  MaterialElastic::printself(stream, indent + 1);
  stream << space << " + Eta : " << eta << std::endl;
  stream << space << " + Ev  : " << Ev << std::endl;
  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */


__END_AKANTU__
