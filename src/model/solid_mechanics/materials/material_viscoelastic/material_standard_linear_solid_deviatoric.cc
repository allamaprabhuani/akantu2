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
  history_integral("history_integral", id),
  dissipated_energy("dissipated_energy", id) {

  AKANTU_DEBUG_IN();

  this->registerParam("Eta",  eta,   1., ParamAccessType(_pat_parsable | _pat_modifiable), "Viscosity");
  this->registerParam("Ev",   Ev,    1., ParamAccessType(_pat_parsable | _pat_modifiable), "Stiffness of the viscous element");
  this->registerParam("Einf", E_inf, 1., ParamAccessType(_pat_readable), "Stiffness of the elastic element");

  UInt stress_size = spatial_dimension * spatial_dimension;

  this->initInternalArray(this->stress_dev, stress_size);
  this->initInternalArray(this->history_integral, stress_size);
  this->initInternalArray(this->dissipated_energy, 1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialStandardLinearSolidDeviatoric<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  updateInternalParameters();

  MaterialElastic<spatial_dimension>::initMaterial();

  this->resizeInternalArray(this->stress_dev);
  this->resizeInternalArray(this->history_integral);
  this->resizeInternalArray(this->dissipated_energy);

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

  Array<Real> & stress_dev_vect  = stress_dev(el_type, ghost_type);
  Array<Real> & history_int_vect = history_integral(el_type, ghost_type);

  Array<Real>::iterator< Matrix<Real> > stress_d = stress_dev_vect.begin(spatial_dimension, spatial_dimension);
  Array<Real>::iterator< Matrix<Real> > history_int = history_int_vect.begin(spatial_dimension, spatial_dimension);

  /// Loop on all quadrature points
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  Matrix<Real> & dev_s = *stress_d;
  Matrix<Real> & h = *history_int;

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

  Array<Real> & stress_dev_vect  = stress_dev(el_type, ghost_type);
  Array<Real> & history_int_vect = history_integral(el_type, ghost_type);

  Array<Real>::iterator< Matrix<Real> > stress_d = stress_dev_vect.begin(spatial_dimension, spatial_dimension);
  Array<Real>::iterator< Matrix<Real> > history_int = history_int_vect.begin(spatial_dimension, spatial_dimension);

  Matrix<Real> s(spatial_dimension, spatial_dimension);

  Real dt = this->model->getTimeStep();
  Real exp_dt_tau = exp( -dt/tau );
  Real exp_dt_tau_2 = exp( -.5*dt/tau );

  Matrix<Real> epsilon_d(spatial_dimension, spatial_dimension);
  Matrix<Real> epsilon_v(spatial_dimension, spatial_dimension);

  /// Loop on all quadrature points
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  Matrix<Real> & dev_s = *stress_d;
  Matrix<Real> & h = *history_int;

  s.clear();
  sigma.clear();

  /// Compute the first invariant of strain
  Real gamma_inf = E_inf / this->E;
  Real gamma_v   = Ev    / this->E;

  this->gradUToEpsilon(grad_u, epsilon_d);
  Real Theta = epsilon_d.trace();
  epsilon_v.eye(1./3. * Theta);
  epsilon_d -= epsilon_v;

  Matrix<Real> U_rond_prim(spatial_dimension, spatial_dimension);

  U_rond_prim.eye(gamma_inf * this->kpa * Theta);

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j) {
      s(i, j) = 2 * this->mu * epsilon_d(i, j);
      h(i, j)  = exp_dt_tau * h(i, j) + exp_dt_tau_2 * (s(i, j) - dev_s(i, j));
      dev_s(i, j) = s(i, j);
      sigma(i, j) =  U_rond_prim(i,j) + gamma_inf * s(i, j) + gamma_v * h(i, j);
    }

  /// Save the deviator of stress
  ++stress_d;
  ++history_int;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  this->updateDissipatedEnergy(el_type, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialStandardLinearSolidDeviatoric<spatial_dimension>::updateDissipatedEnergy(ElementType el_type,
										      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  // if(ghost_type == _ghost) return 0.;
  
  Real tau = 0.;
  tau = eta / Ev;
  Real * dis_energy = dissipated_energy(el_type, ghost_type).storage();

  Array<Real> & stress_dev_vect  = stress_dev(el_type, ghost_type);
  Array<Real> & history_int_vect = history_integral(el_type, ghost_type);

  Array<Real>::iterator< Matrix<Real> > stress_d = stress_dev_vect.begin(spatial_dimension, spatial_dimension);
  Array<Real>::iterator< Matrix<Real> > history_int = history_int_vect.begin(spatial_dimension, spatial_dimension);

  Matrix<Real> q(spatial_dimension, spatial_dimension);
  Matrix<Real> q_rate(spatial_dimension, spatial_dimension);
  Matrix<Real> epsilon_d(spatial_dimension, spatial_dimension);
  Matrix<Real> epsilon_v(spatial_dimension, spatial_dimension);

  Real dt = this->model->getTimeStep();

  Real gamma_v   = Ev / this->E;
  Real alpha = 1. / (2. * this->mu * gamma_v);

  /// Loop on all quadrature points
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  Matrix<Real> & dev_s = *stress_d;
  Matrix<Real> & h = *history_int;

  /// Compute the first invariant of strain
  this->gradUToEpsilon(grad_u, epsilon_d);

  Real Theta = epsilon_d.trace();
  epsilon_v.eye(1./3. * Theta);
  epsilon_d -= epsilon_v;

  q.copy(dev_s);
  q -= h;
  q *= gamma_v;

  q_rate.copy(dev_s);
  q_rate *= gamma_v;
  q_rate -= q;
  q_rate /= tau;

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      *dis_energy +=  (epsilon_d(i, j) - alpha * q(i, j)) * q_rate(i, j) * dt;


  /// Save the deviator of stress
  ++stress_d;
  ++history_int;
  ++dis_energy;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialStandardLinearSolidDeviatoric<spatial_dimension>::getDissipatedEnergy() const {
  AKANTU_DEBUG_IN();

  Real de = 0.;
  const Mesh & mesh = this->model->getFEM().getMesh();

  /// integrate the dissipated energy for each type of elements
  Mesh::type_iterator it  = mesh.firstType(spatial_dimension, _not_ghost);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, _not_ghost);

  for(; it != end; ++it) {
    de += this->model->getFEM().integrate(dissipated_energy(*it, _not_ghost), *it,
					  _not_ghost, this->element_filter(*it, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return de;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialStandardLinearSolidDeviatoric<spatial_dimension>::getDissipatedEnergy(ElementType type, UInt index) const {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = this->model->getFEM().getNbQuadraturePoints(type);
  Array<Real>::const_iterator< Vector<Real> > it = this->dissipated_energy(type, _not_ghost).begin(nb_quadrature_points);
  UInt gindex = (this->element_filter(type, _not_ghost))(index);

  AKANTU_DEBUG_OUT();
  return this->model->getFEM().integrate(it[index], type, gindex);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialStandardLinearSolidDeviatoric<spatial_dimension>::getEnergy(std::string type) {
  if(type == "dissipated") return getDissipatedEnergy();
  else if(type == "dissipated_sls_deviatoric") return getDissipatedEnergy();
  else return MaterialElastic<spatial_dimension>::getEnergy(type);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialStandardLinearSolidDeviatoric<spatial_dimension>::getEnergy(std::string energy_id, ElementType type, UInt index) {
  if(energy_id == "dissipated") return getDissipatedEnergy(type, index);
  else if(energy_id == "dissipated_sls_deviatoric") return getDissipatedEnergy(type, index);
  else return MaterialElastic<spatial_dimension>::getEnergy(energy_id, type, index);
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialStandardLinearSolidDeviatoric);

__END_AKANTU__
