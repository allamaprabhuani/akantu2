/**
 * @file   material_damage_iterative.cc
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Fri Feb 14 14:54:44 2014
 *
 * @brief  Specialization of the class material damage to damage only one gauss
 * point at a time and propagate damage in a linear way. Max principal stress
 * criterion is used as a failure criterion.
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
#include "material_damage_iterative.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialDamageIterative<spatial_dimension>::MaterialDamageIterative(SolidMechanicsModel & model,
						  const ID & id)  :
  Material(model, id),
  MaterialDamage<spatial_dimension>(model, id),
  Sc("Sc", *this),
  equivalent_stress("equivalent_stress", *this),
  norm_max_equivalent_stress(0),
  q_point() {
  AKANTU_DEBUG_IN();

  this->registerParam("Sc",                  Sc,                  _pat_parsable, "critical stress threshold");
  this->registerParam("prescribed_dam",      prescribed_dam, 0.1, _pat_parsable | _pat_modifiable, "increase of damage in every step" );
  this->registerParam("dam_threshold",       dam_threshold,  0.8,  _pat_parsable | _pat_modifiable, "damage threshold at which damage damage will be set to 1" );


  this->use_previous_stress           = true;
  this->use_previous_strain           = true;
  this->Sc.initialize(1);
  this->equivalent_stress.initialize(1);

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialDamageIterative<spatial_dimension>::computeNormalizedEquivalentStress(const Array<Real> & grad_u,
										   ElementType el_type,
										   GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  /// Vector to store eigenvalues of current stress tensor
  Vector<Real> eigenvalues(spatial_dimension);

  Array<Real>::const_iterator<Real> Sc_it = Sc(el_type).begin();
  Array<Real>::iterator<Real> equivalent_stress_it = equivalent_stress(el_type).begin();

  Array<Real>::const_matrix_iterator grad_u_it = grad_u.begin(spatial_dimension,
							      spatial_dimension);
  Array<Real>::const_matrix_iterator grad_u_end = grad_u.end(spatial_dimension,
							     spatial_dimension);
  Real * dam = this->damage(el_type, ghost_type).storage();
  Matrix<Real> sigma(spatial_dimension, spatial_dimension);
  for(;grad_u_it != grad_u_end; ++ grad_u_it) {
    sigma.clear();
    MaterialElastic<spatial_dimension>::computeStressOnQuad(*grad_u_it, sigma, 0.);
    computeDamageAndStressOnQuad(sigma,*dam);

    /// compute eigenvalues
    sigma.eig(eigenvalues);
    /// find max eigenvalue and normalize by tensile strength
    *equivalent_stress_it = *(std::max_element(eigenvalues.storage(),
					       eigenvalues.storage() + spatial_dimension)) / *(Sc_it);
    ++Sc_it;
    ++equivalent_stress_it;
    ++dam;
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialDamageIterative<spatial_dimension>::computeAllStresses(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  /// reset Gauss point to store the information where highest critical stress occurs
  if(ghost_type==_not_ghost) {
    norm_max_equivalent_stress = 0;
    q_point.type = _not_defined;
    q_point.global_num = 0;
    q_point.ghost_type = ghost_type;
  }
  MaterialDamage<spatial_dimension>::computeAllStresses(ghost_type);

  /// find global Gauss point with highest stress
  StaticCommunicator & comm = akantu::StaticCommunicator::getStaticCommunicator();
  Real global_norm_max_equivalent_stress = norm_max_equivalent_stress;
  comm.allReduce(&global_norm_max_equivalent_stress, 1, _so_max);

  if (!Math::are_float_equal(global_norm_max_equivalent_stress,
  			     norm_max_equivalent_stress)) {
    q_point.type = _not_defined;
    q_point.global_num = 0;
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialDamageIterative<spatial_dimension>::findMaxNormalizedEquivalentStress(ElementType el_type,
						      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if(ghost_type==_not_ghost) {

    const Array<Real> & e_stress = equivalent_stress(el_type);

    if (e_stress.begin() != e_stress.end() ) {
      Array<Real>::const_iterator<Real> equivalent_stress_it_max = std::max_element(e_stress.begin(),e_stress.end());
      /// check if max equivalent stress for this element type is greater than the current norm_max_eq_stress
      if (*equivalent_stress_it_max > norm_max_equivalent_stress) {
	norm_max_equivalent_stress = *equivalent_stress_it_max;
	q_point.type = el_type;
	q_point.global_num = equivalent_stress_it_max-e_stress.begin();
      }
    }
  }
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialDamageIterative<spatial_dimension>::computeStress(ElementType el_type,
							       GhostType ghost_type) {
  AKANTU_DEBUG_IN();


  MaterialDamage<spatial_dimension>::computeStress(el_type, ghost_type);

  Real * dam = this->damage(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  computeDamageAndStressOnQuad(sigma,*dam);
  ++dam;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  computeNormalizedEquivalentStress(this->strain(el_type, ghost_type), el_type, ghost_type);
  norm_max_equivalent_stress = 0;
  findMaxNormalizedEquivalentStress(el_type, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
UInt MaterialDamageIterative<spatial_dimension>::updateDamage() {
  UInt nb_damaged_elements = 0;

  AKANTU_DEBUG_ASSERT(prescribed_dam > 0.,
		      "Your prescribed damage must be greater than zero");

  if (q_point.type != _not_defined) {
    //const Array<Real> & Sc_vect = material.getInternal("Sc")(q_point.type);

    const Array<Real> & e_stress_vect = equivalent_stress(q_point.type);


    /// check if damage occurs
    if (e_stress_vect(q_point.global_num) >= 1.) { //Sc_vect(q_point.global_num)) {

      /// update internal field damage and increment # if damaged elements
      Array<Real> & dam_vect = this->damage(q_point.type);

      if (dam_vect(q_point.global_num) < dam_threshold)
	dam_vect(q_point.global_num) += prescribed_dam;
      else dam_vect(q_point.global_num) = 1.;

      nb_damaged_elements += 1;

    }
  }

  StaticCommunicator & comm = akantu::StaticCommunicator::getStaticCommunicator();
  comm.allReduce(&nb_damaged_elements, 1, _so_sum);

  return nb_damaged_elements;

}
/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialDamageIterative<spatial_dimension>::updateEnergiesAfterDamage(ElementType el_type, GhostType ghost_type) {
  MaterialDamage<spatial_dimension>::updateEnergies(el_type, ghost_type);
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialDamageIterative);


__END_AKANTU__
