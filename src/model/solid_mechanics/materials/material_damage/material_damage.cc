/**
 * @file   material_damage.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Chambart <marion.chambart@epfl.ch>
 * @date   Tue Jul 27 11:53:52 2010
 *
 * @brief  Specialization of the material class for the damage material
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
#include "material_damage.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialDamage<spatial_dimension>::MaterialDamage(SolidMechanicsModel & model,
						  const ID & id)  :
  Material(model, id), MaterialElastic<spatial_dimension>(model, id),
  damage("damage", id),
  dissipated_energy("dissipated energy", id),
  strain_prev("previous strain", id),
  stress_prev("previous stress", id),
  int_sigma("integral of sigma", id) {
  AKANTU_DEBUG_IN();

  this->is_non_local = false;
  this->initInternalVector(this->damage, 1);
  this->initInternalVector(this->dissipated_energy, 1);
  this->initInternalVector(this->strain_prev, spatial_dimension * spatial_dimension);
  this->initInternalVector(this->stress_prev, spatial_dimension * spatial_dimension);
  this->initInternalVector(this->int_sigma, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialDamage<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialElastic<spatial_dimension>::initMaterial();

  this->resizeInternalVector(this->damage);
  this->resizeInternalVector(this->dissipated_energy);
  this->resizeInternalVector(this->strain_prev);
  this->resizeInternalVector(this->stress_prev);
  this->resizeInternalVector(this->int_sigma);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Compute the dissipated energy in  each element by a trapezoidal approximation
 * of
 * @f$ Ed = \int_0^{\epsilon}\sigma(\omega)d\omega - \frac{1}{2}\sigma:\epsilon@f$
 */
template<UInt spatial_dimension>
void MaterialDamage<spatial_dimension>::updateDissipatedEnergy(GhostType ghost_type) {
  // compute the dissipated energy per element
  const Mesh & mesh = this->model->getFEM().getMesh();
  Mesh::type_iterator it  = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, ghost_type);

  for(; it != end; ++it) {
    ElementType el_type = *it;
    Vector<Real>::iterator<types::Matrix> sigma =
      this->stress(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);
    Vector<Real>::iterator<types::Matrix> sigma_p =
      stress_prev(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);
    Vector<Real>::iterator<types::Matrix> epsilon =
      this->strain(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);
    Vector<Real>::iterator<types::Matrix> epsilon_p =
      strain_prev(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);


    Vector<Real>::iterator<Real> ints = int_sigma(el_type, ghost_type).begin();
    Vector<Real>::iterator<Real> ed   = dissipated_energy(el_type, ghost_type).begin();
    Vector<Real>::iterator<Real> ed_end  = dissipated_energy(el_type, ghost_type).end();

    for (; ed != ed_end; ++ed, ++ints, ++epsilon, ++sigma, ++epsilon_p, ++sigma_p) {
      Real epot = 0.;
      Real dint = 0.;
      for (UInt i = 0; i < spatial_dimension; ++i) {
	for (UInt j = 0; j < spatial_dimension; ++j) {
	  epot += (*sigma)(i,j) * (*epsilon)(i,j); /// \f$ epot = .5 \sigma : \epsilon \f$
	  dint += .5 * ((*sigma_p)(i,j) + (*sigma)(i,j)) * ((*epsilon)(i,j) - (*epsilon_p)(i,j)); /// \f$ \frac{.5 \sigma(\epsilon(t-h)) + \sigma(\epsilon(t))}{\epsilon(t) - \epsilon(t-h)} \f$

	  (*epsilon_p)(i,j) = (*epsilon)(i,j);
	  (*sigma_p)(i,j) = (*sigma)(i,j);
	}
      }

      epot *= .5;

      *ints += dint;
      *ed = *ints - epot;
    }
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialDamage<spatial_dimension>::computeAllStresses(GhostType ghost_type) {
  Material::computeAllStresses(ghost_type);
  if(!this->is_non_local) this->updateDissipatedEnergy(ghost_type);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialDamage<spatial_dimension>::getDissipatedEnergy() const {
  AKANTU_DEBUG_IN();

  Real de = 0.;
  const Mesh & mesh = this->model->getFEM().getMesh();

  /// integrate the dissipated energy for each type of elements
  Mesh::type_iterator it  = mesh.firstType(spatial_dimension, _not_ghost);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, _not_ghost);

  for(; it != end; ++it) {
    de += this->model->getFEM().integrate(dissipated_energy(*it, _not_ghost), *it,
					  _not_ghost, &this->element_filter(*it, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return de;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialDamage<spatial_dimension>::getEnergy(std::string type) {
  if(type == "dissipated") return getDissipatedEnergy();
  else return Material::getEnergy(type);
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialDamage);


__END_AKANTU__
