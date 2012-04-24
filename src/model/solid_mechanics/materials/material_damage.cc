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
MaterialDamage::MaterialDamage(Model & model, const ID & id)  :
  Material(model, id), MaterialElastic(model, id),
  damage("damage", id),
  dissipated_energy("Dissipated Energy", id),
  strain_prev("Previous Strain", id),
  stress_prev("Previous Stress", id),
  int_sigma("Integral of sigma", id) {
  AKANTU_DEBUG_IN();

  is_non_local = false;
  initInternalVector(this->damage, 1);
  initInternalVector(this->dissipated_energy, 1);
  initInternalVector(this->strain_prev, spatial_dimension * spatial_dimension);
  initInternalVector(this->stress_prev, spatial_dimension * spatial_dimension);
  initInternalVector(this->int_sigma, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialDamage::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialElastic::initMaterial();

  resizeInternalVector(this->damage);
  resizeInternalVector(this->dissipated_energy);
  resizeInternalVector(this->strain_prev);
  resizeInternalVector(this->stress_prev);
  resizeInternalVector(this->int_sigma);

  is_init = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Compute the dissipated energy in  each element by a rectangular approximation
 * of
 * @f$ Ed = \int_0^{\epsilon}\sigma(\omega)d\omega - \frac{1}{2}\sigma:\epsilon@f$
 */
void MaterialDamage::updateDissipatedEnergy(GhostType ghost_type) {
  // compute the dissipated energy per element
  const Mesh & mesh = model->getFEM().getMesh();
  Mesh::type_iterator it  = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, ghost_type);

  for(; it != end; ++it) {
    ElementType el_type = *it;
    Vector<Real>::iterator<types::Matrix> sigma =
      stress(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);
    Vector<Real>::iterator<types::Matrix> sigma_p =
      stress_prev(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);
    Vector<Real>::iterator<types::Matrix> epsilon =
      strain(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);
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
	  epot += (*sigma)(i,j) * (*epsilon)(i,j);
	  dint += .5 * ((*sigma_p)(i,j) + (*sigma)(i,j)) * ((*epsilon)(i,j) - (*epsilon_p)(i,j));
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
Real MaterialDamage::getDissipatedEnergy() const {
  AKANTU_DEBUG_IN();

  Real de = 0.;
  const Mesh & mesh = model->getFEM().getMesh();

  /// integrate the dissipated energy for each type of elements
  Mesh::type_iterator it  = mesh.firstType(spatial_dimension, _not_ghost);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, _not_ghost);

  for(; it != end; ++it) {
    de += model->getFEM().integrate(dissipated_energy(*it, _not_ghost), *it,
				    _not_ghost, &element_filter(*it, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return de;
}

/* -------------------------------------------------------------------------- */
bool MaterialDamage::setParam(const std::string & key, const std::string & value,
			       const ID & id) {
  return MaterialElastic::setParam(key, value, id);
}


/* -------------------------------------------------------------------------- */
void MaterialDamage::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material<_damage> [" << std::endl;
  MaterialElastic::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */

__END_AKANTU__
