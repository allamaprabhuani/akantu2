/**
 * @file   material_mazars_non_local.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marion Chambart <marion.chambart@epfl.ch>
 * @date   Tue Jul 27 11:53:52 2010
 *
 * @brief  Specialization of the material class for the non-local mazars material
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
#include "material_mazars_non_local.hh"
#include "solid_mechanics_model.hh"


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialMazarsNonLocal<spatial_dimension>::MaterialMazarsNonLocal(SolidMechanicsModel & model,
								  const ID & id)  :
  Material(model, id),
  MaterialMazars<spatial_dimension>(model, id),
  MaterialNonLocalParent(model, id),
  Ehat("epsilon_equ", id) {
  AKANTU_DEBUG_IN();

  this->damage_in_compute_stress = false;
  this->is_non_local = true;
  this->initInternalVector(this->Ehat, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialMazarsNonLocal<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialMazars<spatial_dimension>::initMaterial();
  MaterialNonLocalParent::initMaterial();
  this->resizeInternalVector(this->Ehat);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialMazarsNonLocal<spatial_dimension>::computeStress(ElementType el_type,
							      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * damage      = this->damage(el_type, ghost_type).storage();
  Real * epsilon_equ = this->Ehat(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  MaterialMazars<spatial_dimension>::computeStressOnQuad(grad_u, sigma,
							 *damage,
							 *epsilon_equ);
  ++damage;
  ++epsilon_equ;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialMazarsNonLocal<spatial_dimension>::computeNonLocalStresses(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  ByElementTypeReal nl_var("Non local variable", this->id);
  this->initInternalVector(nl_var, 1, true);
  this->resizeInternalVector(nl_var);

  if(this->damage_in_compute_stress)
    this->weightedAvergageOnNeighbours(this->damage, nl_var, 1);
  else
    this->weightedAvergageOnNeighbours(this->Ehat, nl_var, 1);

  Mesh::type_iterator it = this->model->getFEM().getMesh().firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = this->model->getFEM().getMesh().lastType(spatial_dimension, ghost_type);
  for(; it != last_type; ++it) {
    this->computeNonLocalStress(nl_var(*it, ghost_type), *it, ghost_type);
  }

  this->updateDissipatedEnergy(ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialMazarsNonLocal<spatial_dimension>::computeNonLocalStress(Vector<Real> & non_loc_var,
								      ElementType el_type,
								      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * damage;
  Real * epsilon_equ;
  if(this->damage_in_compute_stress){
    damage      = non_loc_var.storage();
    epsilon_equ = this->Ehat(el_type, ghost_type).storage();
  } else {
    damage      = this->damage(el_type, ghost_type).storage();
    epsilon_equ = non_loc_var.storage();
  }

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
  this->computeDamageAndStressOnQuad(grad_u, sigma, *damage, *epsilon_equ);

  ++damage;
  ++epsilon_equ;
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
bool MaterialMazarsNonLocal<spatial_dimension>::parseParam(const std::string & key,
							 const std::string & value,
							 const ID & id) {
  std::stringstream sstr(value);
  if(key == "average_on_damage") { sstr >> this->damage_in_compute_stress; }
  else {
  return MaterialNonLocalParent::parseParam(key, value, id) ||
    MaterialMazars<spatial_dimension>::parseParam(key, value, id);
  }

  return true;
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialMazarsNonLocal<spatial_dimension>::printself(std::ostream & stream,
							  int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "MaterialMazarsNonLocal [" << std::endl;

  if(this->damage_in_compute_stress)
    stream << space << " + Average on damage" << std::endl;
  else 
    stream << space << " + Average on Ehat" << std::endl;

  MaterialMazars<spatial_dimension>::printself(stream, indent + 1);
  MaterialNonLocalParent::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialMazarsNonLocal);


__END_AKANTU__
