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
  MaterialElastic<spatial_dimension>(model, id),
  MaterialMazars<spatial_dimension>(model, id),
  MaterialNonLocalParent(model, id),
  Ehat("Ehat", id) {
  AKANTU_DEBUG_IN();
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

  Real * dam = this->damage(el_type, ghost_type).storage();
  Real * Ehatt = this->Ehat(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;

  MaterialMazars<spatial_dimension>::computeStressOnQuad(grad_u, sigma,
							 *dam, *Ehatt);
  ++dam;
  ++Ehatt;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialMazarsNonLocal<spatial_dimension>::computeNonLocalStress(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  ByElementTypeReal nl_var("Non local variable", this->id);
  this->initInternalVector(nl_var, 1);
  this->resizeInternalVector(nl_var);

#ifdef AKANTU_MAZARS_NON_LOCAL_AVERAGE_DAMAGE
  this->weightedAvergageOnNeighbours(damage, nl_var, 1);
#else
  this->weightedAvergageOnNeighbours(Ehat, nl_var, 1);
#endif

  Mesh::type_iterator it = this->model->getFEM().getMesh().firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = this->model->getFEM().getMesh().lastType(spatial_dimension, ghost_type);

  for(; it != last_type; ++it) {
    this->computeNonLocalStress(nl_var(*it, ghost_type), *it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialMazarsNonLocal<spatial_dimension>::computeNonLocalStress(Vector<Real> & non_loc_var,
								      ElementType el_type,
								      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

#ifdef AKANTU_MAZARS_NON_LOCAL_AVERAGE_DAMAGE
  Real * Ehatt = this->Ehat(el_type, ghost_type).storage();
#else
  Real * dam   = this->damage(el_type, ghost_type).storage();
#endif
  Real * nl_var = non_loc_var.storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;

#ifdef AKANTU_MAZARS_NON_LOCAL_AVERAGE_DAMAGE
  this->computeDamageAndStressOnQuad(grad_u, sigma, *nl_var, *Ehatt);
  ++Ehatt;
#else
  this->computeDamageAndStressOnQuad(grad_u, sigma, *dam, *nl_var);
  ++dam;
#endif
  ++nl_var;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  this->updateDissipatedEnergy(ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
bool MaterialMazarsNonLocal<spatial_dimension>::setParam(const std::string & key,
							 const std::string & value,
							 const ID & id) {
  return MaterialNonLocalParent::setParam(key, value, id) ||
    MaterialMazars<spatial_dimension>::setParam(key, value, id);
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialMazarsNonLocal<spatial_dimension>::printself(std::ostream & stream,
							  int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material<_mazars_non_local> [" << std::endl;
  MaterialMazars<spatial_dimension>::printself(stream, indent + 1);
  MaterialNonLocalParent::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialMazarsNonLocal);


__END_AKANTU__
