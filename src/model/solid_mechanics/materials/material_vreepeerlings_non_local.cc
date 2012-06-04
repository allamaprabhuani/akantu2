/**
 * @file   material_vreepeerlings_non_local.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 * @date   Tue Jul 27 11:53:52 2010
 *
 * @brief  Specialization of the material class for the non-local Vree-Peerlings material
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
#include "material_vreepeerlings_non_local.hh"
#include "solid_mechanics_model.hh"


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialVreePeerlingsNonLocal<spatial_dimension>::MaterialVreePeerlingsNonLocal(SolidMechanicsModel & model,
							     const ID & id)  :
  Material(model, id),
  MaterialElastic<spatial_dimension>(model, id),
  MaterialVreePeerlings<spatial_dimension>(model, id),
  MaterialNonLocalParent(model, id),
  equi_strain("equi_strain", id) {
  AKANTU_DEBUG_IN();

  this->is_non_local = true;

  this->initInternalVector(this->equi_strain, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialVreePeerlingsNonLocal<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialVreePeerlings<spatial_dimension>::initMaterial();
  MaterialNonLocalParent::initMaterial();

  this->resizeInternalVector(this->equi_strain);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialVreePeerlingsNonLocal<spatial_dimension>::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam = this->damage(el_type, ghost_type).storage();
  Real * equi_straint = equi_strain(el_type, ghost_type).storage();
  Real * Kapaq = this->Kapa(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;

  MaterialVreePeerlings<spatial_dimension>::computeStressOnQuad(grad_u, sigma,
								*dam,
								*equi_straint,
								*Kapaq);
  ++dam;
  ++equi_straint;
  ++Kapaq;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialVreePeerlingsNonLocal<spatial_dimension>::computeNonLocalStress(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  ByElementTypeReal nl_var("Non local variable", this->id);
  this->initInternalVector(nl_var, 1);
  this->resizeInternalVector(nl_var);

  this->weightedAvergageOnNeighbours(equi_strain, nl_var, 1);

  Mesh::type_iterator it = this->model->getFEM().getMesh().firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = this->model->getFEM().getMesh().lastType(spatial_dimension, ghost_type);

  for(; it != last_type; ++it) {
    computeNonLocalStress(nl_var(*it, ghost_type), *it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialVreePeerlingsNonLocal<spatial_dimension>::computeNonLocalStress(Vector<Real> & non_loc_var, ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam   = this->damage(el_type, ghost_type).storage();
  Real * Kapaq = this->Kapa(el_type, ghost_type).storage();

  Real * nl_var = non_loc_var.storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;

  this->computeDamageAndStressOnQuad(sigma, *dam, *nl_var, *Kapaq);
  ++dam;
  ++Kapaq;
  ++nl_var;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  this->updateDissipatedEnergy(ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
bool MaterialVreePeerlingsNonLocal<spatial_dimension>::setParam(const std::string & key, const std::string & value,
			       const ID & id) {
  return MaterialNonLocalParent::setParam(key, value, id) ||
    MaterialVreePeerlings<spatial_dimension>::setParam(key, value, id);
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialVreePeerlingsNonLocal<spatial_dimension>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material<_vreepeerlings_non_local> [" << std::endl;
  MaterialVreePeerlings<spatial_dimension>::printself(stream, indent + 1);
  MaterialNonLocalParent::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
INSTANSIATE_MATERIAL(MaterialVreePeerlingsNonLocal);

__END_AKANTU__
