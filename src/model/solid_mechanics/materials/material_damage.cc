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
  damage("damage", id) {
  AKANTU_DEBUG_IN();

  is_non_local = false;
  initInternalVector(this->damage, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialDamage::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialElastic::initMaterial();

  resizeInternalVector(this->damage);

  is_init = true;

  AKANTU_DEBUG_OUT();
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
