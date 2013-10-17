/**
 * @file   material_thermal.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date   Th Oct 17 11:56:37 2013
 *
 * @brief  Specialization of the material class for the thermal material
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
#include "material_thermal.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialThermal<spatial_dimension>::MaterialThermal(SolidMechanicsModel & model, const ID & id)  :
  Material(model, id), 
  delta_t("delta_t", id) {
  AKANTU_DEBUG_IN();

  this->registerParam("alpha"       ,alpha       , 0.   , _pat_parsable | _pat_modifiable, "Thermal expansion coefficient");
  this->registerParam("delta_t_init",delta_t_init, 0.   , _pat_parsable | _pat_modifiable, "Uniform temperature field (for test purposes)");

  initInternalArray(delta_t, 1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialThermal<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  updateInternalParameters();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialThermal<spatial_dimension>::updateInternalParameters() {
  const Mesh & mesh = this->model->getFEM().getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension);
  this->resizeInternalArray(delta_t);
  for (; it != last_type ; ++it) {
    Array<Real> & delta_t_vec = delta_t(*it);
    delta_t_vec.set(delta_t_init);
  }
}

INSTANSIATE_MATERIAL(MaterialThermal);

__END_AKANTU__
