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

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialThermal<spatial_dimension>::MaterialThermal(SolidMechanicsModel & model, const ID & id)  :
  Material(model, id),
  delta_T("delta_T", *this),
  sigma_th_cur("sigma_th_cur", *this),
  sigma_th_prev("sigma_th_prev", *this),
  use_previous_stress_thermal(false) {
  AKANTU_DEBUG_IN();

  this->registerParam("E"      , E      , 0.  , _pat_parsable | _pat_modifiable, "Young's modulus"        );
  this->registerParam("nu"     , nu     , 0.5 , _pat_parsable | _pat_modifiable, "Poisson's ratio"        );
  this->registerParam("alpha"  , alpha  , 0.  , _pat_parsable | _pat_modifiable, "Thermal expansion coefficient");
  this->registerParam("delta_T", delta_T,       _pat_parsable | _pat_modifiable, "Uniform temperature field");

  delta_T.initialize(1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialThermal<dim>::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Array<Real>::iterator<> delta_t_it = this->delta_T(el_type, ghost_type).begin();
  Array<Real>::iterator<> sigma_th_cur_it = this->sigma_th_cur(el_type, ghost_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  if (dim == 1) {
    *sigma_th_cur_it = - this->E * this->alpha * *delta_t_it;
  }
  else {
    *sigma_th_cur_it = - this->E/(1-2*this->nu) * this->alpha * *delta_t_it;
  }

  ++delta_t_it;
  ++sigma_th_cur_it;
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialThermal<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  sigma_th_cur .initialize(1);

  if(use_previous_stress_thermal) {
    sigma_th_prev.initialize(1);
  }

  Material::initMaterial();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void MaterialThermal<dim>::savePreviousState(const GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Material::savePreviousState(ghost_type);
  UInt spatial_dimension = model->getSpatialDimension();
  Mesh::type_iterator it = model->getMesh().firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = model->getMesh().lastType(spatial_dimension, ghost_type);

  if(use_previous_stress_thermal) {
    for (; it != last_type; ++it) {
      sigma_th_prev(*it, ghost_type).copy(sigma_th_cur(*it, ghost_type));
    }
  }

  AKANTU_DEBUG_OUT();
}

INSTANSIATE_MATERIAL(MaterialThermal);

__END_AKANTU__
