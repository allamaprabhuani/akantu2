/**
 * @file   material_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 11:57:43 2010
 *
 * @brief  Implementation of the inline functions of the class material
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique fédérale de Lausanne)
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
inline void Material::addElement(ElementType type, UInt element) {
  element_filter[type]->push_back(element);
  UInt nb_quadrature_points = FEM::getNbQuadraturePoints(type);
  for (UInt i = 0; i < nb_quadrature_points; ++i) {
    strain[type]->push_back(REAL_INIT_VALUE);
    stress[type]->push_back(REAL_INIT_VALUE);
    if(potential_energy[type] != NULL) {
      potential_energy[type]->push_back(REAL_INIT_VALUE);
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void Material::addGhostElement(ElementType type, UInt element) {
  ghost_element_filter[type]->push_back(element);
  UInt nb_quadrature_points = FEM::getNbQuadraturePoints(type);
  for (UInt i = 0; i < nb_quadrature_points; ++i) {
    ghost_strain[type]->push_back(REAL_INIT_VALUE);
    ghost_stress[type]->push_back(REAL_INIT_VALUE);
  }
  // if(potential_energy_vector)
  //   ghost_potential_energy[type]->push_back(REAL_INIT_VALUE);
}

/* -------------------------------------------------------------------------- */
inline UInt Material::getNbDataToPack(__attribute__ ((unused)) const Element & element,
				      __attribute__ ((unused)) GhostSynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
  return 0;
}

/* -------------------------------------------------------------------------- */
inline UInt Material::getNbDataToUnpack(__attribute__ ((unused)) const Element & element,
					__attribute__ ((unused)) GhostSynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
  return 0;
}

/* -------------------------------------------------------------------------- */
inline void Material::packData(__attribute__ ((unused)) Real ** buffer,
			       __attribute__ ((unused)) const Element & element,
			       __attribute__ ((unused)) GhostSynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void Material::unpackData(__attribute__ ((unused)) Real ** buffer,
				 __attribute__ ((unused)) const Element & element,
				 __attribute__ ((unused)) GhostSynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}
