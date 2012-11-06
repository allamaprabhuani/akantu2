/**
 * @file   velocity_weakening_coulomb_inline_impl.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date   Wed May 11 15:14:03 2011
 *
 * @brief  implementation of velocity weakening constant friction coefficient
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

/* -------------------------------------------------------------------------- */
inline Real VelocityWeakeningCoulomb::computeFricCoef(UInt impactor_node_index) {
  AKANTU_DEBUG_IN();
  
  Real friction_coefficient = 0.;

  // find impactor_information for given master
  const ContactRigid::SurfaceToImpactInfoMap & imp_info = this->contact.getImpactorsInformation();
  ContactRigid::SurfaceToImpactInfoMap::const_iterator it;

  it = imp_info.find(this->master_surface);
  AKANTU_DEBUG_ASSERT(it != imp_info.end(), 
		      "Couldn't find impactor information object for master surface " << master_surface);
  ContactRigid::ImpactorInformationPerMaster * impactor_info = it->second;

  bool * is_sticking = impactor_info->node_is_sticking->values;

  if(is_sticking[impactor_node_index*2])
    friction_coefficient = this->static_friction_coefficient;
  else
    friction_coefficient = this->dynamic_friction_coefficient;

  AKANTU_DEBUG_OUT();
  return friction_coefficient;
}
