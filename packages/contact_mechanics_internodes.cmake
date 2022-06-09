#===============================================================================
# @file   contact_mechanics.cmake
#
# @author Mohit Pundir <mohit.pundir@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Wed Jun 23 2021
#
# @brief  package description for contact mechanics
#
#
# @section LICENSE
#
# Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
# 
# Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
# 
# You should have received a copy of the GNU Lesser General Public License along
# with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================


package_declare(contact_mechanics_internodes
  DEPENDS implicit
  DESCRIPTION "Use Contact Mechanics Internodes package of Akantu")

package_declare_sources(contact_mechanics_internodes
  model/contact_mechanics/contact_mechanics_internodes_model.hh
  model/contact_mechanics/contact_mechanics_internodes_model.cc
  model/contact_mechanics/contact_detector_internodes.hh
  model/contact_mechanics/contact_detector_internodes.cc
  )
