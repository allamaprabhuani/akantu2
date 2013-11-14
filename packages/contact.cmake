#===============================================================================
# @file   contact.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Mon Nov 21 18:19:15 2011
#
# @brief  package description for contact
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

option(AKANTU_CONTACT "Use Contact package of Akantu" OFF)

set(AKANTU_CONTACT_FILES
  #cc files
  contact/discretization.cc
  contact/element.cc
  contact/friction.cc
  contact/resolution.cc
  contact/scheme.cc
  contact/search.cc
  contact/surface.cc
  contact/zone.cc
  model/model_manager.cc

  # include files

  contact/contact_common.hh
  contact/contact_manager.hh
  contact/discretization.hh
  contact/element.hh
  contact/friction.hh
  contact/resolution.hh
  contact/scheme.hh
  contact/search.hh
  contact/surface.hh
  contact/zone.hh
  model/model_manager.hh
  )

add_external_package_dependencies(contact cblas)
add_internal_package_dependencies(contact cpparray)

if(AKANTU_CONTACT)
  list(APPEND AKANTU_BOOST_COMPONENTS
    chrono
    system
    )
endif()

#add_internal_package_dependencies(contact optimization)
