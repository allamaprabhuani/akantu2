#===============================================================================
# @file   level_set.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Mon Nov 21 18:19:15 2011
#
# @brief  package description for core
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

option(AKANTU_LEVEL_SET "Add the support for level set computations in Akantu" OFF)

set(AKANTU_LEVEL_SET_FILES
  model/level_set/level_set_model.cc
  model/level_set/sphere.cc
  model/level_set/plane.cc
  model/level_set/container.cc

  model/level_set/level_set_model.hh
  model/level_set/geometry.hh
  model/level_set/sphere.hh
  model/level_set/plane.hh
  model/level_set/container.hh

  model/level_set/level_set_model_inline_impl.cc
  )
