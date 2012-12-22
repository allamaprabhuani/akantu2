#===============================================================================
# @file   iohelper.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Tue Nov 29 15:16:35 2011
#
# @brief  package description for iohelper
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

option(AKANTU_USE_IOHELPER "Add IOHelper support in akantu" ON)
mark_as_advanced(AKANTU_USE_IOHELPER)

if(AKANTU_USE_IOHELPER)
  set(IOHELPER_TARGETS_EXPORT ${AKANTU_TARGETS_EXPORT})
  add_subdirectory(third-party/iohelper)

  list(APPEND AKANTU_EXTERNAL_LIBRARIES iohelper)
  list(APPEND AKANTU_EXTERNAL_LIB_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/third-party/iohelper/src)

  set(AKANTU_IOHELPER_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/third-party/iohelper/src)

  list(APPEND AKANTU_EXPORT_LIST iohelper)
  list(APPEND AKANTU_OPTION_LIST IOHELPER)

  mark_as_advanced(IOHELPER_TESTS)
  set(AKANTU_IOHELPER ON)
else()
  set(AKANTU_IOHELPER OFF)
endif()


set(AKANTU_IOHELPER_FILES
  io/dumper/dumper_iohelper.hh
  io/dumper/dumper_iohelper.cc
  io/dumper/dumper_iohelper_tmpl.hh
  io/dumper/dumper_paraview.hh
  io/dumper/dumper_paraview.cc
  io/dumper/dumper_iohelper_tmpl_elemental_field.hh
  io/dumper/dumper_iohelper_tmpl_homogenizing_field.hh
  io/dumper/dumper_iohelper_tmpl_material_internal_field.hh
  io/dumper/dumper_iohelper_tmpl_nodal_field.hh
  io/dumper/dumper_iohelper_tmpl_quadrature_points_field.hh
  )
