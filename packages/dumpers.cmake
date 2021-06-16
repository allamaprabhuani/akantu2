#===============================================================================
# @file   dumpers.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Tue Jan 16 2018
# @date last modification: Fri Mar 16 2018
#
# @brief  IO dumpers package description
#
#
# @section LICENSE
#
# Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


package_declare(dumpers
  DEFAULT ON
  DESCRIPTION "Dumpers for Akantu"
  DEPENDS INTERFACE IOHelper)

package_declare_sources(dumpers
  io/dumper/dumpable_iohelper.hh
  io/dumper/dumper_compute.hh
  io/dumper/dumper_element_iterator.hh
  io/dumper/dumper_elemental_field.hh
  io/dumper/dumper_generic_elemental_field.hh
  io/dumper/dumper_generic_elemental_field_tmpl.hh
  io/dumper/dumper_homogenizing_field.hh
  io/dumper/dumper_internal_material_field.hh
  io/dumper/dumper_iohelper.cc
  io/dumper/dumper_iohelper.hh
  io/dumper/dumper_iohelper_paraview.cc
  io/dumper/dumper_iohelper_paraview.hh
  io/dumper/dumper_nodal_field.hh
  io/dumper/dumper_padding_helper.hh
  io/dumper/dumper_quadrature_point_iterator.hh
  io/dumper/dumper_text.cc
  io/dumper/dumper_text.hh
  io/dumper/dumper_type_traits.hh
  io/dumper/dumper_variable.hh
  )

package_declare_documentation(dumpers
  "This package activates the IOHelper facilities withing Akantu. This is mandatory if you want to be able to output Paraview files"
  "as well as any Dumper within Akantu."
  )

package_declare_extra_files_to_package(dumpers
  PROJECT
    third-party/cmake/iohelper.cmake
    cmake/Modules/FindIOHelper.cmake
  )
