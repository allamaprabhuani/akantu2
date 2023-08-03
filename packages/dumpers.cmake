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
  DEFAULT OFF
  DESCRIPTION "Dumpers for Akantu"
  DEPENDS HDF5)

package_declare_sources(dumpers
  io/new_dumpers/dumper.hh
  io/new_dumpers/dumper.cc
  io/new_dumpers/dumper_field.hh
  io/new_dumpers/dumper_variable.hh
  io/new_dumpers/support.hh
  io/new_dumpers/support_tmpl.hh



  # io/dumper/dumpable_iohelper.hh
  # io/dumper/dumper_compute.hh
  # io/dumper/dumper_element_iterator.hh
  # io/dumper/dumper_elemental_field.hh
  # io/dumper/dumper_generic_elemental_field.hh
  # io/dumper/dumper_generic_elemental_field_tmpl.hh
  # io/dumper/dumper_homogenizing_field.hh
  # io/dumper/dumper_internal_material_field.hh
  # io/dumper/dumper_iohelper.cc
  # io/dumper/dumper_iohelper.hh
  # io/dumper/dumper_iohelper_paraview.cc
  # io/dumper/dumper_iohelper_paraview.hh
  # io/dumper/dumper_nodal_field.hh
  # io/dumper/dumper_padding_helper.hh
  # io/dumper/dumper_quadrature_point_iterator.hh
  # io/dumper/dumper_text.cc
  # io/dumper/dumper_text.hh
  # io/dumper/dumper_type_traits.hh
  # io/dumper/dumpable.hh
  # io/dumper/dumpable.cc
  # io/dumper/dumpable_dummy.hh
  # io/dumper/dumpable_inline_impl.hh
  # io/dumper/dumper_material_padders.hh
  # io/dumper/dumper_filtered_connectivity.hh
  # io/dumper/dumper_element_partition.hh

  )
