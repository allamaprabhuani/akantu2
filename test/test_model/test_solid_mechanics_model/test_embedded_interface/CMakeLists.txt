#===============================================================================
# Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# This file is part of Akantu
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
# @section DESCRIPTION
#
#===============================================================================


register_test(test_embedded_element_matrix
  SOURCES test_embedded_element_matrix.cc
  FILES_TO_COPY triangle.msh embedded_element.dat
  PACKAGE embedded implicit
  )

register_test(test_embedded_interface_model
  SOURCES test_embedded_interface_model.cc
  FILES_TO_COPY embedded_mesh.msh material.dat matrix
  DIRECTORIES_TO_CREATE paraview
  PACKAGE embedded implicit
  )

register_test(test_embedded_interface_model_prestress
  SOURCES test_embedded_interface_model_prestress.cc
  FILES_TO_COPY embedded_mesh_prestress.msh embedded_mesh_prestress_reinforcement.msh prestress.dat
  PACKAGE embedded implicit
  )
