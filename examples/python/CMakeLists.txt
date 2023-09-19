#===============================================================================
# Copyright (©) 2016-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#===============================================================================


add_example(structural_mechanics_model "structural mechanics example in python"
  PACKAGE structural_mechanics)
add_example(phase_field_model "phase-field example in python"
  PACKAGE phase_field)
add_example(solid_mechanics_cohesive_model "cohesive element examples in python"
  PACKAGE cohesive_element)
add_example(solid_mechanics_model "solid_mechanics_model in python"
  PACKAGE solid_mechanics)
add_example(contact_mechanics_model "contact mechanics example in python"
  PACKAGE contact_mechanics)
package_add_files_to_package(
  examples/python/README.rst
  )
