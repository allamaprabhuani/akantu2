#===============================================================================
# Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


add_akantu_test(test_solid_mechanics_model_igfem "Test solid mechanics model igfem")
add_akantu_test(test_geometry "Test the IGFEM intersection algorithms")
add_akantu_test(test_fe_engine "Test the IGFEM interpolation and integration")
add_akantu_test(patch_tests "Patch tests for IGFEM model")