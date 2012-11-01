#===============================================================================
# @file   CMakeLists.txt
#
# @author David Simon Kammer <david.kammer@epfl.ch>
#
# @date   Wed Oct 27 13:54:26 2010
#
# @brief  configuration for contact tests
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
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================

add_akantu_test(test_regular_grid_neighbor_structure "Regular grid neighbor structure tests" PACKAGE contact)
add_akantu_test(test_contact_2d_neighbor_structure "Contact 2D neighbor structure tests"     PACKAGE contact)