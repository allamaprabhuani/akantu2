#===============================================================================
# @file   contact.cmake
#
# @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
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

option(AKANTU_UNFINISHED "Use to continue work in progress within Akantu" OFF)

set(AKANTU_UNFINISHED_FILES
  #cc files
  #  analysis/analysis.cc
  io/aka_abaqus_parser.hh
  mesh_utils/mesh_io/mesh_io_abaqus.cc

  # include files
  #analysis/analysis.hh
  mesh_utils/mesh_io/mesh_io_abaqus.hh
)
