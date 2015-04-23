#===============================================================================
# @file   40_optimization.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Jan 04 2013
# @date last modification: Wed Jul 30 2014
#
# @brief  Optimization external library interface
#
# @section LICENSE
#
# Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

package_declare(optimization
  DESCRIPTION "Use optimization package in Akantu"
  DEPENDS nlopt
  ADVANCED)


package_declare_sources(optimization
  common/aka_optimize.hh
  common/aka_optimize.cc
  )


package_declare_documentation(optimization
  "This activates the optimization routines of Akantu. This is currently needed by the"
  "contact detection algorithms."
  )
