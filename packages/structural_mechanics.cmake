#===============================================================================
# @file   CMakeLists.txt
# @author Richart Nicolas <nicolas.richart@epfl.ch>
# @date   Fri Sep 29 16:46:30 2010 
#
# @brief package description for structural mechanics
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
option(AKANTU_STRUCTURAL_MECHANICS "Use Structural mechanics model package of Akantu" OFF)
set(AKANTU_STRUCTURAL_MECHANICS_FILES
  model/structural_mechanics/structural_mechanics_model.cc
  model/structural_mechanics/structural_mechanics_model_boundary.cc
  model/structural_mechanics/structural_mechanics_model_inline_impl.cc
  )

set(AKANTU_STRUCTURAL_MECHANICS_DOC
  manual/manual-structuralmechanicsmodel.tex
  )