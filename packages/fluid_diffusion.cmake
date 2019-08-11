#===============================================================================
# @file   heat_transfer.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Nov 21 2011
# @date last modification: Mon Mar 30 2015
#
# @brief  package description for heat transfer
#
# @section LICENSE
#
# Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
# Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
# Solides)
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

package_declare(fluid_diffusion
  DESCRIPTION "Use Fluid Diffusion package of Akantu")

package_declare_sources(fluid_diffusion
  model/fluid_diffusion/fluid_diffusion_model.cc
  model/fluid_diffusion/fluid_diffusion_model.hh
  model/fluid_diffusion/fluid_diffusion_model_inline_impl.cc
  )

#package_declare_documentation_files(fluid_diffusion
#  manual-heattransfermodel.tex
#  )

#package_declare_documentation(heat_transfer
#  "This package activates the heat transfer model within Akantu. "
#  "It has no additional dependencies."
#  )