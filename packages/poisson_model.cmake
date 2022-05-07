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
#
# @section LICENSE
#
# Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


package_declare(poisson_model
  DESCRIPTION "Use Poisson Model package of Akantu")

package_declare_sources(poisson_model
  model/heat_transfer/poisson_model.cc
  model/heat_transfer/poisson_model.hh
  model/heat_transfer/poisson_model_inline_impl.hh

  
  model/heat_transfer/constitutive_law.cc
  model/heat_transfer/constitutive_law.hh
  model/heat_transfer/constitutive_law_selector.hh
  model/heat_transfer/constitutive_law_selector_tmpl.hh

  model/heat_transfer/constitutive_laws/constitutive_law_diffusion.hh
  model/heat_transfer/constitutive_laws/constitutive_law_diffusion.cc
  model/heat_transfer/constitutive_laws/constitutive_law_heat.hh
  model/heat_transfer/constitutive_laws/constitutive_law_heat.cc
  )
