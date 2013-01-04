#===============================================================================
# @file   contact.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
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

option(AKANTU_DEAD_CONTACT "Use Contact package of Akantu" OFF)

set(AKANTU_DEAD_CONTACT_FILES
  #cc files
  model/solid_mechanics/contact.cc
  model/solid_mechanics/contact_search.cc
  model/solid_mechanics/contact_neighbor_structure.cc
  model/solid_mechanics/contact/contact_2d_explicit.cc
  model/solid_mechanics/contact/contact_search_2d_explicit.cc
  model/solid_mechanics/contact/regular_grid_neighbor_structure.cc
  model/solid_mechanics/contact/contact_search_explicit.cc
  model/solid_mechanics/contact/contact_3d_explicit.cc
  model/solid_mechanics/contact/grid_2d_neighbor_structure.cc
  model/solid_mechanics/contact/contact_rigid.cc
  model/solid_mechanics/contact/friction_coefficient.cc
  model/solid_mechanics/contact/friction_coefficient/unique_constant_fric_coef.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_weakening_coulomb.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef/velocity_weakening_exponential.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef/rice_kuwano.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef/rice_kuwano_modified.cc
  model/solid_mechanics/contact/friction_coefficient/simplified_dieterich_fric_coef.cc
  model/solid_mechanics/contact/friction_coefficient/simplified_dieterich_fric_coef/ruina_slowness_fric_coef.cc
  model/solid_mechanics/contact/friction_coefficient/historic_velocity_fric_coef.cc
  # include files

  model/solid_mechanics/contact.hh
  model/solid_mechanics/contact/contact_search_explicit.hh
  model/solid_mechanics/contact/friction_coefficient/unique_constant_fric_coef.hh
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef/velocity_weakening_exponential.hh
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef/rice_kuwano_modified.hh
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef/rice_kuwano.hh
  model/solid_mechanics/contact/friction_coefficient/simplified_dieterich_fric_coef/ruina_slowness_fric_coef.hh
  model/solid_mechanics/contact/friction_coefficient/velocity_weakening_coulomb.hh
  model/solid_mechanics/contact/friction_coefficient/simplified_dieterich_fric_coef.hh
  model/solid_mechanics/contact/friction_coefficient/historic_velocity_fric_coef.hh
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef.hh
  model/solid_mechanics/contact/contact_rigid.hh
  model/solid_mechanics/contact/friction_coefficient.hh
  model/solid_mechanics/contact/contact_search_2d_explicit.hh
  model/solid_mechanics/contact/contact_3d_explicit.hh
  model/solid_mechanics/contact/grid_2d_neighbor_structure.hh
  model/solid_mechanics/contact/regular_grid_neighbor_structure.hh
  model/solid_mechanics/contact/contact_2d_explicit.hh
  model/solid_mechanics/contact_search.hh
  model/solid_mechanics/contact_neighbor_structure.hh

  # inline implementation

  model/solid_mechanics/contact/friction_coefficient/unique_constant_fric_coef_inline_impl.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_weakening_coulomb_inline_impl.cc

  model/solid_mechanics/contact/friction_coefficient/unique_constant_fric_coef_inline_impl.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_weakening_coulomb_inline_impl.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef/rice_kuwano_modified_inline_impl.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef/rice_kuwano_inline_impl.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef/velocity_weakening_exponential_inline_impl.cc
  model/solid_mechanics/contact/friction_coefficient/simplified_dieterich_fric_coef/ruina_slowness_fric_coef_inline_impl.cc
  model/solid_mechanics/contact/friction_coefficient/simplified_dieterich_fric_coef_inline_impl.cc
  model/solid_mechanics/contact/contact_search_explicit_inline_impl.cc
  model/solid_mechanics/contact/contact_search_2d_explicit_inline_impl.cc
  model/solid_mechanics/contact/regular_grid_neighbor_structure_inline_impl.cc
  )



