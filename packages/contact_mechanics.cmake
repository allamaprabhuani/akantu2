#===============================================================================
# @file   contact_mechanics.cmake
#
# @author Mohit Pundir <mohit.pundir@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Wed Jun 23 2021
#
# @brief  package description for contact mechanics
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


package_declare(contact_mechanics
  DEPENDS model_couplers cohesive_element
  DESCRIPTION "Use Contact Mechanics package of Akantu")

package_declare_sources(contact_mechanics
  model/contact_mechanics/contact_mechanics_model.hh
  model/contact_mechanics/contact_mechanics_model.cc
  model/contact_mechanics/contact_detector.hh
  model/contact_mechanics/contact_detector.cc
  model/contact_mechanics/contact_detector_inline_impl.cc
  model/contact_mechanics/contact_element.hh
  model/contact_mechanics/geometry_utils.hh
  model/contact_mechanics/geometry_utils.cc
  model/contact_mechanics/geometry_utils_inline_impl.cc

  model/contact_mechanics/resolution.hh
  model/contact_mechanics/resolution.cc
  model/contact_mechanics/resolution_utils.hh
  model/contact_mechanics/resolution_utils.cc
  model/contact_mechanics/resolutions/resolution_penalty.hh
  model/contact_mechanics/resolutions/resolution_penalty_tmpl.hh
  model/contact_mechanics/resolutions/resolution_penalty.cc
  model/contact_mechanics/resolutions/resolution_penalty_quadratic.hh
  model/contact_mechanics/resolutions/resolution_penalty_quadratic.cc
  
  model/contact_mechanics/surface_selector.hh
  model/contact_mechanics/surface_selector.cc

  model/model_couplers/coupler_solid_contact.hh
  model/model_couplers/coupler_solid_contact_tmpl.hh
  model/model_couplers/coupler_solid_contact.cc
  model/model_couplers/coupler_solid_cohesive_contact.hh
  model/model_couplers/coupler_solid_cohesive_contact.cc
  model/model_couplers/cohesive_contact_solvercallback.hh
  model/model_couplers/cohesive_contact_solvercallback.cc
  )

package_declare_documentation_files(contact_mechanics
  manual-contactmechanicsmodel.tex
  manual-contact-detector.tex
  )

package_declare_documentation(contact_mechanics
  "This package activates the contact mechanics model")
