#===============================================================================
# @file   contact_mechanics.cmake
#
# @author Mohit Pundir <mohit.pundir@epfl.ch>
#
# @date creation: Sun Oct 21 2018
# @date last modification: Sun Oct 21 2018
#
# @brief  package description for contact mechanics
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

package_declare(contact_mechanics
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

  model/contact_mechanics/resolution.hh
  model/contact_mechanics/resolution.cc
  model/contact_mechanics/resolution_utils.hh
  model/contact_mechanics/resolution_utils.cc
  model/contact_mechanics/resolutions/resolution_penalty.hh
  model/contact_mechanics/resolutions/resolution_penalty.cc

  model/contact_mechanics/surface_selector.hh
  model/contact_mechanics/surface_selector.cc)

package_declare_documentation_files(contact_mechanics
  manual-contactmechanicsmodel.tex
  manual-contact-detector.tex
  manual-contact-resolutions.tex
  )

package_declare_documentation(contact_mechanics
  "This package activates the contact mechanics model")
