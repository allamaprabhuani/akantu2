/**
 * @file   integration_scheme.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Oct 23 15:33:36 2015
 *
 * @brief  Common interface to all interface schemes
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "integration_scheme.hh"
#include "dof_manager.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
IntegrationScheme::IntegrationScheme(DOFManager & dof_manager,
                                     const ID & dof_id, UInt order)
    : dof_manager(dof_manager), dof_id(dof_id), order(order) {}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
// /// standard output stream operator for SolutionType
// std::ostream & operator<<(std::ostream & stream,
//                           const IntegrationScheme::SolutionType & type) {
//   switch (type) {
//   case IntegrationScheme::_displacement:
//     stream << "displacement";
//     break;
//   case IntegrationScheme::_temperature:
//     stream << "temperature";
//     break;
//   case IntegrationScheme::_velocity:
//     stream << "velocity";
//     break;
//   case IntegrationScheme::_temperature_rate:
//     stream << "temperature_rate";
//     break;
//   case IntegrationScheme::_acceleration:
//     stream << "acceleration";
//     break;
//   }
//   return stream;
// }

/* -------------------------------------------------------------------------- */
/// standard input stream operator for SolutionType
std::istream & operator>>(std::istream & stream,
                          IntegrationScheme::SolutionType & type) {
  std::string str;
  stream >> str;
  if (str == "displacement")
    type = IntegrationScheme::_displacement;
  else if (str == "temperature")
    type = IntegrationScheme::_temperature;
  else if (str == "velocity")
    type = IntegrationScheme::_velocity;
  else if (str == "temperature_rate")
    type = IntegrationScheme::_temperature_rate;
  else if (str == "acceleration")
    type = IntegrationScheme::_acceleration;
  else {
    stream.setstate(std::ios::failbit);
  }

  return stream;
}

__END_AKANTU__
