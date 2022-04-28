/**
 * Copyright (©) 2022 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
 */

/* -------------------------------------------------------------------------- */
#include "py_integration_scheme.hh"
/* -------------------------------------------------------------------------- */
#include <dof_manager.hh>
#include <integration_scheme.hh>
#include <integration_scheme_1st_order.hh>
#include <integration_scheme_2nd_order.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;

namespace akantu {

/* -------------------------------------------------------------------------- */
void register_integration_schemes(py::module & mod) {
  auto integration_scheme_class =
      py::class_<IntegrationScheme>(mod, "IntegrationScheme");

  py::enum_<IntegrationScheme::SolutionType>(integration_scheme_class,
                                             "SolutionType")
      .value("_not_defined", IntegrationScheme::SolutionType::_not_defined)
      .value("_displacement", IntegrationScheme::SolutionType::_displacement)
      .value("_temperature", IntegrationScheme::SolutionType::_temperature)
      .value("_damage", IntegrationScheme::SolutionType::_damage)
      .value("_velocity", IntegrationScheme::SolutionType::_velocity)
      .value("_temperature_rate",
             IntegrationScheme::SolutionType::_temperature_rate)
      .value("_acceleration", IntegrationScheme::SolutionType::_acceleration)
      .export_values();

  py::class_<IntegrationScheme2ndOrder, IntegrationScheme>(
      mod, "IntegrationScheme2ndOrder");

  py::class_<NewmarkBeta, IntegrationScheme2ndOrder>(mod, "NewmarkBeta")
      .def(py::init<DOFManager &, const ID &, Real, Real>(),
           py::arg("dof_manager"), py::arg("id"), py::arg("alpha"),
           py::arg("beta"));

  py::class_<CentralDifference, NewmarkBeta>(mod, "CentralDifference");
  py::class_<TrapezoidalRule2, NewmarkBeta>(mod, "TrapezoidalRule2");
  py::class_<FoxGoodwin, NewmarkBeta>(mod, "FoxGoodwin");
  py::class_<LinearAceleration, NewmarkBeta>(mod, "LinearAceleration");
}
} // namespace akantu
