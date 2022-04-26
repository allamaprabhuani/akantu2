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
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;

namespace akantu {

/* -------------------------------------------------------------------------- */
void register_integration_schemes(py::module & mod) {
  py::class_<IntegrationScheme>(mod, "IntegrationScheme")
      .def(py::init<DOFManager, const ID &, UInt>());

  py::enum_<IntegrationScheme::SolutionType>()
      .def("_not_defined", _not_defined)
      .def("_displacement", _displacement)
      .def("_temperature", _temperature)
      .def("_damage", _damage)
      .def("_velocity", _velocity)
      .def("_temperature_rate", _temperature_rate)
      .def("_acceleration", _acceleration);

  py::class_<IntegrationScheme2ndOrder, IntegrationScheme>(
      mod, "IntegrationScheme2ndOrder")
      .def(py::init<DOFManager, const ID &>());

  py::class_<NewmarkBeta, IntegrationScheme2ndOrder>(mod, "NewmarkBeta")
      .def(py::init<DOFManager, const ID &, UInt, Real, Real>(),
           py::arg("dof_manager"), py::arg("id"), py::arg("alpha"),
           py::arg("beta"));

  py::class_<CentralDifference, NewmarkBeta>(mod, "CentralDifference");
  py::class_<TrapezoidalRule2, NewmarkBeta>(mod, "TrapezoidalRule2");
  py::class_<FoxGoodwin, NewmarkBeta>(mod, "FoxGoodwin");
  py::class_<LinearAceleration, NewmarkBeta>(mod, "LinearAceleration");
}
} // namespace akantu
