/**
 * @file   py_heat_transfer_model.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Jun 16 2019
 * @date last modification: Sun Jun 16 2019
 *
 * @brief  pybind11 interface to HeatTransferModel
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <heat_transfer_model.hh>
#include <non_linear_solver.hh>
/* -------------------------------------------------------------------------- */
//#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
#define def_deprecated(func_name, mesg)                                        \
  def(func_name, [](py::args, py::kwargs) { AKANTU_ERROR(mesg); })

#define def_function_nocopy(func_name)                                         \
  def(                                                                         \
      #func_name,                                                              \
      [](HeatTransferModel & self) -> decltype(auto) {                         \
        return self.func_name();                                               \
      },                                                                       \
      py::return_value_policy::reference)

#define def_function(func_name)                                                \
  def(#func_name, [](HeatTransferModel & self) -> decltype(auto) {             \
    return self.func_name();                                                   \
  })
/* -------------------------------------------------------------------------- */

void register_heat_transfer_model(py::module & mod) {
  py::class_<HeatTransferModelOptions>(mod, "HeatTransferModelOptions")
      .def(py::init<AnalysisMethod>(),
           py::arg("analysis_method") = _explicit_lumped_mass);

  py::class_<HeatTransferModel, Model>(mod, "HeatTransferModel",
                                       py::multiple_inheritance())
      .def(py::init<Mesh &, UInt, const ID &>(), py::arg("mesh"),
           py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "heat_transfer_model")
      .def(
          "initFull",
          [](HeatTransferModel & self,
             const HeatTransferModelOptions & options) {
            self.initFull(options);
          },
          py::arg("_analysis_method") = HeatTransferModelOptions())
      .def(
          "initFull",
          [](HeatTransferModel & self,
             const AnalysisMethod & _analysis_method) {
            self.initFull(HeatTransferModelOptions(_analysis_method));
          },
          py::arg("_analysis_method"))
      .def("setTimeStep", &HeatTransferModel::setTimeStep, py::arg("time_step"),
           py::arg("solver_id") = "")
      .def_function(getStableTimeStep)
      .def_function_nocopy(getTemperature)
      .def_function_nocopy(getBlockedDOFs)
      .def("getTemperatureGradient", &HeatTransferModel::getTemperatureGradient,
           py::arg("el_type"), py::arg("ghost_type") = _not_ghost,
           py::return_value_policy::reference)
      .def("getKgradT", &HeatTransferModel::getKgradT, py::arg("el_type"),
           py::arg("ghost_type") = _not_ghost,
           py::return_value_policy::reference);
}

} // namespace akantu
