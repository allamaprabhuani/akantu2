/**
 * @file   py_heat_transfer_interface_model.cc
 *
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 *
 * @date creation: Tue May 16 2023
 * @date last modification: Tue May 16 2023
 *
 * @brief  pybind11 interface to HeatTransferInterfaceModel
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
#include <heat_transfer_interface_model.hh>
#include <non_linear_solver.hh>
#include <shape_cohesive_inline_impl.hh>
/* -------------------------------------------------------------------------- */
// #include <pybind11/operators.h>
#include <pybind11/pybind11.h>
// #include <pybind11/stl.h>
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
      [](HeatTransferInterfaceModel & self) -> decltype(auto) {                \
        return self.func_name();                                               \
      },                                                                       \
      py::return_value_policy::reference)

#define def_function(func_name)                                                \
  def(#func_name, [](HeatTransferInterfaceModel & self) -> decltype(auto) {    \
    return self.func_name();                                                   \
  })
/* -------------------------------------------------------------------------- */

void register_heat_transfer_interface_model(py::module & mod) {
  py::class_<HeatTransferInterfaceModel, HeatTransferModel>(
      mod, "HeatTransferInterfaceModel")
      .def(py::init<Mesh &, UInt, const ID &>(), py::arg("mesh"),
           py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "heat_transfer_model")
      .def(
          "initFull",
          [](HeatTransferInterfaceModel & self,
             const HeatTransferModelOptions & options) {
            self.initFull(options);
          },
          py::arg("_analysis_method") = HeatTransferModelOptions())
      .def(
          "initFull",
          [](HeatTransferInterfaceModel & self,
             const AnalysisMethod & _analysis_method) {
            self.initFull(HeatTransferModelOptions(_analysis_method));
          },
          py::arg("_analysis_method"))
      .def("setTimeStep", &HeatTransferInterfaceModel::setTimeStep,
           py::arg("time_step"), py::arg("solver_id") = "")
      .def("getShapeFunctionsCohesive",
           &HeatTransferInterfaceModel::getShapeFunctionsCohesive)
      .def("getTransversalConductivityOnQpoints",
           &HeatTransferInterfaceModel::getTransversalConductivityOnQpoints,
           py::arg("el_type"), py::arg("ghost_type") = _not_ghost,
           py::return_value_policy::reference)
      .def("getLongitudinalConductivityOnQpoints",
           &HeatTransferInterfaceModel::getLongitudinalConductivityOnQpoints,
           py::arg("el_type"), py::arg("ghost_type") = _not_ghost,
           py::return_value_policy::reference)
      .def(
          "getOpening",
          [](HeatTransferInterfaceModel & self, ElementType & el_type,
             GhostType & ghost_type) { self.getOpening(el_type, ghost_type); },
          py::arg("el_type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference);
}
} // namespace akantu
