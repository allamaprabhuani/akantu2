/**
 * @file   py_solid_mechanics_model_cohesive.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Jul 21 2020
 * @date last modification: Tue Sep 29 2020
 *
 * @brief  pybind11 interface to SolidMechanicsModelCohesive
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
#include <non_linear_solver.hh>
#include <solid_mechanics_model_cohesive.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
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
      [](SolidMechanicsModel & self) -> decltype(auto) {                       \
        return self.func_name();                                               \
      },                                                                       \
      py::return_value_policy::reference)

#define def_function(func_name)                                                \
  def(#func_name, [](SolidMechanicsModel & self) -> decltype(auto) {           \
    return self.func_name();                                                   \
  })

void register_solid_mechanics_model_cohesive(py::module & mod) {
  py::class_<CohesiveElementInserter>(mod, "CohesiveElementInserter")
      .def("setLimit", &CohesiveElementInserter::setLimit);

  py::class_<SolidMechanicsModelCohesiveOptions, SolidMechanicsModelOptions>(
      mod, "SolidMechanicsModelCohesiveOptions")
      .def(py::init<AnalysisMethod, bool>(),
           py::arg("analysis_method") = _explicit_lumped_mass,
           py::arg("is_extrinsic") = false);

  py::class_<SolidMechanicsModelCohesive, SolidMechanicsModel>(
      mod, "SolidMechanicsModelCohesive")
      .def(py::init<Mesh &, UInt, const ID &>(), py::arg("mesh"),
           py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "solid_mechanics_model")
      .def(
          "initFull",
          [](SolidMechanicsModel & self, const AnalysisMethod & analysis_method,
             bool is_extrinsic) {
            self.initFull(_analysis_method = analysis_method,
                          _is_extrinsic = is_extrinsic);
          },
          py::arg("_analysis_method"), py::arg("_is_extrinsic") = false)

      .def("checkCohesiveStress",
           &SolidMechanicsModelCohesive::checkCohesiveStress)
      .def("getElementInserter",
           &SolidMechanicsModelCohesive::getElementInserter,
           py::return_value_policy::reference)
      .def("updateAutomaticInsertion",
           &SolidMechanicsModelCohesive::updateAutomaticInsertion);
}

} // namespace akantu
