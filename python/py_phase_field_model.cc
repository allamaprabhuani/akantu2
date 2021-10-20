/**
 * @file   py_phase_field_model.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Sun Jun 16 2019
 * @date last modification: Fri Jun 25 2021
 *
 * @brief  Phase field python binding
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
#include <coupler_solid_phasefield.hh>
#include <non_linear_solver.hh>
#include <phase_field_model.hh>
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
      [](PhaseFieldModel & self) -> decltype(auto) {                           \
        return self.func_name();                                               \
      },                                                                       \
      py::return_value_policy::reference)

#define def_function(func_name)                                                \
  def(#func_name, [](PhaseFieldModel & self) -> decltype(auto) {               \
    return self.func_name();                                                   \
  })
/* -------------------------------------------------------------------------- */

[[gnu::visibility("default")]] void
register_phase_field_model(py::module & mod) {

  py::class_<PhaseFieldModelOptions>(mod, "PhaseFieldModelOptions")
      .def(py::init<AnalysisMethod>(), py::arg("analysis_method") = _static);

  py::class_<PhaseFieldModel, Model>(mod, "PhaseFieldModel",
                                     py::multiple_inheritance())
      .def(py::init<Mesh &, UInt, const ID &, const ModelType>(),
           py::arg("mesh"), py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "phase_field_model",
           py::arg("model_type") = ModelType::_phase_field_model)
      .def(
          "initFull",
          [](PhaseFieldModel & self, const PhaseFieldModelOptions & options) {
            self.initFull(options);
          },
          py::arg("_analysis_method") = PhaseFieldModelOptions())
      .def(
          "initFull",
          [](PhaseFieldModel & self, const AnalysisMethod & analysis_method) {
            self.initFull(_analysis_method = analysis_method);
          },
          py::arg("_analysis_method"))
      .def_deprecated("applyDirichletBC", "Deprecated: use applyBC")
      .def("applyBC",
           [](PhaseFieldModel & self, BC::Dirichlet::DirichletFunctor & func,
              const std::string & element_group) {
             self.applyBC(func, element_group);
           })
      .def("applyBC",
           [](PhaseFieldModel & self, BC::Neumann::NeumannFunctor & func,
              const std::string & element_group) {
             self.applyBC(func, element_group);
           })
      .def("setTimeStep", &PhaseFieldModel::setTimeStep, py::arg("time_step"),
           py::arg("solver_id") = "")
      .def_function(assembleStiffnessMatrix)
      .def_function(assembleInternalForces)
      .def_function_nocopy(getDamage)
      .def_function_nocopy(getInternalForce)
      .def_function_nocopy(getBlockedDOFs)
      .def_function_nocopy(getMesh)
      .def(
          "getPhaseField",
          [](PhaseFieldModel & self, UInt phase_field_id) -> decltype(auto) {
            return self.getPhaseField(phase_field_id);
          },
          py::arg("phase_field_id"), py::return_value_policy::reference)
      .def(
          "getPhaseField",
          [](PhaseFieldModel & self,
             const ID & phase_field_name) -> decltype(auto) {
            return self.getPhaseField(phase_field_name);
          },
          py::arg("phase_field_name"), py::return_value_policy::reference)
      .def("getPhaseFieldIndex", &PhaseFieldModel::getPhaseFieldIndex)
      .def("setPhaseFieldSelector", &PhaseFieldModel::setPhaseFieldSelector);
}

[[gnu::visibility("default")]] void
register_phase_field_coupler(py::module & mod) {

  py::class_<CouplerSolidPhaseField, Model>(mod, "CouplerSolidPhaseField")
      .def(py::init<Mesh &, UInt, const ID &, const ModelType>(),
           py::arg("mesh"), py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "coupler_solid_phasefield",
           py::arg("model_type") = ModelType::_coupler_solid_phasefield)
      .def("solve",
           [](CouplerSolidPhaseField & self, const ID & solid_solver_id,
              const ID & phase_solver_id) {
             self.solve(solid_solver_id, phase_solver_id);
           })
      .def("getSolidMechanicsModel",
           &CouplerSolidPhaseField::getSolidMechanicsModel,
           py::return_value_policy::reference)
      .def("getPhaseFieldModel", &CouplerSolidPhaseField::getPhaseFieldModel,
           py::return_value_policy::reference);
}

} // namespace akantu
