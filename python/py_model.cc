/**
 * @file   py_model.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @author Philip Mueller <philip.paul.mueller@bluemail.ch>
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Jun 16 2019
 * @date last modification: Sat Mar 13 2021
 *
 * @brief  pybind11 interface to Model and parent classes
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
#include <model.hh>
#include <non_linear_solver.hh>
#include <solver_callback.hh>
#include <sparse_matrix_aij.hh>
#include <time_step_solver.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void register_model(py::module & mod) {
  py::class_<ModelSolver, SolverCallback, Parsable>(mod, "ModelSolver",
                                                    py::multiple_inheritance())
      .def(
          "getNonLinearSolver",
          [](ModelSolver & self, const ID & solver_id) -> NonLinearSolver & {
            return self.getNonLinearSolver(solver_id);
          },
          py::arg("solver_id") = "", py::return_value_policy::reference)
      .def(
          "getTimeStepSolver",
          [](ModelSolver & self, const ID & solver_id) -> TimeStepSolver & {
            return self.getTimeStepSolver(solver_id);
          },
          py::arg("solver_id") = "", py::return_value_policy::reference)
      .def(
          "solveStep",
          [](ModelSolver & self, const ID & solver_id) {
            self.solveStep(solver_id);
          },
          py::arg("solver_id") = "")
      .def(
          "solveStep",
          [](ModelSolver & self, SolverCallback & callback,
             const ID & solver_id) { self.solveStep(callback, solver_id); },
          py::arg("callback"), py::arg("solver_id") = "");

  py::class_<Model, ModelSolver>(mod, "Model", py::multiple_inheritance())
      .def("setBaseName", &Model::setBaseName)
      .def("setDirectory", &Model::setDirectory)
      .def("getFEEngine", &Model::getFEEngine, py::arg("name") = "",
           py::return_value_policy::reference)
      .def("getFEEngineBoundary", &Model::getFEEngine, py::arg("name") = "",
           py::return_value_policy::reference)
      .def("addDumpFieldVector", &Model::addDumpFieldVector)
      .def("addDumpField", &Model::addDumpField)
      .def("setBaseNameToDumper", &Model::setBaseNameToDumper)
      .def("addDumpFieldVectorToDumper", &Model::addDumpFieldVectorToDumper)
      .def("addDumpFieldToDumper", &Model::addDumpFieldToDumper)
      .def("dump", [](Model & self) { self.dump(); })
      .def(
          "dump", [](Model & self, UInt step) { self.dump(step); },
          py::arg("step"))
      .def(
          "dump",
          [](Model & self, Real time, UInt step) { self.dump(time, step); },
          py::arg("time"), py::arg("step"))
      .def(
          "dump",
          [](Model & self, const std::string & dumper) { self.dump(dumper); },
          py::arg("dumper_name"))
      .def(
          "dump",
          [](Model & self, const std::string & dumper, UInt step) {
            self.dump(dumper, step);
          },
          py::arg("dumper_name"), py::arg("step"))
      .def(
          "dump",
          [](Model & self, const std::string & dumper, Real time, UInt step) {
            self.dump(dumper, time, step);
          },
          py::arg("dumper_name"), py::arg("time"), py::arg("step"))
      .def("initNewSolver", &Model::initNewSolver)
      .def(
          "getNewSolver",
          [](Model & self, const std::string id,
             const TimeStepSolverType & time,
             const NonLinearSolverType & type) {
            self.getNewSolver(id, time, type);
          },
          py::return_value_policy::reference)
      .def(
          "setIntegrationScheme",
          [](Model & self, const std::string id, const std::string primal,
             const IntegrationSchemeType & scheme_type,
             IntegrationScheme::SolutionType solution_type) {
            self.setIntegrationScheme(id, primal, scheme_type, solution_type);
          },
          py::arg("id"), py::arg("primal"), py::arg("scheme_type"),
          py::arg("solution_type") =
              IntegrationScheme::SolutionType::_not_defined)
      // .def("setIntegrationScheme",
      //      [](Model & self, const std::string id, const std::string primal,
      //         std::unique_ptr<IntegrationScheme> & scheme,
      //         IntegrationScheme::SolutionType solution_type) {
      //        self.setIntegrationScheme(id, primal, scheme, solution_type);
      //      })

      .def("getDOFManager", &Model::getDOFManager,
           py::return_value_policy::reference)
      .def("assembleMatrix", &Model::assembleMatrix);
}

} // namespace akantu
