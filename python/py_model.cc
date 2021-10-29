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
#include <sparse_matrix_aij.hh>
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
  py::class_<DOFManager>(mod, "DOFManager")
      .def("getMatrix", &DOFManager::getMatrix,
           py::return_value_policy::reference)
      .def(
          "getNewMatrix",
          [](DOFManager & self, const std::string & name,
             const std::string & matrix_to_copy_id) -> decltype(auto) {
            return self.getNewMatrix(name, matrix_to_copy_id);
          },
          py::return_value_policy::reference)
      .def(
          "getResidual",
          [](DOFManager & self) -> decltype(auto) {
            return self.getResidual();
          },
          py::return_value_policy::reference)
      .def("getArrayPerDOFs", &DOFManager::getArrayPerDOFs)
      .def(
          "hasMatrix",
          [](DOFManager & self, const ID & name) -> bool {
            return self.hasMatrix(name);
          },
          py::arg("name"))
      .def("assembleToResidual", &DOFManager::assembleToResidual);

  py::class_<NonLinearSolver>(mod, "NonLinearSolver")
      .def(
          "set",
          [](NonLinearSolver & self, const std::string & id, const Real & val) {
            if (id == "max_iterations") {
              self.set(id, int(val));
            } else {
              self.set(id, val);
            }
          })
      .def("set",
           [](NonLinearSolver & self, const std::string & id,
              const SolveConvergenceCriteria & val) { self.set(id, val); });

  py::class_<ModelSolver, Parsable>(mod, "ModelSolver",
                                    py::multiple_inheritance())
      .def("getNonLinearSolver",
           (NonLinearSolver & (ModelSolver::*)(const ID &)) &
               ModelSolver::getNonLinearSolver,
           py::arg("solver_id") = "", py::return_value_policy::reference)
      .def("solveStep", [](ModelSolver & self) { self.solveStep(); })
      .def("solveStep", [](ModelSolver & self, const ID & solver_id) {
        self.solveStep(solver_id);
      });

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
      .def("setIntegrationScheme",
           [](Model & self, const std::string id, const std::string primal,
              const IntegrationSchemeType & scheme) {
             self.setIntegrationScheme(id, primal, scheme);
           })
      .def("getDOFManager", &Model::getDOFManager,
           py::return_value_policy::reference)
      .def("assembleMatrix", &Model::assembleMatrix);
}

} // namespace akantu
