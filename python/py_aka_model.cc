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

[[gnu::visibility("default")]] void register_model(py::module & mod) {

  py::class_<SparseMatrix>(mod, "SparseMatrix")
      .def("getMatrixType", &SparseMatrix::getMatrixType)
      .def("size", &SparseMatrix::size);

  py::class_<SparseMatrixAIJ, SparseMatrix>(mod, "SparseMatrixAIJ")
      .def("getIRN", &SparseMatrixAIJ::getIRN)
      .def("getJCN", &SparseMatrixAIJ::getJCN)
      .def("getA", &SparseMatrixAIJ::getA);

  py::class_<DOFManager>(mod, "DOFManager")
      .def("getMatrix",
           [](DOFManager & self, const std::string & name) {
             return dynamic_cast<akantu::SparseMatrixAIJ &>(
                 self.getMatrix(name));
           },
           py::return_value_policy::reference);

  py::class_<Model>(mod, "Model")
      .def("setBaseName", &Model::setBaseName)
      .def("getFEEngine", &Model::getFEEngine, py::arg("name") = "",
           py::return_value_policy::reference)
      .def("addDumpFieldVector", &Model::addDumpFieldVector)
      .def("addDumpField", &Model::addDumpField)
      .def("setBaseNameToDumper", &Model::setBaseNameToDumper)
      .def("addDumpFieldVectorToDumper", &Model::addDumpFieldVectorToDumper)
      .def("addDumpFieldToDumper", &Model::addDumpFieldToDumper)
      .def("dump", &Model::dump)
      .def("initNewSolver", &Model::initNewSolver)
      .def("getDOFManager", &Model::getDOFManager,
           py::return_value_policy::reference);

  py::class_<NonLinearSolver>(mod, "NonLinearSolver")
      .def(
          "set",
          [](NonLinearSolver & self, const std::string & id, const Real & val) {
            if (id == "max_iterations")
              self.set(id, int(val));
            else if (id == "convergence_type")
              self.set(id, akantu::SolveConvergenceCriteria(UInt(val)));
            else
              self.set(id, val);
          })
      .def("set", &NonLinearSolver::set<SolveConvergenceCriteria>);

  py::class_<ModelSolver>(mod, "ModelSolver")
      .def("getNonLinearSolver",
           (NonLinearSolver & (ModelSolver::*)(const ID &)) &
               ModelSolver::getNonLinearSolver,
           py::arg("solver_id") = "", py::return_value_policy::reference)
      .def("solveStep", &ModelSolver::solveStep, py::arg("solver_id") = "");
}

} // namespace akantu
