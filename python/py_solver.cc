/**
 * @file   py_solver.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 29 2020
 * @date last modification: Sat Mar 06 2021
 *
 * @brief  pybind11 interface to Solver and SparseMatrix
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
#include "py_solver.hh"
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
void register_solvers(py::module & mod) {
  py::class_<SparseMatrix>(mod, "SparseMatrix")
      .def("getMatrixType", &SparseMatrix::getMatrixType)
      .def("size", &SparseMatrix::size)
      .def("zero", &SparseMatrix::zero)
      .def("saveProfile", &SparseMatrix::saveProfile)
      .def("saveMatrix", &SparseMatrix::saveMatrix)
      .def(
          "add", [](SparseMatrix & self, UInt i, UInt j) { self.add(i, j); },
          "Add entry in the profile")
      .def(
          "add",
          [](SparseMatrix & self, UInt i, UInt j, Real value) {
            self.add(i, j, value);
          },
          "Add the value to the matrix")
      .def(
          "add",
          [](SparseMatrix & self, SparseMatrix & A, Real alpha) {
            self.add(A, alpha);
          },
          "Add a matrix to the matrix", py::arg("A"), py::arg("alpha") = 1.)
      .def("__call__", [](const SparseMatrix & self, UInt i, UInt j) {
        return self(i, j);
      });

  py::class_<SparseMatrixAIJ, SparseMatrix>(mod, "SparseMatrixAIJ")
      .def("getIRN", &SparseMatrixAIJ::getIRN)
      .def("getJCN", &SparseMatrixAIJ::getJCN)
      .def("getA", &SparseMatrixAIJ::getA);

  py::class_<SolverVector>(mod, "SolverVector");
}

} // namespace akantu
