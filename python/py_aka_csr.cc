/**
 * Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
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
 */

/* -------------------------------------------------------------------------- */
#include "aka_config.hh"
/* -------------------------------------------------------------------------- */
#include <aka_csr.hh>
#include <element.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

namespace {
  /* ------------------------------------------------------------------------ */
  template <typename T>
  auto register_csr(py::module & mod, const std::string & name) {
    py::class_<CSR<T>, std::shared_ptr<CSR<T>>>(mod, ("CSR" + name).c_str())
        .def(py::init<Int>(), py::arg("nb_rows") = 0)
        .def("begin", [](CSR<T> & self, Idx row) { return self.begin(row); })
        .def("end", [](CSR<T> & self, Idx row) { return self.end(row); });

    py::class_<typename CSR<T>::iterator>(mod, ("CSRIterator" + name).c_str());
  }

} // namespace
/* -------------------------------------------------------------------------- */
void register_aka_csr(py::module & mod) {
  register_csr<Real>(mod, "Real");
  register_csr<Element>(mod, "Element");
}
} // namespace akantu
