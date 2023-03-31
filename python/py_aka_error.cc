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
#include "py_aka_error.hh"
/* -------------------------------------------------------------------------- */
#include <aka_error.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */

void register_error(py::module & mod) {

  py::module mod_debug = mod.def_submodule("debug");

  mod.def("setDebugLevel", [](DebugLevel lvl) {
    debug::setDebugLevel(lvl);
    PyErr_WarnEx(PyExc_DeprecationWarning,
                 "setDebugLevel() is deprecated, it has moved in the "
                 "submodule debug",
                 1);
  });
  mod.def("getDebugLevel", []() {
    PyErr_WarnEx(PyExc_DeprecationWarning,
                 "getDebugLevel() is deprecated, it has moved in the "
                 "submodule debug",
                 1);
    return debug::getDebugLevel();
  });

  mod.def("printBacktrace", [](bool flag) {
    debug::debugger.printBacktrace(flag);
    PyErr_WarnEx(PyExc_DeprecationWarning,
                 "printBacktrace() is deprecated, it has moved in the "
                 "submodule debug",
                 1);
  });

  mod_debug.def("setDebugLevel", &debug::setDebugLevel);
  mod_debug.def("getDebugLevel", &debug::getDebugLevel);
  mod_debug.def("printBacktrace",
                [](bool flag) { debug::debugger.printBacktrace(flag); });

  py::enum_<DebugLevel>(mod, "DebugLevel")
      .value("dblError", dblError)
      .value("dblException", dblException)
      .value("dblCritical", dblCritical)
      .value("dblMajor", dblMajor)
      .value("dblWarning", dblWarning)
      .value("dblInfo", dblInfo)
      .value("dblTrace", dblTrace)
      .value("dblAccessory", dblAccessory)
      .value("dblDebug", dblDebug)
      .value("dblDump", dblDump)
      .value("dblTest", dblTest)
      .export_values();
}

} // namespace akantu
