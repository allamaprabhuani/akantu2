/**
 * @file   py_aka_error.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue May 07 2019
 * @date last modification: Tue Sep 29 2020
 *
 * @brief  pybind11 interface to aka_error
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

  mod.def("setDebugLevel", &debug::setDebugLevel);
  mod.def("getDebugLevel", &debug::getDebugLevel);
  mod.def("printBacktrace",
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
