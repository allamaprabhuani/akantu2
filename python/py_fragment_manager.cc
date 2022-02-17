/**
 * @file   py_fragment_manager.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Mar 29 2021
 * @date last modification: Mon Mar 29 2021
 *
 * @brief  pybind11 interface to FragmentManager
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
#include <fragment_manager.hh>
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
  def(#func_name, [](FragmentManager & self) -> decltype(auto) {               \
    return self.func_name();                                                   \
  })

void register_fragment_manager(py::module & mod) {
  py::class_<FragmentManager, GroupManager>(mod, "FragmentManager")
      .def(py::init<SolidMechanicsModelCohesive &, bool, const ID &>(),
           py::arg("model"), py::arg("dump_data") = true,
           py::arg("ID") = "fragment_manager")
      .def("buildFragments", &FragmentManager::buildFragments,
           py::arg("damage_limit") = 1.)
      .def_function(computeCenterOfMass)
      .def_function(computeVelocity)
      .def_function(computeInertiaMoments)
      .def("computeAllData", &FragmentManager::computeAllData,
           py::arg("damage_limit") = 1.)
      .def_function(computeNbElementsPerFragment)
      .def_function(getNbFragment)
      .def_function(getMass)
      .def_function(getVelocity)
      .def_function(getMomentsOfInertia)
      .def_function(getPrincipalDirections)
      .def_function(getNbElementsPerFragment);
}
} // namespace akantu
