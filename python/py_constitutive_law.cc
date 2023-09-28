/**
 * @file   py_material.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jun 20 2019
 * @date last modification: Fri Apr 09 2021
 *
 * @brief  pybind11 interface to Material
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
#include "py_constitutive_law.hh"
#include "py_aka_array.hh"
#include "py_akantu_pybind11_compatibility.hh"
/* -------------------------------------------------------------------------- */
#include <constitutive_law.hh>
#include <constitutive_laws_handler.hh>
#include <parsable.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

namespace {
  /* ------------------------------------------------------------------------ */
  template <typename T>
  void register_internal_field(py::module & mod, const std::string & name) {
    py::class_<InternalField<T>, ElementTypeMapArray<T>,
               std::shared_ptr<InternalField<T>>>(
        mod, ("InternalField" + name).c_str(), py::multiple_inheritance())
        .def(
            "previous",
            [](InternalField<T> & self, ElementType type, GhostType ghost_type)
                -> Array<T> & { return self.previous(type, ghost_type); },
            py::arg("type"), py::arg("ghost_type") = _not_ghost,
            py::return_value_policy::reference)
        .def(
            "previous",
            [](InternalField<T> & self) -> InternalField<T> & {
              return self.previous();
            },
            py::return_value_policy::reference)
        .def(
            "__call__",
            [](InternalField<T> & self, ElementType type, GhostType ghost_type)
                -> Array<T> & { return self(type, ghost_type); },
            py::arg("type"), py::arg("ghost_type") = _not_ghost,
            py::return_value_policy::reference)
        .def("__repr__",
             [name](const InternalField<T> & self) {
               std::stringstream sstream;
               sstream << std::hex << &self;

               std::string str{"<InternalField" + name + " (0x" +
                               sstream.str() + "): " + self.getRegisterID() +
                               " "};
               for (auto && ghost_type : ghost_types) {
                 for (auto && type : self.elementTypes(ghost_type)) {
                   const auto & array = self(type, ghost_type);
                   if (str.back() == ']') {
                     str += ", ";
                   }
                   str += "[" + std::to_string(type) + ":" +
                          std::to_string(ghost_type) + " - Array<" +
                          debug::demangle<T>() + " " +
                          std::to_string(array.size()) + ", " +
                          std::to_string(array.getNbComponent()) + ">]";
                 }
               }
               str += ">";
               return str;
             })
        .def(
            "getID", [](InternalField<T> & self) { return self.getID(); },
            py::return_value_policy::copy)
        .def("initializeHistory", &InternalField<T>::initializeHistory)
        .def("hasHistory", &InternalField<T>::hasHistory);
  }
} // namespace

/* -------------------------------------------------------------------------- */
void register_constitutive_law_internal_handler(py::module & mod) {
  register_internal_field<Real>(mod, "Real");
  register_internal_field<Int>(mod, "Int");

  py::class_<ConstitutiveLawInternalHandler>(
      mod, "ConstitutiveLawInternalHandler", py::multiple_inheritance())
      .def(py::init<const ID &, Int, const ID &>())
      .def("registerInternalReal",
           [](ConstitutiveLawInternalHandler & self, const std::string & name,
              Int nb_component) -> decltype(auto) {
             self.registerInternal<Real>(name, nb_component);
             return self.getSharedPtrInternal<Real>(name);
           })
      .def("registerInternalInt",
           [](ConstitutiveLawInternalHandler & self, const std::string & name,
              UInt nb_component) -> decltype(auto) {
             self.template registerInternal<Int>(name, nb_component);
             return self.getSharedPtrInternal<Int>(name);
           })
      .def("getInternalReal",
           [](ConstitutiveLawInternalHandler & self,
              const ID & id) -> decltype(auto) {
             return self.getSharedPtrInternal<Real>(id);
           })
      .def("getInternalInt",
           [](ConstitutiveLawInternalHandler & self, const ID & id)
               -> decltype(auto) { return self.getSharedPtrInternal<Int>(id); })
      .def("getElementFilter",
           [](const ConstitutiveLawInternalHandler & self) -> decltype(auto) {
             return self.getElementFilterSharedPtr();
           })
      .def("getSpatialDimension",
           &ConstitutiveLawInternalHandler::getSpatialDimension);
}

} // namespace akantu
