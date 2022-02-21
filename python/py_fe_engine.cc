/**
 * @file   py_fe_engine.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 27 2019
 * @date last modification: Sat Dec 12 2020
 *
 * @brief  pybind11 interface to FEEngine
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
#include "py_aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <element.hh>
#include <fe_engine.hh>
#include <integration_point.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

void register_fe_engine(py::module & mod) {
  py::class_<Element>(mod, "Element")
      .def(py::init([](ElementType type, UInt id) {
        return new Element{type, id, _not_ghost};
      }))
      .def(py::init([](ElementType type, UInt id, GhostType ghost_type) {
        return new Element{type, id, ghost_type};
      }))
      .def("__lt__",
           [](Element & self, const Element & other) { return (self < other); })
      .def("__repr__", [](Element & self) { return std::to_string(self); });

  mod.attr("ElementNull") = ElementNull;

  py::class_<FEEngine>(mod, "FEEngine")
      .def(
          "getNbIntegrationPoints",
          [](FEEngine & fem, const ElementType & type,
             const GhostType & ghost_type) {
            return fem.getNbIntegrationPoints(type, ghost_type);
          },
          py::arg("type"), py::arg("ghost_type") = _not_ghost)
      .def(
          "gradientOnIntegrationPoints",
          [](FEEngine & fem, const Array<Real> & u, Array<Real> & nablauq,
             UInt nb_degree_of_freedom, ElementType type, GhostType ghost_type,
             const Array<UInt> * filter_elements) {
            if (filter_elements == nullptr) {
              // This is due to the ArrayProxy that looses the
              // empty_filter information
              filter_elements = &empty_filter;
            }
            fem.gradientOnIntegrationPoints(u, nablauq, nb_degree_of_freedom,
                                            type, ghost_type, *filter_elements);
          },
          py::arg("u"), py::arg("nablauq"), py::arg("nb_degree_of_freedom"),
          py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::arg("filter_elements") = nullptr)
      .def(
          "interpolateOnIntegrationPoints",
          [](FEEngine & self, const Array<Real> & u, Array<Real> & uq,
             UInt nb_degree_of_freedom, ElementType type, GhostType ghost_type,
             const Array<UInt> * filter_elements) {
            if (filter_elements == nullptr) {
              // This is due to the ArrayProxy that looses the
              // empty_filter information
              filter_elements = &empty_filter;
            }

            self.interpolateOnIntegrationPoints(u, uq, nb_degree_of_freedom,
                                                type, ghost_type,
                                                *filter_elements);
          },
          py::arg("u"), py::arg("uq"), py::arg("nb_degree_of_freedom"),
          py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::arg("filter_elements") = nullptr)
      .def(
          "interpolateOnIntegrationPoints",
          [](FEEngine & self, const Array<Real> & u,
             ElementTypeMapArray<Real> & uq,
             const ElementTypeMapArray<UInt> * filter_elements) {
            self.interpolateOnIntegrationPoints(u, uq, filter_elements);
          },
          py::arg("u"), py::arg("uq"), py::arg("filter_elements") = nullptr)
      .def(
          "computeIntegrationPointsCoordinates",
          [](FEEngine & self, ElementTypeMapArray<Real> & coordinates,
             const ElementTypeMapArray<UInt> * filter_elements)
              -> decltype(auto) {
            return self.computeIntegrationPointsCoordinates(coordinates,
                                                            filter_elements);
          },
          py::arg("coordinates"), py::arg("filter_elements") = nullptr)
      .def(
          "assembleFieldLumped",
          [](FEEngine & fem,
             const std::function<void(Matrix<Real> &, const Element &)> &
                 field_funct,
             const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
             ElementType type, GhostType ghost_type) {
            fem.assembleFieldLumped(field_funct, matrix_id, dof_id, dof_manager,
                                    type, ghost_type);
          },
          py::arg("field_funct"), py::arg("matrix_id"), py::arg("dof_id"),
          py::arg("dof_manager"), py::arg("type"),
          py::arg("ghost_type") = _not_ghost)
      .def(
          "assembleFieldMatrix",
          [](FEEngine & fem,
             const std::function<void(Matrix<Real> &, const Element &)> &
                 field_funct,
             const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
             ElementType type, GhostType ghost_type = _not_ghost) {
            fem.assembleFieldMatrix(field_funct, matrix_id, dof_id, dof_manager,
                                    type, ghost_type);
          },
          py::arg("field_funct"), py::arg("matrix_id"), py::arg("dof_id"),
          py::arg("dof_manager"), py::arg("type"),
          py::arg("ghost_type") = _not_ghost)
      .def("getElementInradius", [](FEEngine & self, const Element & element) {
        return self.getElementInradius(element);
      });

  py::class_<IntegrationPoint>(mod, "IntegrationPoint");
}
} // namespace akantu
