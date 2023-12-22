/**
 * @file   py_mesh.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Philip Mueller <philip.paul.mueller@bluemail.ch>
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Jun 16 2019
 * @date last modification: Mon Mar 15 2021
 *
 * @brief  pybind11 interface to Mesh
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
#include "aka_config.hh"
/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <mesh.hh>
#include <mesh_accessor.hh>
#include <mesh_utils.hh>
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
  decltype(auto) register_element_type_map_array(py::module & mod,
                                                 const std::string & name) {
    return py::class_<ElementTypeMapArray<T>,
                      std::shared_ptr<ElementTypeMapArray<T>>>(
               mod, ("ElementTypeMapArray" + name).c_str())
        .def(py::init<const ID &, const ID &>(),
             py::arg("id") = "by_element_type_array",
             py::arg("parent_id") = "no_parent")
        .def(
            "__call__",
            [](ElementTypeMapArray<T> & self, ElementType type,
               GhostType ghost_type) -> decltype(auto) {
              return self(type, ghost_type);
            },
            py::arg("type"), py::arg("ghost_type") = _not_ghost,
            py::return_value_policy::reference, py::keep_alive<0, 1>())
        .def(
            "elementTypes",
            [](ElementTypeMapArray<T> & self, UInt _dim, GhostType _ghost_type,
               ElementKind _kind) -> std::vector<ElementType> {
              auto types = self.elementTypes(_dim, _ghost_type, _kind);
              std::vector<ElementType> _types;
              for (auto && t : types) {
                _types.push_back(t);
              }
              return _types;
            },
            py::arg("dim") = _all_dimensions,
            py::arg("ghost_type") = _not_ghost, py::arg("kind") = _ek_regular)
        .def(
            "initialize",
            [](ElementTypeMapArray<T> & self, const Mesh & mesh,
               GhostType ghost_type = _casper, UInt nb_component = 1,
               UInt spatial_dimension = UInt(-2),
               ElementKind element_kind = _ek_not_defined,
               bool with_nb_element = false,
               bool with_nb_nodes_per_element = false, T default_value = T(),
               bool do_not_default = false) {
              self.initialize(
                  mesh, _ghost_type = ghost_type, _nb_component = nb_component,
                  _spatial_dimension = (spatial_dimension == UInt(-2)
                                            ? mesh.getSpatialDimension()
                                            : spatial_dimension),
                  _element_kind = element_kind,
                  _with_nb_element = with_nb_element,
                  _with_nb_nodes_per_element = with_nb_nodes_per_element,
                  _default_value = default_value,
                  _do_not_default = do_not_default);
            },
            py::arg("mesh"), py::arg("ghost_type") = _casper,
            py::arg("nb_component") = 1,
            py::arg("spatial_dimension") = UInt(-2),
            py::arg("element_kind") = _ek_not_defined,
            py::arg("with_nb_element") = false,
            py::arg("with_nb_nodes_per_element") = false,
            py::arg("default_value") = T(), py::arg("do_not_default") = false);
  }
} // namespace
/* -------------------------------------------------------------------------- */
void register_mesh(py::module & mod) {
  py::class_<Mesh::PeriodicSlaves>(mod, "PeriodicSlaves")
      .def(
          "__iter__",
          [](Mesh::PeriodicSlaves & _this) {
            py::make_iterator(_this.begin(), _this.end());
          },
          py::keep_alive<0, 1>());

  py::class_<MeshData>(mod, "MeshData")
      .def(
          "getElementalDataUInt",
          [](MeshData & _this, const ID & name) -> decltype(auto) {
            return _this.getElementalData<UInt>(name);
          },
          py::return_value_policy::reference)
      .def(
          "getElementalDataReal",
          [](MeshData & _this, const ID & name) -> decltype(auto) {
            return _this.getElementalData<Real>(name);
          },
          py::return_value_policy::reference);

  py::class_<Mesh, GroupManager, Dumpable, MeshData>(mod, "Mesh",
                                                     py::multiple_inheritance())
      .def(py::init<UInt, const ID &>(), py::arg("spatial_dimension"),
           py::arg("id") = "mesh")
      .def("read", &Mesh::read, py::arg("filename"),
           py::arg("mesh_io_type") = _miot_auto, "read the mesh from a file")
      .def(
          "getNodes",
          [](Mesh & self) -> decltype(auto) { return self.getNodes(); },
          py::return_value_policy::reference)
      .def("isLocalOrMasterNode", &Mesh::isLocalOrMasterNode)
      .def("getNbNodes", &Mesh::getNbNodes)
      .def("getGlobalNodesIds", &Mesh::getGlobalNodesIds)
      .def("getLowerBounds", &Mesh::getLowerBounds)
      .def("getUpperBounds", &Mesh::getUpperBounds)
      .def(
          "getNodesFlags",
          [](const Mesh & self) -> decltype(auto) {
            return self.getNodesFlags();
          },
          py::return_value_policy::reference)
      .def(
          "getConnectivity",
          [](Mesh & self, ElementType type) -> decltype(auto) {
            return self.getConnectivity(type);
          },
          py::return_value_policy::reference)
      .def(
          "getConnectivity",
          [](Mesh & self, ElementType type,
             GhostType ghost_type) -> decltype(auto) {
            return self.getConnectivity(type, ghost_type);
          },
          py::return_value_policy::reference)
      .def(
          "getConnectivities",
          [](Mesh & self) -> decltype(auto) {
            return self.getConnectivities();
          },
          py::return_value_policy::reference)
      .def(
          "addConnectivityType",
          [](Mesh & self, ElementType type, GhostType ghost_type) -> void {
            self.addConnectivityType(type, ghost_type);
          },
          py::arg("type"), py::arg("ghost_type") = _not_ghost)
      .def("distribute", [](Mesh & self) { self.distribute(); })
      .def("isDistributed",
           [](const Mesh & self) { return self.isDistributed(); })
      .def("fillNodesToElements", &Mesh::fillNodesToElements,
           py::arg("dimension") = _all_dimensions)
      .def("getAssociatedElements",
           [](Mesh & self, const UInt & node, py::list list) {
             Array<Element> elements;
             self.getAssociatedElements(node, elements);
             for (auto && element : elements) {
               list.append(element);
             }
           })
      .def("makePeriodic",
           [](Mesh & self, const SpatialDirection & direction) {
             self.makePeriodic(direction);
           })
      .def(
          "getNbElement",
          [](Mesh & self, const UInt spatial_dimension, GhostType ghost_type,
             ElementKind kind) {
            return self.getNbElement(spatial_dimension, ghost_type, kind);
          },
          py::arg("spatial_dimension") = _all_dimensions,
          py::arg("ghost_type") = _not_ghost, py::arg("kind") = _ek_not_defined)
      .def(
          "getNbElement",
          [](Mesh & self, ElementType type, GhostType ghost_type) {
            return self.getNbElement(type, ghost_type);
          },
          py::arg("type"), py::arg("ghost_type") = _not_ghost)
      .def_static(
          "getSpatialDimension",
          [](ElementType & type) { return Mesh::getSpatialDimension(type); })
      .def_static("getNaturalSpaceDimension",
                  [](ElementType & type) {
                    return Mesh::getNaturalSpaceDimension(type);
                  })
      .def_static("getKind",
                  [](ElementType & type) { return Mesh::getKind(type); })
      .def_static(
          "getNbNodesPerElement",
          [](ElementType & type) { return Mesh::getNbNodesPerElement(type); })
      .def(
          "getDataReal",
          [](Mesh & _this, const ID & name, ElementType type,
             GhostType ghost_type) -> decltype(auto) {
            return _this.getData<Real>(name, type, ghost_type);
          },
          py::arg("name"), py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)
      .def(
          "hasDataReal",
          [](Mesh & _this, const ID & name, ElementType type,
             GhostType ghost_type) -> bool {
            return _this.hasData<Real>(name, type, ghost_type);
          },
          py::arg("name"), py::arg("type"), py::arg("ghost_type") = _not_ghost)
      .def("isPeriodic", [](const Mesh & _this) { return _this.isPeriodic(); })
      .def("getPeriodicMaster", &Mesh::getPeriodicMaster)
      .def("getPeriodicSlaves", &Mesh::getPeriodicSlaves)
      .def("isPeriodicSlave", &Mesh::isPeriodicSlave)
      .def("isPeriodicMaster", &Mesh::isPeriodicMaster)
      .def(
          "getMeshFacets",
          [](const Mesh & self) -> const Mesh & {
            return self.getMeshFacets();
          },
          py::return_value_policy::reference)
      .def("initMeshFacets", &Mesh::initMeshFacets,
           py::arg("id") = "mesh_facets", py::return_value_policy::reference)
      .def(
          "elementTypes",
          [](Mesh & self, UInt spatial_dimension, GhostType ghost_type,
             ElementKind kind) -> std::vector<ElementType> {
            auto types = self.elementTypes(spatial_dimension, ghost_type, kind);
            std::vector<ElementType> _types;
            for (auto && t : types) {
              _types.push_back(t);
            }
            return _types;
          },
          py::arg("spatial_dimension") = _all_dimensions,
          py::arg("ghost_type") = _not_ghost, py::arg("kind") = _ek_regular)
      .def(
          "getDataPointerUInt",
          [](Mesh & self, const std::string & data_name, ElementType el_type,
             GhostType ghost_type, UInt nb_component, bool size_to_nb_element,
             bool resize_with_parent) -> decltype(auto) {
            return self.getDataPointer<UInt>(data_name, el_type, ghost_type,
                                             nb_component, size_to_nb_element,
                                             resize_with_parent);
          },
          py::return_value_policy::reference, py::arg("data_name"),
          py::arg("el_type"), py::arg("ghost_type") = _not_ghost,
          py::arg("nb_component") = 1, py::arg("size_to_nb_element") = true,
          py::arg("resize_with_parent") = false);

  /* ------------------------------------------------------------------------
   */
  py::class_<MeshUtils>(mod, "MeshUtils")
      .def_static("buildFacets", &MeshUtils::buildFacets);

  py::class_<MeshAccessor>(mod, "MeshAccessor")
      .def(py::init<Mesh &>(), py::arg("mesh"))
      .def(
          "resizeConnectivity",
          [](MeshAccessor & self, UInt new_size, ElementType type, GhostType gt)
              -> void { self.resizeConnectivity(new_size, type, gt); },
          py::arg("new_size"), py::arg("type"),
          py::arg("ghost_type") = _not_ghost)
      .def(
          "resizeNodes",
          [](MeshAccessor & self, UInt new_size) -> void {
            self.resizeNodes(new_size);
          },
          py::arg("new_size"))
      .def("makeReady", &MeshAccessor::makeReady);

  register_element_type_map_array<Real>(mod, "Real");
  register_element_type_map_array<UInt>(mod, "UInt");
  register_element_type_map_array<bool>(mod, "bool");
  // register_element_type_map_array<std::string>(mod, "String");
}
} // namespace akantu
