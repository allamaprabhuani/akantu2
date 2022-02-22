/**
 * @file   py_structural_mechanics_model.cc
 *
 * @author Philip Mueller <philip.paul.mueller@bluemail.ch>
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Feb 03 2021
 * @date last modification: Thu Apr 01 2021
 *
 * @brief  pybind11 interface to StructuralMechanicsModel
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
#include <structural_mechanics_model.hh>
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
      [](StructuralMechanicsModel & self) -> decltype(auto) {                  \
        return self.func_name();                                               \
      },                                                                       \
      py::return_value_policy::reference)

#define def_function_(func_name)                                               \
  def(#func_name, [](StructuralMechanicsModel & self) -> decltype(auto) {      \
    return self.func_name();                                                   \
  })

#define def_plainmember(M) def_readwrite(#M, &StructuralMaterial::M)
/* -------------------------------------------------------------------------- */

void register_structural_mechanics_model(pybind11::module & mod) {
  /* First we have to register the material class
   * The wrapper aims to mimic the behaviour of the real material.
   */
  py::class_<StructuralMaterial>(mod, "StructuralMaterial")
      .def(py::init<>())
      .def(py::init<const StructuralMaterial &>())
      .def_plainmember(E)
      .def_plainmember(A)
      .def_plainmember(I)
      .def_plainmember(Iz)
      .def_plainmember(Iy)
      .def_plainmember(GJ)
      .def_plainmember(rho)
      .def_plainmember(t)
      .def_plainmember(nu);

  /* Now we create the structural model wrapper
   * Note that this is basically a port from the solid mechanic part.
   */
  py::class_<StructuralMechanicsModel, Model>(mod, "StructuralMechanicsModel")
      .def(py::init<Mesh &, UInt, const ID &>(), py::arg("mesh"),
           py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "structural_mechanics_model")
      .def(
          "initFull",
          [](StructuralMechanicsModel & self,
             const AnalysisMethod & analysis_method) -> void {
            self.initFull(_analysis_method = analysis_method);
          },
          py::arg("_analysis_method"))
      .def("initFull",
           [](StructuralMechanicsModel & self) -> void { self.initFull(); })
      .def_function_nocopy(getExternalForce)
      .def_function_nocopy(getDisplacement)
      .def_function_nocopy(getInternalForce)
      .def_function_nocopy(getVelocity)
      .def_function_nocopy(getAcceleration)
      .def_function_nocopy(getInternalForce)
      .def_function_nocopy(getBlockedDOFs)
      .def_function_nocopy(getMesh)

      .def("setTimeStep", &StructuralMechanicsModel::setTimeStep,
           py::arg("time_step"), py::arg("solver_id") = "")
      .def(
          "getElementMaterial",
          [](StructuralMechanicsModel & self, const ElementType & type,
             GhostType ghost_type) -> decltype(auto) {
            return self.getElementMaterial(type, ghost_type);
          },
          "This function returns the map that maps elements to materials.",
          py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)
      .def(
          "getMaterialByElement",
          [](StructuralMechanicsModel & self, Element element)
              -> decltype(auto) { return self.getMaterialByElement(element); },
          "This function returns the `StructuralMaterial` instance that is "
          "associated with element `element`.",
          py::arg("element"), py::return_value_policy::reference)
      .def(
          "addMaterial",
          [](StructuralMechanicsModel & self, StructuralMaterial & mat,
             const ID & name) -> UInt { return self.addMaterial(mat, name); },
          "This function adds the `StructuralMaterial` `mat` to `self`."
          " The function returns the ID of the new material.",
          py::arg("mat"), py::arg("name") = "")
      .def(
          "getMaterial",
          [](StructuralMechanicsModel & self, UInt material_index)
              -> decltype(auto) { return self.getMaterial(material_index); },
          "This function returns the `i`th material of `self`."
          " Note a reference is returned which allows to modify the material "
          "inside `self`.",
          py::arg("i"), py::return_value_policy::reference)
      .def(
          "getMaterial",
          [](StructuralMechanicsModel & self, const ID & name)
              -> decltype(auto) { return self.getMaterial(name); },
          "This function returns the material with name `i` of `self`."
          " Note a reference is returned which allows to modify the material "
          "inside `self`.",
          py::arg("i"), py::return_value_policy::reference)
      .def(
          "getNbMaterials",
          [](StructuralMechanicsModel & self) { return self.getNbMaterials(); },
          "Returns the number of different materials inside `self`.")
      .def("getKineticEnergy", &StructuralMechanicsModel::getKineticEnergy,
           "Compute kinetic energy")
      .def("getPotentialEnergy", &StructuralMechanicsModel::getPotentialEnergy,
           "Compute potential energy")
      .def("getEnergy", &StructuralMechanicsModel::getEnergy,
           "Compute the specified energy")
      .def("getLumpedMass",
       	   [](StructuralMechanicsModel& self) -> decltype(auto) {
     	   	return self.getLumpedMass(); },
       	   py::return_value_policy::reference_internal)
      .def("getLumpedMassMutable",
       	   [](StructuralMechanicsModel& self) -> decltype(auto) {
      	   	return self.getLumpedMassMutable(); },
      	   py::return_value_policy::reference_internal)

      .def("computeShadyLumpedMass",
      	   [](StructuralMechanicsModel& self,  GhostType ghost_type) -> void {
      	   	return self.computeShadyLumpedMass(ghost_type); },
           py::arg("ghost_type") = _not_ghost )
      ;

} // End: register structural mechanical model

} // namespace akantu
