/**
 * @file   py_contact_mechanics_internodes_model.cc
 *
 * @author Moritz Waldleben <moritz.waldleben@epfl.ch>
 *
 * @date creation: Thu Jul 09 2022
 * @date last modification: Thu Jul 09 2022
 *
 * @brief  Contact mechanics internodes python binding
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
#include <contact_detector_internodes.hh>
#include <contact_mechanics_internodes_model.hh>
#include <parsable.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
#define def_function_nocopy(func_name)                                         \
  def(                                                                         \
      #func_name,                                                              \
      [](ContactMechanicsInternodesModel & self) -> decltype(auto) {           \
        return self.func_name();                                               \
      },                                                                       \
      py::return_value_policy::reference)

#define def_function_nocopy_detector(func_name)                                \
  def(                                                                         \
      #func_name,                                                              \
      [](ContactDetectorInternodes & self) -> decltype(auto) {                 \
        return self.func_name();                                               \
      },                                                                       \
      py::return_value_policy::reference)

#define def_function(func_name)                                                \
  def(#func_name, [](ContactMechanicsInternodesModel & self) -> decltype(auto) { \
    return self.func_name();                                                   \
  })

/* -------------------------------------------------------------------------- */
void register_contact_mechanics_internodes_model(py::module & mod) {
    py::class_<ContactDetectorInternodes>(mod, "ContactDetectorInternodes",
                      py::multiple_inheritance())
      .def(py::init<Mesh &, const ID &>(), py::arg("mesh"),
           py::arg("id") = "contact_detector_internodes")
      .def_function_nocopy_detector(findContactNodes)
      .def_function_nocopy_detector(getInitialMasterNodeGroup)
      .def_function_nocopy_detector(getInitialSlaveNodeGroup)
      .def_function_nocopy_detector(getMasterNodeGroup)
      .def_function_nocopy_detector(getSlaveNodeGroup)
      .def_function_nocopy_detector(getMasterRadiuses)
      .def_function_nocopy_detector(getSlaveRadiuses)

      .def(
          "constructInterpolationMatrix",
          [](ContactDetectorInternodes & self, NodeGroup & ref_node_group,
            NodeGroup & eval_node_group, Array<Real> eval_radiuses)
              -> decltype(auto) { return self.constructInterpolationMatrix(ref_node_group, eval_node_group, eval_radiuses); },
          py::arg("ref_node_group"), py::arg("eval_node_group"), py::arg("eval_radiuses"), py::return_value_policy::reference);

  /* ------------------------------------------------------------------------ */
    py::class_<ContactMechanicsInternodesModel, Model>(mod,
        "ContactMechanicsInternodesModel", py::multiple_inheritance())
        .def(py::init<Mesh &, UInt, const ID &, std::shared_ptr<DOFManager>,
                      const ModelType>(),
             py::arg("mesh"), py::arg("spatial_dimension") = _all_dimensions,
             py::arg("id") = "contact_mechanics_internodes_model",
             py::arg("dof_manager") = nullptr,
             py::arg("model_type") = ModelType::_solid_mechanics_model)
        .def(
            "initFull",
            [](ContactMechanicsInternodesModel & self,
               const AnalysisMethod & analysis_method) {
              self.initFull(_analysis_method = analysis_method);
            },
            py::arg("_analysis_method"))
        .def_function_nocopy(assembleInternodesMatrix)
        .def_function_nocopy(getSolidMechanicsModel)
        .def_function_nocopy(getContactDetectorInternodes)
        .def_function_nocopy(getLambdas)
        .def("assembleInterfaceMass",
            [](ContactMechanicsInternodesModel & self, NodeGroup & contact_node_group)
                -> decltype(auto) { return self.assembleInterfaceMass(contact_node_group); },
            py::arg("contact_node_group"), py::return_value_policy::copy);
}

} // namespace akantu
