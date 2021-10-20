/**
 * @file   py_contact_mechanics_model.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jun 20 2019
 * @date last modification: Thu Jun 24 2021
 *
 * @brief  Contact mechanics python binding
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
#include <contact_detector.hh>
#include <contact_element.hh>
#include <contact_mechanics_model.hh>
#include <geometry_utils.hh>
#include <mesh_events.hh>
#include <parsable.hh>
#include <surface_selector.hh>
/* -------------------------------------------------------------------------- */
#include <algorithm>
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
      [](ContactMechanicsModel & self) -> decltype(auto) {                     \
        return self.func_name();                                               \
      },                                                                       \
      py::return_value_policy::reference)

#define def_function(func_name)                                                \
  def(#func_name, [](ContactMechanicsModel & self) -> decltype(auto) {         \
    return self.func_name();                                                   \
  })

namespace {
  class ContactElementsView {
  public:
    ContactElementsView(const Array<ContactElement> & contact_elements)
        : contact_elements(contact_elements) {}

    auto begin() const { return contact_elements.begin(); }
    auto end() const { return contact_elements.end(); }

    auto size() const { return contact_elements.size(); }
    auto operator[](size_t i) const { return contact_elements(i); }

    auto contains(const ContactElement & contact_element) const {
      return std::find(contact_elements.begin(), contact_elements.end(),
                       contact_element) != contact_elements.end();
    }

  private:
    const Array<ContactElement> & contact_elements;
  };
} // namespace

/* -------------------------------------------------------------------------- */
void register_contact_mechanics_model(py::module & mod) {
  py::class_<ContactDetector>(mod, "ContactDetector",
                              py::multiple_inheritance())
      .def(py::init<Mesh &, const ID &>(), py::arg("mesh"),
           py::arg("id") = "contact_detector")
      .def(py::init<Mesh &, Array<Real>, const ID &>(), py::arg("mesh"),
           py::arg("positions"), py::arg("id") = "contact_detector")
      .def("setSurfaceSelector", &ContactDetector::setSurfaceSelector);

  py::class_<SurfaceSelector, std::shared_ptr<SurfaceSelector>>(
      mod, "SurfaceSelector", py::multiple_inheritance())
      .def(py::init<Mesh &>(), py::arg("mesh"));

  py::class_<PhysicalSurfaceSelector, SurfaceSelector,
             std::shared_ptr<PhysicalSurfaceSelector>>(
      mod, "PhysicalSurfaceSelector")
      .def(py::init<Mesh &>(), py::arg("mesh"));

  py::class_<CohesiveSurfaceSelector, SurfaceSelector,
             std::shared_ptr<CohesiveSurfaceSelector>>(
      mod, "CohesiveSurfaceSelector")
      .def(py::init<Mesh &>(), py::arg("mesh"));

  py::class_<AllSurfaceSelector, SurfaceSelector,
             std::shared_ptr<AllSurfaceSelector>>(mod, "AllSurfaceSelector")
      .def(py::init<Mesh &>(), py::arg("mesh"));

  py::class_<ContactMechanicsModelOptions>(mod, "ContactMechanicsModelOptions")
      .def(py::init<AnalysisMethod>(),
           py::arg("analysis_method") = _explicit_contact);

  /* ------------------------------------------------------------------------ */
  py::class_<ContactElementsView>(mod, "ContactElementsView")
      .def("__iter__",
           [](const ContactElementsView & self) {
             return py::make_iterator(self.begin(), self.end());
           })
      .def("__size__",
           [](const ContactElementsView & self) { return self.size(); })
      .def(
          "__contains__",
          [](const ContactElementsView & self, const ContactElement & element) {
            return self.contains(element);
          })
      .def("__getitem__",
           [](const ContactElementsView & self, size_t i) { return self[i]; });

  /* ------------------------------------------------------------------------ */
  py::class_<ContactMechanicsModel, Model>(mod, "ContactMechanicsModel",
                                           py::multiple_inheritance())
      .def(py::init<Mesh &, UInt, const ID &, std::shared_ptr<DOFManager>,
                    const ModelType>(),
           py::arg("mesh"), py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "contact_mechanics_model",
           py::arg("dof_manager") = nullptr,
           py::arg("model_type") = ModelType::_contact_mechanics_model)
      .def(
          "initFull",
          [](ContactMechanicsModel & self,
             const ContactMechanicsModelOptions & options) {
            self.initFull(options);
          },
          py::arg("options") = ContactMechanicsModelOptions())
      .def(
          "initFull",
          [](ContactMechanicsModel & self,
             const AnalysisMethod & analysis_method) {
            self.initFull(_analysis_method = analysis_method);
          },
          py::arg("_analysis_method"))
      .def_function(search)
      .def_function(assembleStiffnessMatrix)
      .def_function(assembleInternalForces)
      .def_function_nocopy(getExternalForce)
      .def_function_nocopy(getNormalForce)
      .def_function_nocopy(getTangentialForce)
      .def_function_nocopy(getInternalForce)
      .def_function_nocopy(getGaps)
      .def_function_nocopy(getNormals)
      .def_function_nocopy(getNodalArea)
      .def_function_nocopy(getContactDetector)
      .def("getContactElements", [](ContactMechanicsModel & self) {
        return ContactElementsView(self.getContactElements());
      });

  py::class_<ContactElement>(mod, "ContactElement")
      .def(py::init<>())
      .def_readwrite("master", &ContactElement::master)
      .def_readwrite("slave", &ContactElement::slave)
      .def("__repr__", [](ContactElement & self) {
        return "{master: " + std::to_string(self.master) +
               ", slave: " + std::to_string(self.slave) + "}";
      });

  py::class_<GeometryUtils>(mod, "GeometryUtils")
      .def_static(
          "normal",
          py::overload_cast<const Mesh &, const Array<Real> &, const Element &,
                            Vector<Real> &, bool>(&GeometryUtils::normal),
          py::arg("mesh"), py::arg("positions"), py::arg("element"),
          py::arg("normal"), py::arg("outward") = true)
      .def_static(
          "covariantBasis",
          py::overload_cast<const Mesh &, const Array<Real> &, const Element &,
                            const Vector<Real> &, Vector<Real> &,
                            Matrix<Real> &>(&GeometryUtils::covariantBasis),
          py::arg("mesh"), py::arg("positions"), py::arg("element"),
          py::arg("normal"), py::arg("natural_projection"), py::arg("basis"))
      .def_static("curvature", &GeometryUtils::curvature)
      .def_static("contravariantBasis", &GeometryUtils::contravariantBasis,
                  py::arg("covariant_basis"), py::arg("basis"))
      .def_static("realProjection",
                  py::overload_cast<const Mesh &, const Array<Real> &,
                                    const Vector<Real> &, const Element &,
                                    const Vector<Real> &, Vector<Real> &>(
                      &GeometryUtils::realProjection),
                  py::arg("mesh"), py::arg("positions"), py::arg("slave"),
                  py::arg("element"), py::arg("normal"), py::arg("projection"))
      .def_static("isBoundaryElement", &GeometryUtils::isBoundaryElement);
}

} // namespace akantu
