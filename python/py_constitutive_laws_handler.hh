/* -------------------------------------------------------------------------- */
#include <constitutive_laws_handler.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_PY_CONSTITUTIVE_LAWS_HANDLER_HH_
#define AKANTU_PY_CONSTITUTIVE_LAWS_HANDLER_HH_

namespace akantu {

template <class ConstitutiveLawType, class Model_>
void register_constitutive_laws_handler(pybind11::module & mod) {
  namespace py = pybind11;
  using CLH = ConstitutiveLawsHandler<ConstitutiveLawType, Model_>;
  std::string name = "ConstitutiveLawsHandler" +
                     debug::demangle<ConstitutiveLawType>() +
                     debug::demangle<Model_>();

  py::class_<CLH, Model_>(mod, name.c_str(), py::multiple_inheritance())
      .def(py::init<Mesh &, ModelType &, Int, const ID &>())
      .def(
          "registerNewConstitutiveLaw",
          [](CLH & self, const ID & name, const ID & type,
             const ID & opt_param) -> decltype(auto) {
            return self.registerNewConstitutiveLaw(name, type, opt_param);
          },
          py::arg("name"), py::arg("type"), py::arg("opt_param") = "",
          py::return_value_policy::reference)
      .def("reassignConstitutiveLaw", &CLH::reassignConstitutiveLaw)
      .def("isInternal", &CLH::isInternal, py::arg("field_name"),
           py::arg("element_kind") = _ek_regular)
      .def(
          "flattenInternalReal",
          [](CLH & self, const std::string & field_name, ElementKind kind,
             GhostType ghost_type) -> decltype(auto) {
            return self.template flattenInternal<Real>(field_name, kind,
                                                       ghost_type);
          },
          py::arg("field_name"), py::arg("kind") = _ek_regular,
          py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)
      .def(
          "inflateInternalReal",
          [](CLH & self, const std::string & field_name,
             const std::shared_ptr<ElementTypeMapArray<Real>> & field,
             ElementKind kind, GhostType ghost_type) {
            return self.template inflateInternal<Real>(field_name, *field,
                                                       ghost_type, kind);
          },
          py::arg("field_name"), py::arg("field"),
          py::arg("kind") = _ek_regular, py::arg("ghost_type") = _not_ghost)
      .def(
          "getConstitutiveLaws",
          [](CLH & self) -> decltype(auto) {
            return self.getConstitutiveLaws();
          },
          py::return_value_policy::reference)
      .def(
          "getConstitutiveLaw",
          [](CLH & self, Idx cl_index) -> decltype(auto) {
            return self.getConstitutiveLaw(cl_index);
          },
          py::return_value_policy::reference)
      .def(
          "getConstitutiveLaw",
          [](CLH & self, const std::string & name) -> decltype(auto) {
            return self.getConstitutiveLaw(name);
          },
          py::return_value_policy::reference)
      .def(
          "getConstitutiveLaw",
          [](CLH & self, const Element & element) -> decltype(auto) {
            return self.getConstitutiveLaw(element);
          },
          py::return_value_policy::reference)
      .def(
          "getConstitutiveLawIndex",
          [](CLH & self, const std::string & name) -> decltype(auto) {
            return self.getConstitutiveLawIndex(name);
          },
          py::return_value_policy::reference)
      .def(
          "getNbConstitutiveLaws",
          [](CLH & self) -> decltype(auto) {
            return self.getNbConstitutiveLaws();
          },
          py::return_value_policy::reference)
      .def(
          "getConstitutiveLawByElement",
          [](const CLH & self) -> decltype(auto) {
            return self.getConstitutiveLawByElement();
          },
          py::return_value_policy::reference)
      .def(
          "getConstitutiveLawLocalNumbering",
          [](const CLH & self) -> decltype(auto) {
            return self.getConstitutiveLawLocalNumbering();
          },
          py::return_value_policy::reference);
}
} // namespace akantu

#endif /* AKANTU_PY_CONSTITUTIVE_LAWS_HANDLER_HH_ */
