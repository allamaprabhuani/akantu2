/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <contact_mechanics_model.hh>
#include <surface_selector.hh>
#include <non_linear_solver.hh>
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
  def(#func_name,                                                              \
      [](ContactMechanicsModel & self) -> decltype(auto) {                     \
        return self.func_name();                                               \
      },                                                                       \
      py::return_value_policy::reference)

#define def_function(func_name)                                                \
  def(#func_name, [](ContactMechanicsModel & self) -> decltype(auto) {         \
    return self.func_name();                                                   \
  })
/* -------------------------------------------------------------------------- */

void register_contact_mechanics_model(py::module & mod) {
  py::class_<ContactMechanicsModelOptions>(mod, "ContactMechanicsModelOptions")
      .def(py::init<AnalysisMethod>(),
           py::arg("analysis_method") = _explicit_contact);

  py::class_<ContactMechanicsModel, Model>(mod, "ContactMechanicsModel",
                                           py::multiple_inheritance())
      .def(py::init<Mesh &, UInt, const ID &, const MemoryID &,
                    std::shared_ptr<DOFManager>, const ModelType>(),
           py::arg("mesh"), py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "contact_mechanics_model", py::arg("memory_id") = 0,
           py::arg("dof_manager") = nullptr,
           py::arg("model_type") = ModelType::_contact_mechanics_model)
      .def("initFull",
           [](ContactMechanicsModel & self,
              const ContactMechanicsModelOptions & options) {
             self.initFull(options);
           },
           py::arg("_analysis_method") = ContactMechanicsModelOptions())
      .def("initFull",
           [](ContactMechanicsModel & self,
              const AnalysisMethod & analysis_method) {
             self.initFull(_analysis_method = analysis_method);
           },
           py::arg("_analysis_method"))
      .def_function(search)
      .def_function(assembleStiffnessMatrix)
      .def_function(assembleInternalForces)
      .def_function_nocopy(getExternalForce)
      .def_function_nocopy(getInternalForce)
      .def_function_nocopy(getGaps)
      .def_function_nocopy(getNormals)
      .def_function_nocopy(getNodalArea)
      .def("dump", py::overload_cast<>(&ContactMechanicsModel::dump))
      .def("dump",
           py::overload_cast<const std::string &>(&ContactMechanicsModel::dump))
      .def("dump", py::overload_cast<const std::string &, UInt>(
                       &ContactMechanicsModel::dump))
      .def("dump", py::overload_cast<const std::string &, Real, UInt>(
                       &ContactMechanicsModel::dump));

  py::class_<ContactDetector>(mod, "ContactDetector")
    .def(py::init<Mesh &, const ID &, const MemoryID &>(),
	 py::arg("mesh"), py::arg("id") = "contact_detector",
	 py::arg("memory_id") = 0)
    .def("search", [](ContactDetector & self, std::map<UInt, ContactElement> & contact_map) {
	self.search(contact_map);
      });
  
  py::class_<ContactElement>(mod, "ContactElement")
    .def(py::init<>());
}

} // namespace akantu
