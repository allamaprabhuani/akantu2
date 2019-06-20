/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <coupler_solid_contact.hh>
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
      [](CouplerSolidContact & self) -> decltype(auto) {                       \
        return self.func_name();                                               \
      },                                                                       \
      py::return_value_policy::reference)

#define def_function(func_name)                                                \
  def(#func_name, [](CouplerSolidContact & self) -> decltype(auto) {           \
    return self.func_name();                                                   \
  })
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
void register_model_couplers(py::module & mod) {
  py::class_<CouplerSolidContactOptions>(mod, "CouplerSolidContactOptions")
      .def(py::init<AnalysisMethod>(),
           py::arg("analysis_method") = _explicit_dynamic_contact);

  py::class_<CouplerSolidContact, Model>(mod, "CouplerSolidContact")
      .def(py::init<Mesh &, UInt, const ID &, const MemoryID &,
                    std::shared_ptr<DOFManager>, const ModelType>(),
           py::arg("mesh"), py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "coupler_solid_contact", py::arg("memory_id") = 0,
           py::arg("dof_manager") = nullptr,
           py::arg("model_type") = ModelType::_coupler_solid_contact)
      .def("initFull",
           [](CouplerSolidContact & self,
              const CouplerSolidContactOptions & options) {
             self.initFull(options);
           },
           py::arg("_analysis_method") = CouplerSolidContactOptions())
      .def("initFull",
           [](CouplerSolidContact & self,
              const AnalysisMethod & analysis_method) {
             self.initFull(_analysis_method = analysis_method);
           },
           py::arg("_analysis_method"))
      .def("setTimeStep", &CouplerSolidContact::setTimeStep,
           py::arg("time_step"), py::arg("solver_id") = "")
      .def("getSolidMechanicsModel",
           &CouplerSolidContact::getSolidMechanicsModel,
           py::return_value_policy::reference)
      .def("getContactMechanicsModel",
           &CouplerSolidContact::getContactMechanicsModel,
           py::return_value_policy::reference)
      .def("dump", py::overload_cast<>(&CouplerSolidContact::dump))
      .def("dump",
           py::overload_cast<const std::string &>(&CouplerSolidContact::dump))
      .def("dump", py::overload_cast<const std::string &, UInt>(
                       &CouplerSolidContact::dump))
      .def("dump", py::overload_cast<const std::string &, Real, UInt>(
                       &CouplerSolidContact::dump));
}

} // namespace akantu
