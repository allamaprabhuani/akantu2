/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <non_linear_solver.hh>
#include <solid_mechanics_model_cohesive.hh>
/* -------------------------------------------------------------------------- */
//#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
#define def_deprecated(func_name, mesg)                                        \
  def(func_name, [](py::args, py::kwargs) { AKANTU_ERROR(mesg); })

#define def_function_nocopy(func_name)                                         \
  def(#func_name,                                                              \
      [](SolidMechanicsModel & self) -> decltype(auto) {                       \
        return self.func_name();                                               \
      },                                                                       \
      py::return_value_policy::reference)

#define def_function(func_name)                                                \
  def(#func_name, [](SolidMechanicsModel & self) -> decltype(auto) {           \
    return self.func_name();                                                   \
  })
/* -------------------------------------------------------------------------- */

[[gnu::visibility("default")]] void
register_solid_mechanics_model_cohesive(py::module & mod) {

  py::class_<SolidMechanicsModelCohesiveOptions, SolidMechanicsModelOptions>(
      mod, "SolidMechanicsModelCohesiveOptions")
      .def(py::init<AnalysisMethod, bool>(),
           py::arg("analysis_method") = _explicit_lumped_mass,
           py::arg("extrinsic") = false);

  py::class_<SolidMechanicsModelCohesive, SolidMechanicsModel>(
      mod, "SolidMechanicsModelCohesive")
      .def(py::init<Mesh &, UInt, const ID &, const MemoryID &>(),
           py::arg("mesh"), py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "solid_mechanics_model", py::arg("memory_id") = 0)
      .def("checkCohesiveStress",
           &SolidMechanicsModelCohesive::checkCohesiveStress);
}

} // namespace akantu
