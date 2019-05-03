/* -------------------------------------------------------------------------- */
#include "py_aka_array.cc"
/* -------------------------------------------------------------------------- */
#include <non_linear_solver.hh>
#include <solid_mechanics_model.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
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
register_solid_mechanics_model(py::module & mod) {

  py::class_<SolidMechanicsModelOptions>(mod, "SolidMechanicsModelOptions")
      .def(py::init<AnalysisMethod>(),
           py::arg("analysis_method") = _explicit_lumped_mass);

  py::class_<Model>(mod, "Model")
      .def("setBaseName", &Model::setBaseName)
      .def("getFEEngine", &Model::getFEEngine, py::arg("name") = "",
           py::return_value_policy::reference)
      .def("addDumpFieldVector", &Model::addDumpFieldVector)
      .def("addDumpField", &Model::addDumpField)
      .def("dump", &Model::dump);

  py::class_<NonLinearSolver>(mod, "NonLinearSolver")
      .def(
          "set",
          [](NonLinearSolver & self, const std::string & id, const Real & val) {
            if (id == "max_iterations")
              self.set(id, int(val));
            else if (id == "convergence_type")
              self.set(id, akantu::SolveConvergenceCriteria(UInt(val)));
            else
              self.set(id, val);
          })
      .def("set", &NonLinearSolver::set<SolveConvergenceCriteria>);

  py::class_<ModelSolver>(mod, "ModelSolver")
      .def("getNonLinearSolver",
           (NonLinearSolver & (ModelSolver::*)(const ID &)) &
               ModelSolver::getNonLinearSolver,
           py::arg("solver_id") = "", py::return_value_policy::reference)
      .def("solveStep", &ModelSolver::solveStep, py::arg("solver_id") = "");

  py::class_<SolidMechanicsModel, Model, ModelSolver>(mod,
                                                      "SolidMechanicsModel")
      .def(py::init<Mesh &, UInt, const ID &, const MemoryID &,
                    const ModelType>(),
           py::arg("mesh"), py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "solid_mechanics_model", py::arg("memory_id") = 0,
           py::arg("model_type") = ModelType::_solid_mechanics_model)
      .def("initFull",
           [](SolidMechanicsModel & self,
              const SolidMechanicsModelOptions & options) {
             self.initFull(options);
           },
           py::arg("_analysis_method") = SolidMechanicsModelOptions())
      .def("initFull",
           [](SolidMechanicsModel & self,
              const AnalysisMethod & _analysis_method) {
             self.initFull(SolidMechanicsModelOptions(_analysis_method));
           },
           py::arg("_analysis_method"))
      .def_deprecated("applyDirichletBC", "Deprecated: use applyBC")
      .def("applyBC",
           [](SolidMechanicsModel & self, BC::DirichletFunctor & func,
              const std::string & element_group) {
             self.applyBC(func, element_group);
           })
      .def("applyBC",
           [](SolidMechanicsModel & self, BC::NeumannFunctor & func,
              const std::string & element_group) {
             self.applyBC(func, element_group);
           })
      .def("setTimeStep", &SolidMechanicsModel::setTimeStep,
           py::arg("time_step"), py::arg("solver_id") = "")
      .def("getEnergy",
           py::overload_cast<const std::string &>(
               &SolidMechanicsModel::getEnergy),
           py::arg("energy_id"))
      .def_function(assembleStiffnessMatrix)
      .def_function(assembleInternalForces)
      .def_function(getStableTimeStep)
      .def_function_nocopy(getExternalForce)
      .def_function_nocopy(getDisplacement)
      .def_function_nocopy(getPreviousDisplacement)
      .def_function_nocopy(getIncrement)
      .def_function_nocopy(getMass)
      .def_function_nocopy(getVelocity)
      .def_function_nocopy(getAcceleration)
      .def_function_nocopy(getInternalForce)
      .def_function_nocopy(getBlockedDOFs)
      .def_function_nocopy(getIncrementFlag)
      .def_function_nocopy(getMesh);
}

} // namespace akantu
