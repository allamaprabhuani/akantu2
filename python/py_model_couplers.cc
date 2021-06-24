/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <cohesive_contact_solvercallback.hh>
#include <coupler_solid_cohesive_contact.hh>
#include <coupler_solid_contact.hh>
#include <non_linear_solver.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

namespace {
  template <class CouplerSolidContact_>
  decltype(auto) register_coupler_solid_contact(py::module & mod,
                                                const std::string & name) {
    auto py_class_ =
        py::class_<CouplerSolidContact_, Model>(mod, name.c_str(),
                                                py::multiple_inheritance())
            .def(py::init<Mesh &, UInt, const ID &, std::shared_ptr<DOFManager>,
                          const ModelType>(),
                 py::arg("mesh"),
                 py::arg("spatial_dimension") = _all_dimensions,
                 py::arg("id") = "coupler_solid_contact",
                 py::arg("dof_manager") = nullptr,
                 py::arg("model_type") = ModelType::_coupler_solid_contact)
            .def("applyBC",
                 [](CouplerSolidContact_ & self,
                    BC::Dirichlet::DirichletFunctor & func,
                    const std::string & element_group) {
                   self.applyBC(func, element_group);
                 })
            .def("applyBC",
                 [](CouplerSolidContact_ & self,
                    BC::Neumann::NeumannFunctor & func,
                    const std::string & element_group) {
                   self.applyBC(func, element_group);
                 })

            .def("setTimeStep", &CouplerSolidContact_::setTimeStep,
                 py::arg("time_step"), py::arg("solver_id") = "")
            .def("getContactMechanicsModel",
                 &CouplerSolidContact_::getContactMechanicsModel,
                 py::return_value_policy::reference);
    return py_class_;
  }
} // namespace

/* -------------------------------------------------------------------------- */
void register_model_couplers(py::module & mod) {
  register_coupler_solid_contact<CouplerSolidContact>(mod,
                                                      "CouplerSolidContact")
      .def(
          "initFull",
          [](SolidMechanicsModel & self,
             const AnalysisMethod & analysis_method) {
            self.initFull(_analysis_method = analysis_method);
          },
          py::arg("_analysis_method") = _explicit_lumped_mass)
      .def("getSolidMechanicsModel",
           &CouplerSolidContact::getSolidMechanicsModel,
           py::return_value_policy::reference);

  register_coupler_solid_contact<CouplerSolidCohesiveContact>(
      mod, "CouplerSolidCohesiveContact")
      .def(
          "initFull",
          [](SolidMechanicsModel & self, const AnalysisMethod & analysis_method,
             bool is_extrinsic) {
            self.initFull(_analysis_method = analysis_method,
                          _is_extrinsic = is_extrinsic);
          },
          py::arg("_analysis_method") = _explicit_lumped_mass,
          py::arg("_is_extrinsic") = false)
      .def("checkCohesiveStress",
           &CouplerSolidCohesiveContact::checkCohesiveStress)
      .def("getSolidMechanicsModelCohesive",
           &CouplerSolidCohesiveContact::getSolidMechanicsModelCohesive,
           py::return_value_policy::reference);
}

} // namespace akantu
