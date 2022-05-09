/**
 * @file   py_poisson_model.cc
 *
 * @author Mohit Pundir <mpundir@ethz.ch>
 *
 * @date creation: Mon May 09 2022
 * @date last modification: Mon May 09 2022
 *
 * @brief  pybind11 interface to Poisson Model
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
#include <non_linear_solver.hh>
#include <poisson_model.hh>
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
      [](PoissonModel & self) -> decltype(auto) {                       \
        return self.func_name();                                               \
      },                                                                       \
      py::return_value_policy::reference)

#define def_function(func_name)                                                \
  def(#func_name, [](PoissonModel & self) -> decltype(auto) {           \
    return self.func_name();                                                   \
  })

/* -------------------------------------------------------------------------- */
void register_poisson_model(py::module & mod) {

  py::class_<PoissonModelOptions>(mod, "PoissonModelOptions")
      .def(py::init<AnalysisMethod>(),
           py::arg("_analysis_method") = _explicit_lumped_mass);

  py::class_<PoissonModel, Model>(mod, "PoissonModel",
				  py::multiple_inheritance())
      .def(py::init<Mesh &, UInt, const ID &, std::shared_ptr<DOFManager>,
                    const ModelType>(),
           py::arg("mesh"), py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "poisson_model",
           py::arg("dof_manager") = nullptr,
           py::arg("model_type") = ModelType::_poisson_model)
      .def(
          "initFull",
          [](PoissonModel & self,
             const PoissonModelOptions & options) {
            self.initFull(options);
          },
          py::arg("option") = PoissonModelOptions())
      .def(
          "initFull",
          [](PoissonModel & self,
             const AnalysisMethod & analysis_method) {
            self.initFull(_analysis_method = analysis_method);
          },
          py::arg("_analysis_method"))
      .def_deprecated("applyDirichletBC", "Deprecated: use applyBC")
      .def("applyBC",
           [](PoissonModel & self,
              BC::Dirichlet::DirichletFunctor & func,
              const std::string & element_group) {
             self.applyBC(func, element_group);
           })
      .def("applyBC",
           [](PoissonModel & self, BC::Neumann::NeumannFunctor & func,
              const std::string & element_group) {
             self.applyBC(func, element_group);
           })
      .def("setTimeStep", &PoissonModel::setTimeStep,
           py::arg("time_step"), py::arg("solver_id") = "")
  

      .def_function(assembleStiffnessMatrix)
      .def_function(assembleInternalDofRate)
      .def_function(assembleCapacity)
      .def_function(assembleCapacityLumped)
      .def_function(getStableTimeStep)
      .def_function_nocopy(getDof)
      .def_function_nocopy(getInternalDofRate)
      .def_function_nocopy(getExternalDofRate)
      .def_function_nocopy(getBlockedDOFs)
      .def_function_nocopy(getMesh)
      .def(
          "getMaterial",
          [](PoissonModel & self, UInt constitutive_law_id) -> decltype(auto) {
            return self.getConstitutiveLaw(constitutive_law_id);
          },
          py::arg("constitutive_law_id"), py::return_value_policy::reference)
      .def(
          "getConstitutiveLaw",
          [](PoissonModel & self, const ID & constitutive_law_name)
              -> decltype(auto) { return self.getConstitutiveLaw(constitutive_law_name); },
          py::arg("constitutive_law_name"), py::return_value_policy::reference)
      .def("getNbConstitutiveLaws", &PoissonModel::getNbConstitutiveLaws)
      .def("getConstitutiveLawIndex", &PoissonModel::getConstitutiveLawIndex)
      .def("setConstitutiveLawSelector",
           [](PoissonModel & self,
              std::shared_ptr<ConstitutiveLawSelector> constitutive_law_selector) {
             std::cout << (*constitutive_law_selector)(ElementNull) << std::endl;
             self.setConstitutiveLawSelector(constitutive_law_selector);
           })
      .def("getConstitutiveLawSelector", &PoissonModel::getConstitutiveLawSelector)
      .def(
          "getConstitutiveLawByElement",
          [](const PoissonModel & self) -> decltype(auto) {
            return self.getConstitutiveLawByElement();
          },
          py::return_value_policy::reference, py::keep_alive<0, 1>())
      .def("reassignConstitutiveLaw", &PoissonModel::reassignConstitutiveLaw)
      .def(
          "registerNewConstitutiveLaw",
          [](PoissonModel & self, const ID & law_name,
             const ID & law_type, const ID & opt_param) -> decltype(auto) {
            return self.registerNewConstitutiveLaw(law_name, law_type, opt_param);
          },
          py::arg("constitutive_law_name"), py::arg("constitutive_law_type"),
          py::arg("option") = "", py::return_value_policy::reference)
      .def("initConstitutiveLaws", &PoissonModel::initConstitutiveLaws);
}

} // namespace akantu
