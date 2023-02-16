/**
 * @file   py_fluid_diffusion_model.cc
 *
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 *
 * @date creation: Thu Dec 08 2022
 * @date last modification: Thu Dec 08 2022
 *
 * @brief  pybind11 interface to FluidDiffusionModel
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
#include <fluid_diffusion_model.hh>
#include <non_linear_solver.hh>
/* -------------------------------------------------------------------------- */
// #include <pybind11/operators.h>
#include <pybind11/pybind11.h>
// #include <pybind11/stl.h>
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
      [](FluidDiffusionModel & self) -> decltype(auto) {                       \
        return self.func_name();                                               \
      },                                                                       \
      py::return_value_policy::reference)

#define def_function(func_name)                                                \
  def(#func_name, [](FluidDiffusionModel & self) -> decltype(auto) {           \
    return self.func_name();                                                   \
  })
/* -------------------------------------------------------------------------- */

void register_fluid_diffusion_model(py::module & mod) {
  py::class_<FluidDiffusionModelOptions>(mod, "FluidDiffusionModelOptions")
      .def(py::init<AnalysisMethod>(),
           py::arg("analysis_method") = _explicit_lumped_mass);

  py::class_<FluidDiffusionModel, Model>(mod, "FluidDiffusionModel",
                                         py::multiple_inheritance())
      .def(py::init<Mesh &, UInt, const ID &>(), py::arg("mesh"),
           py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "fluid_diffusion_model")
      .def(
          "initFull",
          [](FluidDiffusionModel & self,
             const FluidDiffusionModelOptions & options) {
            self.initFull(options);
          },
          py::arg("_analysis_method") = FluidDiffusionModelOptions())
      .def(
          "initFull",
          [](FluidDiffusionModel & self,
             const AnalysisMethod & _analysis_method) {
            self.initFull(FluidDiffusionModelOptions(_analysis_method));
          },
          py::arg("_analysis_method"))
      .def("setTimeStep", &FluidDiffusionModel::setTimeStep,
           py::arg("time_step"), py::arg("solver_id") = "")
      .def_function(getStableTimeStep)
      .def_function_nocopy(getPressure)
      .def_function_nocopy(getBlockedDOFs)
      .def_function_nocopy(getExternalFlux)
      .def_function_nocopy(getMesh)
      .def("getPressureGradient", &FluidDiffusionModel::getPressureGradient,
           py::arg("el_type"), py::arg("ghost_type") = _not_ghost,
           py::return_value_policy::reference)
      .def("getKgradP", &FluidDiffusionModel::getKgradP, py::arg("el_type"),
           py::arg("ghost_type") = _not_ghost,
           py::return_value_policy::reference)
      .def("applyBC",
           [](FluidDiffusionModel & self,
              BC::Dirichlet::DirichletFunctor & func,
              const std::string & element_group) {
             self.applyBC(func, element_group);
           })
      .def("applyBC",
           [](FluidDiffusionModel & self, BC::Neumann::NeumannFunctor & func,
              const std::string & element_group) {
             self.applyBC(func, element_group);
           });
}

} // namespace akantu
