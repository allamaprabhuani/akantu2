/**
 * @file   py_material_selector.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed May 26 2021
 * @date last modification: Wed May 26 2021
 *
 * @brief  Material selector python binding
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
#include "py_constitutive_law_selector.hh"
/* -------------------------------------------------------------------------- */
#include <constitutive_law_selector.hh>
#if defined(AKANTU_SOLID_MECHANICS)
#include <material_selector.hh>
#include <solid_mechanics_model.hh>
#endif
#if defined(AKANTU_COHESIVE_ELEMENT)
#include <material_selector_cohesive.hh>
#include <solid_mechanics_model_cohesive.hh>
#endif
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace akantu {
namespace {
  template <class Base = ConstitutiveLawSelector>
  class PyConstitutiveLawSelector : public Base {
  public:
    /* Inherit the constructors */
    using Base::Base;

    Idx operator()(const Element & element) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE_NAME(Idx, ConstitutiveLawSelector,
                             "__call__", operator(), element);
    }
  };

  template <class ConstitutiveLawSelectorDaughter>
  decltype(auto) register_material_selectors(py::module & mod,
                                             const std::string & class_name) {
    return py::class_<
        ConstitutiveLawSelectorDaughter, ConstitutiveLawSelector,
        PyConstitutiveLawSelector<ConstitutiveLawSelectorDaughter>,
        std::shared_ptr<ConstitutiveLawSelectorDaughter>>(mod,
                                                          class_name.c_str());
  }
} // namespace

void register_constitutive_law_selector(py::module & mod) {
  py::class_<ConstitutiveLawSelector, PyConstitutiveLawSelector<>,
             std::shared_ptr<ConstitutiveLawSelector>>(
      mod, "ConstitutiveLawSelector")
      .def(py::init())
      .def("setFallback",
           [](ConstitutiveLawSelector & self, Int f) { self.setFallback(f); })
      .def("setFallback",
           [](ConstitutiveLawSelector & self,
              const std::shared_ptr<ConstitutiveLawSelector> &
                  fallback_selector) { self.setFallback(fallback_selector); })
      .def("__call__", &ConstitutiveLawSelector::operator());

#if defined(AKANTU_SOLID_MECHANICS)
  register_material_selectors<DefaultMaterialSelector>(
      mod, "DefaultMaterialSelector")
      .def(py::init<const ElementTypeMapArray<Idx> &>());

  register_material_selectors<MeshDataMaterialSelector<std::string>>(
      mod, "MeshDataMaterialSelectorString")
      .def(py::init<const std::string &, const SolidMechanicsModel &, Int>(),
           py::arg("name"), py::arg("model"), py::arg("first_index") = 1);
#endif

#if defined(AKANTU_COHESIVE_ELEMENT)
  register_material_selectors<DefaultMaterialCohesiveSelector>(
      mod, "DefaultMaterialCohesiveSelector")
      .def(py::init<const SolidMechanicsModelCohesive &>());

  register_material_selectors<MeshDataMaterialCohesiveSelector>(
      mod, "MeshDataMaterialCohesiveSelector")
      .def(py::init<const SolidMechanicsModelCohesive &>());

  register_material_selectors<MaterialCohesiveRulesSelector>(
      mod, "MaterialCohesiveRulesSelector")
      .def(py::init<const SolidMechanicsModelCohesive &,
                    const MaterialCohesiveRules &, const ID &>(),
           py::arg("model"), py::arg("rules"),
           py::arg("mesh_data_id") = "physical_names");
#endif
}
} // namespace akantu
