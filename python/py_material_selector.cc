/* -------------------------------------------------------------------------- */
#include "py_material_selector.hh"
#include "py_akantu_pybind11_compatibility.hh"
/* -------------------------------------------------------------------------- */
#include <material_selector.hh>
#include <solid_mechanics_model.hh>
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
  template <class Base = MaterialSelector>
  class PyMaterialSelector : public Base {
  public:
    /* Inherit the constructors */
    using Base::Base;

    ~PyMaterialSelector() override = default;

    UInt operator()(const Element & element) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE_NAME(UInt, MaterialSelector, "__call__", operator(),
                             element);
    }
  };

  template <class MaterialSelectorDaughter>
  decltype(auto) register_material_selectors(py::module & mod,
                                             const std::string & class_name) {
    return py::class_<MaterialSelectorDaughter, MaterialSelector,
                      PyMaterialSelector<MaterialSelectorDaughter>,
                      std::shared_ptr<MaterialSelectorDaughter>>(
        mod, class_name.c_str());
  }
} // namespace

void register_material_selector(py::module & mod) {
  py::class_<MaterialSelector, PyMaterialSelector<>,
             std::shared_ptr<MaterialSelector>>(mod, "MaterialSelector")
      .def(py::init())
      .def("setFallback",
           [](MaterialSelector & self, UInt f) { self.setFallback(f); })
      .def("setFallback",
           [](MaterialSelector & self,
              const std::shared_ptr<MaterialSelector> & fallback_selector) {
             self.setFallback(fallback_selector);
           })
      .def("__call__", &MaterialSelector::operator());

  register_material_selectors<DefaultMaterialSelector>(
      mod, "DefaultMaterialSelector")
      .def(py::init<const ElementTypeMapArray<UInt>>());

  register_material_selectors<MeshDataMaterialSelector<std::string>>(
      mod, "MeshDataMaterialSelectorString")
      .def(py::init<const std::string &, const SolidMechanicsModel &, UInt>(),
           py::arg("name"), py::arg("model"), py::arg("first_index") = 1);

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
