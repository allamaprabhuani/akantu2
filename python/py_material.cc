/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <material_selector.hh>
#include <solid_mechanics_model.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

template <typename _Material> class PyMaterial : public _Material {

public:
  /* Inherit the constructors */
  using _Material::_Material;

  virtual ~PyMaterial(){};
  void initMaterial() override {
    PYBIND11_OVERLOAD(void, _Material, initMaterial);
  };
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override {
    PYBIND11_OVERLOAD_PURE(void, _Material, computeStress, el_type, ghost_type);
  }
  void computeTangentModuli(const ElementType & el_type,
                            Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override {
    PYBIND11_OVERLOAD(void, _Material, computeTangentModuli, el_type,
                      tangent_matrix, ghost_type);
  }

  void computePotentialEnergy(ElementType el_type) override {
    PYBIND11_OVERLOAD(void, _Material, computePotentialEnergy, el_type);
  }

  Real getPushWaveSpeed(const Element & element) const override {
    PYBIND11_OVERLOAD(Real, _Material, getPushWaveSpeed, element);
  }

  Real getShearWaveSpeed(const Element & element) const override {
    PYBIND11_OVERLOAD(Real, _Material, getShearWaveSpeed, element);
  }

  void registerInternal(const std::string & name, UInt nb_component) {
    this->internals[name] = std::make_shared<InternalField<Real>>(name, *this);
    AKANTU_DEBUG_INFO("alloc internal " << name << " "
                                        << &this->internals[name]);

    this->internals[name]->initialize(nb_component);
  }

  auto & getInternals() { return this->internals; }

protected:
  std::map<std::string, std::shared_ptr<InternalField<Real>>> internals;
};

/* -------------------------------------------------------------------------- */

template <typename T>
void register_element_type_map_array(py::module & mod,
                                     const std::string & name) {

  py::class_<ElementTypeMapArray<T>, std::shared_ptr<ElementTypeMapArray<T>>>(
      mod, ("ElementTypeMapArray" + name).c_str())
      .def("__call__",
           [](ElementTypeMapArray<T> & self, ElementType & type,
              const GhostType & ghost_type) -> decltype(auto) {
             return self(type, ghost_type);
           },
           py::arg("type"), py::arg("ghost_type") = _not_ghost,
           py::return_value_policy::reference)
      .def("elementTypes",
           [](ElementTypeMapArray<T> & self, UInt _dim, GhostType _ghost_type,
              ElementKind _kind) -> decltype(auto) {
             auto types = self.elementTypes(_dim, _ghost_type, _kind);
             std::vector<ElementType> _types;
             for (auto && t : types) {
               _types.push_back(t);
             }
             return _types;
           },
           py::arg("dim") = _all_dimensions, py::arg("ghost_type") = _not_ghost,
           py::arg("kind") = _ek_regular);

  py::class_<InternalField<T>, ElementTypeMapArray<T>,
             std::shared_ptr<InternalField<T>>>(
      mod, ("InternalField" + name).c_str());
}

/* -------------------------------------------------------------------------- */
template <typename _Material>
void define_material(py::module & mod, const std::string & name) {

  auto mat = py::class_<_Material, PyMaterial<_Material>, Parsable>(
      mod, name.c_str(), py::multiple_inheritance());

  mat.def(py::init<SolidMechanicsModel &, const ID &>())
      .def("getGradU",
           [](Material & self, ElementType el_type,
              GhostType ghost_type = _not_ghost) -> decltype(auto) {
             return self.getGradU(el_type, ghost_type);
           },
           py::arg("el_type"), py::arg("ghost_type") = _not_ghost,
           py::return_value_policy::reference)
      .def("getStress",
           [](Material & self, ElementType el_type,
              GhostType ghost_type = _not_ghost) -> decltype(auto) {
             return self.getStress(el_type, ghost_type);
           },
           py::arg("el_type"), py::arg("ghost_type") = _not_ghost,
           py::return_value_policy::reference)
      .def("getPotentialEnergy",
           [](Material & self, ElementType el_type) -> decltype(auto) {
             return self.getPotentialEnergy(el_type);
           },
           py::return_value_policy::reference)
      .def("initMaterial", &Material::initMaterial)
      .def("getModel", &Material::getModel)
      .def("registerInternal",
           [](Material & self, const std::string & name, UInt nb_component) {
             return dynamic_cast<PyMaterial<Material> &>(self).registerInternal(
                 name, nb_component);
           })
      .def_property_readonly(
          "internals",
          [](Material & self) {
            return dynamic_cast<PyMaterial<Material> &>(self).getInternals();
          })
      .def_property_readonly("element_filter",
                             [](Material & self) -> decltype(auto) {
                               return self.getElementFilter();
                             },
                             py::return_value_policy::reference);
}

/* -------------------------------------------------------------------------- */

[[gnu::visibility("default")]] void register_material(py::module & mod) {
  py::class_<MaterialFactory>(mod, "MaterialFactory")
      .def_static("getInstance",
                  []() -> MaterialFactory & { return Material::getFactory(); },
                  py::return_value_policy::reference)
      .def("registerAllocator",
           [](MaterialFactory & self, const std::string id, py::function func) {
             self.registerAllocator(
                 id,
                 [func, id](UInt dim, const ID &, SolidMechanicsModel & model,
                            const ID & option) -> std::unique_ptr<Material> {
                   py::object obj = func(dim, id, model, option);
                   auto & ptr = py::cast<Material &>(obj);

                   obj.release();
                   return std::unique_ptr<Material>(&ptr);
                 });
           });

  register_element_type_map_array<Real>(mod, "Real");
  register_element_type_map_array<UInt>(mod, "UInt");

  define_material<Material>(mod, "Material");
}

/* -------------------------------------------------------------------------- */
template <typename T>
void register_data_material_selector(py::module & mod,
                                          const std::string & name) {
  py::class_<ElementDataMaterialSelector<T>, MaterialSelector,
             std::shared_ptr<ElementDataMaterialSelector<T>>>(
      mod, ("ElementDataMaterialSelector" + name).c_str());

  py::class_<MeshDataMaterialSelector<T>, ElementDataMaterialSelector<T>,
             std::shared_ptr<MeshDataMaterialSelector<T>>>(
      mod, ("MeshDataMaterialSelector" + name).c_str())
      .def(py::init<const std::string &, SolidMechanicsModel &, UInt>(),
           py::arg("name"), py::arg("model"), py::arg("first_index") = 1);
}

/* -------------------------------------------------------------------------- */
void register_material_selector(py::module & mod) {
  py::class_<MaterialSelector, std::shared_ptr<MaterialSelector>>(
      mod, "MaterialSelector");

  register_data_material_selector<std::string>(mod, "String");
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
