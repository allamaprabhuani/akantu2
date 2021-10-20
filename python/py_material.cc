/**
 * @file   py_material.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jun 20 2019
 * @date last modification: Fri Apr 09 2021
 *
 * @brief  pybind11 interface to Material
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
#include "py_akantu_pybind11_compatibility.hh"
/* -------------------------------------------------------------------------- */
#include <material_selector.hh>
#include <solid_mechanics_model.hh>
#if defined(AKANTU_COHESIVE_ELEMENT)
#include <solid_mechanics_model_cohesive.hh>
#endif
#include <material_elastic.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

namespace {
  template <typename _Material> class PyMaterial : public _Material {
  public:
    /* Inherit the constructors */
    using _Material::_Material;

    ~PyMaterial() override = default;

    void initMaterial() override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, _Material, initMaterial, );
    };
    void computeStress(ElementType el_type,
                       GhostType ghost_type = _not_ghost) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE_PURE(void, _Material, computeStress, el_type,
                             ghost_type);
    }
    void computeTangentModuli(ElementType el_type, Array<Real> & tangent_matrix,
                              GhostType ghost_type = _not_ghost) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, _Material, computeTangentModuli, el_type,
                        tangent_matrix, ghost_type);
    }

    void computePotentialEnergy(ElementType el_type) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, _Material, computePotentialEnergy, el_type);
    }

    Real getPushWaveSpeed(const Element & element) const override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(Real, _Material, getPushWaveSpeed, element);
    }

    Real getShearWaveSpeed(const Element & element) const override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(Real, _Material, getShearWaveSpeed, element);
    }

    template <typename T>
    void registerInternal(const std::string & name, UInt nb_component) {
      auto && internal = std::make_shared<InternalField<T>>(name, *this);
      AKANTU_DEBUG_INFO("alloc internal " << name << " "
                                          << &this->internals[name]);

      internal->initialize(nb_component);
      this->internals[name] = internal;
    }

  protected:
    std::map<std::string, std::shared_ptr<ElementTypeMapBase>> internals;
  };

  /* ------------------------------------------------------------------------ */
  template <typename T>
  void register_internal_field(py::module & mod, const std::string & name) {
    py::class_<InternalField<T>, ElementTypeMapArray<T>,
               std::shared_ptr<InternalField<T>>>(
        mod, ("InternalField" + name).c_str());
  }

  /* ------------------------------------------------------------------------ */
  template <typename _Material>
  void register_material_classes(py::module & mod, const std::string & name) {
    py::class_<_Material, Material, Parsable, PyMaterial<_Material>>(
        mod, name.c_str(), py::multiple_inheritance())
        .def(py::init<SolidMechanicsModel &, const ID &>());
  }
} // namespace

/* -------------------------------------------------------------------------- */
void register_material(py::module & mod) {
  py::class_<MaterialFactory>(mod, "MaterialFactory")
      .def_static(
          "getInstance",
          []() -> MaterialFactory & { return Material::getFactory(); },
          py::return_value_policy::reference)
      .def("registerAllocator",
           [](MaterialFactory & self, const std::string id, py::function func) {
             self.registerAllocator(
                 id,
                 [func, id](UInt dim, const ID & /*unused*/,
                            SolidMechanicsModel & model,
                            const ID & option) -> std::unique_ptr<Material> {
                   py::object obj = func(dim, id, model, option);
                   auto & ptr = py::cast<Material &>(obj);

                   obj.release();
                   return std::unique_ptr<Material>(&ptr);
                 });
           })
      .def("getPossibleAllocators", &MaterialFactory::getPossibleAllocators);

  register_internal_field<Real>(mod, "Real");
  register_internal_field<UInt>(mod, "UInt");

  py::class_<Material, Parsable, PyMaterial<Material>>(
      mod, "Material", py::multiple_inheritance())
      .def(py::init<SolidMechanicsModel &, const ID &>())
      .def(
          "getGradU",
          [](Material & self, ElementType el_type,
             GhostType ghost_type = _not_ghost) -> decltype(auto) {
            return self.getGradU(el_type, ghost_type);
          },
          py::arg("el_type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)
      .def(
          "getStress",
          [](Material & self, ElementType el_type,
             GhostType ghost_type = _not_ghost) -> decltype(auto) {
            return self.getStress(el_type, ghost_type);
          },
          py::arg("el_type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)
      .def(
          "getPotentialEnergy",
          [](Material & self, ElementType el_type) -> decltype(auto) {
            return self.getPotentialEnergy(el_type);
          },
          py::return_value_policy::reference)
      .def("initMaterial", &Material::initMaterial)
      .def("getModel", &Material::getModel)
      .def("registerInternalReal",
           [](Material & self, const std::string & name, UInt nb_component) {
             return dynamic_cast<PyMaterial<Material> &>(self)
                 .registerInternal<Real>(name, nb_component);
           })
      .def("registerInternalUInt",
           [](Material & self, const std::string & name, UInt nb_component) {
             return dynamic_cast<PyMaterial<Material> &>(self)
                 .registerInternal<UInt>(name, nb_component);
           })
      .def(
          "getInternalReal",
          [](Material & self, const ID & id) -> decltype(auto) {
            return self.getInternal<Real>(id);
          },
          py::arg("id"), py::return_value_policy::reference)
      .def(
          "getInternalUInt",
          [](Material & self, const ID & id) -> decltype(auto) {
            return self.getInternal<UInt>(id);
          },
          py::arg("id"), py::return_value_policy::reference)
      .def(
          "getElementFilter",
          [](Material & self) -> decltype(auto) {
            return self.getElementFilter();
          },
          py::return_value_policy::reference)
      .def("getPushWaveSpeed", &Material::getPushWaveSpeed)
      .def("getShearWaveSpeed", &Material::getShearWaveSpeed);

  register_material_classes<MaterialElastic<2>>(mod, "MaterialElastic2D");
  register_material_classes<MaterialElastic<3>>(mod, "MaterialElastic3D");
}

} // namespace akantu
