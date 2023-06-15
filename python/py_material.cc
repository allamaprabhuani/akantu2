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
#include <material_cohesive.hh>
#include <material_cohesive_bilinear.hh>
#include <material_cohesive_linear.hh>
#include <material_cohesive_linear_friction.hh>
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

    template <typename T>
    void setDefaultValueToInternal(const ID & int_id, const T value) {
      auto & internal_field = this->template getInternal<T>(int_id);
      internal_field.setDefaultValue(value);
    }

  protected:
    std::map<std::string, std::shared_ptr<ElementTypeMapBase>> internals;
  };

  /* ------------------------------------------------------------------------ */
  template <typename _Material>
  void register_material_classes(py::module & mod, const std::string & name) {
    py::class_<_Material, Material, Parsable, PyMaterial<_Material>>(
        mod, name.c_str(), py::multiple_inheritance())
        .def(py::init<SolidMechanicsModel &, const ID &>())
        .def("registerInternalReal",
             [](Material & self, const std::string & name, Int nb_component) {
               return dynamic_cast<PyMaterial<_Material> &>(self)
                   .template registerInternal<Real>(name, nb_component);
             })
        .def("registerInternalUInt",
             [](Material & self, const std::string & name, Int nb_component) {
               return dynamic_cast<PyMaterial<_Material> &>(self)
                   .template registerInternal<UInt>(name, nb_component);
             })
        .def("setDefaultValueToInternalReal",
             [](Material & self, const ID & int_id, const Real value) {
               return dynamic_cast<PyMaterial<_Material> &>(self)
                   .template setDefaultValueToInternal<Real>(int_id, value);
             })
        .def("setDefaultValueToInternalUInt",
             [](Material & self, const ID & int_id, const UInt value) {
               return dynamic_cast<PyMaterial<_Material> &>(self)
                   .template setDefaultValueToInternal<UInt>(int_id, value);
             });
  }

#if defined(AKANTU_COHESIVE_ELEMENT)
  // trampoline for the cohesive materials
  template <typename _Material>
  class PyMaterialCohesive : public PyMaterial<_Material> {
    using Parent = PyMaterial<_Material>;

  public:
    using Parent::Parent;

    void checkInsertion(bool check_only) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, _Material, checkInsertion, check_only);
    }

    void computeTraction(const Array<Real> & normal, ElementType el_type,
                         GhostType ghost_type = _not_ghost) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE_PURE(void, _Material, computeTraction, normal, el_type,
                             ghost_type);
    }

    void computeTangentTraction(ElementType el_type,
                                Array<Real> & tangent_matrix,
                                const Array<Real> & normal,
                                GhostType ghost_type = _not_ghost) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, _Material, computeTangentTraction, el_type,
                        tangent_matrix, normal, ghost_type);
    }

    void computeStress(ElementType /*el_type*/,
                       GhostType /*ghost_type*/ = _not_ghost) final {}

    void computeTangentModuli(ElementType /*el_type*/,
                              Array<Real> & /*tangent_matrix*/,
                              GhostType /*ghost_type*/ = _not_ghost) final {}

    template <typename T>
    void registerInternal(const std::string & name, Int nb_component) {
      auto && internal =
          std::make_shared<CohesiveInternalField<T>>(name, *this);
      AKANTU_DEBUG_INFO("alloc internal " << name << " "
                                          << &this->internals[name]);

      internal->initialize(nb_component);
      this->internals[name] = internal;
    }
  };

  // trampoline for the cohesive material inheritance where computeTraction is
  // not pure virtual
  template <typename _Material>
  class PyMaterialCohesiveDaughters : public PyMaterialCohesive<_Material> {
    using Parent = PyMaterialCohesive<_Material>;

  public:
    using Parent::Parent;

    void computeTraction(const Array<Real> & normal, ElementType el_type,
                         GhostType ghost_type = _not_ghost) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, _Material, computeTraction, normal, el_type,
                        ghost_type);
    }

    void computeTangentTraction(ElementType el_type,
                                Array<Real> & tangent_matrix,
                                const Array<Real> & normal,
                                GhostType ghost_type = _not_ghost) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, _Material, computeTangentTraction, el_type,
                        tangent_matrix, normal, ghost_type);
    }
  };

  template <typename _Material>
  void register_material_cohesive_classes(py::module & mod,
                                          const std::string & name) {
    py::class_<_Material, MaterialCohesive,
               PyMaterialCohesiveDaughters<_Material>>(
        mod, name.c_str(), py::multiple_inheritance())
        .def(py::init<SolidMechanicsModelCohesive &, const ID &>())
        .def("registerInternalReal",
             [](_Material & self, const std::string & name, UInt nb_component) {
               return dynamic_cast<PyMaterialCohesive<_Material> &>(self)
                   .template registerInternal<Real>(name, nb_component);
             })
        .def("registerInternalUInt",
             [](_Material & self, const std::string & name, UInt nb_component) {
               return dynamic_cast<PyMaterialCohesive<_Material> &>(self)
                   .template registerInternal<UInt>(name, nb_component);
             })
        .def("setDefaultValueToInternalReal",
             [](_Material & self, const ID & int_id, const Real value) {
               return dynamic_cast<PyMaterialCohesive<_Material> &>(self)
                   .template setDefaultValueToInternal<Real>(int_id, value);
             })
        .def("setDefaultValueToInternalUInt",
             [](_Material & self, const ID & int_id, const UInt value) {
               return dynamic_cast<PyMaterialCohesive<_Material> &>(self)
                   .template setDefaultValueToInternal<UInt>(int_id, value);
             })
        .def("computeTraction",
             [](_Material & self, const Array<Real> & normal,
                ElementType el_type, GhostType ghost_type) {
               return dynamic_cast<PyMaterialCohesiveDaughters<_Material> &>(
                          self)
                   .computeTraction(normal, el_type, ghost_type);
             })
        .def("computeTangentTraction", [](_Material & self, ElementType el_type,
                                          Array<Real> & tangent_matrix,
                                          const Array<Real> & normal,
                                          GhostType ghost_type) {
          return dynamic_cast<PyMaterialCohesiveDaughters<_Material> &>(self)
              .computeTangentTraction(el_type, tangent_matrix, normal,
                                      ghost_type);
        });

    ;
  }

#endif

  /* ------------------------------------------------------------------------ */
  template <typename T>
  void register_internal_field(py::module & mod, const std::string & name) {
    py::class_<InternalField<T>, ElementTypeMapArray<T>,
               std::shared_ptr<InternalField<T>>>(
        mod, ("InternalField" + name).c_str());
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
          py::arg("el_type"), py::return_value_policy::reference)
      .def(
          "getPotentialEnergy",
          [](Material & self, ElementType el_type, UInt index) -> Real {
            return self.getPotentialEnergy(el_type, index);
          },
          py::arg("el_type"), py::arg("index"))
      .def("getPotentialEnergy",
           [](Material & self) -> Real { return self.getPotentialEnergy(); })
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
      .def("setDefaultValueToInternalReal",
           [](Material & self, const ID & int_id, const Real value) {
             return dynamic_cast<PyMaterial<Material> &>(self)
                 .setDefaultValueToInternal<Real>(int_id, value);
           })
      .def("setDefaultValueToInternalUInt",
           [](Material & self, const ID & int_id, const UInt value) {
             return dynamic_cast<PyMaterial<Material> &>(self)
                 .setDefaultValueToInternal<UInt>(int_id, value);
           })
      .def(
          "getInternalReal",
          [](Material & self, const ID & id) -> decltype(auto) {
            return self.getInternal<Real>(id);
          },
          py::arg("id"), py::return_value_policy::reference)
      .def(
          "getInternalRealPrevious",
          [](Material & self, const ID & id) -> decltype(auto) {
            return self.getInternal<Real>(id).previous();
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
      /*
       * These functions override the `Parsable` interface.
       * This ensure that the `updateInternalParameters()` function is called.
       */
      .def(
          "setReal",
          [](Material & self, const ID & name, const Real value) -> void {
            self.setParam(name, value);
            return;
          },
          py::arg("name"), py::arg("value"))
      .def(
          "setBool",
          [](Material & self, const ID & name, const bool value) -> void {
            self.setParam(name, value);
            return;
          },
          py::arg("name"), py::arg("value"))
      .def(
          "setString",
          [](Material & self, const ID & name,
             const std::string & value) -> void {
            self.setParam(name, value);
            return;
          },
          py::arg("name"), py::arg("value"))
      .def(
          "setInt",
          [](Material & self, const ID & name, const int value) -> void {
            self.setParam(name, value);
            return;
          },
          py::arg("name"), py::arg("value"))
      .def("getPushWaveSpeed", &Material::getPushWaveSpeed)
      .def("getShearWaveSpeed", &Material::getShearWaveSpeed)
      .def("__repr__", [](Material & self) {
        std::stringstream sstr;
        sstr << self;
        return sstr.str();
      });

  register_material_classes<MaterialElastic<2>>(mod, "MaterialElastic2D");
  register_material_classes<MaterialElastic<3>>(mod, "MaterialElastic3D");

#if defined(AKANTU_COHESIVE_ELEMENT)
  /* ------------------------------------------------------------------------ */
  py::class_<MaterialCohesive, Material, Parsable,
             PyMaterialCohesive<MaterialCohesive>>(mod, "MaterialCohesive",
                                                   py::multiple_inheritance())
      .def(py::init<SolidMechanicsModelCohesive &, const ID &>())
      .def("registerInternalReal",
           [](MaterialCohesive & self, const std::string & name,
              UInt nb_component) {
             return dynamic_cast<PyMaterialCohesive<MaterialCohesive> &>(self)
                 .registerInternal<Real>(name, nb_component);
           })
      .def("registerInternalUInt",
           [](MaterialCohesive & self, const std::string & name,
              UInt nb_component) {
             return dynamic_cast<PyMaterialCohesive<MaterialCohesive> &>(self)
                 .registerInternal<UInt>(name, nb_component);
           })
      .def(
          "getFacetFilter",
          [](MaterialCohesive & self) -> decltype(auto) {
            return self.getFacetFilter();
          },
          py::return_value_policy::reference)
      .def(
          "getFacetFilter",
          [](MaterialCohesive & self, ElementType type,
             GhostType ghost_type) -> decltype(auto) {
            return self.getFacetFilter(type, ghost_type);
          },
          py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)
      .def(
          "getOpening",
          [](MaterialCohesive & self, ElementType type, GhostType ghost_type)
              -> decltype(auto) { return self.getOpening(type, ghost_type); },
          py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)
      .def(
          "getTraction",
          [](MaterialCohesive & self, ElementType type, GhostType ghost_type)
              -> decltype(auto) { return self.getTraction(type, ghost_type); },
          py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)
      .def(
          "computeTraction",
          [](MaterialCohesive & self, GhostType ghost_type) {
            return self.computeTraction(ghost_type);
          },
          py::arg("ghost_type") = _not_ghost)
      .def(
          "computeTraction",
          [](MaterialCohesive & self, const Array<Real> & normal,
             ElementType el_type, GhostType ghost_type) {
            return dynamic_cast<PyMaterialCohesive<MaterialCohesive> &>(self)
                .computeTraction(normal, el_type, ghost_type);
          },
          py::arg("normal"), py::arg("el_type"),
          py::arg("ghost_type") = _not_ghost)
      .def(
          "computeTangentTraction",
          [](MaterialCohesive & self, ElementType el_type,
             Array<Real> & tangent_matrix, const Array<Real> & normal,
             GhostType ghost_type) {
            return dynamic_cast<PyMaterialCohesive<MaterialCohesive> &>(self)
                .computeTangentTraction(el_type, tangent_matrix, normal,
                                        ghost_type);
          },
          py::arg("el_type"), py::arg("tangent_matrix"), py::arg("normal"),
          py::arg("ghost_type") = _not_ghost)
      .def(
          "computeNormal",
          [](MaterialCohesive & self, const Array<Real> & position,
             Array<Real> & normal, ElementType type, GhostType ghost_type) {
            self.computeNormal(position, normal, type, ghost_type);
          },
          py::arg("position"), py::arg("normal"), py::arg("type"),
          py::arg("ghost_type") = _not_ghost)
      .def("getNormalsAtQuads", &MaterialCohesive::getNormalsAtQuads);

  register_material_cohesive_classes<MaterialCohesiveLinear<2>>(
      mod, "MaterialCohesiveLinear2D");
  register_material_cohesive_classes<MaterialCohesiveLinear<3>>(
      mod, "MaterialCohesiveLinear3D");
  register_material_cohesive_classes<MaterialCohesiveBilinear<2>>(
      mod, "MaterialCohesiveBilinear2D");
  register_material_cohesive_classes<MaterialCohesiveBilinear<3>>(
      mod, "MaterialCohesiveBilinear3D");
  register_material_cohesive_classes<MaterialCohesiveLinearFriction<2>>(
      mod, "MaterialCohesiveLinearFriction2D");
  register_material_cohesive_classes<MaterialCohesiveLinearFriction<3>>(
      mod, "MaterialCohesiveLinearFriction3D");
#endif
}

} // namespace akantu
