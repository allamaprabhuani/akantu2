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
 * @brief  pybind11 interface to ConstitutiveLaw
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
#include <constitutive_law_selector.hh>
#include <poisson_model.hh>
#include <constitutive_law.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

namespace {
  template <typename _ConstitutiveLaw> class PyConstitutiveLaw : public _ConstitutiveLaw {
  public:
    /* Inherit the constructors */
    using _ConstitutiveLaw::_ConstitutiveLaw;

    ~PyConstitutiveLaw() override = default;

    void initConstitutiveLaw() override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, _ConstitutiveLaw, initConstitutiveLaw, );
    };
    void computeFlux(ElementType el_type,
                       GhostType ghost_type = _not_ghost) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE_PURE(void, _ConstitutiveLaw, computeFlux, el_type,
                             ghost_type);
    }
    void computeTangentModuli(ElementType el_type, Array<Real> & tangent_matrix,
                              GhostType ghost_type = _not_ghost) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, _ConstitutiveLaw, computeTangentModuli, el_type,
                        tangent_matrix, ghost_type);
    }

    Real getCelerity() const override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(Real, _ConstitutiveLaw, getCelerity);
    }

    Real getEffectiveCapacity() const override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(Real, _ConstitutiveLaw, getEffectiveCapacity);
    }

    template <typename T>
    void registerInternal(const std::string & name, UInt nb_component) {
      auto && internal = std::make_shared<InternalConstitutiveLaw<T>>(name, *this);
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
  void register_internal_constitutive_law(py::module & mod, const std::string & name) {
    py::class_<InternalConstitutiveLaw<T>, ElementTypeMapArray<T>,
               std::shared_ptr<InternalConstitutiveLaw<T>>>(
        mod, ("InternalConstitutiveLaw" + name).c_str());
  }

  /* ------------------------------------------------------------------------ */
  //template <typename _ConstitutiveLaw>
  //void register_constitutive_law_classes(py::module & mod, const std::string & name) {
  //  py::class_<_ConstitutiveLaw, ConstitutiveLaw, Parsable, PyConstitutiveLaw<_ConstitutiveLaw>>(
  //      mod, name.c_str(), py::multiple_inheritance())
  //      .def(py::init<PoissonModel &, const ID &>());
  //}
} // namespace

/* -------------------------------------------------------------------------- */
void register_constitutive_law(py::module & mod) {
  py::class_<ConstitutiveLawFactory>(mod, "ConstitutiveLawFactory")
      .def_static(
          "getInstance",
          []() -> ConstitutiveLawFactory & { return ConstitutiveLaw::getFactory(); },
          py::return_value_policy::reference)
      .def("registerAllocator",
           [](ConstitutiveLawFactory & self, const std::string id, py::function func) {
             self.registerAllocator(
                 id,
                 [func, id](const ID & /*unused*/,
                            PoissonModel & model,
                            const ID & option) -> std::unique_ptr<ConstitutiveLaw> {
                   py::object obj = func(id, model, option);
                   auto & ptr = py::cast<ConstitutiveLaw &>(obj);

                   obj.release();
                   return std::unique_ptr<ConstitutiveLaw>(&ptr);
                 });
           })
      .def("getPossibleAllocators", &ConstitutiveLawFactory::getPossibleAllocators);

  register_internal_constitutive_law<Real>(mod, "Real");
  register_internal_constitutive_law<UInt>(mod, "UInt");

  py::class_<ConstitutiveLaw, Parsable, PyConstitutiveLaw<ConstitutiveLaw>>(
      mod, "ConstitutiveLaw", py::multiple_inheritance())
      .def(py::init<PoissonModel &, const ID &>())
      .def(
          "getGradientDof",
          [](ConstitutiveLaw & self, ElementType el_type,
             GhostType ghost_type = _not_ghost) -> decltype(auto) {
            return self.getGradientDof(el_type, ghost_type);
          },
          py::arg("el_type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)
      .def(
          "getFluxDof",
          [](ConstitutiveLaw & self, ElementType el_type,
             GhostType ghost_type = _not_ghost) -> decltype(auto) {
            return self.getFluxDof(el_type, ghost_type);
          },
          py::arg("el_type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)
      .def("initConstitutiveLaw", &ConstitutiveLaw::initConstitutiveLaw)
      .def("getModel", &ConstitutiveLaw::getModel)
      .def("registerInternalReal",
           [](ConstitutiveLaw & self, const std::string & name, UInt nb_component) {
             return dynamic_cast<PyConstitutiveLaw<ConstitutiveLaw> &>(self)
                 .registerInternal<Real>(name, nb_component);
           })
      .def("registerInternalUInt",
           [](ConstitutiveLaw & self, const std::string & name, UInt nb_component) {
             return dynamic_cast<PyConstitutiveLaw<ConstitutiveLaw> &>(self)
                 .registerInternal<UInt>(name, nb_component);
           })
      .def(
          "getInternalReal",
          [](ConstitutiveLaw & self, const ID & id) -> decltype(auto) {
            return self.getInternal<Real>(id);
          },
          py::arg("id"), py::return_value_policy::reference)
      .def(
          "getInternalUInt",
          [](ConstitutiveLaw & self, const ID & id) -> decltype(auto) {
            return self.getInternal<UInt>(id);
          },
          py::arg("id"), py::return_value_policy::reference)
      .def(
          "getElementFilter",
          [](ConstitutiveLaw & self) -> decltype(auto) {
            return self.getElementFilter();
          },
          py::return_value_policy::reference)

      /*
       * These functions override the `Parsable` interface.
       * This ensure that the `updateInternalParameters()` function is called.
       */
      .def(
          "setReal",
          [](ConstitutiveLaw & self, const ID & name, const Real value) -> void {
            self.setParam(name, value);
            return;
          },
          py::arg("name"), py::arg("value"))
      .def(
          "setBool",
          [](ConstitutiveLaw & self, const ID & name, const bool value) -> void {
            self.setParam(name, value);
            return;
          },
          py::arg("name"), py::arg("value"))
      .def(
          "setString",
          [](ConstitutiveLaw & self, const ID & name,
             const std::string & value) -> void {
            self.setParam(name, value);
            return;
          },
          py::arg("name"), py::arg("value"))
      .def(
          "setInt",
          [](ConstitutiveLaw & self, const ID & name, const int value) -> void {
            self.setParam(name, value);
            return;
          },
          py::arg("name"), py::arg("value"))

      .def("getCelerity", &ConstitutiveLaw::getCelerity)
      .def("getEffectiveCapacity", &ConstitutiveLaw::getEffectiveCapacity)
      .def("__repr__", [](ConstitutiveLaw & self) {
        std::stringstream sstr;
        sstr << self;
        return sstr.str();
      });

  //register_constitutive_law_classes<ConstitutiveLaw>(mod, "ConstitutiveLaw");
}

} // namespace akantu
