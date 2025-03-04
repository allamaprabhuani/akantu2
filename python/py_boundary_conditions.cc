/**
 * Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
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
 */

/* -------------------------------------------------------------------------- */
#include "py_boundary_conditions.hh"
#include "py_akantu_pybind11_compatibility.hh"
/* -------------------------------------------------------------------------- */
#include <boundary_condition_functor.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;

namespace akantu {

namespace {
  /* ------------------------------------------------------------------------ */
  template <typename daughter = BC::Dirichlet::DirichletFunctor>
  class PyDirichletFunctor : public daughter {
  public:
    /* Inherit the constructors */
    using daughter::daughter;

    /* Trampoline (need one for each virtual function) */
    void operator()(Int node, VectorProxy<bool> & flags,
                    VectorProxy<Real> & primal,
                    const VectorProxy<const Real> & coord) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE_NAME(void, daughter, "__call__", operator(), node,
                             flags, primal, coord);
    }
  };

  /* ------------------------------------------------------------------------ */
  template <typename daughter = BC::Neumann::NeumannFunctor>
  class PyNeumannFunctor : public daughter {
  public:
    /* Inherit the constructors */
    using daughter::daughter;

    /* Trampoline (need one for each virtual function) */
    void operator()(const IntegrationPoint & quad_point,
                    VectorProxy<Real> & dual,
                    const VectorProxy<const Real> & coord,
                    const VectorProxy<const Real> & normals) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE_PURE_NAME(void, daughter, "__call__", operator(),
                                  quad_point, dual, coord, normals);
    }
  };

  /* ------------------------------------------------------------------------ */
  template <typename Functor, typename Constructor>
  decltype(auto) register_dirichlet_functor(py::module mod, const char * name,
                                            Constructor && cons) {
    py::class_<Functor, PyDirichletFunctor<Functor>,
               BC::Dirichlet::DirichletFunctor>(mod, name)
        .def(cons);
  }

  /* ------------------------------------------------------------------------ */
  template <typename Functor, typename Constructor>
  decltype(auto) register_neumann_functor(py::module mod, const char * name,
                                          Constructor && cons) {
    py::class_<Functor, PyNeumannFunctor<Functor>, BC::Neumann::NeumannFunctor>(
        mod, name)
        .def(cons);
  }
} // namespace

/* -------------------------------------------------------------------------- */
void register_boundary_conditions(py::module & mod) {

  py::class_<BC::Functor>(mod, "BCFunctor");
  py::class_<BC::Dirichlet::DirichletFunctor, PyDirichletFunctor<>,
             BC::Functor>(mod, "DirichletFunctor")
      .def(py::init())
      .def(py::init<SpatialDirection>());

  py::class_<BC::Neumann::NeumannFunctor, PyNeumannFunctor<>, BC::Functor>(
      mod, "NeumannFunctor")
      .def(py::init());

  register_dirichlet_functor<BC::Dirichlet::FixedValue>(
      mod, "FixedValue", py::init<Real, BC::Axis>());

  register_dirichlet_functor<BC::Dirichlet::IncrementValue>(
      mod, "IncrementValue", py::init<Real, BC::Axis>());

  register_dirichlet_functor<BC::Dirichlet::Increment>(
      mod, "Increment", py::init<Vector<Real> &>());

  register_neumann_functor<BC::Neumann::FromHigherDim>(
      mod, "FromHigherDim", py::init<Matrix<Real> &>());

  register_neumann_functor<BC::Neumann::FromSameDim>(
      mod, "FromSameDim", py::init<Vector<Real> &>());

  register_neumann_functor<BC::Neumann::FreeBoundary>(mod, "FreeBoundary",
                                                      py::init());
}
} // namespace akantu
