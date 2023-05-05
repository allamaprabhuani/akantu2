/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "fe_engine.hh"
#include "integration_point.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_BOUNDARY_CONDITION_FUNCTOR_HH_
#define AKANTU_BOUNDARY_CONDITION_FUNCTOR_HH_

/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
namespace BC {
  using Axis = ::akantu::SpatialDirection;

  /* ---------------------------------------------------------------------- */
  struct Functor {
    enum Type { _dirichlet, _neumann };
    virtual ~Functor() = default;
  };

  /* ---------------------------------------------------------------------- */
  namespace Dirichlet {

    class DirichletFunctor : public Functor {
    public:
      DirichletFunctor() = default;
      explicit DirichletFunctor(Axis ax) : axis(ax) {}

      virtual void operator()(Idx node, VectorProxy<bool> & flags,
                              VectorProxy<Real> & primal,
                              const VectorProxy<const Real> & coord) {
        // Implementation with copies for backward compatibility
        Vector<bool> flags_ = flags;
        Vector<Real> primal_ = primal;
        Vector<Real> coord_ = coord;

        (*this)(node, flags_, primal_, coord_);

        flags = flags_;
        primal = primal_;
      }

      virtual void operator()(Idx /*node*/, Vector<bool> & /*flags*/,
                              Vector<Real> & /*primal*/,
                              const Vector<Real> & /*coord*/) {
        AKANTU_TO_IMPLEMENT();
      }

    public:
      static const Type type = _dirichlet;

    protected:
      Axis axis{_x};
    };

    /* ---------------------------------------------------------------------- */
    class FlagOnly : public DirichletFunctor {
    public:
      explicit FlagOnly(Axis ax = _x) : DirichletFunctor(ax) {}

    public:
      inline void operator()(Idx node, VectorProxy<bool> & flags,
                             VectorProxy<Real> & primal,
                             const VectorProxy<const Real> & coord) override;
    };

    /* ---------------------------------------------------------------------- */
    class FixedValue : public DirichletFunctor {
    public:
      FixedValue(Real val, Axis ax = _x) : DirichletFunctor(ax), value(val) {}

    public:
      inline void operator()(Idx node, VectorProxy<bool> & flags,
                             VectorProxy<Real> & primal,
                             const VectorProxy<const Real> & coord) override;

    protected:
      Real value;
    };

    /* ---------------------------------------------------------------------- */
    class IncrementValue : public DirichletFunctor {
    public:
      IncrementValue(Real val, Axis ax = _x)
          : DirichletFunctor(ax), value(val) {}

    public:
      inline void operator()(Idx node, VectorProxy<bool> & flags,
                             VectorProxy<Real> & primal,
                             const VectorProxy<const Real> & coord) override;

      inline void setIncrement(Real val) { this->value = val; }

    protected:
      Real value;
    };

    /* ---------------------------------------------------------------------- */
    class Increment : public DirichletFunctor {
    public:
      explicit Increment(const Vector<Real> & val)
          : DirichletFunctor(_x), value(val) {}

    public:
      inline void operator()(Idx node, VectorProxy<bool> & flags,
                             VectorProxy<Real> & primal,
                             const VectorProxy<const Real> & coord) override;

      inline void setIncrement(const Vector<Real> & val) { this->value = val; }

    protected:
      Vector<Real> value;
    };
  } // namespace Dirichlet

  /* ------------------------------------------------------------------------ */
  /* Neumann                                                                  */
  /* ------------------------------------------------------------------------ */
  namespace Neumann {

    class NeumannFunctor : public Functor {

    protected:
      NeumannFunctor() = default;

    public:
      virtual void operator()(const IntegrationPoint & quad_point,
                              VectorProxy<Real> & dual,
                              const VectorProxy<const Real> & coord,
                              const VectorProxy<const Real> & normals) = 0;

      ~NeumannFunctor() override = default;

    public:
      static const Type type = _neumann;
    };

    /* ---------------------------------------------------------------------- */
    class FromHigherDim : public NeumannFunctor {
    public:
      explicit FromHigherDim(const Matrix<Real> & mat) : bc_data(mat) {}
      ~FromHigherDim() override = default;

    public:
      inline void operator()(const IntegrationPoint & quad_point,
                             VectorProxy<Real> & dual,
                             const VectorProxy<const Real> & coord,
                             const VectorProxy<const Real> & normals) override;

    protected:
      Matrix<Real> bc_data;
    };

    /* ---------------------------------------------------------------------- */
    class FromSameDim : public NeumannFunctor {
    public:
      explicit FromSameDim(const Vector<Real> & vec) : bc_data(vec) {}
      ~FromSameDim() override = default;

    public:
      inline void operator()(const IntegrationPoint & quad_point,
                             VectorProxy<Real> & dual,
                             const VectorProxy<const Real> & coord,
                             const VectorProxy<const Real> & normals) override;

    protected:
      Vector<Real> bc_data;
    };

    /* ---------------------------------------------------------------------- */
    class FreeBoundary : public NeumannFunctor {
    public:
      inline void operator()(const IntegrationPoint & quad_point,
                             VectorProxy<Real> & dual,
                             const VectorProxy<const Real> & coord,
                             const VectorProxy<const Real> & normals) override;
    };
  } // namespace Neumann
} // namespace BC
} // namespace akantu

#include "boundary_condition_functor_inline_impl.hh"

#endif /* AKANTU_BOUNDARY_CONDITION_FUNCTOR_HH_ */
