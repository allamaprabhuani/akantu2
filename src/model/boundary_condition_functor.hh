/**
 * @file   boundary_condition_functor.hh
 *
 * @author Dana Christen <dana.christen@gmail.com>
 *
 * @date   Wed Mar 18 11:30:00 2013
 *
 * @brief  XXX
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_BOUNDARY_CONDITION_FUNCTOR_HH__
#define __AKANTU_BOUNDARY_CONDITION_FUNCTOR_HH__

#include "aka_common.hh"
#include "fem.hh"


/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__


namespace BC {

enum Axis {
  _x = 0,
  _y = 1,
  _z = 2,
};

class Functor {

protected:
  Functor() {}

public:
  enum Type {
    _dirichlet,
    _neumann
  };


};


namespace Dirichlet {

/* -------------------------------------------------------------------------- */
class DirichletFunctor : public Functor {

protected:
  DirichletFunctor() : axis(_x) {}
  DirichletFunctor(Axis ax) : axis(ax) {}

public:
  void operator()(UInt node, Vector<bool> & flags, Vector<Real> & primal, const Vector<Real> & coord) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

public:
  static const Type type = _dirichlet;

protected:
  Axis axis;
};

/* -------------------------------------------------------------------------- */
class FlagOnly : public DirichletFunctor {

public:
  FlagOnly(Axis ax = _x) : DirichletFunctor(ax) {}

public:
  inline void operator()(UInt node,
			 Vector<bool> & flags,
			 Vector<Real> & primal,
			 const Vector<Real> & coord) const;

};

/* -------------------------------------------------------------------------- */
class FixedValue : public DirichletFunctor {

public:
  FixedValue(Real val, Axis ax = _x) : DirichletFunctor(ax), value(val) {}

public:
  inline void operator()(UInt node, 
			 Vector<bool> & flags, 
			 Vector<Real> & primal,
			 const Vector<Real> & coord) const;

private:
  Real value;
};

/* -------------------------------------------------------------------------- */
class IncrementValue : public DirichletFunctor {

public:
  IncrementValue(Real val, Axis ax = _x) : DirichletFunctor(ax), value(val) {}

public:
  inline void operator()(UInt node, 
			 Vector<bool> & flags, 
			 Vector<Real> & primal, 
			 const Vector<Real> & coord) const;

  inline void setIncrement(Real val) { this->value = val; }
  
protected:
  Real value;
};
  

} //end namespace Dirichlet

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

namespace Neumann {

/* -------------------------------------------------------------------------- */
class NeumannFunctor : public Functor {

protected:
  NeumannFunctor() {}
public:
  void operator()(UInt node, Vector<bool> & flags, Vector<Real> & primal, const Vector<Real> & coord) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

public:
  static const Type type = _neumann;

};
/* -------------------------------------------------------------------------- */
class FromHigherDim : public NeumannFunctor {

public:
  FromHigherDim(Matrix<Real> mat) : bc_data(mat) {}

public:
  inline void operator()(QuadraturePoint quad_point, Vector<Real> & dual, const Vector<Real> & coord, const Vector<Real> & normals) const;

protected:
  Matrix<Real> bc_data;
};

/* -------------------------------------------------------------------------- */
class FromSameDim : public NeumannFunctor {

public:
  FromSameDim(Vector<Real> vec) : bc_data(vec) {}

public:
  inline void operator()(QuadraturePoint quad_point, Vector<Real> & dual, const Vector<Real> & coord, const Vector<Real> & normals) const;

protected:
  Vector<Real> bc_data;
};

/* -------------------------------------------------------------------------- */
class FreeBoundary : public NeumannFunctor {

public:
  inline void operator()(QuadraturePoint quad_point, Vector<Real> & dual, const Vector<Real> & coord, const Vector<Real> & normals) const;
};

typedef FromHigherDim FromStress;
typedef FromSameDim   FromTraction;

} //end namespace Neumann

} //end namespace BC

__END_AKANTU__

#include "boundary_condition_functor_inline_impl.cc"

#endif /* __AKANTU_BOUNDARY_CONDITION_FUNCTOR_HH__ */

