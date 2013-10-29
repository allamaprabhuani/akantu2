/**
 * @file   boundary_condition.hh
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

#ifndef __AKANTU_BOUNDARY_CONDITION_HH__
#define __AKANTU_BOUNDARY_CONDITION_HH__

#include "aka_common.hh"
#include "boundary_condition_functor.hh"
#include "mesh.hh"
#include "fem.hh"
#include "sub_boundary.hh"
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

template <class ModelType>
class BoundaryCondition {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
private:

  /* ------------------------------------------------------------------------ */
  /* Constructors / Destructors / Initializers                                */
  /* ------------------------------------------------------------------------ */
public:
  BoundaryCondition() : model(NULL), primal(NULL), dual(NULL), primal_increment(NULL) {}

  void initBC(ModelType & ptr, Array<Real> & primal, Array<Real> & dual);
  void initBC(ModelType & ptr, Array<Real> & primal,
	      Array<Real> & primal_increment, Array<Real> & dual);
  /* ------------------------------------------------------------------------ */
  /* Methods and accessors                                                    */
  /* ------------------------------------------------------------------------ */
public:

  //inline void initBoundaryCondition();
  template<typename FunctorType>
  inline void applyBC(const FunctorType & func);

  template<class FunctorType>
  inline void applyBC(const FunctorType & func, const std::string & boundary_name);

  template<class FunctorType>
  inline void applyBC(const FunctorType & func, const SubBoundary & boundary_ref);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

public:

  template<class FunctorType, BC::Functor::Type type = FunctorType::type>
  struct TemplateFunctionWrapper;

private:

  ModelType * model;
  Array<Real> * primal;
  Array<Real> * dual;
  Array<Real> * primal_increment;
};

#include "boundary_condition_tmpl.hh"

__END_AKANTU__


#endif /* __AKANTU_BOUNDARY_CONDITION_HH__ */

