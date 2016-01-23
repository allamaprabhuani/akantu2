/**
 * @file   boundary_condition_python_functor.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Wed Aug 31 2011
 * @date last modification: Fri Nov 13 2015
 *
 * @brief  Interface for BC::Functor written in python
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "boundary_condition_python_functor.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__


namespace BC {



  void PythonFunctorDirichlet::operator ()(UInt node,
					   Vector<bool> & flags,
					   Vector<Real> & primal,
					   const Vector<Real> & coord) const{

    this->callFunctor<void>("operator",node,flags,primal,coord);
  }


  
  void PythonFunctorNeumann::operator()(const IntegrationPoint & quad_point,
					Vector<Real> & dual,
					const Vector<Real> & coord,
					const Vector<Real> & normals) const{

    this->callFunctor<void>("operator",quad_point,dual,coord,normals);
  }

  
}//end namespace BC


__END_AKANTU__

