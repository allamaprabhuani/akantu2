/**
 * @file   aka_contact_common.hh
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @date   Wed Mar 14 11:16:00 2012
 *
 * @brief  Forward declarations for contact classes
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

#ifndef __AKANTU_CONTACT_COMMON_HH__
#define __AKANTU_CONTACT_COMMON_HH__

#include <iostream>

using std::is_arithmetic;
#include <stdexcept>
#include <array/expr.hpp>

#include "aka_visitor.hh"

//#define DEBUG_CONTACT 1


__BEGIN_AKANTU__

using std::cout;
using std::endl;


typedef array::Array<1,Real> vector_type;
typedef array::Array<2,Real> matrix_type;

using array::transpose;


//! Discretization types
enum Discretization_type { N2N_t, N2S_t };

//! Scheme types
enum Scheme_type { Implicit_t, Explicit_t };

//! Contact type
enum Contact_type { Self_contact_t, No_self_contact_t };

struct EmptyType {};
class NullType {};

struct Command_args {
  static int argc_;
  static char** argv_;
};




//! This functor is called when the visitor is not implemented for a particular object.
template <class R>
struct Discretization_visitor_default {
    template <class T> R operator()(T& t)
    { cout<<"*** WARNING *** No implementation for discretization visitor"<<endl; }
};

template <int>
class Contact_discretization;

//! Discretization type list
typedef MakeTypelist <Contact_discretization<N2N_t> >::Result Contact_discretization_list;

//! Discretization visitor.
typedef Visitor<void, Contact_discretization_list, Mutable, Discretization_visitor_default> Contact_discretization_visitor;



template <typename value_type = Real>
array::Array<2,value_type> eye(UInt d) {
  array::Array<2,value_type> I(d);
  for (UInt i=0; i<d; ++i)
    I(i,i) = 1.;
  return I;
}


__END_AKANTU__


#endif /* __AKANTU_CONTACT_COMMON_HH__ */
