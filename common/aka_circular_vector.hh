/**
 * @file   aka_circular_vector.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Oct 17 13:39:39 2011
 *
 * @brief  class of circular vector
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

#ifndef __AKANTU_AKA_CIRCULAR_VECTOR_HH__
#define __AKANTU_AKA_CIRCULAR_VECTOR_HH__

/* -------------------------------------------------------------------------- */
#include <typeinfo>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_vector.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

template<class T>
class CircularVector : protected Vector<T> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef typename Vector<T>::value_type      value_type;
  typedef typename Vector<T>::reference       reference;
  typedef typename Vector<T>::pointer_type    pointer_type;
  typedef typename Vector<T>::const_reference const_reference;

  /// Allocation of a new vector with a default value
  CircularVector(UInt size, UInt nb_component = 1,
		 const_reference value = value_type(), const ID & id = "") : 
    Vector<T>(size, nb_component, value, id), 
    start_position(0), 
    end_position(size-1) {
    AKANTU_DEBUG_IN();
    
    AKANTU_DEBUG_OUT();
  };

  virtual ~CircularVector() {
    AKANTU_DEBUG_IN();
    
    AKANTU_DEBUG_OUT();
  };
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// advance start and end position by one
  __aka_inline__ void makeStep();

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

private:

  /* ------------------------------------------------------------------------ */
  /* Operators                                                                */
  /* ------------------------------------------------------------------------ */
public:
  __aka_inline__ reference operator()(UInt i, UInt j = 0);
  __aka_inline__ const_reference operator()(UInt i, UInt j = 0) const;
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  UInt getSize() const{ return this->size; };
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// indice of first element in this circular vector
  UInt start_position;

  /// indice of last element in this circular vector
  UInt end_position;
};


/* -------------------------------------------------------------------------- */
/* __aka_inline__ functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "aka_circular_vector_inline_impl.cc"
#endif

/// standard output stream operator
template <typename T>
inline std::ostream & operator <<(std::ostream & stream, const CircularVector<T> & _this)
{
  _this.printself(stream);
  return stream;
}



__END_AKANTU__



#endif /* __AKANTU_AKA_CIRCULAR_VECTOR_HH__ */


