/**
 * @file   common.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jun 11 09:48:06 2010
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 * @section DESCRIPTION
 *
 * All common things to be included in the projects files
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COMMON_HH__
#define __AKANTU_COMMON_HH__

/* -------------------------------------------------------------------------- */
#include <iostream>
#include <string>
#include <exception>
#include <map>

/* -------------------------------------------------------------------------- */
#include <cstdlib>

/* -------------------------------------------------------------------------- */
#define __BEGIN_AKANTU__ namespace akantu {
#define __END_AKANTU__ };

/* -------------------------------------------------------------------------- */
#include "error.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* Memory types                                                               */
/* -------------------------------------------------------------------------- */

typedef unsigned int MemoryID;

typedef std::string VectorID;

/* -------------------------------------------------------------------------- */
/* Mesh types                                                                 */
/* -------------------------------------------------------------------------- */

typedef std::string MeshID;

enum ElementType {
  _not_defined  = 0,
  _triangle_1   = 1,
  _triangle_2   = 2,
  _tetrahedra_1 = 3,
  _tetrahedra_2 = 4,
  _max_element_type
};

/* -------------------------------------------------------------------------- */
//! standard output stream operator for ElementType
inline std::ostream & operator <<(std::ostream & stream, ElementType type)
{
  switch(type)
    {
    case _triangle_1   : stream << "triangle 1st order"  ; break;
    case _triangle_2   : stream << "triangle 2nd order"  ; break;
    case _tetrahedra_1 : stream << "tetrahedra 1st order"; break;
    case _tetrahedra_2 : stream << "tetrahedra 2nd order"; break;
    default : stream << "unknown ElementType (" << type << ")"; break;
    }
  return stream;
}

__END_AKANTU__

/* -------------------------------------------------------------------------- */
/* Global defines                                                             */
/* -------------------------------------------------------------------------- */

#define AKANTU_MIN_ALLOCATION 2000

#define AKANTU_INDENT " "

#define MAX_NUMBER_OF_NODE_PER_ELEMENT 10 // tetrahedron of second order

/* -------------------------------------------------------------------------- */
#define AKANTU_SET_MACRO(name, variable, type)	\
  inline void set##name (type variable) {	\
    this->variable = variable;			\
  }

#define AKANTU_GET_MACRO(name, variable, type)	\
  inline const type get##name () const {	\
    return this->variable;			\
  }

#define AKANTU_GET_MACRO_SCALAR(name, variable, type)	\
  inline type get##name () const {	\
    return this->variable;			\
  }


#endif /* __AKANTU_COMMON_HH__ */
