/**
 * @file   aka_common.hh
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
#include <cmath>
#include <cstdlib>
#include <cstring>

/* -------------------------------------------------------------------------- */
#include <iostream>
#include <iomanip>
#include <string>
#include <exception>
#include <vector>
#include <map>
#include <set>

/* -------------------------------------------------------------------------- */
#define __BEGIN_AKANTU__ namespace akantu {
#define __END_AKANTU__ };

/* -------------------------------------------------------------------------- */
#include "aka_error.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* Common types                                                               */
/* -------------------------------------------------------------------------- */

typedef double Real;
typedef unsigned int UInt;
typedef int Int;

typedef std::string ID;

/* -------------------------------------------------------------------------- */
/* Memory types                                                               */
/* -------------------------------------------------------------------------- */

typedef UInt MemoryID;

typedef ID VectorID;

/* -------------------------------------------------------------------------- */
/* Mesh/FEM/Model types                                                       */
/* -------------------------------------------------------------------------- */

typedef ID MeshID;
typedef ID FEMID;
typedef ID ModelID;
typedef ID MaterialID;

enum ElementType {
  _not_defined  = 0,
  _line_1       = 1,
  _line_2       = 2,
  _triangle_1   = 3,
  _triangle_2   = 4,
  _tetrahedra_1 = 5,
  _tetrahedra_2 = 6,
  _max_element_type
};

enum MaterialType {
  _elastic = 0,
  _max_material_type
};

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
  inline type get##name () const {		\
    return variable;				\
  }

#define AKANTU_GET_MACRO_NOT_CONST(name, variable, type)	\
  inline type get##name () {					\
    return variable;						\
  }


/* -------------------------------------------------------------------------- */
//! standard output stream operator for ElementType
inline std::ostream & operator <<(std::ostream & stream, ElementType type)
{
  switch(type)
    {
    case _line_1       : stream << "line_1"  ; break;
    case _line_2       : stream << "line_2"  ; break;
    case _triangle_1   : stream << "triangle_1"  ; break;
    case _triangle_2   : stream << "triangle_2"  ; break;
    case _tetrahedra_1 : stream << "tetrahedra_1"; break;
    case _tetrahedra_2 : stream << "tetrahedra_2"; break;
    case _not_defined  :
    case _max_element_type :  stream << "unknown ElementType (" << type << ")"; break;
    }
  return stream;
}

/* -------------------------------------------------------------------------- */
void initiakize();

/* -------------------------------------------------------------------------- */
void finalize ();

__END_AKANTU__

#endif /* __AKANTU_COMMON_HH__ */
