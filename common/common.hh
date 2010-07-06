/**
 * @file   common.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jun 11 09:48:06 2010
 *
 * @section LICENSE
 *
 * <insert lisence here>
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
