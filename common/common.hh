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

#ifndef __MYFEM_COMMON_HH__
#define __MYFEM_COMMON_HH__

/* -------------------------------------------------------------------------- */
#include <iostream>
#include <string>
#include <exception>
#include <map>

/* -------------------------------------------------------------------------- */
#include <cstdlib>

/* -------------------------------------------------------------------------- */
#define __BEGIN_MYFEM__ namespace myfem {
#define __END_MYFEM__ };

/* -------------------------------------------------------------------------- */
#include "error.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_MYFEM__

/* -------------------------------------------------------------------------- */
/* Memory types                                                               */
/* -------------------------------------------------------------------------- */

typedef unsigned int MemoryID;

typedef std::string VectorID;

__END_MYFEM__

/* -------------------------------------------------------------------------- */
/* Global defines                                                             */
/* -------------------------------------------------------------------------- */

#define MYFEM_MIN_ALLOCATION 2000

#define MYFEM_INDENT " "

#define MAX_NUMBER_OF_NODE_PER_ELEMENT 10 // tetrahedron of second order

/* -------------------------------------------------------------------------- */
#define MYFEM_SET_MACRO(name, variable, type)	\
  inline void set##name (type variable) {	\
    this->variable = variable;			\
  }

#define MYFEM_GET_MACRO(name, variable, type)	\
  inline const type get##name () const {	\
    return this->variable;			\
  }


#endif /* __MYFEM_COMMON_HH__ */
