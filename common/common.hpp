/**
 * @file   common.hpp
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

#ifndef __MYFEM_COMMON__
#define __MYFEM_COMMON__

/* -------------------------------------------------------------------------- */
#include <iostream>
#include <string>
#include <exception>

/* -------------------------------------------------------------------------- */
#include <cstdlib>

/* -------------------------------------------------------------------------- */
#define __BEGIN_MYFEM__ namespace myfem {
#define __END_MYFEM__ };

/* -------------------------------------------------------------------------- */
#include "error.hpp"

/* -------------------------------------------------------------------------- */

__BEGIN_MYFEM__

/// Type Codes
enum TypeCode {
  _integer = 0,
  _unsigned_integer = 1,
  _float = 2,
  _double = 3
};

/// to size of type table (in bytes)
extern int SizeOfType[];

/// to string table
extern const char* TypeCodeName[];

/// to size of type macro
#define MYFEM_SIZEOF(x) (myfem::SizeOfType[(x)])

/// to string macro
#define MYFEM_TYPECODE(x) (myfem::TypeCodeName[(x)])

__END_MYFEM__

/* -------------------------------------------------------------------------- */
#define MYFEM_MIN_ALLOCATION 200

/* -------------------------------------------------------------------------- */
#define MYFEM_SET_MACRO(name, variable, type)	\
  inline void set##name (type variable) {	\
    this->variable = variable;			\
  }

#define MYFEM_GET_MACRO(name, variable, type)	\
  inline type get##name () {			\
    return this->variable;			\
  }


#endif //__MYFEM_COMMON__
