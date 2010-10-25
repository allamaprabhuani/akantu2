/**
 * @file   aka_common.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jun 11 09:48:06 2010
 *
 * @namespace akantu
 *
 * @section LICENSE
 *
 * \<insert license here\>
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
#include <list>
#include <limits>
#include <algorithm>

/* -------------------------------------------------------------------------- */
#define __BEGIN_AKANTU__ namespace akantu {
#define __END_AKANTU__ }

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

#ifdef AKANTU_NDEBUG
  static const Real REAL_INIT_VALUE = 0;
#else
  static const Real REAL_INIT_VALUE = std::numeric_limits<Real>::quiet_NaN();
#endif


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

typedef UInt Surface;

enum ElementType {
  _not_defined  = 0,
  _line_1       = 1, // implemented
  _line_2       = 2, // implemented
  _triangle_1   = 3, // implemented
  _triangle_2   = 4, // implemented
  _tetrahedra_1 = 5, // implemented
  _tetrahedra_2 = 6,
  _max_element_type,
  _point
};

enum MaterialType {
  _elastic = 0,
  _max_material_type
};

/* -------------------------------------------------------------------------- */
/* Contact                                                                    */
/* -------------------------------------------------------------------------- */
typedef ID ContactID;
typedef ID ContactSearchID;
typedef ID ContactNeighborStructureID;

enum ContactType {
  _ct_not_defined = 0,
  _ct_2d_expli = 1
};

enum ContactSearchType {
  _cst_not_defined = 0,
  _cst_2d_expli = 1
};

enum ContactNeighborStructureType {
  _cnst_not_defined = 0,
  _cnst_regular_grid
};

/* -------------------------------------------------------------------------- */
/* Ghosts handling                                                            */
/* -------------------------------------------------------------------------- */

typedef ID SynchronizerID;

/// @CommunicatorType type of communication method to use
enum CommunicatorType {
  _communicator_mpi,
  _communicator_dummy
};

/// @enum GhostSynchronizationTag type of synchronizations
enum GhostSynchronizationTag {
  /// SolidMechanicsModel tags
  _gst_smm_mass,      /// synchronization of the SolidMechanicsModel.mass
  _gst_smm_residual,  /// synchronization of the SolidMechanicsModel.current_position
  _gst_smm_boundary,  /// synchronization of the boundary, forces, velocities and displacement
  /// Test tag
  _gst_test
};

/// @enum GhostType type of ghost
enum GhostType {
  _not_ghost,
  _ghost,
  _casper  // not used but a real cute ghost
};

/// @enum SynchronizerOperation reduce operation that the synchronizer can perform
enum SynchronizerOperation {
  _so_sum,
  _so_min,
  _so_max,
  _so_null
};

/* -------------------------------------------------------------------------- */
/* Global defines                                                             */
/* -------------------------------------------------------------------------- */

#define AKANTU_MIN_ALLOCATION 2000

#define AKANTU_INDENT " "

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

#define AKANTU_GET_MACRO_BY_ELEMENT_TYPE(name, variable, type)	\
  inline type get##name (const ::akantu::ElementType & el_type) const {	\
    AKANTU_DEBUG_IN();							\
    AKANTU_DEBUG_ASSERT(variable[el_type] != NULL,			\
			"No " << #variable << " of type "		\
			<< el_type << " in " << id);			\
    AKANTU_DEBUG_OUT();							\
    return *variable[el_type];						\
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
    case _not_defined  : stream << "undefined" ; break;
    case _max_element_type :  stream << "ElementType(" << (int) type << ")"; break;
    case _point        : stream << "point"; break;
    }
  return stream;
}

/* -------------------------------------------------------------------------- */
void initialize(int * argc, char *** argv);

/* -------------------------------------------------------------------------- */
void finalize ();

/* -------------------------------------------------------------------------- */
/*
 * For intel compiler annoying remark
 */
#if defined(__INTEL_COMPILER)
/// remark #981: operands are evaluated in unspecified order
#pragma warning ( disable : 981 )

/// remark #383: value copied to temporary, reference to temporary used
#pragma warning ( disable : 383 )

#endif //defined(__INTEL_COMPILER)

/* -------------------------------------------------------------------------- */
/* string manipulation                                                        */
/* -------------------------------------------------------------------------- */
inline void to_lower(std::string & str) {
  std::transform(str.begin(),
		 str.end(),
		 str.begin(),
		 (int(*)(int))std::tolower);
}

/* -------------------------------------------------------------------------- */
inline void trim(std::string & to_trim) {
  size_t first = to_trim.find_first_not_of(" \t");
  if (first != std::string::npos) {
    size_t last = to_trim.find_last_not_of(" \t");
    to_trim = to_trim.substr(first, last - first + 1);
  } else to_trim = "";
}

__END_AKANTU__

#endif /* __AKANTU_COMMON_HH__ */
