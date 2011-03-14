/**
 * @file   aka_common.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jun 11 09:48:06 2010
 *
 * @namespace akantu
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

static const Real UINT_INIT_VALUE = 0;
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
typedef ID SparseMatrixID;
typedef ID SolverID;
typedef ID ShapeID;
typedef ID IntegratorID;

typedef UInt Surface;


/* -------------------------------------------------------------------------- */
// BOOST PART: TOUCH ONLY IF YOU KNOW WHAT YOU ARE DOING
#include <boost/preprocessor.hpp>

#define AKANTU_BOOST_CASE_MACRO(r,macro,type)	\
  case type : { macro(type); break;}

#define AKANTU_BOOST_ELEMENT_SWITCH(macro)				\
  switch(type) {							\
    BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_CASE_MACRO,macro,AKANTU_ELEMENT_TYPE) \
  case _not_defined:							\
  case _max_element_type:  {						\
    AKANTU_DEBUG_ERROR("Wrong type : " << type);			\
    break;								\
  }									\
  }

#define AKANTU_BOOST_LIST_MACRO(r,macro,type)	\
  macro(type)

#define AKANTU_BOOST_ELEMENT_LIST(macro)				\
  BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_LIST_MACRO,macro,AKANTU_ELEMENT_TYPE)

/* -------------------------------------------------------------------------- */

/// @boost sequence of element to loop on in global tasks
#define AKANTU_ELEMENT_TYPE			\
  (_segment_2)					\
  (_segment_3)					\
  (_triangle_3)					\
  (_triangle_6)					\
  (_tetrahedron_4)				\
  (_tetrahedron_10)				\
  (_quadrangle_4)				\
  (_hexahedron_8)                               \
  (_point)


/// @enum ElementType type of element potentially contained in a Mesh
enum ElementType {
  _not_defined     = 0,
  _segment_2       = 1, /// first  order segment
  _segment_3       = 2, /// second order segment
  _triangle_3      = 3, /// first  order triangle
  _triangle_6      = 4, /// second order triangle
  _tetrahedron_4   = 5, /// first  order tetrahedron
  _tetrahedron_10  = 6, /// second order tetrahedron @remark not implemented yet
  _quadrangle_4,        /// first  order quadrangle
  _hexahedron_8,        /// first  order hexahedron
  _point,               /// point only for some algorithm to be generic like mesh partitioning
  _max_element_type
};

/// @enum MaterialType different materials implemented
enum MaterialType {
  _elastic = 0,
  _max_material_type
};


typedef void (*BoundaryFunction)(double *,double *);

/// @enum BoundaryFunctionType type of function passed for boundary conditions
enum BoundaryFunctionType {
  _bft_stress,
  _bft_forces
};

/// @enum SparseMatrixType type of sparse matrix used
enum SparseMatrixType {
  _unsymmetric,
  _symmetric
};

/* -------------------------------------------------------------------------- */
/* Contact                                                                    */
/* -------------------------------------------------------------------------- */
typedef ID ContactID;
typedef ID ContactSearchID;
typedef ID ContactNeighborStructureID;

enum ContactType {
  _ct_not_defined = 0,
  _ct_2d_expli    = 1,
  _ct_3d_expli    = 2,
  _ct_rigid       = 3
};

enum ContactSearchType {
  _cst_not_defined = 0,
  _cst_2d_expli    = 1,
  _cst_expli       = 2
};

enum ContactNeighborStructureType {
  _cnst_not_defined  = 0,
  _cnst_regular_grid = 1,
  _cnst_2d_grid = 2
};

/* -------------------------------------------------------------------------- */
/* Ghosts handling                                                            */
/* -------------------------------------------------------------------------- */

typedef ID SynchronizerID;

/// @enum CommunicatorType type of communication method to use
enum CommunicatorType {
  _communicator_mpi,
  _communicator_dummy
};

/// @enum GhostSynchronizationTag type of synchronizations
enum GhostSynchronizationTag {
  /// SolidMechanicsModel tags
  _gst_smm_mass,      /// synchronization of the SolidMechanicsModel.mass
  _gst_smm_for_strain,  /// synchronization of the SolidMechanicsModel.current_position
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
    case _segment_2        : stream << "segment_2"     ; break;
    case _segment_3        : stream << "segment_3"     ; break;
    case _triangle_3       : stream << "triangle_3"    ; break;
    case _triangle_6       : stream << "triangle_6"    ; break;
    case _tetrahedron_4    : stream << "tetrahedron_4" ; break;
    case _tetrahedron_10   : stream << "tetrahedron_10"; break;
    case _quadrangle_4     : stream << "quadrangle_4"  ; break;
    case _hexahedron_8     : stream << "hexahedron_8"  ; break;
    case _not_defined      : stream << "undefined"     ; break;
    case _max_element_type : stream << "ElementType(" << (int) type << ")"; break;
    case _point            : stream << "point"; break;
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
