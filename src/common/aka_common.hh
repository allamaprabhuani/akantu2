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
#include <list>
#include <limits>

/* -------------------------------------------------------------------------- */
#define __BEGIN_AKANTU__ namespace akantu {
#define __END_AKANTU__ }

/* -------------------------------------------------------------------------- */
#include "aka_config.hh"
#include "aka_error.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* Common types                                                               */
/* -------------------------------------------------------------------------- */

typedef double Real;
typedef unsigned int UInt;
typedef unsigned long long UInt64;
typedef signed int Int;

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

/* -------------------------------------------------------------------------- */
/* Mesh/FEM/Model types                                                       */
/* -------------------------------------------------------------------------- */

typedef UInt Surface;
typedef std::pair<Surface, Surface> SurfacePair;
typedef std::list< SurfacePair > SurfacePairList;

/* -------------------------------------------------------------------------- */

/// @boost sequence of element to loop on in global tasks
#define AKANTU_REGULAR_ELEMENT_TYPE		\
  (_not_defined)				\
  (_segment_2)					\
  (_segment_3)					\
  (_triangle_3)					\
  (_triangle_6)					\
  (_tetrahedron_4)				\
  (_tetrahedron_10)				\
  (_quadrangle_4)				\
  (_quadrangle_8)				\
  (_hexahedron_8)                               \
  (_point)					\
  (_bernoulli_beam_2)

#define AKANTU_COHESIVE_ELEMENT_TYPE		\
  (_cohesive_2d_4)				\
  (_cohesive_2d_6)



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
  _quadrangle_8,        /// second order quadrangle
  _hexahedron_8,        /// first  order hexahedron
  _point,               /// point only for some algorithm to be generic like mesh partitioning
  _bernoulli_beam_2,    /// bernoulli beam 2D
  _cohesive_2d_4,       /// first order 2D cohesive
  _cohesive_2d_6,       /// second order 2D cohesive
  _max_element_type
};

enum ElementKind {
  _ek_not_defined,
  _ek_regular,
  _ek_cohesive
};


/// enum AnalysisMethod type of solving method used to solve the equation of motion
enum AnalysisMethod {
  _static,
  _implicit_dynamic,
  _explicit_dynamic
};

/// myfunction(double * position, double * stress/force, double * normal, unsigned int material_id)
typedef void (*BoundaryFunction)(double *,double *, double*, unsigned int);

/// @enum BoundaryFunctionType type of function passed for boundary conditions
enum BoundaryFunctionType {
  _bft_stress,
  _bft_traction
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
/* Friction                                                                   */
/* -------------------------------------------------------------------------- */
typedef ID FrictionID;

/* -------------------------------------------------------------------------- */
/* Ghosts handling                                                            */
/* -------------------------------------------------------------------------- */

typedef ID SynchronizerID;

/// @enum CommunicatorType type of communication method to use
enum CommunicatorType {
  _communicator_mpi,
  _communicator_dummy
};

/// @enum SynchronizationTag type of synchronizations
enum SynchronizationTag {
  /// SolidMechanicsModel tags
  _gst_smm_mass,                  /// synchronization of the SolidMechanicsModel.mass
  _gst_smm_for_strain,            /// synchronization of the SolidMechanicsModel.current_position
  _gst_smm_boundary,              /// synchronization of the boundary, forces, velocities and displacement
  _gst_smm_uv,                    /// synchronization of the nodal velocities and displacement
  _gst_smm_res,                   /// synchronization of the nodal residual
  _gst_smm_init_mat,              /// synchronization of the data to initialize materials
  /// HeatTransfer tags
  _gst_htm_capacity,              /// synchronization of the nodal heat capacity
  _gst_htm_temperature,           /// synchronization of the nodal temperature
  _gst_htm_gradient_temperature,  /// synchronization of the element gradient temperature
  /// Material non local
  _gst_mnl_damage,                /// synchronization of data to average in non local material
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
#define AKANTU_INCLUDE_INLINE_IMPL

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

#define AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(name, variable, type)	\
  inline const Vector<type> &						\
  get##name (const ::akantu::ElementType & el_type,			\
	     const GhostType & ghost_type = _not_ghost) const {		\
    AKANTU_DEBUG_IN();							\
									\
    AKANTU_DEBUG_OUT();							\
    return variable(el_type, ghost_type);				\
  }

#define AKANTU_GET_MACRO_BY_ELEMENT_TYPE(name, variable, type)		\
  inline Vector<type> &							\
  get##name (const ::akantu::ElementType & el_type,			\
	     const GhostType & ghost_type = _not_ghost) {		\
    AKANTU_DEBUG_IN();							\
									\
    AKANTU_DEBUG_OUT();							\
    return variable(el_type, ghost_type);				\
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
    case _quadrangle_8     : stream << "quadrangle_8"  ; break;
    case _hexahedron_8     : stream << "hexahedron_8"  ; break;
    case _bernoulli_beam_2 : stream << "bernoulli_beam_2"; break;
    case _cohesive_2d_4    : stream << "cohesive_2d_4" ; break;
    case _cohesive_2d_6    : stream << "cohesive_2d_6" ; break;
    case _not_defined      : stream << "undefined"     ; break;
    case _max_element_type : stream << "ElementType(" << (int) type << ")"; break;
    case _point            : stream << "point"; break;
    }
  return stream;
}

/// standard output stream operator for GhostType
inline std::ostream & operator <<(std::ostream & stream, GhostType type)
{
  switch(type)
    {
    case _not_ghost : stream << "not_ghost"; break;
    case _ghost     : stream << "ghost"    ; break;
    case _casper    : stream << "Casper the friendly ghost"; break;
    }
  return stream;
}

/* -------------------------------------------------------------------------- */
void initialize(int & argc, char ** & argv);

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
inline void to_lower(std::string & str);
/* -------------------------------------------------------------------------- */
inline void trim(std::string & to_trim);

__END_AKANTU__

#include "aka_common_inline_impl.cc"

/* -------------------------------------------------------------------------- */
// BOOST PART: TOUCH ONLY IF YOU KNOW WHAT YOU ARE DOING
#include <boost/preprocessor.hpp>

#define AKANTU_EXCLUDE_ELEMENT_TYPE(type)			\
  AKANTU_DEBUG_ERROR("Type (" << type << ") not handled by this function")

#define AKANTU_BOOST_CASE_MACRO(r,macro,type)	\
  case type : { macro(type); break; }

#define AKANTU_BOOST_ELEMENT_SWITCH(macro1, list1, macro2, list2)	\
  do {									\
    switch(type) {							\
      BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_CASE_MACRO, macro1, list1)	\
      BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_CASE_MACRO, macro2, list2)	\
    case _max_element_type:  {						\
      AKANTU_DEBUG_ERROR("Wrong type : " << type);			\
      break;								\
    }									\
    }									\
  } while(0)

#define AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(macro)			\
  AKANTU_BOOST_ELEMENT_SWITCH(macro,					\
			      AKANTU_REGULAR_ELEMENT_TYPE,		\
			      AKANTU_EXCLUDE_ELEMENT_TYPE,		\
			      AKANTU_COHESIVE_ELEMENT_TYPE)

#define AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(macro)			\
  AKANTU_BOOST_ELEMENT_SWITCH(macro,					\
			      AKANTU_COHESIVE_ELEMENT_TYPE,		\
			      AKANTU_EXCLUDE_ELEMENT_TYPE,		\
			      AKANTU_REGULAR_ELEMENT_TYPE)

#define AKANTU_BOOST_LIST_MACRO(r,macro,type)	\
  macro(type)

#define AKANTU_BOOST_ELEMENT_LIST(macro,list)			\
  BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_LIST_MACRO,macro,list)

#define AKANTU_BOOST_REGULAR_ELEMENT_LIST(macro)	\
  AKANTU_BOOST_ELEMENT_LIST(macro, AKANTU_REGULAR_ELEMENT_TYPE)

#define AKANTU_BOOST_COHESIVE_ELEMENT_LIST(macro)	\
  AKANTU_BOOST_ELEMENT_LIST(macro, AKANTU_COHESIVE_ELEMENT_TYPE)


#endif /* __AKANTU_COMMON_HH__ */
