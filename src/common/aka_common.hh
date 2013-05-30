/**
 * @file   aka_common.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Jun 14 19:12:20 2010
 *
 * @brief  common type descriptions for akantu
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
#include "aka_safe_enum.hh"
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

// int 2 type construct
template <int d>
struct Int2Type {
  enum { result = d };
};

// type 2 type construct
template <class T>
class Type2Type {
  typedef T OriginalType;
};


/* -------------------------------------------------------------------------- */
/* Memory types                                                               */
/* -------------------------------------------------------------------------- */

typedef UInt MemoryID;

/* -------------------------------------------------------------------------- */
/* Mesh/FEM/Model types                                                       */
/* -------------------------------------------------------------------------- */

typedef std::string Surface;
typedef std::pair<Surface, Surface> SurfacePair;
typedef std::list< SurfacePair > SurfacePairList;

/* -------------------------------------------------------------------------- */

extern const UInt _all_dimensions;

/// @boost sequence of element to loop on in global tasks
#define AKANTU_REGULAR_ELEMENT_TYPE		\
  (_point_1)					\
  (_segment_2)					\
  (_segment_3)					\
  (_triangle_3)					\
  (_triangle_6)					\
  (_quadrangle_4)				\
  (_quadrangle_8)				\
  (_tetrahedron_4)				\
  (_tetrahedron_10)				\
  (_hexahedron_8)

#if defined(AKANTU_STRUCTURAL_MECHANICS)
#define AKANTU_STRUCTURAL_ELEMENT_TYPE		\
  (_bernoulli_beam_2)				\
  (_bernoulli_beam_3)
#else
#define AKANTU_STRUCTURAL_ELEMENT_TYPE
#endif

#if defined(AKANTU_COHESIVE_ELEMENT)
#  define AKANTU_COHESIVE_ELEMENT_TYPE		\
  (_cohesive_2d_4)				\
  (_cohesive_2d_6)				\
  (_cohesive_1d_2)
#else
#  define AKANTU_COHESIVE_ELEMENT_TYPE
#endif

#define AKANTU_ALL_ELEMENT_TYPE					\
  AKANTU_REGULAR_ELEMENT_TYPE					\
  AKANTU_COHESIVE_ELEMENT_TYPE					\
  AKANTU_STRUCTURAL_ELEMENT_TYPE

#define AKANTU_NOT_STRUCTURAL_ELEMENT_TYPE			\
  AKANTU_REGULAR_ELEMENT_TYPE					\
  AKANTU_COHESIVE_ELEMENT_TYPE

/// @enum ElementType type of elements
enum ElementType {
  _not_defined,
  _point_1,
  _segment_2,         ///< first order segment
  _segment_3,         ///< second order segment
  _triangle_3,        ///< first order triangle
  _triangle_6,        ///< second order triangle
  _tetrahedron_4,     ///< first order tetrahedron
  _tetrahedron_10,    ///< second order tetrahedron
  _quadrangle_4,      ///< first order quadrangle
  _quadrangle_8,      ///< second order quadrangle
  _hexahedron_8,      ///< first order hexahedron
  _bernoulli_beam_2,  ///< Bernoulli beam 2D
  _bernoulli_beam_3,  ///< Bernoulli beam 3D
#if defined(AKANTU_COHESIVE_ELEMENT)
  _cohesive_2d_4,     ///< first order 2D cohesive
  _cohesive_2d_6,     ///< second order 2D cohesive
  _cohesive_1d_2,     ///< first order 1D cohesive
#endif
  _max_element_type
};

/// @enum GeometricalType type of element potentially contained in a Mesh
enum GeometricalType {
  _gt_point,             ///< point @remark only for some algorithm to be generic like mesh partitioning */
  _gt_segment_2,         ///< 2 nodes segment
  _gt_segment_3,         ///< 3 nodes segment
  _gt_triangle_3,        ///< 3 nodes triangle
  _gt_triangle_6,        ///< 6 nodes triangle
  _gt_quadrangle_4,      ///< 4 nodes quadrangle
  _gt_quadrangle_8,      ///< 8 nodes quadrangle
  _gt_tetrahedron_4,     ///< 4 nodes tetrahedron
  _gt_tetrahedron_10,    ///< 10 nodes tetrahedron
  _gt_hexahedron_8,      ///< 8 nodes hexahedron
#if defined(AKANTU_COHESIVE_ELEMENT)
  _gt_cohesive_2d_4,     ///< 4 nodes 2D cohesive
  _gt_cohesive_2d_6,     ///< 6 nodes 2D cohesive
  _gt_cohesive_1d_2,     ///< 2 nodes 1D cohesive
#endif
  _gt_not_defined
};

/// @enum InterpolationType type of elements
enum InterpolationType {
  _itp_lagrange_point_1,           ///< zeroth (!) order lagrangian point (for compatibility purposes)
  _itp_lagrange_segment_2,         ///< first order lagrangian segment
  _itp_lagrange_segment_3,         ///< second order lagrangian segment
  _itp_lagrange_triangle_3,        ///< first order lagrangian triangle
  _itp_lagrange_triangle_6,        ///< second order lagrangian triangle
  _itp_lagrange_quadrangle_4,      ///< first order lagrangian quadrangle
  _itp_serendip_quadrangle_8,      /**< second order serendipian quadrangle
				      @remark used insted of the 9 node
				      lagrangian element */
  _itp_lagrange_tetrahedron_4,     ///< first order lagrangian tetrahedron
  _itp_lagrange_tetrahedron_10,    ///< second order lagrangian tetrahedron
  _itp_lagrange_hexahedron_8,      ///< first order lagrangian hexahedron
  _itp_bernoulli_beam,             ///< Bernoulli beam
  _itp_not_defined
};

/// @enum InterpolationKind the familly of interpolation types
enum InterpolationKind {
  _itk_lagrangian,
  _itk_structural
};


//! standard output stream operator for ElementType
inline std::ostream & operator <<(std::ostream & stream, ElementType type);


enum ElementKind {
  _ek_regular,
  _ek_cohesive,
  _ek_structural,
  _ek_not_defined
};


/// enum MeshIOType type of mesh reader/writer
enum MeshIOType {
  _miot_auto,
  _miot_gmsh,
  _miot_diana
};

/// enum AnalysisMethod type of solving method used to solve the equation of motion
enum AnalysisMethod {
  _static,
  _implicit_dynamic,
  _explicit_lumped_mass,
  _explicit_consistent_mass
};

/// myfunction(double * position, double * stress/force, double * normal, unsigned int material_id)
typedef void (*BoundaryFunction)(double *,double *, double*, unsigned int);

/// @enum BoundaryFunctionType type of function passed for boundary conditions
enum BoundaryFunctionType {
  _bft_stress,
  _bft_traction,
  _bft_traction_local
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
  //--- SolidMechanicsModel tags ---
  _gst_smm_mass,         //< synchronization of the SolidMechanicsModel.mass
  _gst_smm_for_strain,   //< synchronization of the SolidMechanicsModel.current_position
  _gst_smm_boundary,     //< synchronization of the boundary, forces, velocities and displacement
  _gst_smm_uv,           //< synchronization of the nodal velocities and displacement
  _gst_smm_res,          //< synchronization of the nodal residual
  _gst_smm_init_mat,     //< synchronization of the data to initialize materials
  _gst_smm_stress,       //< synchronization of the stresses to compute the internal forces
  _gst_smmc_facets,      //< synchronization of facet data to setup facet synch
  _gst_smmc_normals,     //< synchronization of facet normals to setup facet synch
  //--- HeatTransfer tags ---
  _gst_htm_capacity,     //< synchronization of the nodal heat capacity
  _gst_htm_temperature,  //< synchronization of the nodal temperature
  _gst_htm_gradient_temperature,  //< synchronization of the element gradient temperature
  //--- LevelSet tags ---
  /// synchronization of the nodal level set value phi
  _gst_htm_phi,
  /// synchronization of the element gradient phi
  _gst_htm_gradient_phi,
  //--- Material non local ---
  _gst_mnl_for_average,  //< synchronization of data to average in non local material
  _gst_mnl_weight,       //< synchronization of data for the weight computations
  //--- General tags ---
  _gst_test,             //< Test tag
  _gst_material_id       //< synchronization of the material ids
};

/// standard output stream operator for SynchronizationTag
inline std::ostream & operator <<(std::ostream & stream, SynchronizationTag type);

/// @enum GhostType type of ghost
enum GhostType {
  _not_ghost,
  _ghost,
  _casper  // not used but a real cute ghost
};

/* -------------------------------------------------------------------------- */
struct GhostType_def {
  typedef GhostType type;
  static const type _begin_ = _not_ghost;
  static const type _end_   = _casper;
};

typedef safe_enum<GhostType_def> ghost_type_t;

/// standard output stream operator for GhostType
inline std::ostream & operator <<(std::ostream & stream, GhostType type);


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
template<class T>
struct is_scalar {
  enum{ value = false };
};

#define AKANTU_SPECIFY_IS_SCALAR(type)		\
  template<>					\
  struct is_scalar<type> {			\
    enum { value = true };			\
  }


AKANTU_SPECIFY_IS_SCALAR(Real);
AKANTU_SPECIFY_IS_SCALAR(UInt);
AKANTU_SPECIFY_IS_SCALAR(Int);
AKANTU_SPECIFY_IS_SCALAR(bool);

template < typename T1, typename T2 >
struct is_same {
  enum { value = false };      // is_same represents a bool.
};

template < typename T >
struct is_same<T, T> {
  enum { value = true };
};

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

#define AKANTU_GET_MACRO_BY_SUPPORT_TYPE(name, variable, type,		\
					 support, con)			\
  inline con Array<type> &						\
  get##name (const support & el_type,					\
	     const GhostType & ghost_type = _not_ghost) con {		\
    return variable(el_type, ghost_type);				\
  }

#define AKANTU_GET_MACRO_BY_ELEMENT_TYPE(name, variable, type)		\
  AKANTU_GET_MACRO_BY_SUPPORT_TYPE(name, variable, type, ElementType,)
#define AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(name, variable, type)	\
  AKANTU_GET_MACRO_BY_SUPPORT_TYPE(name, variable, type, ElementType, const)

#define AKANTU_GET_MACRO_BY_GEOMETRIE_TYPE(name, variable, type)	\
  AKANTU_GET_MACRO_BY_SUPPORT_TYPE(name, variable, type, GeometricalType,)
#define AKANTU_GET_MACRO_BY_GEOMETRIE_TYPE_CONST(name, variable, type)	\
  AKANTU_GET_MACRO_BY_SUPPORT_TYPE(name, variable, type, GeometricalType, const)

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
inline std::string to_lower(const std::string & str);
/* -------------------------------------------------------------------------- */
inline std::string trim(const std::string & to_trim);

__END_AKANTU__

/* -------------------------------------------------------------------------- */
// BOOST PART: TOUCH ONLY IF YOU KNOW WHAT YOU ARE DOING
#include <boost/preprocessor.hpp>

#define AKANTU_BOOST_CASE_MACRO(r,macro,type)	\
  case type : { macro(type); break; }

#define AKANTU_BOOST_ELEMENT_SWITCH(macro1, list1)			\
  do {									\
    switch(type) {							\
      BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_CASE_MACRO, macro1, list1)	\
    default: {								\
      AKANTU_DEBUG_ERROR("Type (" << UInt(type) << ") not handled by this function"); \
    }									\
    }									\
  } while(0)


#define AKANTU_BOOST_ALL_ELEMENT_SWITCH(macro)				\
  AKANTU_BOOST_ELEMENT_SWITCH(macro,					\
			      AKANTU_ALL_ELEMENT_TYPE)

#define AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(macro)			\
  AKANTU_BOOST_ELEMENT_SWITCH(macro,					\
			      AKANTU_REGULAR_ELEMENT_TYPE)

#define AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(macro)			\
  AKANTU_BOOST_ELEMENT_SWITCH(macro,					\
			      AKANTU_COHESIVE_ELEMENT_TYPE)

#define AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(macro)			\
  AKANTU_BOOST_ELEMENT_SWITCH(macro,					\
			      AKANTU_STRUCTURAL_ELEMENT_TYPE)

#define AKANTU_BOOST_LIST_MACRO(r,macro,type)	\
  macro(type)

#define AKANTU_BOOST_ELEMENT_LIST(macro, list)			\
  BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_LIST_MACRO,macro,list)

#define AKANTU_BOOST_ALL_ELEMENT_LIST(macro)			\
  AKANTU_BOOST_ELEMENT_LIST(macro, AKANTU_ALL_ELEMENT_TYPE)

#define AKANTU_BOOST_REGULAR_ELEMENT_LIST(macro)		\
  AKANTU_BOOST_ELEMENT_LIST(macro, AKANTU_REGULAR_ELEMENT_TYPE)

#define AKANTU_BOOST_STRUCTURAL_ELEMENT_LIST(macro)			\
  AKANTU_BOOST_ELEMENT_LIST(macro, AKANTU_STRUCTURAL_ELEMENT_TYPE)

#define AKANTU_BOOST_COHESIVE_ELEMENT_LIST(macro)			\
  AKANTU_BOOST_ELEMENT_LIST(macro, AKANTU_COHESIVE_ELEMENT_TYPE)

#include "aka_common_inline_impl.cc"


#include "aka_fwd.hh"


#endif /* __AKANTU_COMMON_HH__ */
