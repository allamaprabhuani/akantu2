/**
 * @file   element_class.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Jul 04 2014
 *
 * @brief  Declaration of the ElementClass main class and the
 * Integration and Interpolation elements
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "aka_types.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ELEMENT_CLASS_HH__
#define __AKANTU_ELEMENT_CLASS_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<ElementType element_type>
struct ElementClassProperty {
  static const GeometricalType geometrical_type = _gt_not_defined;
  static const InterpolationType interpolation_type = _itp_not_defined;
  static const ElementKind element_kind = _ek_regular;
  static const UInt spatial_dimension = 0;
  static const GaussIntergrationType gauss_integration_type = _git_not_defined;
  static const UInt minimal_integration_order = 0;
};

#define AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(elem_type, geom_type,	\
					     interp_type,		\
					     elem_kind,			\
					     sp,			\
					     gauss_int_type,		\
					     min_int_order)		\
  template<>								\
  struct ElementClassProperty<elem_type> {				\
    static const GeometricalType geometrical_type = geom_type;		\
    static const InterpolationType interpolation_type = interp_type;	\
    static const ElementKind element_kind = elem_kind;			\
    static const UInt spatial_dimension = sp;				\
    static const GaussIntergrationType gauss_integration_type = gauss_int_type;	\
    static const UInt minimal_integration_order = min_int_order;	\
  }

/* -------------------------------------------------------------------------- */
/* Geometry                                                                   */
/* -------------------------------------------------------------------------- */
template<GeometricalType geometrical_type>
struct GeometricalShape {
  static const GeometricalShapeType shape = _gst_point;
};

template<GeometricalShapeType shape>
struct GeometricalShapeContains {
  template <class vector_type>
  static inline bool contains(const vector_type & coord);
};

#define AKANTU_DEFINE_SHAPE(geom_type, geom_shape)			\
  template<>								\
  struct GeometricalShape<geom_type> {					\
    static const GeometricalShapeType shape = geom_shape;		\
  }


/* -------------------------------------------------------------------------- */
template< GeometricalType geometrical_type,
	  GeometricalShapeType shape = GeometricalShape<geometrical_type>::shape >
class GeometricalElement {
public:
  /// compute the in-radius
  static inline Real getInradius(__attribute__((unused)) const Matrix<Real> & coord) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// true if the natural coordinates are in the element
  template <class vector_type>
  static inline bool contains(const vector_type & coord);
public:
  static AKANTU_GET_MACRO_NOT_CONST(SpatialDimension,   spatial_dimension,    UInt);
  static AKANTU_GET_MACRO_NOT_CONST(NbNodesPerElement,  nb_nodes_per_element, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(NbFacetTypes,  nb_facet_types, UInt);
  static inline UInt getNbFacetsPerElement(UInt t);
  static inline UInt getNbFacetsPerElement();
  static inline const MatrixProxy<UInt> getFacetLocalConnectivityPerElement(UInt t = 0);
protected:
  /// Number of nodes per element
  static UInt nb_nodes_per_element;
  /// spatial dimension of the element
  static UInt spatial_dimension;
  /// number of different facet types
  static UInt nb_facet_types;
  /// number of facets for element
  static UInt nb_facets[];
  /// storage of the facet local connectivity
  static UInt facet_connectivity_vect[];
  /// local connectivity of facets
  static UInt * facet_connectivity[];
private:
  /// Type of the facet elements
  static UInt nb_nodes_per_facet[];
};

/* -------------------------------------------------------------------------- */
/* Interpolation                                                              */
/* -------------------------------------------------------------------------- */
template<InterpolationType interpolation_type>
struct InterpolationPorperty {
  static const InterpolationKind kind = _itk_not_defined;
  static const UInt nb_nodes_per_element = 0;
  static const UInt natural_space_dimension = 0;
};

#define AKANTU_DEFINE_INTERPOLATION_TYPE_PROPERTY(itp_type,		\
						  itp_kind,		\
						  nb_nodes,		\
						  ndim)			\
  template<>								\
  struct InterpolationPorperty<itp_type> {				\
    static const InterpolationKind kind = itp_kind;			\
    static const UInt nb_nodes_per_element = nb_nodes;			\
    static const UInt natural_space_dimension = ndim;			\
  };

#include "interpolation_element_tmpl.hh"

/* -------------------------------------------------------------------------- */
template<InterpolationType interpolation_type,
	 InterpolationKind kind = InterpolationPorperty<interpolation_type>::kind>
class InterpolationElement {
public:
  typedef InterpolationPorperty<interpolation_type> interpolation_property;

  /// compute the shape values for a given set of points in natural coordinates
  static inline void computeShapes(const Matrix<Real> & natural_coord,
				   Matrix<Real> & N);

  /// compute the shape values for a given point in natural coordinates
  template <class vector_type>
  static inline void computeShapes(__attribute__((unused)) const vector_type & natural_coord,
				   __attribute__((unused)) vector_type & N) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /**
   * compute @f$ B_{ij} = \frac{\partial N_j}{\partial S_i} @f$ the variation of
   * shape functions along with variation of natural coordinates on a given set
   * of points in natural coordinates
   */
  static inline void computeDNDS(const Matrix<Real> & natural_coord,
				 Tensor3<Real> & dnds);
  /**
   * compute @f$ B_{ij} = \frac{\partial N_j}{\partial S_i} @f$ the variation of shape functions along with
   * variation of natural coordinates on a given point in natural
   * coordinates
   */
  template <class vector_type, class matrix_type>
  static inline void computeDNDS(__attribute__((unused)) const vector_type & natural_coord,
				 __attribute__((unused)) matrix_type & dnds) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// compute jacobian (or integration variable change factor) for a given point
  /// in the case of spatial_dimension != natural_space_dimension
  static inline void computeSpecialJacobian(__attribute__((unused)) const Matrix<Real> & J,
					    __attribute__((unused)) Real & jacobians) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// interpolate a field given (arbitrary) natural coordinates
  static inline void interpolateOnNaturalCoordinates(const Vector<Real> & natural_coords,
						     const Matrix<Real> & nodal_values,
						     Vector<Real> & interpolated);

  /// interpolate a field given the shape functions on the interpolation point
  static inline void interpolate(const Matrix<Real> & nodal_values,
				 const Vector<Real> & shapes,
				 Vector<Real> & interpolated);

  /// interpolate a field given the shape functions on the interpolations points
  static inline void interpolate(const Matrix<Real> & nodal_values,
				 const Matrix<Real> & shapes,
				 Matrix<Real> & interpolated);


  /// compute the gradient of a given field on the given natural coordinates
  static inline void gradientOnNaturalCoordinates(const Vector<Real> & natural_coords,
						  const Matrix<Real> & f,
						  Matrix<Real> & gradient);

public:
  static AKANTU_GET_MACRO_NOT_CONST(ShapeSize, InterpolationPorperty<interpolation_type>::nb_nodes_per_element, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(ShapeDerivativesSize, (InterpolationPorperty<interpolation_type>::nb_nodes_per_element * InterpolationPorperty<interpolation_type>::natural_space_dimension), UInt);
  static AKANTU_GET_MACRO_NOT_CONST(NaturalSpaceDimension, InterpolationPorperty<interpolation_type>::natural_space_dimension, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(NbNodesPerInterpolationElement,  InterpolationPorperty<interpolation_type>::nb_nodes_per_element, UInt);
};

/* -------------------------------------------------------------------------- */
/* Integration                                                                */
/* -------------------------------------------------------------------------- */
template<GaussIntergrationType git_class, UInt max_order>
struct GaussIntegrationTypeData {
  /// quadrature points in natural coordinates
  static Real quad_positions[];
  /// weights for the Gauss integration
  static Real quad_weights[];
  /// Number of quadrature points per element
  static UInt nb_quadrature_points;
};

template<ElementType type, UInt order = ElementClassProperty<type>::minimal_integration_order>
class GaussIntegrationElement {
public:
  static UInt getNbQuadraturePoints();
  static const Matrix<Real> getQuadraturePoints();
  static const Vector<Real> getWeights();
};

/* -------------------------------------------------------------------------- */
/* ElementClass                                                               */
/* -------------------------------------------------------------------------- */
template<ElementType element_type,
	 ElementKind element_kind = ElementClassProperty<element_type>::element_kind>
class ElementClass :
  public GeometricalElement<ElementClassProperty<element_type>::geometrical_type>,
  public InterpolationElement<ElementClassProperty<element_type>::interpolation_type> {
protected:
  typedef GeometricalElement<ElementClassProperty<element_type>::geometrical_type> geometrical_element;
  typedef InterpolationElement<ElementClassProperty<element_type>::interpolation_type> interpolation_element;

  typedef ElementClassProperty<element_type> element_property;
  typedef typename interpolation_element::interpolation_property interpolation_property;
public:
  /**
   * compute @f$ J = \frac{\partial x_j}{\partial s_i} @f$ the variation of real
   * coordinates along with variation of natural coordinates on a given point in
   * natural coordinates
   */
  static inline void computeJMat(const Matrix<Real> & dnds,
				 const Matrix<Real> & node_coords,
				 Matrix<Real> & J);

  /**
   * compute the Jacobian matrix by computing the variation of real coordinates
   * along with variation of natural coordinates on a given set of points in
   * natural coordinates
   */
  static inline void computeJMat(const Tensor3<Real> & dnds,
				 const Matrix<Real> & node_coords,
				 Tensor3<Real> & J);

  /// compute the jacobians of a serie of natural coordinates
  static inline void computeJacobian(const Matrix<Real> & natural_coords,
				     const Matrix<Real> & node_coords,
				     Vector<Real> & jacobians);

  /// compute jacobian (or integration variable change factor) for a set of points
  static inline void computeJacobian(const Tensor3<Real> & J,
				     Vector<Real> & jacobians);

  /// compute jacobian (or integration variable change factor) for a given point
  static inline void computeJacobian(const Matrix<Real> & J,
				     Real & jacobians);

  /// compute shape derivatives (input is dxds) for a set of points
  static inline void computeShapeDerivatives(const Tensor3<Real> & J,
					     const Tensor3<Real> & dnds,
					     Tensor3<Real> & shape_deriv);

  /// compute shape derivatives (input is dxds) for a given point
  static inline void computeShapeDerivatives(const Matrix<Real> & J,
					     const Matrix<Real> & dnds,
					     Matrix<Real> & shape_deriv);

  /// compute the normal of a surface defined by the function f
  static inline void computeNormalsOnNaturalCoordinates(const Matrix<Real> & coord,
							Matrix<Real> & f,
							Matrix<Real> & normals);

  /// get natural coordinates from real coordinates
  static inline void inverseMap(const Vector<Real> & real_coords,
				const Matrix<Real> & node_coords,
				Vector<Real> & natural_coords,
				Real tolerance = 1e-8);

  /// get natural coordinates from real coordinates
  static inline void inverseMap(const Matrix<Real> & real_coords,
				const Matrix<Real> & node_coords,
				Matrix<Real> & natural_coords,
				Real tolerance = 1e-8);
public:
  static AKANTU_GET_MACRO_NOT_CONST(Kind, element_kind, ElementKind);
  static AKANTU_GET_MACRO_NOT_CONST(SpatialDimension, ElementClassProperty<element_type>::spatial_dimension, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(P1ElementType, p1_type,    const ElementType &);
  static const ElementType & getFacetType(UInt t = 0) { return facet_type[t]; }
  static ElementType * getFacetTypeInternal() { return facet_type; }
protected:
  /// Type of the facet elements
  static ElementType facet_type[];
  /// type of element P1 associated
  static ElementType p1_type;
};


/* -------------------------------------------------------------------------- */
#include "element_class_tmpl.hh"

/* -------------------------------------------------------------------------- */
#include "element_classes/element_class_point_1_inline_impl.cc"
#include "element_classes/element_class_segment_2_inline_impl.cc"
#include "element_classes/element_class_segment_3_inline_impl.cc"
#include "element_classes/element_class_triangle_3_inline_impl.cc"
#include "element_classes/element_class_triangle_6_inline_impl.cc"
#include "element_classes/element_class_tetrahedron_4_inline_impl.cc"
#include "element_classes/element_class_tetrahedron_10_inline_impl.cc"
#include "element_classes/element_class_quadrangle_4_inline_impl.cc"
#include "element_classes/element_class_quadrangle_8_inline_impl.cc"
#include "element_classes/element_class_hexahedron_8_inline_impl.cc"
#include "element_classes/element_class_hexahedron_20_inline_impl.cc"
#include "element_classes/element_class_pentahedron_6_inline_impl.cc"
#include "element_classes/element_class_pentahedron_15_inline_impl.cc"

__END_AKANTU__


#if defined(AKANTU_STRUCTURAL_MECHANICS)
#  include "element_class_structural.hh"
#endif

#if defined(AKANTU_IGFEM)
#  include "element_class_igfem.hh"
#endif



#endif /* __AKANTU_ELEMENT_CLASS_HH__ */
