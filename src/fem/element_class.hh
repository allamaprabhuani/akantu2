/**
  * @file   element_class.hh
  *
  * @author Nicolas Richart <nicolas.richart@epfl.ch>
  *
  * @date   Wed Nov 14 11:00:21 2012
  *
  * @brief  Desciption of the element classes
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
};

#define AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(elem_type, geom_type, interp_type, elem_kind, sp) \
  template<>								\
  struct ElementClassProperty<elem_type> {				\
    static const GeometricalType geometrical_type = geom_type;		\
    static const InterpolationType interpolation_type = interp_type;	\
    static const ElementKind element_kind = elem_kind;			\
    static const UInt spatial_dimension = sp;				\
  }

/* -------------------------------------------------------------------------- */
/* Geometrie                                                                  */
/* -------------------------------------------------------------------------- */
enum GeometricalShapeType {
  _gst_point,
  _gst_triangle,
  _gst_square
};

template<GeometricalType geometrical_type>
struct GeometricalShape {
  static const GeometricalShapeType shape = _gst_point;
};

template<GeometricalShapeType shape>
struct GeometricalShapeContains {
  static inline bool contains(const types::Vector<Real> & coord);
};

#define AKANTU_DEFINE_SHAPE(geom_type, geom_shape)			\
  template<>								\
  struct GeometricalShape<geom_type> {					\
    static const GeometricalShapeType shape = geom_shape;		\
  }

AKANTU_DEFINE_SHAPE(_gt_point, _gst_point);
AKANTU_DEFINE_SHAPE(_gt_segment_2, _gst_square);
AKANTU_DEFINE_SHAPE(_gt_segment_3, _gst_square);
AKANTU_DEFINE_SHAPE(_gt_quadrangle_4, _gst_square);
AKANTU_DEFINE_SHAPE(_gt_quadrangle_8, _gst_square);
AKANTU_DEFINE_SHAPE(_gt_hexahedron_8, _gst_square);
AKANTU_DEFINE_SHAPE(_gt_triangle_3, _gst_triangle);
AKANTU_DEFINE_SHAPE(_gt_triangle_6, _gst_triangle);
AKANTU_DEFINE_SHAPE(_gt_tetrahedron_4,  _gst_triangle);
AKANTU_DEFINE_SHAPE(_gt_tetrahedron_10, _gst_triangle);

/* -------------------------------------------------------------------------- */
template< GeometricalType geometrical_type,
	  GeometricalShapeType shape = GeometricalShape<geometrical_type>::shape >
class GeometricalElement {
public:
  /// compute the in-radius
  static inline Real getInradius(__attribute__((unused)) const types::Matrix<Real> & coord) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// true if the natural coordinates are in the element
  static inline bool contains(const types::Vector<Real> & coord);
public:
  static AKANTU_GET_MACRO_NOT_CONST(SpatialDimension,   spatial_dimension,    UInt);
  static AKANTU_GET_MACRO_NOT_CONST(NbNodesPerElement,  nb_nodes_per_element, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(NbFacetsPerElement, nb_facets,            UInt);
  static inline const types::Matrix<UInt> getFacetLocalConnectivityPerElement();
protected:
  /// Number of nodes per element
  static UInt nb_nodes_per_element;
  /// spatial dimension of the element
  static UInt spatial_dimension;
  /// number of facets for element
  static UInt nb_facets;
  /// local connectivity of facets
  static UInt facet_connectivity[];
private:
  /// Type of the facet elements
  static UInt nb_nodes_per_facet;
};

/* -------------------------------------------------------------------------- */
/* Interpolation                                                              */
/* -------------------------------------------------------------------------- */
template<InterpolationType interpolation_type>
struct InterpolationPorperty {
  static const InterpolationKind kind = _itk_lagrangian;
};

template<>
struct InterpolationPorperty<_itp_bernoulli_beam> {
  static const InterpolationKind kind = _itk_structural;
};

/* -------------------------------------------------------------------------- */
template<InterpolationType interpolation_type,
	 InterpolationKind kind = InterpolationPorperty<interpolation_type>::kind>
class InterpolationElement {
public:
  /// compute the shape values for a given set of points in natural coordinates
  static inline void computeShapes(const types::Matrix<Real> & natural_coord,
				   types::Matrix<Real> & N);

  /// compute the shape values for a given point in natural coordinates
  static inline void computeShapes(__attribute__((unused)) const types::Vector<Real> & natural_coord,
				   __attribute__((unused)) types::Vector<Real> & N) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /**
   * compute @f$ B_{ij} = \frac{\partial N_j}{\partial S_i} @f$ the variation of
   * shape functions along with variation of natural coordinates on a given set
   * of points in natural coordinates
   */
  static inline void computeDNDS(const types::Matrix<Real> & natural_coord,
				 types::Tensor3<Real> & dnds);
  /**
   * compute @f$ B_{ij} = \frac{\partial N_j}{\partial S_i} @f$ the variation of shape functions along with
   * variation of natural coordinates on a given point in natural
   * coordinates
   */
  static inline void computeDNDS(__attribute__((unused)) const types::Vector<Real> & natural_coord,
				 __attribute__((unused)) types::Matrix<Real> & dnds) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// compute jacobian (or integration variable change factor) for a given point
  /// in the case of spatial_dimension != natural_space_dimension
  static inline void computeSpecialJacobian(__attribute__((unused)) const types::Matrix<Real> & J,
					    __attribute__((unused)) Real & jacobians) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// interpolate a field given (arbitrary) natural coordinates
  static inline void interpolateOnNaturalCoordinates(const types::Vector<Real> & natural_coords,
						     const types::Matrix<Real> & nodal_values,
						     types::Vector<Real> & interpolated);


  /// compute the gradient of a given field on the given natural coordinates
  static inline void gradientOnNaturalCoordinates(const types::Vector<Real> & natural_coords,
						  const types::Matrix<Real> & f,
						  types::Matrix<Real> & gradient);

public:
  static AKANTU_GET_MACRO_NOT_CONST(ShapeSize, nb_nodes_per_element, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(ShapeDerivativesSize, (nb_nodes_per_element * natural_space_dimension), UInt);
  static AKANTU_GET_MACRO_NOT_CONST(NaturalSpaceDimension, natural_space_dimension, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(NbNodesPerInterpolationElement,  nb_nodes_per_element, UInt);
protected:
  /// number of nodes per element
  static UInt nb_nodes_per_element;
  /// dimension of the natural space of the element
  static UInt natural_space_dimension;
};

/* -------------------------------------------------------------------------- */
/* Integration                                                                */
/* -------------------------------------------------------------------------- */
template<ElementType element_type>
class GaussIntegrationElement {
public:
  static AKANTU_GET_MACRO_NOT_CONST(NbQuadraturePoints, nb_quadrature_points, UInt);
  static const types::Matrix<Real> getQuadraturePoints();
  static const types::Vector<Real> getWeights();
private:
  /// quadrature points in natural coordinates
  static Real quad[];
  /// weights for the Gauss integration
  static Real weights[];
  /// Number of quadrature points per element
  static UInt nb_quadrature_points;
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
public:
  /**
   * compute @f$ J = \frac{\partial x_j}{\partial s_i} @f$ the variation of real
   * coordinates along with variation of natural coordinates on a given point in
   * natural coordinates
   */
  static inline void computeJMat(const types::Matrix<Real> & dnds,
				 const types::Matrix<Real> & node_coords,
				 types::Matrix<Real> & J);

  /**
   * compute the Jacobian matrix by computing the variation of real coordinates
   * along with variation of natural coordinates on a given set of points in
   * natural coordinates
   */
  static inline void computeJMat(const types::Tensor3<Real> & dnds,
				 const types::Matrix<Real> & node_coords,
				 types::Tensor3<Real> & J);

  /// compute the jacobians of a serie of natural coordinates
  static inline void computeJacobian(const types::Matrix<Real> & natural_coords,
				     const types::Matrix<Real> & node_coords,
				     types::Vector<Real> & jacobians);

  /// compute jacobian (or integration variable change factor) for a set of points
  static inline void computeJacobian(const types::Tensor3<Real> & J,
				     types::Vector<Real> & jacobians);

  /// compute jacobian (or integration variable change factor) for a given point
  static inline void computeJacobian(const types::Matrix<Real> & J,
				     Real & jacobians);

  /// compute shape derivatives (input is dxds) for a set of points
  static inline void computeShapeDerivatives(const types::Tensor3<Real> & J,
					     const types::Tensor3<Real> & dnds,
					     types::Tensor3<Real> & shape_deriv);

  /// compute shape derivatives (input is dxds) for a given point
  static inline void computeShapeDerivatives(const types::Matrix<Real> & J,
					     const types::Matrix<Real> & dnds,
					     types::Matrix<Real> & shape_deriv);

  /// compute the normal of a surface defined by the function f
  static inline void computeNormalsOnNaturalCoordinates(const types::Matrix<Real> & coord,
							types::Matrix<Real> & f,
							types::Matrix<Real> & normals);

  /// get natural coordinates from real coordinates
  static inline void inverseMap(const types::Vector<Real> & real_coords,
				const types::Matrix<Real> & node_coords,
				types::Vector<Real> & natural_coords,
				Real tolerance = 1e-8);
public:
  static AKANTU_GET_MACRO_NOT_CONST(Kind, element_kind, ElementKind);
  static AKANTU_GET_MACRO_NOT_CONST(SpatialDimension, ElementClassProperty<element_type>::spatial_dimension, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(P1ElementType, p1_type,    const ElementType &);
  static AKANTU_GET_MACRO_NOT_CONST(FacetType,     facet_type, const ElementType &);
private:
  /// Type of the facet elements
  static ElementType facet_type;
  /// type of element P1 associated
  static ElementType p1_type;
};


/* -------------------------------------------------------------------------- */
#include "element_class_tmpl.hh"

/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_point_1, _gt_point, _itp_not_defined, _ek_regular, 0);
#include "element_classes/element_class_segment_2_inline_impl.cc"
#include "element_classes/element_class_segment_3_inline_impl.cc"
#include "element_classes/element_class_triangle_3_inline_impl.cc"
#include "element_classes/element_class_triangle_6_inline_impl.cc"
#include "element_classes/element_class_tetrahedron_4_inline_impl.cc"
#include "element_classes/element_class_tetrahedron_10_inline_impl.cc"
#include "element_classes/element_class_quadrangle_4_inline_impl.cc"
#include "element_classes/element_class_quadrangle_8_inline_impl.cc"
#include "element_classes/element_class_hexahedron_8_inline_impl.cc"

#if defined(AKANTU_STRUCTURAL_MECHANICS)
#  include "element_class_structural.hh"
#  include "element_classes/element_class_bernoulli_beam_inline_impl.cc"
#endif

__END_AKANTU__

#endif /* __AKANTU_ELEMENT_CLASS_HH__ */
