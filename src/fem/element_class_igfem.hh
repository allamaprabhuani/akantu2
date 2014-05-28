/**
 * @file   element_class_igfem.hh
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu May 22 13:50:46 2014
 *
 * @brief  Specialization for interface-enriched finite elements
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
/* Interpolation                                                              */
/* -------------------------------------------------------------------------- */
template<InterpolationType interpolation_type>
class InterpolationElement<interpolation_type, _itk_igfem> {
public:
  typedef InterpolationPorperty<interpolation_type> interpolation_property;

  /// compute the shape values for a given set of points in natural coordinates
  static inline void computeShapes(const Matrix<Real> & natural_coord,
				   Matrix<Real> & N) {
    ElementClass<ElementClassProperty<elem_type>::regular_el_type>::computeShapes(natural_coord, N);
  }

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
				 Tensor3<Real> & dnds) {
    ElementClass<ElementClassProperty<elem_type>::regular_el_type>::computeShapes(natural_coord, N);
  }
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
						     Vector<Real> & interpolated) {
    ElementClass<ElementClassProperty<elem_type>::regular_el_type>::computeShapes(natural_coord, N);
  }


  /// compute the gradient of a given field on the given natural coordinates
  static inline void gradientOnNaturalCoordinates(const Vector<Real> & natural_coords,
						  const Matrix<Real> & f,
						  Matrix<Real> & gradient) {
    ElementClass<ElementClassProperty<elem_type>::regular_el_type>::computeShapes(natural_coord, N);
  }

public:
  static AKANTU_GET_MACRO_NOT_CONST(ShapeSize, InterpolationPorperty<interpolation_type>::nb_nodes_per_element, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(ShapeDerivativesSize, (InterpolationPorperty<interpolation_type>::nb_nodes_per_element * InterpolationPorperty<interpolation_type>::natural_space_dimension), UInt);
  static AKANTU_GET_MACRO_NOT_CONST(NaturalSpaceDimension, InterpolationPorperty<interpolation_type>::natural_space_dimension, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(NbNodesPerInterpolationElement,  InterpolationPorperty<interpolation_type>::nb_nodes_per_element, UInt);
};

/* -------------------------------------------------------------------------- */
/// define ElementClassProperty for parent elements
#define AKANTU_DEFINE_IGFEM_PARENT_ELEMENT_CLASS_PROPERTY(elem_type,		\
						   geom_type,		\
						   interp_type,		\
						   regular_el_type,	\
						   parent_el_type,	\
						   elem_kind,		\
						   sp,			\
						   gauss_int_type,	\
						   min_int_order)	\
  template<>								\
  struct ElementClassProperty<elem_type> {				\
    static const GeometricalType geometrical_type = geom_type;		\
    static const InterpolationType interpolation_type = interp_type;	\
    static const ElementType regular_element_type = regular_el_type;	\
    static const ElementKind element_kind = elem_kind;			\
    static const UInt spatial_dimension = sp;				\
    static const GaussIntergrationType gauss_integration_type = gauss_int_type;	\
    static const UInt minimal_integration_order = min_int_order;	\
    static const bool is_subelement = false;				\
  }


/// define ElementClassProperty for subelements
#define AKANTU_DEFINE_IGFEM_SUBELEMENT_CLASS_PROPERTY(elem_type,	\
						      geom_type,	\
						      interp_type,	\
						      regular_el_type,	\
						      parent_el_type,	\
						      elem_kind,	\
						      sp,		\
						      gauss_int_type,	\
						      min_int_order)	\
  template<>								\
  struct ElementClassProperty<elem_type> {				\
    static const GeometricalType geometrical_type = geom_type;		\
    static const InterpolationType interpolation_type = interp_type;	\
    static const ElementType regular_element_type = regular_el_type;	\
    static const ElementType parent_element_type = parent_el_type;      \
    static const ElementKind element_kind = elem_kind;			\
    static const UInt spatial_dimension = sp;				\
    static const GaussIntergrationType gauss_integration_type = gauss_int_type;	\
    static const UInt minimal_integration_order = min_int_order;	\
    static const bool is_subelement = true;				\
  }

/* -------------------------------------------------------------------------- */
/* ElementClass                                                               */
/* -------------------------------------------------------------------------- */
template<ElementType element_type>
class ElementClass<element_type, _ek_igfem> :
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
				 Matrix<Real> & J) {
    ElementClass<ElementClassProperty<elem_type>::regular_el_type>::computeShapes(natural_coord, N);
  }

  /**
   * compute the Jacobian matrix by computing the variation of real coordinates
   * along with variation of natural coordinates on a given set of points in
   * natural coordinates
   */
  static inline void computeJMat(const Tensor3<Real> & dnds,
				 const Matrix<Real> & node_coords,
				 Tensor3<Real> & J) {
    ElementClass<ElementClassProperty<elem_type>::regular_el_type>::computeShapes(natural_coord, N);
  }

  /// compute the jacobians of a serie of natural coordinates
  static inline void computeJacobian(const Matrix<Real> & natural_coords,
				     const Matrix<Real> & node_coords,
				     Vector<Real> & jacobians) {
    ElementClass<ElementClassProperty<elem_type>::regular_el_type>::computeShapes(natural_coord, N);
  }

  /// compute jacobian (or integration variable change factor) for a set of points
  static inline void computeJacobian(const Tensor3<Real> & J,
				     Vector<Real> & jacobians) {
    ElementClass<ElementClassProperty<elem_type>::regular_el_type>::computeShapes(natural_coord, N);
  }

  /// compute jacobian (or integration variable change factor) for a given point
  static inline void computeJacobian(const Matrix<Real> & J,
				     Real & jacobians) {
    ElementClass<ElementClassProperty<elem_type>::regular_el_type>::computeShapes(natural_coord, N);
  }

  /// compute shape derivatives (input is dxds) for a set of points
  static inline void computeShapeDerivatives(const Tensor3<Real> & J,
					     const Tensor3<Real> & dnds,
					     Tensor3<Real> & shape_deriv) {
    ElementClass<ElementClassProperty<elem_type>::regular_el_type>::computeShapes(natural_coord, N);
  }

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
				Real tolerance = 1e-8) {
    ElementClass<ElementClassProperty<elem_type>::regular_el_type>::computeShapes(natural_coord, N);
  }

public:
  static AKANTU_GET_MACRO_NOT_CONST(Kind, element_kind, ElementKind);
  static AKANTU_GET_MACRO_NOT_CONST(SpatialDimension, ElementClassProperty<element_type>::spatial_dimension, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(P1ElementType, p1_type,    const ElementType &);
  static const ElementType & getFacetType(UInt t = 0) { return facet_type[t]; }
private:
  /// Type of the facet elements
  static ElementType facet_type[];
  /// type of element P1 associated
  static ElementType p1_type;
};


/* -------------------------------------------------------------------------- */
