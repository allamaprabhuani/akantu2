/**
 * @file   element_class_structural.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Sat Nov 17 19:32:51 2012
 *
 * @brief  Specialization for structural elements
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
template<InterpolationType interpolation_type>
class InterpolationElement<interpolation_type, _itk_structural> {
public:
  /// compute the shape values for a given set of points in natural coordinates
  static inline void computeShapes(const Matrix<Real> & natural_coord,
				   Matrix<Real> & N,
				   const Matrix<Real> & real_nodal_coord,
				   UInt n = 0) {
    UInt nb_points = natural_coord.cols();
    for (UInt p = 0; p < nb_points; ++p) {
      Vector<Real> Np = N(p);
      computeShapes(natural_coord(p), Np, real_nodal_coord, n);
    }
  }

  /// compute the shape values for a given point in natural coordinates
  static inline void computeShapes(const Vector<Real> & natural_coord,
				   Vector<Real> & N,
				   const Matrix<Real> & real_nodal_coord,
				   UInt n = 0);
  /**
   * compute @f$ B_{ij} = \frac{\partial N_j}{\partial S_i} @f$ the variation of
   * shape functions along with variation of natural coordinates on a given set
   * of points in natural coordinates
   */
  static inline void computeDNDS(const Matrix<Real> & natural_coord,
				 Tensor3<Real> & dnds,
				 const Matrix<Real> & real_nodal_coord,
				 UInt n = 0) {
    for (UInt i = 0; i < natural_coord.cols(); ++i) {
      Matrix<Real> dnds_t = dnds(i);
      computeDNDS(natural_coord(i), dnds_t, real_nodal_coord, n);
    }
  }

  /**
   * compute @f$ B_{ij} = \frac{\partial N_j}{\partial S_i} @f$ the variation of shape functions along with
   * variation of natural coordinates on a given point in natural
   * coordinates
   */
  static inline void computeDNDS(const Vector<Real> & natural_coord,
				 Matrix<Real> & dnds,
				 const Matrix<Real> & real_nodal_coord,
				 UInt n = 0);

public:
  static AKANTU_GET_MACRO_NOT_CONST(NbShapeFunctions, nb_shape_functions, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(ShapeSize, nb_nodes_per_element, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(ShapeDerivativesSize, (nb_nodes_per_element * natural_space_dimension), UInt);
  static AKANTU_GET_MACRO_NOT_CONST(NaturalSpaceDimension, natural_space_dimension, UInt);
protected:
  /// nb shape functions
  static UInt nb_shape_functions;
  /// number of nodes per element
  static UInt nb_nodes_per_element;
  /// dimension of the natural space of the element
  static UInt natural_space_dimension;
};

/* -------------------------------------------------------------------------- */
template<ElementType element_type>
class ElementClass<element_type, _ek_structural> :
  public GeometricalElement<ElementClassProperty<element_type>::geometrical_type>,
  public InterpolationElement<ElementClassProperty<element_type>::interpolation_type> {
protected:
  typedef GeometricalElement<ElementClassProperty<element_type>::geometrical_type> geometrical_element;
  typedef InterpolationElement<ElementClassProperty<element_type>::interpolation_type> interpolation_element;
public:
  /// compute shape derivatives (input is dxds) for a set of points
  static inline void computeShapeDerivatives(const Matrix<Real> & natural_coord,
					     Tensor3<Real> & shape_deriv,
					     const Matrix<Real> & real_nodal_coord,
					     UInt n = 0) {
    UInt nb_points = natural_coord.cols();
    for (UInt p = 0; p < nb_points; ++p) {
      Matrix<Real> shape_deriv_p = shape_deriv(p);
      interpolation_element::computeDNDS(natural_coord(p), shape_deriv_p, real_nodal_coord, n);
    }
  }

  /// compute jacobian (or integration variable change factor) for a given point
  static inline void computeJacobian(const Matrix<Real> & natural_coords,
				     const Matrix<Real> & nodal_coords,
				     Vector<Real> & jacobians);
public:
  static AKANTU_GET_MACRO_NOT_CONST(Kind, _ek_structural, ElementKind);
  static AKANTU_GET_MACRO_NOT_CONST(P1ElementType, _not_defined, const ElementType);
  static AKANTU_GET_MACRO_NOT_CONST(FacetType,     _not_defined, const ElementType);
  static AKANTU_GET_MACRO_NOT_CONST(SpatialDimension, ElementClassProperty<element_type>::spatial_dimension, UInt);
};
