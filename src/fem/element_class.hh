/**
 * @file   element_class.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jun 16 18:48:25 2010
 *
 * @brief
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

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ELEMENT_CLASS_HH__
#define __AKANTU_ELEMENT_CLASS_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_math.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__


/**
 * Class describing  the different  type of element  for mesh or  finite element
 * purpose
 *
 * @tparam type the element type for the specialization of the element class
 */
template<ElementType type> class ElementClass {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /**
   * compute  the  shape  functions,  the shape  functions derivatives  and  the
   * jacobians
   * @param[in] coord coordinates of the nodes
   * @param[out] shape shape functions [nb_quad*node_per_elem]
   * @param[out] shape_deriv shape functions derivatives [nb_quad*node_per_elem*spatial_dim]
   * @param[out] jacobian  jacobians * integration weights [nb_quad]
   */
  inline static void preComputeStandards(const Real * coord,
					 const UInt dimension,
					 Real * shape,
					 Real * shape_deriv,
					 Real * jacobian);
  /// compute the shape values for a point given in natural coordinates
  inline static void computeShapes(const Real * natural_coords, Real * shapes);


  /// compute the shape values for a set of points given in natural coordinates
  inline static void computeShapes(const Real * natural_coords,
				   const UInt nb_points,
				   Real * shapes);

  inline static void computeShapes(const Real * natural_coords,
				   const UInt nb_points,
				   Real * shapes,
				   const Real * local_coord,
				   UInt id = 0);
  inline static void computeShapes(const Real * natural_coords,
				   Real * shapes,
				   const Real * local_coord,
				   UInt id = 0);

  /**
   * compute dxds the variation of real coordinates along with
   * variation of natural coordinates on a given point in natural
   * coordinates
   */
  inline static void computeDXDS(const Real * dnds,
				 const Real * node_coords,
				 const UInt dimension,
				 Real * dxds);

  /**
   * compute dxds the variation of real coordinates along with
   * variation of natural coordinates on a given set of points in
   * natural coordinates
   */
  inline static void computeDXDS(const Real * dnds,
				 const UInt nb_points,
				 const Real * node_coords,
				 const UInt dimension, Real * dxds);

  /**
   * compute dnds the variation of real shape functions along with
   * variation of natural coordinates on a given point in natural
   * coordinates
   */
  inline static void computeDNDS(const Real * natural_coords,
				 Real * dnds);

  /**
   * compute dnds the variation of shape functions along with
   * variation of natural coordinates on a given set of points in
   * natural coordinates
   */
  inline static void computeDNDS(const Real * natural_coords,
				 const UInt nb_points,
				 Real * dnds);


  /// compute jacobian (or integration variable change factor) for a set of points
  inline static void computeJacobian(const Real * dxds,
				     const UInt nb_points,
				     const UInt dimension,
				     Real * jac);

  /// compute jacobian (or integration variable change factor) for a given point
  inline static void computeJacobian(const Real * dxds,
				     const UInt dimension,
				     Real & jac);

  /// compute shape derivatives (input is dxds) for a set of points
  inline static void computeShapeDerivatives(const Real * dxds,
					     const Real * dnds,
					     const UInt nb_points,
					     const UInt dimension,
					     Real * shape_deriv);
  /// compute shape derivatives (input is dxds) for a given point
  inline static void computeShapeDerivatives(const Real * dxds,
					     const Real * dnds,
					     Real * shape_deriv);

  inline static void computeShapeDerivatives(const Real * natural_coords,
					     const UInt  nb_points,
					     const UInt dimension,
					     Real * shape_deriv,
					     const Real * local_coord,
					     UInt id = 0);

  inline static void computeShapeDerivatives(const Real * natural_coords,
					     Real * shape_deriv,
					     const Real * local_coord,
					     UInt id);

  /// compute normals on quad points
  inline static void computeNormalsOnQuadPoint(const Real * dxds,
					       const UInt dimension,
					       Real * normals);


  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {};

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                 */
  /* ------------------------------------------------------------------------ */
public:

  static AKANTU_GET_MACRO_NOT_CONST(NbNodesPerElement, nb_nodes_per_element, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(P1ElementType, p1_element_type, ElementType);
  static AKANTU_GET_MACRO_NOT_CONST(NbQuadraturePoints, nb_quadrature_points, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(SpatialDimension, spatial_dimension, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(FacetElementType, facet_type, const ElementType &);
  static AKANTU_GET_MACRO_NOT_CONST(NbFacetsPerElement, nb_facets, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(FacetLocalConnectivityPerElement, facet_connectivity, UInt**);
  static AKANTU_GET_MACRO_NOT_CONST(NbShapeFunctions, nb_shape_functions, UInt);

  static inline Real * getQuadraturePoints();
  static inline UInt getShapeSize();
  static inline UInt getShapeDerivativesSize();

  /// compute the in-radius
  static inline Real getInradius(const Real * coord);

  static inline Real * getGaussIntegrationWeights();

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// Number of nodes per element
  static UInt nb_nodes_per_element;

  /// Number of quadrature points per element
  static UInt nb_quadrature_points;

  /// Dimension of the element
  static UInt spatial_dimension;

  /// Type of the facet elements
  static ElementType facet_type;

  /// number of facets for element
  static UInt nb_facets;

  /// local connectivity of facets
  static UInt * facet_connectivity[];

  /// vectorial connectivity of facets
  static UInt vec_facet_connectivity[];

  /// type of element P1 associated
  static ElementType p1_element_type;

  /// quadrature points in natural coordinates
  static Real quad[];

  /// Number of shape functions
  static UInt nb_shape_functions;

  /// weights for the Gauss integration
  static Real gauss_integration_weights[];
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "element_class_inline_impl.cc"
#endif

/* -------------------------------------------------------------------------- */

__END_AKANTU__


#endif /* __AKANTU_ELEMENT_CLASS_HH__ */
