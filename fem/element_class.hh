/**
 * @file   element_class.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jun 16 18:48:25 2010
 *
 * @brief
 *
 * @section LICENSE
 *
 * <insert license here>
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
  inline static void shapeFunctions(const Real * coord,
			     Real * shape,
			     Real * shape_deriv,
			     Real * jacobian);

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {};

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                 */
  /* ------------------------------------------------------------------------ */
public:

  static AKANTU_GET_MACRO_NOT_CONST(NbNodesPerElement, nb_nodes_per_element, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(NbNodesPerElementP1, nb_nodes_per_element_p1, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(NbQuadraturePoints, nb_quadrature_points, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(SpatialDimension, spatial_dimension, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(SurfaceElementType, surface_type, const ElementType &);

  static inline UInt getShapeSize();
  static inline UInt getShapeDerivatiesSize();
  static inline UInt getJacobiansSize();

  /// compute the in-radius
  static inline Real getInradius(const Real * coord);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// Number of nodes per element
  static UInt nb_nodes_per_element;

  /// Number of nodes per element
  static UInt nb_nodes_per_element_p1;

  /// Number of quadrature points per element
  static UInt nb_quadrature_points;

  /// Dimension of the element
  static UInt spatial_dimension;

  /// Type of the facet elements
  static ElementType surface_type;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "element_class_inline_impl.cc"

/* -------------------------------------------------------------------------- */

__END_AKANTU__


#endif /* __AKANTU_ELEMENT_CLASS_HH__ */
