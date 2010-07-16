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

#ifndef __AKANTU_ELEMENT_CLASS_HH__
#define __AKANTU_ELEMENT_CLASS_HH__

/* -------------------------------------------------------------------------- */
#include "common.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

template<ElementType type> class ElementClass {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  ElementClass();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /** 
   * compute  the  shape  functions,  the shape  functions derivatives  and  the
   * jacobians
   * @param[in] coord coordinates of the nodes
   * @param[out] shape shape functions [nb_quad * node_per_elem]
   * @param[out] dshape shape functions derivatives []
   * @param[out] jacobian  jacobians * integration weights [nb_quad]
   */
  void shapeFunctions(const Real * coord,
		      Real * shape,
		      Real * dshape,
		      Real * jacobian);


  /// compute the volume of an element
  Real volume(const double * coord);

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {};

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                 */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(NbNodesPerElement, nb_nodes_per_element, UInt);

  AKANTU_GET_MACRO(NbQuadraturePoInts, nb_quadrature_points, UInt);

  AKANTU_GET_MACRO(SpatialDimention, spatial_dimension, UInt);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  UInt nb_nodes_per_element;

  UInt nb_quadrature_points;

  UInt spatial_dimension;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "element_class_inline_impl.cc"

/* -------------------------------------------------------------------------- */

__END_AKANTU__


#endif /* __AKANTU_ELEMENT_CLASS_HH__ */
