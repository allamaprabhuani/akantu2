/**
 * @file   shape_lagrange.hh
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @date   Thu Feb 10 11:44:29 2011
 *
 * @brief  lagrangian shape functions class
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

#ifndef __AKANTU_SHAPE_LAGRANGE_HH__
#define __AKANTU_SHAPE_LAGRANGE_HH__

#include "shape_functions.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */

class ShapeLagrange : public ShapeFunctions{
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  ShapeLagrange(Mesh & mesh,ShapeID id="shapeLagrange");
  virtual ~ShapeLagrange(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// pre compute all shapes on the element control points from natural coordinates
  template <ElementType type>
  void precomputeShapesOnControlPoints(GhostType ghost_type);

  /// pre compute all shapes on the element control points from natural coordinates
  template <ElementType type>
  void precomputeShapeDerivativesOnControlPoints(GhostType ghost_type);

  /// interpolate nodal values on the control points
  template <ElementType type>
  void interpolateOnControlPoints(const Vector<Real> &u,
				  Vector<Real> &uq,
				  UInt nb_degre_of_freedom,
				  GhostType ghost_type = _not_ghost,
				  const Vector<UInt> * filter_elements = NULL) const;

  /// compute the gradient of u on the control points
  template <ElementType type>
  void gradientOnControlPoints(const Vector<Real> &u,
			       Vector<Real> &nablauq,
			       UInt nb_degre_of_freedom,
			       GhostType ghost_type = _not_ghost,
			       const Vector<UInt> * filter_elements = NULL) const;

  /// set the control points for a given element
  template <ElementType type>
  void setControlPointsByType(Vector<Real> & control_points);

  /// multiply a field by shape functions
  template <ElementType type>
  void fieldTimesShapes(const Vector<Real> & field,
			Vector<Real> & fiedl_times_shapes,
			GhostType ghost_type);


  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  /// compute the shape derivatives on control points for a given element
  template <ElementType type>
  inline void computeShapeDerivativesOnCPointsByElement(UInt spatial_dimension, 
							Real * node_coords,
							UInt nb_nodes_per_element,
							Real * natural_coords,
							UInt nb_points,
							Real * shapesd);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// get a the shapes vector
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Shapes, shapes, const Vector<Real> &);

  /// get a the shapes derivatives vector
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(ShapesDerivatives, shapes_derivatives, const Vector<Real> &);

  /// get a the ghost shapes vector
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(GhostShapes, ghost_shapes,
  				   const Vector<Real> &);

  /// get a the ghost shapes derivatives vector
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(GhostShapesDerivatives, ghost_shapes_derivatives,
  				   const Vector<Real> &);

  /// get the size of the shapes returned by the element class
  static inline UInt getShapeSize(const ElementType & type);

  /// get the size of the shapes derivatives returned by the element class
  static inline UInt getShapeDerivativesSize(const ElementType & type);


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// shape functions for all elements
  ByElementTypeReal shapes;

  /// shape functions for all elements
  ByElementTypeReal ghost_shapes;

  /// shape derivatives for all elements
  ByElementTypeReal shapes_derivatives;

  /// shape derivatives for all elements
  ByElementTypeReal ghost_shapes_derivatives;

  /// shape functions for all elements
  ByElementTypeReal control_points;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "shape_lagrange_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const ShapeLagrange & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_SHAPE_LAGRANGE_HH__ */
