/**
 * @file   shape_linked.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Apr  1 11:46:44 2011
 *
 * @brief  shape class for element with different set of shapes functions
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

#include "shape_functions.hh"

#ifndef __AKANTU_SHAPE_LINKED_HH__
#define __AKANTU_SHAPE_LINKED_HH__

__BEGIN_AKANTU__

class ShapeLinked : public ShapeFunctions {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  typedef ByElementType<Vector<Real> **> ByElementTypeMultiReal;

  ShapeLinked(Mesh & mesh, const ShapeID & id = "shape_linked");
  virtual ~ShapeLinked(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// pre compute all shapes on the element control points from natural coordinates
  template <ElementType type>
  void precomputeShapesOnControlPoints(const GhostType & ghost_type);

  /// pre compute all shapes on the element control points from natural coordinates
  template <ElementType type>
  void precomputeShapeDerivativesOnControlPoints(const GhostType & ghost_type);

  /// interpolate nodal values on the control points
  template <ElementType type>
  void interpolateOnControlPoints(const Vector<Real> &u,
				  Vector<Real> &uq,
				  UInt nb_degre_of_freedom,
				  const GhostType & ghost_type = _not_ghost,
				  const Vector<UInt> * filter_elements = NULL,
				  bool accumulate = false,
				  UInt id_shape = 0,
				  UInt num_degre_of_freedom_to_interpolate = 0,
				  UInt num_degre_of_freedom_interpolated = 0) const;


  /// compute the gradient of u on the control points
  template <ElementType type>
  void gradientOnControlPoints(const Vector<Real> &u,
			       Vector<Real> &nablauq,
			       UInt nb_degre_of_freedom,
			       const GhostType & ghost_type = _not_ghost,
			       const Vector<UInt> * filter_elements = NULL,
			       bool accumulate = false,
			       UInt id_shape = 0,
			       UInt num_degre_of_freedom_to_interpolate = 0,
			       UInt num_degre_of_freedom_interpolated = 0) const;

  /// multiply a field by shape functions
  template <ElementType type>
  void fieldTimesShapes(const Vector<Real> & field,
			Vector<Real> & fiedl_times_shapes,
			const GhostType & ghost_type) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// get a the shapes vector
  inline const Vector<Real> & getShapes(const ElementType & type,
					const GhostType & ghost_type,
					UInt id = 0) const;

  /// get a the shapes derivatives vector
  inline const Vector<Real> & getShapesDerivatives(const ElementType & type,
						   const GhostType & ghost_type,
						   UInt id = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// shape functions for all elements
  ByElementTypeMultiReal shapes;

  /// shape derivatives for all elements
  ByElementTypeMultiReal shapes_derivatives;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "shape_linked_inline_impl.cc"

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const ShapeLinked & _this)
// {
//   _this.printself(stream);
//   return stream;
// }


__END_AKANTU__

#endif /* __AKANTU_SHAPE_LINKED_HH__ */
