/**
 * @file   shape_functions.hh
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @date   Thu Feb 10 11:35:29 2011
 *
 * @brief shape function class
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

#ifndef __AKANTU_SHAPE_FUNCTIONS_HH__
#define __AKANTU_SHAPE_FUNCTIONS_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
class ShapeFunctions : public Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  ShapeFunctions(Mesh & mesh,
		 const ID & id = "shape",
		 const MemoryID & memory_id = 0) :
    Memory(memory_id), mesh(&mesh),
    control_points("control_points", id, memory_id) {
  };
  virtual ~ShapeFunctions(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// pre compute all shapes on the element control points from natural coordinates
  template <ElementType type>
  void precomputeShapesOnControlPoints(const Real * natural_coords,
				       const UInt nb_points,
				       GhostType ghost_type);

  /// pre compute all shapes on the element control points from natural coordinates
  template <ElementType type>
  void precomputeShapesDerivativesOnControlPoints(const Real * natural_coords,
						  const UInt nb_points,
						  const UInt dimension,
						  GhostType ghost_type);

  /// interpolate nodal values on the control points
  template <ElementType type>
  void interpolateOnControlPoints(const Vector<Real> &u,
				  Vector<Real> &uq,
				  UInt nb_degre_of_freedom,
				  GhostType ghost_type = _not_ghost,
				  const Vector<UInt> * filter_elements = NULL) const;


  /// multiply a field by shape functions
  template <ElementType type>
  void fieldTimesShapes(const Vector<Real> & field,
			Vector<Real> & fieal_times_shapes);

  // /// compute the gradient of u on the constrol points
  // template <ElementType type>
  // void gradientOnControlPoints(const Vector<Real> &u,
  // 			       Vector<Real> &nablauq,
  // 			       UInt nb_degre_of_freedom,
  // 			       GhostType ghost_type = _not_ghost,
  // 			       const Vector<UInt> * filter_elements = NULL) const;


  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {
    std::string space;
    for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
    stream << space << "Shapes [" << std::endl;
    control_points.printself(stream, indent + 1);
    stream << space << "]" << std::endl;
  };

  /// set the control points for a given element
  template <ElementType type>
  void setControlPointsByType(const Vector<Real> & control_points,
			      const GhostType & ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the size of the shapes returned by the element class
  static inline UInt getShapeSize(const ElementType & type);

  /// get the size of the shapes derivatives returned by the element class
  static inline UInt getShapeDerivativesSize(const ElementType & type);


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

protected:
  Mesh * mesh;

  ID id;

  /// shape functions for all elements
  ByElementTypeReal control_points;
};


#include "shape_functions_inline_impl.cc"


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const ShapeFunctions & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_SHAPE_FUNCTIONS_HH__ */
