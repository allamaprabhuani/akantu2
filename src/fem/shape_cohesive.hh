/**
 * @file   shape_cohesive.hh
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Feb 03 16:28:03 2012
 *
 * @brief  shape functions for cohesive elements description
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

#include "aka_vector.hh"
#include "shape_functions.hh"
#include "shape_lagrange.hh"

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SHAPE_COHESIVE_HH__
#define __AKANTU_SHAPE_COHESIVE_HH__

__BEGIN_AKANTU__

struct CohesiveReduceFunctionMean {
  inline Real operator()(Real u_plus, Real u_minus) {
    return .5*(u_plus + u_minus);
  }
};

struct CohesiveReduceFunctionOpening {
  inline Real operator()(Real u_plus, Real u_minus) {
    return (u_plus - u_minus);
  }
};



template <class ShapeFunction>
class ShapeCohesive : public ShapeFunctions {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  ShapeCohesive(const Mesh & mesh,
		const ID & id = "shape_cohesive",
		const MemoryID & memory_id = 0);

  virtual ~ShapeCohesive() { if(sub_type_shape_function) delete sub_type_shape_function; };

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
  template <ElementType type, class ReduceFunction>
  void interpolateOnControlPoints(const Vector<Real> &u,
				  Vector<Real> &uq,
				  UInt nb_degree_of_freedom,
				  const GhostType ghost_type = _not_ghost,
				  const Vector<UInt> * filter_elements = NULL) const;

  template <ElementType type>
  void interpolateOnControlPoints(const Vector<Real> &u,
				  Vector<Real> &uq,
				  UInt nb_degree_of_freedom,
				  const GhostType ghost_type = _not_ghost,
				  const Vector<UInt> * filter_elements = NULL) const {
    interpolateOnControlPoints<type, CohesiveReduceFunctionMean>(u, uq, nb_degree_of_freedom, ghost_type, filter_elements);
  }


  /// compute the gradient of u on the control points
  template <ElementType type, class ReduceFunction>
  void variationOnControlPoints(const Vector<Real> &u,
				Vector<Real> &nablauq,
				UInt nb_degree_of_freedom,
				GhostType ghost_type = _not_ghost,
				const Vector<UInt> * filter_elements = NULL) const;

  template <ElementType type, class ReduceFunction>
  void computeNormalsOnControlPoints(const Vector<Real> &u,
				     Vector<Real> &normals_u,
				     GhostType ghost_type = _not_ghost,
				     const Vector<UInt> * filter_elements = NULL) const;

  /// multiply a field by shape functions
  template <ElementType type>
  void fieldTimesShapes(__attribute__((unused)) const Vector<Real> & field,
			__attribute__((unused)) Vector<Real> & fiedl_times_shapes,
			__attribute__((unused)) GhostType ghost_type) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// function to print the contain of the class
  //  virtual void printself(std::ostream & stream, int indent = 0) const {};

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  const Vector<Real> & getShapes(const ElementType & type,
				 GhostType ghost_type = _not_ghost) const;

  const Vector<Real> & getShapesDerivatives(const ElementType & type,
					    GhostType ghost_type = _not_ghost) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// real shape function implementation
  //\todo QUESTION!!!!!!! why having a member of the class you derive ? what do you gain ? containing the storage of the shapes 
  ShapeFunction * sub_type_shape_function;
};


typedef ShapeCohesive<ShapeLagrange> ShapeCohesiveLagrange;

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "shape_cohesive_inline_impl.cc"

/// standard output stream operator
template <class ShapeFunction>
inline std::ostream & operator <<(std::ostream & stream, const ShapeCohesive<ShapeFunction> & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_SHAPE_COHESIVE_HH__ */
