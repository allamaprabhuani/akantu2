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

template<>
class ShapeLagrange<_ek_cohesive> : public ShapeFunctions {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ShapeLagrange(const Mesh & mesh,
		const ID & id = "shape_cohesive",
		const MemoryID & memory_id = 0);

  virtual ~ShapeLagrange() { }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  inline void initShapeFunctions(const Vector<Real> & nodes,
				 const types::Matrix<Real> & control_points,
				 const ElementType & type,
				 const GhostType & ghost_type);

  /// extract the nodal values and store them per element
  template <ElementType type, class ReduceFunction>
  void extractNodalToElementField(const Vector<Real> & nodal_f,
				  Vector<Real> & elemental_f,
				  const GhostType & ghost_type = _not_ghost,
				  const Vector<UInt> * filter_elements = NULL) const;

  /// pre compute all shapes on the element control points from natural coordinates
  template <ElementType type>
  void precomputeShapesOnControlPoints(const Vector<Real> & nodes,
				       GhostType ghost_type);

  /// pre compute all shapes on the element control points from natural coordinates
  template <ElementType type>
  void precomputeShapeDerivativesOnControlPoints(const Vector<Real> & nodes,
						 GhostType ghost_type);

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

  /// compute the gradient of u on the control points in the natural coordinates
  template <ElementType type>
  void gradientOnControlPoints(const Vector<Real> &u,
			       Vector<Real> &nablauq,
			       UInt nb_degree_of_freedom,
			       GhostType ghost_type = _not_ghost,
			       const Vector<UInt> * filter_elements = NULL) const {
    variationOnControlPoints<type, CohesiveReduceFunctionMean>(u, nablauq, nb_degree_of_freedom, ghost_type, filter_elements);
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
  void fieldTimesShapes(const Vector<Real> & field,
			Vector<Real> & fiedl_times_shapes,
			GhostType ghost_type) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get a the shapes vector
  inline const Vector<Real> & getShapes(const ElementType & el_type,
					const GhostType & ghost_type = _not_ghost) const;

  /// get a the shapes derivatives vector
  inline const Vector<Real> & getShapesDerivatives(const ElementType & el_type,
						   const GhostType & ghost_type = _not_ghost) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// shape functions for all elements
  ByElementTypeVector<Real, InterpolationType> shapes;

  /// shape functions derivatives for all elements
  ByElementTypeVector<Real, InterpolationType> shapes_derivatives;

};

// __END_AKANTU__
// #include "shape_lagrange.hh"
// __BEGIN_AKANTU__

// template<>
// class ShapeLagrange<_ek_cohesive> : public ShapeCohesive< ShapeLagrange<_ek_regular> > {
// public:
//   ShapeLagrange(const Mesh & mesh,
// 		const ID & id = "shape_cohesive",
// 		const MemoryID & memory_id = 0) :
//     ShapeCohesive< ShapeLagrange<_ek_regular> >(mesh, id, memory_id) { }

//   virtual ~ShapeLagrange() { };
// };


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
