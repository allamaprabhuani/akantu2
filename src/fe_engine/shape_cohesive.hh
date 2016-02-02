/**
 * @file   shape_cohesive.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue Feb 15 2011
 * @date last modification: Thu Oct 22 2015
 *
 * @brief  shape functions for cohesive elements description
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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

#include "aka_array.hh"
#include "shape_functions.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SHAPE_COHESIVE_HH__
#define __AKANTU_SHAPE_COHESIVE_HH__

__BEGIN_AKANTU__

struct CohesiveReduceFunctionMean {
  inline Real operator()(Real u_plus, Real u_minus) {
    return .5 * (u_plus + u_minus);
  }
};

struct CohesiveReduceFunctionOpening {
  inline Real operator()(Real u_plus, Real u_minus) {
    return (u_plus - u_minus);
  }
};

template <> class ShapeLagrange<_ek_cohesive> : public ShapeFunctions {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ShapeLagrange(const Mesh & mesh, const ID & id = "shape_cohesive",
                const MemoryID & memory_id = 0);

  virtual ~ShapeLagrange() {}

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  inline void initShapeFunctions(const Array<Real> & nodes,
                                 const Matrix<Real> & integration_points,
                                 const ElementType & type,
                                 const GhostType & ghost_type);

  /// extract the nodal values and store them per element
  template <ElementType type, class ReduceFunction>
  void extractNodalToElementField(
      const Array<Real> & nodal_f, Array<Real> & elemental_f,
      const GhostType & ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const;

  /// pre compute all shapes on the element integration points from natural
  /// coordinates
  template <ElementType type>
  void precomputeShapesOnIntegrationPoints(const Array<Real> & nodes,
                                           GhostType ghost_type);

  /// pre compute all shape derivatives on the element integration points from
  /// natural coordinates
  template <ElementType type>
  void precomputeShapeDerivativesOnIntegrationPoints(const Array<Real> & nodes,
                                                     GhostType ghost_type);

  /// interpolate nodal values on the integration points
  template <ElementType type, class ReduceFunction>
  void interpolateOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & uq, UInt nb_degree_of_freedom,
      const GhostType ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const;

  template <ElementType type>
  void interpolateOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & uq, UInt nb_degree_of_freedom,
      const GhostType ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const {
    interpolateOnIntegrationPoints<type, CohesiveReduceFunctionMean>(
        u, uq, nb_degree_of_freedom, ghost_type, filter_elements);
  }

  /// compute the gradient of u on the integration points in the natural
  /// coordinates
  template <ElementType type>
  void gradientOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & nablauq, UInt nb_degree_of_freedom,
      GhostType ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const {
    variationOnIntegrationPoints<type, CohesiveReduceFunctionMean>(
        u, nablauq, nb_degree_of_freedom, ghost_type, filter_elements);
  }

  /// compute the gradient of u on the integration points
  template <ElementType type, class ReduceFunction>
  void variationOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & nablauq, UInt nb_degree_of_freedom,
      GhostType ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const;

  /// compute the normals to the field u on integration points
  template <ElementType type, class ReduceFunction>
  void computeNormalsOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & normals_u,
      GhostType ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const;

  /// multiply a field by shape functions
  template <ElementType type>
  void fieldTimesShapes(__attribute__((unused)) const Array<Real> & field,
                        __attribute__((unused))
                        Array<Real> & fiedl_times_shapes,
                        __attribute__((unused)) GhostType ghost_type) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get a the shapes vector
  inline const Array<Real> &
  getShapes(const ElementType & el_type,
            const GhostType & ghost_type = _not_ghost) const;

  /// get a the shapes derivatives vector
  inline const Array<Real> &
  getShapesDerivatives(const ElementType & el_type,
                       const GhostType & ghost_type = _not_ghost) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// shape functions for all elements
  ElementTypeMapArray<Real, InterpolationType> shapes;

  /// shape functions derivatives for all elements
  ElementTypeMapArray<Real, InterpolationType> shapes_derivatives;
};

// __END_AKANTU__
// #include "shape_lagrange.hh"
// __BEGIN_AKANTU__

// template<>
// class ShapeLagrange<_ek_cohesive> : public ShapeCohesive<
// ShapeLagrange<_ek_regular> > {
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
inline std::ostream & operator<<(std::ostream & stream,
                                 const ShapeCohesive<ShapeFunction> & _this) {
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AKANTU_SHAPE_COHESIVE_HH__ */
