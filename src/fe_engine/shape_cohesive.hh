/**
 * @file   shape_cohesive.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue Feb 15 2011
 * @date last modification: Fri May 14 2021
 *
 * @brief  shape functions for cohesive elements description
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
#include "shape_lagrange.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SHAPE_COHESIVE_HH_
#define AKANTU_SHAPE_COHESIVE_HH_

namespace akantu {

struct CohesiveReduceFunctionMean {
  inline Real operator()(Real u_plus, Real u_minus) {
    return .5 * (u_plus + u_minus);
  }
};

struct CohesiveReduceFunctionDifference {
  inline Real operator()(Real u_plus, Real u_minus) {
    return (u_plus - u_minus);
  }
};

struct ExtendingOperators {
  static Matrix<Real>
  getAveragingOperator(const ElementType & type,
                       const UInt & nb_degree_of_freedom = 1) {
    AKANTU_DEBUG_IN();
    auto kind = Mesh::getKind(type);
    AKANTU_DEBUG_ASSERT(kind == _ek_cohesive,
                        "Extending operators work only for cohesive elements");
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    // averaging operator
    Matrix<Real> A(nb_degree_of_freedom * nb_nodes_per_element / 2,
                   nb_degree_of_freedom * nb_nodes_per_element);

    for (UInt i = 0; i < nb_degree_of_freedom * nb_nodes_per_element / 2; ++i) {
      A(i, i) = 0.5;
      A(i, i + nb_degree_of_freedom * nb_nodes_per_element / 2) = 0.5;
    }
    AKANTU_DEBUG_OUT();
    return A;
  };

  static Matrix<Real>
  getDifferencingOperator(const ElementType & type,
                          const UInt & nb_degree_of_freedom = 1) {
    AKANTU_DEBUG_IN();
    auto kind = Mesh::getKind(type);
    AKANTU_DEBUG_ASSERT(kind == _ek_cohesive,
                        "Extending operators work only for cohesive elements");
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    // averaging operator
    Matrix<Real> A(nb_degree_of_freedom * nb_nodes_per_element / 2,
                   nb_degree_of_freedom * nb_nodes_per_element);

    for (UInt i = 0; i < nb_degree_of_freedom * nb_nodes_per_element / 2; ++i) {
      A(i, i) = 1;
      A(i, i + nb_degree_of_freedom * nb_nodes_per_element / 2) = -1;
    }
    AKANTU_DEBUG_OUT();
    return A;
  }
};

template <> class ShapeLagrange<_ek_cohesive> : public ShapeLagrangeBase {
  /* ------------------------------------------------------------------------
   */
  /* Constructors/Destructors */
  /* ------------------------------------------------------------------------
   */
public:
  ShapeLagrange(const Mesh & mesh, UInt spatial_dimension,
                const ID & id = "shape_cohesive");

  ~ShapeLagrange() override = default;

  /* ------------------------------------------------------------------------
   */
  /* Methods */
  /* ------------------------------------------------------------------------
   */
public:
  inline void initShapeFunctions(const Array<Real> & nodes,
                                 const Matrix<Real> & integration_points,
                                 ElementType type, GhostType ghost_type);

  /// extract the nodal values and store them per element
  template <ElementType type, class ReduceFunction>
  void extractNodalToElementField(
      const Array<Real> & nodal_f, Array<Real> & elemental_f,
      GhostType ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const;

  /// computes the shape functions derivatives for given interpolation points
  template <ElementType type>
  void computeShapeDerivativesOnIntegrationPoints(
      const Array<Real> & nodes, const Matrix<Real> & integration_points,
      Array<Real> & shape_derivatives, GhostType ghost_type,
      const Array<UInt> & filter_elements = empty_filter) const;

  void computeShapeDerivativesOnIntegrationPoints(
      const Array<Real> & nodes, const Matrix<Real> & integration_points,
      Array<Real> & shape_derivative, ElementType type, GhostType ghost_type,
      const Array<UInt> & filter_elements) const override;

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

  template <ElementType type>
  void computeShapeDerivativesOnIntegrationPointsLowerDimension(
      const Array<Real> & nodes, const Matrix<Real> & integration_points,
      Array<Real> & shape_derivatives, const GhostType & ghost_type,
      const Array<UInt> & filter_elements = empty_filter) const;

  /// interpolate nodal values on the integration points
  template <ElementType type, class ReduceFunction>
  void interpolateOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & uq, UInt nb_degree_of_freedom,
      GhostType ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const;

  template <ElementType type>
  void interpolateOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & uq, UInt nb_degree_of_freedom,
      GhostType ghost_type = _not_ghost,
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
    gradientOnIntegrationPoints<type, CohesiveReduceFunctionMean>(
        u, nablauq, nb_degree_of_freedom, ghost_type, filter_elements);
  }

  /// compute gradient on the facet element placed in the base of dim - 1
  template <ElementType type, class ReduceFunction>
  void gradientOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & nablauq, UInt nb_degree_of_freedom,
      GhostType ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const;

  /* ------------------------------------------------------------------------
   */
  template <ElementType type>
  void computeBtD(const Array<Real> & Ds, Array<Real> & BtDs,
                  GhostType ghost_type = _not_ghost,
                  const Array<UInt> & filter_elements = empty_filter) const;

  template <ElementType type>
  void
  computeExtendedBtD(const Array<Real> & Ds, Array<Real> & AtBtDs,
                     GhostType ghost_type = _not_ghost,
                     const Array<UInt> & filter_elements = empty_filter) const;

  template <ElementType type>
  void computeBtDB(const Array<Real> & Ds, Array<Real> & BtDBs, UInt order_d,
                   GhostType ghost_type,
                   const Array<UInt> & filter_elements) const {
    AKANTU_TO_IMPLEMENT();
  };

  /// multiply a field by shape functions
  template <ElementType type>
  void computeNtb(const Array<Real> & bs, Array<Real> & Ntbs,
                  GhostType ghost_type = _not_ghost,
                  const Array<UInt> & filter_elements = empty_filter) const;

  /// extend shape functions to the number of nodes in cohesive element
  template <ElementType type>
  void
  computeExtendedNtb(const Array<Real> & bs, Array<Real> & Ntbs,
                     GhostType ghost_type = _not_ghost,
                     const Array<UInt> & filter_elements = empty_filter) const;

  template <ElementType type>
  void computeNtbN(const Array<Real> & /*bs*/, Array<Real> & /*NtbNs*/,
                   GhostType /*ghost_type*/,
                   const Array<UInt> & /*filter_elements*/) const {
    AKANTU_TO_IMPLEMENT();
  }

  /* ------------------------------------------------------------------------
   */
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

  /// computes bases alligned with cohesive elements
  template <ElementType type, class ReduceFunction>
  void computeAllignedBasisOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & basis_u, GhostType ghost_type,
      const Array<UInt> & filter_elements = empty_filter) const;
};

/// standard output stream operator
template <class ShapeFunction>
inline std::ostream & operator<<(std::ostream & stream,
                                 const ShapeCohesive<ShapeFunction> & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "shape_cohesive_inline_impl.hh"

#endif /* AKANTU_SHAPE_COHESIVE_HH_ */
