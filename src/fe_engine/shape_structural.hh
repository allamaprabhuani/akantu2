/**
 * Copyright (©) 2011-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
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
 */

/* -------------------------------------------------------------------------- */
#include "shape_functions.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SHAPE_STRUCTURAL_HH_
#define AKANTU_SHAPE_STRUCTURAL_HH_

namespace akantu {

template <ElementKind kind> class ShapeStructural : public ShapeFunctions {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  // Ctors/Dtors should be explicitely implemented for _ek_structural
public:
  ShapeStructural(Mesh & mesh, Int spatial_dimension,
                  const ID & id = "shape_structural");
  ~ShapeStructural() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override {
    std::string space(indent, AKANTU_INDENT);
    stream << space << "ShapesStructural [" << std::endl;
    rotation_matrices.printself(stream, indent + 1);
    ShapeFunctions::printself(stream, indent + 1);
    stream << space << "]" << std::endl;
  }

private:
  template <ElementType type>
  void computeShapesOnIntegrationPointsInternal(
      const Array<Real> & nodes, const Matrix<Real> & integration_points,
      Array<Real> & shapes, GhostType ghost_type,
      const Array<Idx> & filter_elements = empty_filter,
      bool mass = false) const;

public:
  /// compute shape functions on given integration points
  template <ElementType type>
  void computeShapesOnIntegrationPoints(
      const Array<Real> & nodes, const Matrix<Real> & integration_points,
      Array<Real> & shapes, GhostType ghost_type,
      const Array<Idx> & filter_elements = empty_filter) const {
    this->template computeShapesOnIntegrationPointsInternal<type>(
        nodes, integration_points, shapes, ghost_type, filter_elements, false);
  }

  template <ElementType type>
  void computeShapesMassOnIntegrationPoints(
      const Array<Real> & nodes, const Matrix<Real> & integration_points,
      Array<Real> & shapes, GhostType ghost_type,
      const Array<Idx> & filter_elements = empty_filter) const {
    this->template computeShapesOnIntegrationPointsInternal<type>(
        nodes, integration_points, shapes, ghost_type, filter_elements, true);
  }

  /// initialization function for structural elements
  inline void initShapeFunctions(const Array<Real> & nodes,
                                 const Matrix<Real> & integration_points,
                                 ElementType type, GhostType ghost_type);

  /// precompute the rotation matrices for the elements dofs
  template <ElementType type>
  void precomputeRotationMatrices(const Array<Real> & nodes,
                                  GhostType ghost_type);

  /// pre compute all shapes on the element integration points from natural
  /// coordinates
  template <ElementType type>
  void precomputeShapesOnIntegrationPoints(const Array<Real> & nodes,
                                           GhostType ghost_type);

  /// pre compute all shapes on the element integration points from natural
  /// coordinates
  template <ElementType type>
  void precomputeShapeDerivativesOnIntegrationPoints(const Array<Real> & nodes,
                                                     GhostType ghost_type);

  /// interpolate nodal values on the integration points
  template <ElementType type>
  void interpolateOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & uq, Int nb_dof,
      GhostType ghost_type = _not_ghost,
      const Array<Idx> & filter_elements = empty_filter) const;

  /// compute the gradient of u on the integration points
  template <ElementType type>
  void gradientOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & nablauq, Int nb_dof,
      GhostType ghost_type = _not_ghost,
      const Array<Idx> & filter_elements = empty_filter) const;

  /// interpolate on physical point
  template <ElementType type, typename D1, typename D2, typename D3,
            aka::enable_if_t<aka::are_vectors<D1, D3>::value> * = nullptr>
  void interpolate(const Eigen::MatrixBase<D1> & /*real_coords*/, Idx /*elem*/,
                   const Eigen::MatrixBase<D2> & /*nodal_values*/,
                   Eigen::MatrixBase<D3> & /*interpolated*/,
                   GhostType /*ghost_type*/) const {
    AKANTU_TO_IMPLEMENT();
  }

  /// compute the shapes on a provided point
  template <ElementType type, typename D1, typename D2,
            aka::enable_if_t<aka::are_vectors<D1, D2>::value> * = nullptr>
  void computeShapes(const Eigen::MatrixBase<D1> & /*real_coords*/,
                     Idx /*elem*/, Eigen::MatrixBase<D2> & /*shapes*/,
                     GhostType /*ghost_type*/) const {
    AKANTU_TO_IMPLEMENT();
  }

  /// compute the shape derivatives on a provided point
  template <ElementType type, typename D1>
  void computeShapeDerivatives(const Eigen::MatrixBase<D1> & /*real_coords*/,
                               Idx /*elem*/, Tensor3Base<Real> & /*shapes*/,
                               GhostType /*ghost_type*/) const {
    AKANTU_TO_IMPLEMENT();
  }

  /// get the rotations vector
  inline const Array<Real> &
  getRotations(ElementType el_type,
               GhostType /*ghost_type*/ = _not_ghost) const {
    return rotation_matrices(el_type);
  }

  /* ------------------------------------------------------------------------ */
  template <ElementType type>
  void computeBtD(const Array<Real> & /*Ds*/, Array<Real> & /*BtDs*/,
                  GhostType /*ghost_type*/,
                  const Array<Idx> & /*filter_elements*/) const {
    AKANTU_TO_IMPLEMENT();
  }

  template <ElementType type>
  void computeBtDB(const Array<Real> & /*Ds*/, Array<Real> & /*BtDBs*/,
                   Int /*order_d*/, GhostType /*ghost_type*/,
                   const Array<Idx> & /*filter_elements*/) const {
    AKANTU_TO_IMPLEMENT();
  }

  template <ElementType type>
  void computeNtbN(const Array<Real> & /*bs*/, Array<Real> & /*NtbNs*/,
                   GhostType /*ghost_type*/,
                   const Array<Idx> & /*filter_elements*/) const {
    AKANTU_TO_IMPLEMENT();
  }

  /// multiply a field by shape functions
  template <ElementType type>
  void computeNtb(const Array<Real> & /*bs*/, Array<Real> & /*Ntbs*/,
                  GhostType /*ghost_type*/,
                  const Array<Idx> & /*filter_elements*/ = empty_filter) const {
    AKANTU_TO_IMPLEMENT();
  }

protected:
  ElementTypeMapArray<Real> rotation_matrices;
};

} // namespace akantu

#include "shape_structural_inline_impl.hh"

#endif /* AKANTU_SHAPE_STRUCTURAL_HH_ */
