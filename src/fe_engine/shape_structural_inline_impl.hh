/**
 * Copyright (©) 2017-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "mesh_iterators.hh"
#include "shape_structural.hh"
/* -------------------------------------------------------------------------- */

//#ifndef AKANTU_SHAPE_STRUCTURAL_INLINE_IMPL_HH_
//#define AKANTU_SHAPE_STRUCTURAL_INLINE_IMPL_HH_

namespace akantu {

namespace {
  /// Extract nodal coordinates per elements
  template <ElementType type>
  std::unique_ptr<Array<Real>> getNodesPerElement(const Mesh & mesh,
                                                  const Array<Real> & nodes,
                                                  GhostType ghost_type) {
    const auto dim = ElementClass<type>::getSpatialDimension();
    const auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    auto nodes_per_element =
        std::make_unique<Array<Real>>(0, dim * nb_nodes_per_element);
    FEEngine::extractNodalToElementField(mesh, nodes, *nodes_per_element, type,
                                         ghost_type);
    return nodes_per_element;
  }
} // namespace

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeStructural<kind>::computeShapesOnIntegrationPointsInternal(
    const Array<Real> & nodes, const Matrix<Real> & integration_points,
    Array<Real> & shapes, GhostType ghost_type,
    const Array<Idx> & filter_elements, bool mass) const {

  auto nb_points = integration_points.cols();
  auto nb_element = mesh.getConnectivity(type, ghost_type).size();
  const auto nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement();

  shapes.resize(nb_element * nb_points);

  const auto nb_dofs = ElementClass<type>::getNbDegreeOfFreedom();
  auto nb_rows = nb_dofs;

  if (mass) {
    nb_rows = ElementClass<type>::getNbStressComponents();
  }

#if !defined(AKANTU_NDEBUG)
  Int size_of_shapes = nb_rows * nb_dofs * nb_nodes_per_element;
  AKANTU_DEBUG_ASSERT(shapes.getNbComponent() == size_of_shapes,
                      "The shapes array does not have the correct "
                          << "number of component");
#endif

  const auto nb_cols_shapes =
      ElementClass<type>::getNbNodesPerInterpolationElement() * nb_dofs;
  auto nodes_per_element = getNodesPerElement<type>(mesh, nodes, ghost_type);

  auto shapes_view = make_view(shapes, nb_rows, nb_cols_shapes, nb_points);
  auto nodes_view = make_view<Eigen::Dynamic, nb_nodes_per_element>(
      *nodes_per_element, mesh.getSpatialDimension(), nb_nodes_per_element);
  auto R_view =
      make_view<nb_dofs, nb_dofs>(rotation_matrices(type, ghost_type));

  auto loop_core = [&](auto && data) {
    auto & N = std::get<0>(data);
    auto & X = std::get<1>(data);
    auto & R = std::get<2>(data);

    Matrix<Real, nb_cols_shapes, nb_cols_shapes> T;
    T.zero();

    for (Int i = 0; i < nb_nodes_per_element; ++i) {
      T.template block<nb_dofs, nb_dofs>(i * nb_dofs, i * nb_dofs) = R;
    }

    if (not mass) {
      ElementClass<type>::computeShapes(integration_points, X, T, N);
    } else {
      ElementClass<type>::computeShapesMass(integration_points, X, T, N);
    }
  };

  if (filter_elements == empty_filter) {
    for (auto && data : zip(shapes_view, nodes_view, R_view)) {
      loop_core(data);
    }

  } else {
    for (auto && data :
         filter(filter_elements, zip(shapes_view, nodes_view, R_view))) {
      loop_core(data);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeStructural<kind>::precomputeRotationMatrices(
    const Array<Real> & nodes, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const auto spatial_dimension = mesh.getSpatialDimension();
  const auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  const auto nb_element = mesh.getNbElement(type, ghost_type);
  const auto nb_dof = ElementClass<type>::getNbDegreeOfFreedom();

  if (not this->rotation_matrices.exists(type, ghost_type)) {
    this->rotation_matrices.alloc(0, nb_dof * nb_dof, type, ghost_type);
  }

  auto & rot_matrices = this->rotation_matrices(type, ghost_type);
  rot_matrices.resize(nb_element);

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type);

  bool has_extra_normal = mesh.hasData<Real>("extra_normal", type, ghost_type);
  Array<Real>::const_vector_iterator extra_normal;
  if (has_extra_normal) {
    extra_normal = mesh.getData<Real>("extra_normal", type, ghost_type)
                       .begin(spatial_dimension);
  }

  for (auto && tuple :
       zip(make_view(x_el, spatial_dimension, nb_nodes_per_element),
           make_view(rot_matrices, nb_dof, nb_dof))) {
    // compute shape derivatives
    auto & X = std::get<0>(tuple);
    auto & R = std::get<1>(tuple);

    if (has_extra_normal) {
      ElementClass<type>::computeRotationMatrix(R, X, *extra_normal);
      ++extra_normal;
    } else {
      ElementClass<type>::computeRotationMatrix(
          R, X, Vector<Real>(spatial_dimension));
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <ElementKind kind>
template <ElementType type>
void ShapeStructural<kind>::precomputeShapesOnIntegrationPoints(
    const Array<Real> & nodes, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const auto & natural_coords = integration_points(type, ghost_type);
  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_points = integration_points(type, ghost_type).cols();
  auto nb_element = mesh.getNbElement(type, ghost_type);
  auto nb_dof = ElementClass<type>::getNbDegreeOfFreedom();
  const auto dim = ElementClass<type>::getSpatialDimension();
  const auto spatial_dimension = mesh.getSpatialDimension();
  const auto natural_spatial_dimension =
      ElementClass<type>::getNaturalSpaceDimension();

  auto itp_type = FEEngine::getInterpolationType(type);
  if (not shapes.exists(itp_type, ghost_type)) {
    auto size_of_shapes = this->getShapeSize(type);
    this->shapes.alloc(0, size_of_shapes, itp_type, ghost_type);
  }

  auto & rot_matrices = this->rotation_matrices(type, ghost_type);
  auto & shapes_ = this->shapes(itp_type, ghost_type);
  shapes_.resize(nb_element * nb_points);

  auto nodes_per_element = getNodesPerElement<type>(mesh, nodes, ghost_type);

  for (auto && tuple :
       zip(make_view(shapes_, nb_dof, nb_dof * nb_nodes_per_element, nb_points),
           make_view(*nodes_per_element, dim, nb_nodes_per_element),
           make_view(rot_matrices, nb_dof, nb_dof))) {
    auto && N = std::get<0>(tuple);
    auto && X = std::get<1>(tuple);
    auto && RDOFs = std::get<2>(tuple);

    Matrix<Real> T(N.size(1), N.size(1));
    T.zero();

    for (Idx i = 0; i < nb_nodes_per_element; ++i) {
      T.block(i * RDOFs.rows(), i * RDOFs.cols(), RDOFs.rows(), RDOFs.cols()) =
          RDOFs;
    }

    auto R = RDOFs.block(0, 0, spatial_dimension, spatial_dimension);
    // Rotate to local basis
    auto x =
        (R * X).block(0, 0, natural_spatial_dimension, nb_nodes_per_element);

    ElementClass<type>::computeShapes(natural_coords, x, T, N);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeStructural<kind>::precomputeShapeDerivativesOnIntegrationPoints(
    const Array<Real> & nodes, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const auto & natural_coords = integration_points(type, ghost_type);
  const auto spatial_dimension = mesh.getSpatialDimension();
  const auto natural_spatial_dimension =
      ElementClass<type>::getNaturalSpaceDimension();
  const auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  const auto nb_points = natural_coords.cols();
  const auto nb_dof = ElementClass<type>::getNbDegreeOfFreedom();
  const auto nb_element = mesh.getNbElement(type, ghost_type);
  const auto nb_stress_components = ElementClass<type>::getNbStressComponents();
  const auto nb_cols_shaped = nb_dof * nb_nodes_per_element;

  auto itp_type = FEEngine::getInterpolationType(type);
  if (not this->shapes_derivatives.exists(itp_type, ghost_type)) {
    auto size_of_shapesd = this->getShapeDerivativesSize(type);
    this->shapes_derivatives.alloc(0, size_of_shapesd, itp_type, ghost_type);
  }

  auto & rot_matrices = this->rotation_matrices(type, ghost_type);

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type);

  auto & shapesd = this->shapes_derivatives(itp_type, ghost_type);
  shapesd.resize(nb_element * nb_points);

  for (auto && tuple :
       zip(make_view(x_el, spatial_dimension, nb_nodes_per_element),
           make_view(shapesd, nb_stress_components, nb_cols_shaped, nb_points),
           make_view<nb_dof, nb_dof>(rot_matrices))) {
    // compute shape derivatives
    auto & X = std::get<0>(tuple);
    auto & B = std::get<1>(tuple);
    auto & RDOFs = std::get<2>(tuple);

    Tensor3<Real> dnds(natural_spatial_dimension,
                       ElementClass<type>::interpolation_property::dnds_columns,
                       B.size(2));
    ElementClass<type>::computeDNDS(natural_coords, X, dnds);

    Tensor3<Real> J(natural_spatial_dimension, natural_spatial_dimension,
                    natural_coords.cols());

    // Computing the coordinates of the element in the natural space
    Matrix<Real, nb_cols_shaped, nb_cols_shaped> T;
    T.zero();

    for (Int i = 0; i < nb_nodes_per_element; ++i) {
      T.template block<nb_dof, nb_dof>(i * nb_dof, i * nb_dof) = RDOFs;
    }

    auto && R = RDOFs.block(0, 0, spatial_dimension, spatial_dimension);
    // Rotate to local basis
    auto x =
        (R * X).template block<natural_spatial_dimension, nb_nodes_per_element>(
            0, 0);

    ElementClass<type>::computeJMat(natural_coords, x, J);
    ElementClass<type>::computeShapeDerivatives(J, dnds, T, B);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeStructural<kind>::interpolateOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & out_uq, Int nb_dof,
    GhostType ghost_type, const Array<Idx> & filter_elements) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(out_uq.getNbComponent() == nb_dof,
                      "The output array shape is not correct");

  auto itp_type = FEEngine::getInterpolationType(type);
  const auto & shapes_ = shapes(itp_type, ghost_type);

  auto nb_element = mesh.getNbElement(type, ghost_type);
  auto nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement();
  auto nb_quad_points_per_element = integration_points(type, ghost_type).cols();

  Array<Real> u_el(0, nb_nodes_per_element * nb_dof);
  FEEngine::extractNodalToElementField(mesh, in_u, u_el, type, ghost_type,
                                       filter_elements);

  auto nb_quad_points = nb_quad_points_per_element * u_el.size();
  out_uq.resize(nb_quad_points);

  auto out_it =
      make_view(out_uq, nb_dof, 1, nb_quad_points_per_element).begin();
  auto shapes_it = make_view(shapes_, nb_dof, nb_dof * nb_nodes_per_element,
                             nb_quad_points_per_element)
                       .begin();
  auto u_it = make_view(u_el, nb_dof * nb_nodes_per_element, 1,
                        nb_quad_points_per_element)
                  .begin();

  for_each_element(nb_element, filter_elements, [&](auto && el) {
    auto & uq = *out_it;
    const auto & u = *u_it;
    auto && N = shapes_it[el];

    for (auto && q : arange(uq.size(2))) {
      uq(q) = N(q) * u(q);
    }

    ++out_it;
    ++u_it;
  });
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeStructural<kind>::gradientOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & out_nablauq, Int nb_dof,
    GhostType ghost_type, const Array<Idx> & filter_elements) const {
  AKANTU_DEBUG_IN();

  auto itp_type = FEEngine::getInterpolationType(type);
  const auto & shapesd = shapes_derivatives(itp_type, ghost_type);

  auto nb_element = mesh.getNbElement(type, ghost_type);
  auto element_dimension = ElementClass<type>::getSpatialDimension();
  auto nb_quad_points_per_element = integration_points(type, ghost_type).cols();
  auto nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement();

  Array<Real> u_el(0, nb_nodes_per_element * nb_dof);
  FEEngine::extractNodalToElementField(mesh, in_u, u_el, type, ghost_type,
                                       filter_elements);

  auto nb_quad_points = nb_quad_points_per_element * u_el.size();
  out_nablauq.resize(nb_quad_points);

  auto out_it =
      make_view(out_nablauq, element_dimension, 1, nb_quad_points_per_element)
          .begin();
  auto shapesd_it =
      make_view(shapesd, element_dimension, nb_dof * nb_nodes_per_element,
                nb_quad_points_per_element)
          .begin();
  auto u_it = make_view(u_el, nb_dof * nb_nodes_per_element, 1,
                        nb_quad_points_per_element)
                  .begin();

  for_each_element(nb_element, filter_elements, [&](auto && el) {
    auto & nablau = *out_it;
    const auto & u = *u_it;
    auto B = shapesd_it[el];

    for (auto && q : arange(nablau.size(2))) {
      nablau(q) = B(q) * u(q);
    }

    ++out_it;
    ++u_it;
  });

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <>
template <ElementType type>
void ShapeStructural<_ek_structural>::computeBtD(
    const Array<Real> & Ds, Array<Real> & BtDs, GhostType ghost_type,
    const Array<Idx> & filter_elements) const {
  auto itp_type = ElementClassProperty<type>::interpolation_type;

  auto nb_stress = ElementClass<type>::getNbStressComponents();
  auto nb_dof_per_element = ElementClass<type>::getNbDegreeOfFreedom() *
                            mesh.getNbNodesPerElement(type);

  const auto & shapes_derivatives =
      this->shapes_derivatives(itp_type, ghost_type);

  Array<Real> shapes_derivatives_filtered(0,
                                          shapes_derivatives.getNbComponent());
  auto && view = make_view(shapes_derivatives, nb_stress, nb_dof_per_element);
  auto B_it = view.begin();
  auto B_end = view.end();

  if (filter_elements != empty_filter) {
    FEEngine::filterElementalData(this->mesh, shapes_derivatives,
                                  shapes_derivatives_filtered, type, ghost_type,
                                  filter_elements);
    auto && view =
        make_view(shapes_derivatives_filtered, nb_stress, nb_dof_per_element);
    B_it = view.begin();
    B_end = view.end();
  }

  for (auto && values : zip(range(B_it, B_end), make_view(Ds, nb_stress),
                            make_view(BtDs, BtDs.getNbComponent()))) {
    const auto & B = std::get<0>(values);
    const auto & D = std::get<1>(values);
    auto & Bt_D = std::get<2>(values);
    Bt_D = B.transpose() * D;
  }
}

/* -------------------------------------------------------------------------- */
template <>
template <ElementType type>
void ShapeStructural<_ek_structural>::computeNtb(
    const Array<Real> & bs, Array<Real> & Ntbs, GhostType ghost_type,
    const Array<Idx> & filter_elements) const {
  constexpr auto itp_type = ElementClassProperty<type>::interpolation_type;
  constexpr auto nb_dof = ElementClass<type>::getNbDegreeOfFreedom();
  constexpr auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

  const auto & shapes = this->shapes(itp_type, ghost_type);

  Array<Real> shapes_filtered(0, shapes.getNbComponent());
  auto && view = make_view<nb_dof, nb_dof * nb_nodes_per_element>(shapes);
  auto N_it = view.begin();
  auto N_end = view.end();

  if (filter_elements != empty_filter) {
    FEEngine::filterElementalData(this->mesh, shapes, shapes_filtered, type,
                                  ghost_type, filter_elements);
    auto && view =
        make_const_view<nb_dof, nb_dof * nb_nodes_per_element>(shapes_filtered);
    N_it = view.begin();
    N_end = view.end();
  }

  for (auto && [N, b, Nt_b] :
       zip(range(N_it, N_end), make_view<nb_dof>(bs),
           make_view<nb_dof * nb_nodes_per_element>(Ntbs))) {
    Nt_b = N.transpose() * b;
  }
}
/* -------------------------------------------------------------------------- */
template <ElementKind kind>
inline void ShapeStructural<kind>::initShapeFunctions(
    const Array<Real> & /* unused */, const Matrix<Real> & /* unused */,
    ElementType /* unused */, GhostType /* unused */) {
  AKANTU_TO_IMPLEMENT();
}

template <>
inline void ShapeStructural<_ek_structural>::initShapeFunctions(
    const Array<Real> & nodes, const Matrix<Real> & integration_points,
    ElementType type, GhostType ghost_type) {
  tuple_dispatch<ElementTypes_t<_ek_structural>>(
      [&](auto && enum_type) {
        constexpr ElementType type = aka::decay_v<decltype(enum_type)>;
        this->setIntegrationPointsByType<type>(integration_points, ghost_type);
        this->precomputeRotationMatrices<type>(nodes, ghost_type);
        this->precomputeShapesOnIntegrationPoints<type>(nodes, ghost_type);
        this->precomputeShapeDerivativesOnIntegrationPoints<type>(nodes,
                                                                  ghost_type);
      },
      type);
}
/* -------------------------------------------------------------------------- */

} // namespace akantu

//#endif /* AKANTU_SHAPE_STRUCTURAL_INLINE_IMPL_HH_ */
