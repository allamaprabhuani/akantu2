/**
 * @file   shape_cohesive_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Feb 03 2012
 * @date last modification: Tue Sep 29 2020
 *
 * @brief  ShapeCohesive inline implementation
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
#include "mesh_iterators.hh"
#include "shape_cohesive.hh"

/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SHAPE_COHESIVE_INLINE_IMPL_HH_
#define AKANTU_SHAPE_COHESIVE_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
inline ShapeLagrange<_ek_cohesive>::ShapeLagrange(const Mesh & mesh,
                                                  UInt spatial_dimension,
                                                  const ID & id)
    : ShapeLagrangeBase(mesh, spatial_dimension, _ek_cohesive, id) {}

#define INIT_SHAPE_FUNCTIONS(type)                                             \
  setIntegrationPointsByType<type>(integration_points, ghost_type);            \
  precomputeShapesOnIntegrationPoints<type>(nodes, ghost_type);                \
  precomputeShapeDerivativesOnIntegrationPoints<type>(nodes, ghost_type);

/* -------------------------------------------------------------------------- */
inline void ShapeLagrange<_ek_cohesive>::initShapeFunctions(
    const Array<Real> & nodes, const Matrix<Real> & integration_points,
    ElementType type, GhostType ghost_type) {
  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(INIT_SHAPE_FUNCTIONS);
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_cohesive>::computeShapeDerivativesOnIntegrationPoints(
    const Array<Real> & nodes, const Matrix<Real> & integration_points,
    Array<Real> & shape_derivatives, GhostType ghost_type,
    const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  this->computeShapeDerivativesOnIntegrationPointsLowerDimension<type>(
      nodes, integration_points, shape_derivatives, ghost_type,
      filter_elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void
ShapeLagrange<_ek_cohesive>::computeShapeDerivativesOnIntegrationPoints(
    const Array<Real> & nodes, const Matrix<Real> & integration_points,
    Array<Real> & shape_derivatives, ElementType type, GhostType ghost_type,
    const Array<UInt> & filter_elements) const {
#define AKANTU_COMPUTE_SHAPES(type)                                            \
  computeShapeDerivativesOnIntegrationPoints<type>(                            \
      nodes, integration_points, shape_derivatives, ghost_type,                \
      filter_elements);

  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(AKANTU_COMPUTE_SHAPES);

#undef AKANTU_COMPUTE_SHAPES
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_cohesive>::precomputeShapesOnIntegrationPoints(
    const Array<Real> & nodes, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  Matrix<Real> & natural_coords = integration_points(type, ghost_type);
  UInt size_of_shapes = ElementClass<type>::getShapeSize();

  Array<Real> & shapes_tmp =
      shapes.alloc(0, size_of_shapes, itp_type, ghost_type);

  this->computeShapesOnIntegrationPoints<type>(nodes, natural_coords,
                                               shapes_tmp, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_cohesive>::precomputeShapeDerivativesOnIntegrationPoints(
    const Array<Real> & nodes, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  Matrix<Real> & natural_coords = integration_points(type, ghost_type);
  UInt size_of_shapesd = ElementClass<type>::getShapeDerivativesSize();
  UInt spatial_dimension = mesh.getSpatialDimension();

  Array<Real> & shapes_derivatives_tmp =
      shapes_derivatives.alloc(0, size_of_shapesd, itp_type, ghost_type);

  this->computeShapeDerivativesOnIntegrationPointsLowerDimension<type>(
      nodes, natural_coords, shapes_derivatives_tmp, ghost_type);

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_cohesive>::
    computeShapeDerivativesOnIntegrationPointsLowerDimension(
        const Array<Real> & nodes, const Matrix<Real> & integration_points,
        Array<Real> & shape_derivatives, const GhostType & ghost_type,
        const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt natural_dimension = ElementClass<type>::getNaturalSpaceDimension();
  const auto facet_type = Mesh::getFacetType(type);

  UInt nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();

  UInt nb_points = integration_points.cols();
  UInt nb_element = mesh.getConnectivity(type, ghost_type).size();

  UInt size_of_shapesd = ElementClass<type>::getShapeDerivativesSize();
  AKANTU_DEBUG_ASSERT(shape_derivatives.getNbComponent() == size_of_shapesd,
                      "The shapes_derivatives array does not have the correct "
                          << "number of component");
  shape_derivatives.resize(nb_element * nb_points);

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  this->extractNodalToElementField<type, CohesiveReduceFunctionMean>(
      nodes, x_el, ghost_type, filter_elements);

  Real * shapesd_val = shape_derivatives.storage();
  Array<Real>::matrix_iterator x_it =
      x_el.begin(spatial_dimension, nb_nodes_per_element);

  Array<Real> x_el_equiv(nb_element, natural_dimension * nb_nodes_per_element,
                         0.);

  if (spatial_dimension == 3) {
    AKANTU_DEBUG_ASSERT(
        facet_type == _triangle_3,
        "The shape derivatives are calculated only for first order cohesives");
    for (auto && data :
         zip(make_view(x_el, spatial_dimension, nb_nodes_per_element),
             make_view(x_el_equiv, natural_dimension, nb_nodes_per_element))) {
      const Matrix<Real> & x = std::get<0>(data);
      auto & x_eq = std::get<1>(data);
      // compute triangle sides
      Vector<Real> AB = Vector<Real>(x(1)) - Vector<Real>(x(0));
      Vector<Real> AC = Vector<Real>(x(2)) - Vector<Real>(x(0));
      Vector<Real> BC = Vector<Real>(x(2)) - Vector<Real>(x(1));
      Real a = AB.norm();
      Real b = BC.norm();
      Real c = AC.norm();
      // assign projected coordinates to the nodes of equivalent element
      // A: [0,0]; B: [a, 0]; C: [xc, yc];
      // https://math.stackexchange.com/questions/50227/how-do-i-map-a-3d-triangle-into-2d
      x_eq(0, 1) = a;
      auto yc =
          sqrt((a + b - c) * (a - b + c) * (-a + b + c) * (a + b + c)) / 2 / a;
      auto xc = sqrt(c * c - yc * yc);
      x_eq(0, 2) = xc;
      x_eq(1, 2) = yc;
    }
  } else if (spatial_dimension == 2) {
    /// array of equivalent coordinates (works for 1D element in 2D only)
    AKANTU_DEBUG_ASSERT(
        facet_type == _segment_2,
        "The shape derivatives are calculated only for first order cohesives");
    for (auto && data :
         zip(make_view(x_el, spatial_dimension, nb_nodes_per_element),
             make_view(x_el_equiv, natural_dimension, nb_nodes_per_element))) {
      const Matrix<Real> & x = std::get<0>(data);
      auto & x_eq = std::get<1>(data);
      for (auto node_nb : arange(1, nb_nodes_per_element)) {
        Vector<Real> delta_x = Vector<Real>(x(node_nb)) - Vector<Real>(x(0));
        x_eq(natural_dimension - 1, node_nb) = delta_x.norm();
      }
    }
  }

  if (filter_elements != empty_filter)
    nb_element = filter_elements.size();

  for (auto && data : enumerate(
           make_view(x_el_equiv, natural_dimension, nb_nodes_per_element))) {
    auto & elem = std::get<0>(data);
    const auto & X = std::get<1>(data);
    if (filter_elements != empty_filter)
      shapesd_val = shape_derivatives.storage() +
                    filter_elements(elem) * size_of_shapesd * nb_points;

    Tensor3<Real> B(shapesd_val, natural_dimension, nb_nodes_per_element,
                    nb_points);
    computeShapeDerivativesOnCPointsByElement<facet_type>(X, integration_points,
                                                          B);

    if (filter_elements == empty_filter)
      shapesd_val += size_of_shapesd * nb_points;
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_cohesive>::computeBtD(
    const Array<Real> & Ds, Array<Real> & BtDs, GhostType ghost_type,
    const Array<UInt> & filter_elements) const {
  this->computeExtendedBtD<type>(Ds, BtDs, ghost_type, filter_elements);
}
/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_cohesive>::computeExtendedBtD(
    const Array<Real> & Ds, Array<Real> & AtBtDs, GhostType ghost_type,
    const Array<UInt> & filter_elements) const {
  auto itp_type = ElementClassProperty<type>::interpolation_type;
  const auto & shapes_derivatives =
      this->shapes_derivatives(itp_type, ghost_type);

  auto natural_dimension = ElementClass<type>::getNaturalSpaceDimension();
  auto nb_nodes_per_element = mesh.getNbNodesPerElement(type);

  Array<Real> shapes_derivatives_filtered(0,
                                          shapes_derivatives.getNbComponent());
  auto && view = make_view(shapes_derivatives, natural_dimension,
                           nb_nodes_per_element / 2);
  auto B_it = view.begin();
  auto B_end = view.end();

  if (filter_elements != empty_filter) {
    FEEngine::filterElementalData(this->mesh, shapes_derivatives,
                                  shapes_derivatives_filtered, type, ghost_type,
                                  filter_elements);
    auto && view = make_view(shapes_derivatives_filtered, natural_dimension,
                             nb_nodes_per_element / 2);
    B_it = view.begin();
    B_end = view.end();
  }

  auto A = ExtendingOperators::getAveragingOperator(type);
  Matrix<Real> B_A(natural_dimension, nb_nodes_per_element);

  for (auto && values :
       zip(range(B_it, B_end),
           make_view(Ds, Ds.getNbComponent() / natural_dimension,
                     natural_dimension),
           make_view(AtBtDs, AtBtDs.getNbComponent() / nb_nodes_per_element,
                     nb_nodes_per_element))) {
    const auto & B = std::get<0>(values);
    const auto & D = std::get<1>(values);
    B_A.mul<false, false>(B, A);
    auto & At_Bt_D = std::get<2>(values);
    // transposed due to the storage layout of B
    At_Bt_D.template mul<false, false>(D, B_A);
  }
}
/* -------------------------------------------------------------------------- */
template <ElementType type, class ReduceFunction>
void ShapeLagrange<_ek_cohesive>::extractNodalToElementField(
    const Array<Real> & nodal_f, Array<Real> & elemental_f,
    GhostType ghost_type, const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_itp_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();
  UInt nb_degree_of_freedom = nodal_f.getNbComponent();
  UInt nb_element = this->mesh.getNbElement(type, ghost_type);

  const auto & conn_array = this->mesh.getConnectivity(type, ghost_type);
  auto conn = conn_array.begin(conn_array.getNbComponent() / 2, 2);

  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
  }

  elemental_f.resize(nb_element);

  Array<Real>::matrix_iterator u_it =
      elemental_f.begin(nb_degree_of_freedom, nb_nodes_per_itp_element);

  ReduceFunction reduce_function;

  auto compute = [&](const auto & el) {
    Matrix<Real> & u = *u_it;
    Matrix<UInt> el_conn(conn[el]);

    // compute the average/difference of the nodal field loaded from cohesive
    // element
    for (UInt n = 0; n < el_conn.rows(); ++n) {
      UInt node_plus = el_conn(n, 0);
      UInt node_minus = el_conn(n, 1);
      for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
        Real u_plus = nodal_f(node_plus, d);
        Real u_minus = nodal_f(node_minus, d);
        u(d, n) = reduce_function(u_plus, u_minus);
      }
    }

    ++u_it;
  };

  for_each_element(nb_element, filter_elements, compute);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type, class ReduceFunction>
void ShapeLagrange<_ek_cohesive>::interpolateOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & out_uq, UInt nb_degree_of_freedom,
    GhostType ghost_type, const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;

  AKANTU_DEBUG_ASSERT(this->shapes.exists(itp_type, ghost_type),
                      "No shapes for the type "
                          << this->shapes.printType(itp_type, ghost_type));

  UInt nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();
  Array<Real> u_el(0, nb_degree_of_freedom * nb_nodes_per_element);
  this->extractNodalToElementField<type, ReduceFunction>(in_u, u_el, ghost_type,
                                                         filter_elements);

  this->template interpolateElementalFieldOnIntegrationPoints<type>(
      u_el, out_uq, ghost_type, shapes(itp_type, ghost_type), filter_elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type, class ReduceFunction>
void ShapeLagrange<_ek_cohesive>::variationOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & nablauq, UInt nb_degree_of_freedom,
    GhostType ghost_type, const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;

  AKANTU_DEBUG_ASSERT(
      this->shapes_derivatives.exists(itp_type, ghost_type),
      "No shapes for the type "
          << this->shapes_derivatives.printType(itp_type, ghost_type));

  UInt nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();
  Array<Real> u_el(0, nb_degree_of_freedom * nb_nodes_per_element);
  this->extractNodalToElementField<type, ReduceFunction>(in_u, u_el, ghost_type,
                                                         filter_elements);

  this->template gradientElementalFieldOnIntegrationPoints<type>(
      u_el, nablauq, ghost_type, shapes_derivatives(itp_type, ghost_type),
      filter_elements);

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <ElementType type, class ReduceFunction>
void ShapeLagrange<_ek_cohesive>::gradientOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & out_nablauq,
    UInt nb_degree_of_freedom, GhostType ghost_type,
    const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  AKANTU_DEBUG_ASSERT(
      shapes_derivatives.exists(itp_type, ghost_type),
      "No shapes derivatives for the type "
          << shapes_derivatives.printType(itp_type, ghost_type));

  UInt nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();

  Array<Real> u_el(0, nb_degree_of_freedom * nb_nodes_per_element);
  this->extractNodalToElementField<type, ReduceFunction>(in_u, u_el, ghost_type,
                                                         filter_elements);

  this->gradientElementalFieldOnIntegrationPoints<type>(
      u_el, out_nablauq, ghost_type, shapes_derivatives(itp_type, ghost_type),
      filter_elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type, class ReduceFunction>
void ShapeLagrange<_ek_cohesive>::computeNormalsOnIntegrationPoints(
    const Array<Real> & u, Array<Real> & normals_u, GhostType ghost_type,
    const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt nb_element = this->mesh.getNbElement(type, ghost_type);
  UInt nb_points = this->integration_points(type, ghost_type).cols();
  UInt spatial_dimension = this->mesh.getSpatialDimension();

  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
  }

  normals_u.resize(nb_points * nb_element);

  Array<Real> tangents_u(0, (spatial_dimension * (spatial_dimension - 1)));

  if (spatial_dimension > 1) {
    tangents_u.resize(nb_element * nb_points);
    this->template variationOnIntegrationPoints<type, ReduceFunction>(
        u, tangents_u, spatial_dimension, ghost_type, filter_elements);
  }

  Real * tangent = tangents_u.storage();

  if (spatial_dimension == 3) {
    for (auto & normal : make_view(normals_u, spatial_dimension)) {
      Math::vectorProduct3(tangent, tangent + spatial_dimension,
                           normal.storage());

      normal /= normal.norm();
      tangent += spatial_dimension * 2;
    }
  } else if (spatial_dimension == 2) {
    for (auto & normal : make_view(normals_u, spatial_dimension)) {
      Vector<Real> a1(tangent, spatial_dimension);

      normal(0) = -a1(1);
      normal(1) = a1(0);
      normal.normalize();

      tangent += spatial_dimension;
    }
  } else if (spatial_dimension == 1) {
    const auto facet_type = Mesh::getFacetType(type);
    const auto & mesh_facets = mesh.getMeshFacets();
    const auto & facets = mesh_facets.getSubelementToElement(type, ghost_type);
    const auto & segments =
        mesh_facets.getElementToSubelement(facet_type, ghost_type);

    Real values[2];

    for (auto el : arange(nb_element)) {
      if (filter_elements != empty_filter) {
        el = filter_elements(el);
      }

      for (UInt p = 0; p < 2; ++p) {
        Element facet = facets(el, p);
        Element segment = segments(facet.element)[0];
        Vector<Real> barycenter(values + p, 1);
        mesh.getBarycenter(segment, barycenter);
      }

      Real difference = values[0] - values[1];

      AKANTU_DEBUG_ASSERT(difference != 0.,
                          "Error in normal computation for cohesive elements");

      normals_u(el) = difference / std::abs(difference);
    }
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <ElementType type, class ReduceFunction>
void ShapeLagrange<_ek_cohesive>::computeAllignedBasisOnIntegrationPoints(
    const Array<Real> & u, Array<Real> & basis_u, GhostType ghost_type,
    const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt nb_element = this->mesh.getNbElement(type, ghost_type);
  UInt nb_points = this->integration_points(type, ghost_type).cols();
  UInt spatial_dimension = this->mesh.getSpatialDimension();

  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
  }

  basis_u.resize(nb_points * nb_element);

  Array<Real> tangents_u(0, (spatial_dimension * (spatial_dimension - 1)));

  if (spatial_dimension > 1) {
    tangents_u.resize(nb_element * nb_points);
    this->template variationOnIntegrationPoints<type, ReduceFunction>(
        u, tangents_u, spatial_dimension, ghost_type, filter_elements);
  }

  Real * tangent = tangents_u.storage();
  Real * basis_storage = basis_u.storage();
  for (auto && data :
       zip(make_view(basis_u, spatial_dimension, spatial_dimension),
           make_view(tangents_u, spatial_dimension, (spatial_dimension - 1)))) {
    auto && basis = std::get<0>(data);
    auto && tangents = std::get<1>(data);
    for (auto i : arange(spatial_dimension - 1)) {
      basis(i) = tangents(i);
      Vector<Real> basis_vec(basis(i), false);
      basis_vec /= basis_vec.norm();
    }

    if (spatial_dimension == 3) {
      Vector<Real> normal(spatial_dimension);
      Math::vectorProduct3(tangent, tangent + spatial_dimension,
                           basis_storage + 2 * spatial_dimension);
      Vector<Real> basis_vec(basis(2), false);
      basis_vec /= basis_vec.norm();
      tangent += spatial_dimension * 2;
      basis_storage += spatial_dimension * 2;
    } else if (spatial_dimension == 2) {
      Vector<Real> a1(tangent, spatial_dimension);
      Vector<Real> basis_vec(basis(1), false);

      basis_vec(0) = -a1(1);
      basis_vec(1) = a1(0);
      basis_vec /= basis_vec.norm();

      tangent += spatial_dimension;
    }

    else if (spatial_dimension == 1) {
      AKANTU_ERROR("1D mesh is currently not supported");
    }

    AKANTU_DEBUG_OUT();
  }
}
} // namespace akantu

#endif /* AKANTU_SHAPE_COHESIVE_INLINE_IMPL_HH_ */
