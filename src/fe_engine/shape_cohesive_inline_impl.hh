/**
 * Copyright (©) 2012-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
//#include "shape_cohesive.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SHAPE_COHESIVE_INLINE_IMPL_HH_
#define AKANTU_SHAPE_COHESIVE_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
inline ShapeLagrange<_ek_cohesive>::ShapeLagrange(const Mesh & mesh,
                                                  Int spatial_dimension,
                                                  const ID & id)
    : ShapeLagrangeBase(mesh, spatial_dimension, _ek_cohesive, id) {}

/* -------------------------------------------------------------------------- */
inline void ShapeLagrange<_ek_cohesive>::initShapeFunctions(
    const Array<Real> & nodes, const Matrix<Real> & integration_points,
    ElementType type, GhostType ghost_type) {
  tuple_dispatch<ElementTypes_t<_ek_cohesive>>(
      [&](auto && enum_type) {
        constexpr ElementType type = aka::decay_v<decltype(enum_type)>;
        this->setIntegrationPointsByType<type>(integration_points, ghost_type);
        this->precomputeShapesOnIntegrationPoints<type>(nodes, ghost_type);
        this->precomputeShapeDerivativesOnIntegrationPoints<type>(nodes,
                                                                  ghost_type);
      },
      type);
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_cohesive>::computeShapeDerivativesOnIntegrationPoints(
    const Array<Real> &, const Matrix<Real> & integration_points,
    Array<Real> & shape_derivatives, GhostType ghost_type,
    const Array<Idx> & filter_elements) const {
  AKANTU_DEBUG_IN();

  auto size_of_shapesd = ElementClass<type>::getShapeDerivativesSize();
  auto spatial_dimension = ElementClass<type>::getNaturalSpaceDimension();
  auto nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();

  auto nb_points = integration_points.cols();
  auto nb_element = mesh.getConnectivity(type, ghost_type).size();

  AKANTU_DEBUG_ASSERT(shape_derivatives.getNbComponent() == size_of_shapesd,
                      "The shapes_derivatives array does not have the correct "
                          << "number of component");

  shape_derivatives.resize(nb_element * nb_points);
  auto * shapesd_val = shape_derivatives.data();

  auto compute = [&](const auto & el) {
    auto ptr = shapesd_val + el * nb_points * size_of_shapesd;
    Tensor3Proxy<Real> B(ptr, spatial_dimension, nb_nodes_per_element,
                         nb_points);
    ElementClass<type>::computeDNDS(integration_points, B);
  };

  for_each_element(nb_element, filter_elements, compute);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void
ShapeLagrange<_ek_cohesive>::computeShapeDerivativesOnIntegrationPoints(
    const Array<Real> & nodes, const Ref<const MatrixXr> integration_points,
    Array<Real> & shape_derivatives, ElementType type, GhostType ghost_type,
    const Array<Idx> & filter_elements) const {
  tuple_dispatch<ElementTypes_t<_ek_cohesive>>(
      [&](auto && enum_type) {
        constexpr ElementType type = aka::decay_v<decltype(enum_type)>;
        this->computeShapeDerivativesOnIntegrationPoints<type>(
            nodes, integration_points, shape_derivatives, ghost_type,
            filter_elements);
      },
      type);
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_cohesive>::precomputeShapesOnIntegrationPoints(
    const Array<Real> & nodes, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto itp_type = ElementClassProperty<type>::interpolation_type;
  auto & natural_coords = integration_points(type, ghost_type);
  auto size_of_shapes = ElementClass<type>::getShapeSize();

  auto & shapes_tmp = shapes.alloc(0, size_of_shapes, itp_type, ghost_type);

  this->computeShapesOnIntegrationPoints<type>(nodes, natural_coords,
                                               shapes_tmp, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_cohesive>::precomputeShapeDerivativesOnIntegrationPoints(
    const Array<Real> & nodes, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto itp_type = ElementClassProperty<type>::interpolation_type;
  auto & natural_coords = integration_points(type, ghost_type);
  auto size_of_shapesd = ElementClass<type>::getShapeDerivativesSize();

  auto & shapes_derivatives_tmp =
      shapes_derivatives.alloc(0, size_of_shapesd, itp_type, ghost_type);

  this->computeShapeDerivativesOnIntegrationPoints<type>(
      nodes, natural_coords, shapes_derivatives_tmp, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type, class ReduceFunction>
void ShapeLagrange<_ek_cohesive>::extractNodalToElementField(
    const Array<Real> & nodal_f, Array<Real> & elemental_f,
    GhostType ghost_type, const Array<Idx> & filter_elements) const {
  AKANTU_DEBUG_IN();

  auto nb_nodes_per_itp_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();
  auto nb_degree_of_freedom = nodal_f.getNbComponent();
  auto nb_element = this->mesh.getNbElement(type, ghost_type);

  const auto & conn_array = this->mesh.getConnectivity(type, ghost_type);
  auto conn = make_view(conn_array, conn_array.getNbComponent() / 2, 2).begin();

  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
  }

  elemental_f.resize(nb_element);

  auto u_it =
      make_view(elemental_f, nb_degree_of_freedom, nb_nodes_per_itp_element)
          .begin();

  auto nodal_f_it = make_view(nodal_f, nb_degree_of_freedom).begin();

  ReduceFunction reduce_function;

  auto compute = [&](const auto & el) {
    auto && u = *u_it;
    auto && el_conn = conn[el];

    // compute the average/difference of the nodal field loaded from cohesive
    // element
    for (Int n = 0; n < el_conn.rows(); ++n) {
      auto node_plus = el_conn(n, 0);
      auto node_minus = el_conn(n, 1);

      u(n) = reduce_function(nodal_f_it[node_plus], nodal_f_it[node_minus]);
    }

    ++u_it;
  };

  for_each_element(nb_element, filter_elements, compute);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type, class ReduceFunction>
void ShapeLagrange<_ek_cohesive>::interpolateOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & out_uq, Int nb_degree_of_freedom,
    GhostType ghost_type, const Array<Idx> & filter_elements) const {
  AKANTU_DEBUG_IN();

  auto itp_type = ElementClassProperty<type>::interpolation_type;

  AKANTU_DEBUG_ASSERT(this->shapes.exists(itp_type, ghost_type),
                      "No shapes for the type "
                          << this->shapes.printType(itp_type, ghost_type));

  auto nb_nodes_per_element =
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
    const Array<Real> & in_u, Array<Real> & nablauq, Int nb_degree_of_freedom,
    GhostType ghost_type, const Array<Idx> & filter_elements) const {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;

  AKANTU_DEBUG_ASSERT(
      this->shapes_derivatives.exists(itp_type, ghost_type),
      "No shapes for the type "
          << this->shapes_derivatives.printType(itp_type, ghost_type));

  auto nb_nodes_per_element =
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
void ShapeLagrange<_ek_cohesive>::computeNormalsOnIntegrationPoints(
    const Array<Real> & u, Array<Real> & normals_u, GhostType ghost_type,
    const Array<Idx> & filter_elements) const {
  AKANTU_DEBUG_IN();

  auto nb_element = this->mesh.getNbElement(type, ghost_type);
  auto nb_points = this->integration_points(type, ghost_type).cols();
  auto spatial_dimension = this->mesh.getSpatialDimension();

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

  if (spatial_dimension == 3) {
    for (auto && data :
         zip(make_view<3>(normals_u), make_view<3, 2>(tangents_u))) {
      auto && n = std::get<0>(data);
      auto && ts = std::get<1>(data);

      n = (ts(0).cross(ts(1))).normalized();
    }
  } else if (spatial_dimension == 2) {
    for (auto && data :
         zip(make_view<2>(normals_u), make_view<2>(tangents_u))) {
      auto && n = std::get<0>(data);
      auto && t = std::get<1>(data);

      n(0) = -t(1);
      n(1) = t(0);
      n.normalize();
    }
  } else if (spatial_dimension == 1) {
    const auto facet_type = Mesh::getFacetType(type);
    const auto & mesh_facets = mesh.getMeshFacets();
    const auto & facets = mesh_facets.getSubelementToElement(type, ghost_type);
    const auto & segments =
        mesh_facets.getElementToSubelement(facet_type, ghost_type);

    Matrix<Real> barycenter(1, 2);

    for (auto el : arange(nb_element)) {
      if (filter_elements != empty_filter) {
        el = filter_elements(el);
      }

      for (Int p = 0; p < 2; ++p) {
        const Element & facet = facets(el, p);
        const Element & segment = segments(facet.element)[0];
        mesh.getBarycenter(segment, barycenter(p));
      }

      Real difference = barycenter(0, 0) - barycenter(0, 1);

      AKANTU_DEBUG_ASSERT(difference != 0.,
                          "Error in normal computation for cohesive elements");

      normals_u(el) = difference / std::abs(difference);
    }
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu

#endif /* AKANTU_SHAPE_COHESIVE_INLINE_IMPL_HH_ */
