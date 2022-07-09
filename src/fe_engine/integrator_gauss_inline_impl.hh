/**
 * @file   integrator_gauss_inline_impl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Feb 15 2011
 * @date last modification: Tue Oct 27 2020
 *
 * @brief  inline function of gauss integrator
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
#include "fe_engine.hh"
#include "mesh_iterators.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
namespace debug {
  struct IntegratorGaussException : public Exception {};
} // namespace debug
/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
inline void IntegratorGauss<kind, IntegrationOrderFunctor>::integrateOnElement(
    const Array<Real> & f, Real * intf, Int nb_degree_of_freedom,
    const Idx elem, GhostType ghost_type) const {
  auto & jac_loc = jacobians(type, ghost_type);

  auto nb_quadrature_points = ElementClass<type>::getNbQuadraturePoints();
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == nb_degree_of_freedom,
                      "The vector f do not have the good number of component.");

  auto * f_val = f.data() + elem * f.getNbComponent();
  auto * jac_val = jac_loc.data() + elem * nb_quadrature_points;

  integrate(f_val, jac_val, intf, nb_degree_of_freedom, nb_quadrature_points);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
inline Real IntegratorGauss<kind, IntegrationOrderFunctor>::integrate(
    const Vector<Real> & in_f, Idx index, GhostType ghost_type) const {
  const Array<Real> & jac_loc = jacobians(type, ghost_type);

  auto nb_quadrature_points =
      GaussIntegrationElement<type>::getNbQuadraturePoints();
  AKANTU_DEBUG_ASSERT(in_f.size() == nb_quadrature_points,
                      "The vector f do not have nb_quadrature_points entries.");

  auto * jac_val = jac_loc.data() + index * nb_quadrature_points;
  Real intf;

  integrate(in_f.data(), jac_val, &intf, 1, nb_quadrature_points);

  return intf;
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
inline void IntegratorGauss<kind, IntegrationOrderFunctor>::integrate(
    const Real * f, const Real * jac, Real * inte, Int nb_degree_of_freedom,
    Int nb_quadrature_points) const {
  Eigen::Map<VectorXr> inte_v(inte, nb_quadrature_points);
  Eigen::Map<const VectorXr> cjac(jac, nb_quadrature_points);

  Eigen::Map<const MatrixXr> fq(f, nb_degree_of_freedom, nb_quadrature_points);

  inte_v.zero();
  inte_v = fq * cjac;
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
inline const Matrix<Real> &
IntegratorGauss<kind, IntegrationOrderFunctor>::getIntegrationPoints(
    GhostType ghost_type) const {
  AKANTU_DEBUG_ASSERT(
      quadrature_points.exists(type, ghost_type),
      "Quadrature points for type "
          << quadrature_points.printType(type, ghost_type)
          << " have not been initialized."
          << " Did you use 'computeQuadraturePoints' function ?");
  return (quadrature_points(type, ghost_type));
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
inline Int
IntegratorGauss<kind, IntegrationOrderFunctor>::getNbIntegrationPoints(
    GhostType ghost_type) const {
  AKANTU_DEBUG_ASSERT(
      quadrature_points.exists(type, ghost_type),
      "Quadrature points for type "
          << quadrature_points.printType(type, ghost_type)
          << " have not been initialized."
          << " Did you use 'computeQuadraturePoints' function ?");
  return quadrature_points(type, ghost_type).cols();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type, Int polynomial_degree>
inline Matrix<Real>
IntegratorGauss<kind, IntegrationOrderFunctor>::getIntegrationPoints() const {
  return GaussIntegrationElement<type,
                                 polynomial_degree>::getQuadraturePoints();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type, Int polynomial_degree>
inline Vector<Real>
IntegratorGauss<kind, IntegrationOrderFunctor>::getIntegrationWeights() const {
  return GaussIntegrationElement<type, polynomial_degree>::getWeights();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
inline void
IntegratorGauss<kind, IntegrationOrderFunctor>::computeQuadraturePoints(
    GhostType ghost_type) {
  auto & quads = quadrature_points(type, ghost_type);
  constexpr auto polynomial_degree =
      IntegrationOrderFunctor::template getOrder<type>();
  quads =
      GaussIntegrationElement<type, polynomial_degree>::getQuadraturePoints();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type, class D1, class D2, class D3>
inline void IntegratorGauss<kind, IntegrationOrderFunctor>::
    computeJacobianOnQuadPointsByElement(
        const Eigen::MatrixBase<D1> & node_coords,
        const Eigen::MatrixBase<D2> & quad,
        Eigen::MatrixBase<D3> & jacobians) const {
  // jacobian
  ElementClass<type>::computeJacobian(quad, node_coords, jacobians);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
IntegratorGauss<kind, IntegrationOrderFunctor>::IntegratorGauss(
    const Mesh & mesh, Int spatial_dimension, const ID & id)
    : Integrator(mesh, spatial_dimension, id) {}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void IntegratorGauss<kind, IntegrationOrderFunctor>::checkJacobians(
    GhostType ghost_type) const {
  auto nb_quadrature_points = this->quadrature_points(type, ghost_type).cols();
  auto nb_element = mesh.getConnectivity(type, ghost_type).size();
  auto * jacobians_val = jacobians(type, ghost_type).data();

  for (Idx i = 0; i < nb_element * nb_quadrature_points; ++i, ++jacobians_val) {
    if (*jacobians_val < 0) {
      AKANTU_CUSTOM_EXCEPTION_INFO(debug::IntegratorGaussException{},
                                   "Negative jacobian computed,"
                                       << " possible problem in the element "
                                          "node ordering (Quadrature Point "
                                       << i % nb_quadrature_points << ":"
                                       << i / nb_quadrature_points << ":"
                                       << type << ":" << ghost_type << ")");
    }
  }
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void IntegratorGauss<kind, IntegrationOrderFunctor>::
    computeJacobiansOnIntegrationPoints(
        const Array<Real> & nodes, const Matrix<Real> & quad_points,
        Array<Real> & jacobians, GhostType ghost_type,
        const Array<Idx> & filter_elements) const {
  auto spatial_dimension = mesh.getSpatialDimension();
  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_quadrature_points = quad_points.cols();
  auto nb_element = mesh.getNbElement(type, ghost_type);

  jacobians.resize(nb_element * nb_quadrature_points);

  auto jacobians_it = make_view(jacobians, nb_quadrature_points).begin();
  auto jacobians_begin = jacobians_it;

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type,
                                       filter_elements);

  auto x_it = make_view(x_el, spatial_dimension, nb_nodes_per_element).begin();

  nb_element = x_el.size();

  //  Matrix<Real> local_coord(spatial_dimension, nb_nodes_per_element);
  for (Idx elem = 0; elem < nb_element; ++elem, ++x_it) {
    const auto & x = *x_it;
    if (filter_elements != empty_filter) {
      jacobians_it = jacobians_begin + filter_elements(elem);
    }

    computeJacobianOnQuadPointsByElement<type>(x, quad_points, *jacobians_it);

    if (filter_elements == empty_filter) {
      ++jacobians_it;
    }
  }
}

/* -------------------------------------------------------------------------- */
#if defined(AKANTU_STRUCTURAL_MECHANICS)
template <>
template <ElementType type>
void IntegratorGauss<_ek_structural, DefaultIntegrationOrderFunctor>::
    computeJacobiansOnIntegrationPoints(
        const Array<Real> & nodes, const Matrix<Real> & quad_points,
        Array<Real> & jacobians, GhostType ghost_type,
        const Array<Int> & filter_elements) const {
  const auto spatial_dimension = mesh.getSpatialDimension();
  const auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  const auto nb_quadrature_points = quad_points.cols();
  const auto nb_dofs = ElementClass<type>::getNbDegreeOfFreedom();

  auto nb_element = mesh.getNbElement(type, ghost_type);
  jacobians.resize(nb_element * nb_quadrature_points);

  auto jacobians_it = make_view(jacobians, nb_quadrature_points).begin();
  auto jacobians_begin = jacobians_it;

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type,
                                       filter_elements);

  auto x_it = make_view(x_el, spatial_dimension, nb_nodes_per_element).begin();

  nb_element = x_el.size();

  const auto has_extra_normal =
      mesh.hasData<Real>("extra_normal", type, ghost_type);
  Array<Real>::const_vector_iterator extra_normal [[gnu::unused]];
  Array<Real>::const_vector_iterator extra_normal_begin [[gnu::unused]];
  if (has_extra_normal) {
    extra_normal = mesh.getData<Real>("extra_normal", type, ghost_type)
                       .begin(spatial_dimension);
    extra_normal_begin = extra_normal;
  }

  //  Matrix<Real> local_coord(spatial_dimension, nb_nodes_per_element);
  for (Idx elem = 0; elem < nb_element; ++elem, ++x_it) {
    if (filter_elements != empty_filter) {
      jacobians_it = jacobians_begin + filter_elements(elem);
      if (has_extra_normal) {
        extra_normal = extra_normal_begin + filter_elements(elem);
      }
    }

    const auto & X = *x_it;
    auto & J = *jacobians_it;
    Matrix<Real> R(nb_dofs, nb_dofs);

    if (has_extra_normal) {
      ElementClass<type>::computeRotationMatrix(R, X, *extra_normal);
    } else {
      ElementClass<type>::computeRotationMatrix(R, X, Vector<Real>(X.rows()));
    }

    const Int natural_space = ElementClass<type>::getNaturalSpaceDimension();
    const Int nb_nodes = ElementClass<type>::getNbNodesPerElement();
    // Extracting relevant lines
    auto x = (R.block(0, 0, spatial_dimension, spatial_dimension) * X)
                 .template block<natural_space, nb_nodes>(0, 0);

    computeJacobianOnQuadPointsByElement<type>(x, quad_points, J);

    if (filter_elements == empty_filter) {
      ++jacobians_it;
      if (has_extra_normal) {
        ++extra_normal;
      }
    }
  }
}
#endif

/* -------------------------------------------------------------------------- */
#if defined(AKANTU_COHESIVE_ELEMENT)
template <>
template <ElementType type>
void IntegratorGauss<_ek_cohesive, DefaultIntegrationOrderFunctor>::
    computeJacobiansOnIntegrationPoints(
        const Array<Real> & nodes, const Matrix<Real> & quad_points,
        Array<Real> & jacobians, GhostType ghost_type,
        const Array<Int> & filter_elements) const {
  auto spatial_dimension = mesh.getSpatialDimension();
  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_quadrature_points = quad_points.cols();
  auto nb_element = mesh.getNbElement(type, ghost_type);

  jacobians.resize(nb_element * nb_quadrature_points);

  auto jacobians_begin = make_view(jacobians, nb_quadrature_points).begin();

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type,
                                       filter_elements);

  auto x_it = make_view(x_el, spatial_dimension, nb_nodes_per_element).begin();

  auto nb_nodes_per_subelement = nb_nodes_per_element / 2;
  Matrix<Real> x(spatial_dimension, nb_nodes_per_subelement);

  nb_element = x_el.size();
  Idx l_el = 0;

  auto compute = [&](const auto & el) {
    auto && J = jacobians_begin[el];
    auto && X = x_it[l_el];
    ++l_el;

    for (Int n = 0; n < nb_nodes_per_subelement; ++n)
      x(n) = (X(n) + X(n + nb_nodes_per_subelement)) / 2.;

    if (type == _cohesive_1d_2) {
      J(0) = 1;
    } else {
      this->computeJacobianOnQuadPointsByElement<type>(x, quad_points, J);
    }
  };

  for_each_element(nb_element, filter_elements, compute);
}
#endif
/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void IntegratorGauss<kind, IntegrationOrderFunctor>::
    precomputeJacobiansOnQuadraturePoints(const Array<Real> & nodes,
                                          GhostType ghost_type) {
  auto & jacobians_tmp = jacobians.alloc(0, 1, type, ghost_type);

  this->computeJacobiansOnIntegrationPoints<type>(
      nodes, quadrature_points(type, ghost_type), jacobians_tmp, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type, Int polynomial_degree>
void IntegratorGauss<kind, IntegrationOrderFunctor>::multiplyJacobiansByWeights(
    Array<Real> & jacobians, const Array<Int> & filter_elements) const {
  constexpr auto nb_quadrature_points =
      GaussIntegrationElement<type, polynomial_degree>::getNbQuadraturePoints();
  auto && weights =
      GaussIntegrationElement<type, polynomial_degree>::getWeights();

  auto && view = make_view<nb_quadrature_points>(jacobians);

  if (filter_elements != empty_filter) {
    for (auto && J : filter(filter_elements, view))
      J.array() *= weights.array();
  } else {
    for (auto && J : view)
      J.array() *= weights.array();
  }
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
void IntegratorGauss<kind, IntegrationOrderFunctor>::integrate(
    const Array<Real> & in_f, Array<Real> & intf, Int nb_degree_of_freedom,
    const Array<Real> & jacobians, Int nb_element) const {
  intf.resize(nb_element);
  if (nb_element == 0) {
    return;
  }

  auto nb_points = jacobians.size() / nb_element;

  for (auto && data : zip(make_view(in_f, nb_degree_of_freedom, nb_points),
                          make_view(intf, nb_degree_of_freedom),
                          make_view(jacobians, nb_points))) {
    auto && f = std::get<0>(data);
    auto && int_f = std::get<1>(data);
    auto && J = std::get<2>(data);

    int_f = f * J;
  }
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void IntegratorGauss<kind, IntegrationOrderFunctor>::integrate(
    const Array<Real> & in_f, Array<Real> & intf, Int nb_degree_of_freedom,
    GhostType ghost_type, const Array<Int> & filter_elements) const {
  AKANTU_DEBUG_ASSERT(jacobians.exists(type, ghost_type),
                      "No jacobians for the type "
                          << jacobians.printType(type, ghost_type));

  const auto & jac_loc = jacobians(type, ghost_type);
  if (filter_elements != empty_filter) {
    auto nb_element = filter_elements.size();
    Array<Real> filtered_J(0, jac_loc.getNbComponent());
    FEEngine::filterElementalData(mesh, jac_loc, filtered_J, type, ghost_type,
                                  filter_elements);
    this->integrate(in_f, intf, nb_degree_of_freedom, filtered_J, nb_element);
  } else {
    auto nb_element = mesh.getNbElement(type, ghost_type);
    this->integrate(in_f, intf, nb_degree_of_freedom, jac_loc, nb_element);
  }
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type, Int polynomial_degree>
void IntegratorGauss<kind, IntegrationOrderFunctor>::integrate(
    const Array<Real> & in_f, Array<Real> & intf, Int nb_degree_of_freedom,
    GhostType ghost_type) const {
  auto quads = this->getIntegrationPoints<type, polynomial_degree>();

  Array<Real> jacobians;
  this->computeJacobiansOnIntegrationPoints<type>(mesh.getNodes(), quads,
                                                  jacobians, ghost_type);
  this->multiplyJacobiansByWeights<type, polynomial_degree>(jacobians);

  this->integrate(in_f, intf, nb_degree_of_freedom, jacobians,
                  mesh.getNbElement(type, ghost_type));
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type, Int polynomial_degree>
Real IntegratorGauss<kind, IntegrationOrderFunctor>::integrate(
    const Array<Real> & in_f, GhostType ghost_type) const {
  Array<Real> intfv(0, 1);
  integrate<type, polynomial_degree>(in_f, intfv, 1, ghost_type);

  auto res = Math::reduce(intfv);
  return res;
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
Real IntegratorGauss<kind, IntegrationOrderFunctor>::integrate(
    const Array<Real> & in_f, GhostType ghost_type,
    const Array<Int> & filter_elements) const {
  AKANTU_DEBUG_ASSERT(jacobians.exists(type, ghost_type),
                      "No jacobians for the type "
                          << jacobians.printType(type, ghost_type));

  Array<Real> intfv(0, 1);
  integrate<type>(in_f, intfv, 1, ghost_type, filter_elements);

  auto res = Math::reduce(intfv);
  return res;
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
void IntegratorGauss<kind, IntegrationOrderFunctor>::
    integrateOnIntegrationPoints(const Array<Real> & in_f, Array<Real> & intf,
                                 Int nb_degree_of_freedom,
                                 const Array<Real> & jacobians,
                                 Int nb_element) const {
  auto nb_points = jacobians.size() / nb_element;

  intf.resize(nb_element * nb_points);

  auto J_it = jacobians.begin();
  auto f_it = in_f.begin(nb_degree_of_freedom);
  auto inte_it = intf.begin(nb_degree_of_freedom);

  for (Idx el = 0; el < nb_element; ++el, ++J_it, ++f_it, ++inte_it) {
    const auto & J = *J_it;
    const auto & f = *f_it;
    auto & inte_f = *inte_it;

    inte_f = f;
    inte_f *= J;
  }
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void IntegratorGauss<kind, IntegrationOrderFunctor>::
    integrateOnIntegrationPoints(const Array<Real> & in_f, Array<Real> & intf,
                                 Int nb_degree_of_freedom, GhostType ghost_type,
                                 const Array<Int> & filter_elements) const {
  AKANTU_DEBUG_ASSERT(jacobians.exists(type, ghost_type),
                      "No jacobians for the type "
                          << jacobians.printType(type, ghost_type));

  const auto & jac_loc = this->jacobians(type, ghost_type);

  if (filter_elements != empty_filter) {

    auto nb_element = filter_elements.size();
    auto filtered_J =
        std::make_shared<Array<Real>>(0, jac_loc.getNbComponent());
    FEEngine::filterElementalData(mesh, jac_loc, *filtered_J, type, ghost_type,
                                  filter_elements);

    this->integrateOnIntegrationPoints(in_f, intf, nb_degree_of_freedom,
                                       *filtered_J, nb_element);
  } else {
    auto nb_element = mesh.getNbElement(type, ghost_type);
    this->integrateOnIntegrationPoints(in_f, intf, nb_degree_of_freedom,
                                       jac_loc, nb_element);
  }
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
inline void
IntegratorGauss<kind, IntegrationOrderFunctor>::onElementsAddedByType(
    const Array<Idx> & elements, GhostType ghost_type) {
  const auto & nodes = mesh.getNodes();

  computeQuadraturePoints<type>(ghost_type);

  if (not jacobians.exists(type, ghost_type)) {
    jacobians.alloc(0, 1, type, ghost_type);
  }

  this->computeJacobiansOnIntegrationPoints(
      nodes, quadrature_points(type, ghost_type), jacobians(type, ghost_type),
      type, ghost_type, elements);

  constexpr auto polynomial_degree =
      IntegrationOrderFunctor::template getOrder<type>();

  this->multiplyJacobiansByWeights<type, polynomial_degree>(
      this->jacobians(type, ghost_type), elements);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
void IntegratorGauss<kind, IntegrationOrderFunctor>::onElementsAdded(
    const Array<Element> & new_elements) {
  for (auto elements_range : MeshElementsByTypes(new_elements)) {
    auto type = elements_range.getType();
    auto ghost_type = elements_range.getGhostType();

    if (mesh.getSpatialDimension(type) != _spatial_dimension) {
      continue;
    }

    if (mesh.getKind(type) != kind) {
      continue;
    }

    tuple_dispatch<ElementTypes_t<kind>>(
        [&](auto && enum_type) {
          constexpr auto type = std::decay_t<decltype(enum_type)>::value;

          this->template onElementsAddedByType<type>(
              elements_range.getElements(), ghost_type);
        },
        type);
  }
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
inline void IntegratorGauss<kind, IntegrationOrderFunctor>::initIntegrator(
    const Array<Real> & nodes, GhostType ghost_type) {
  computeQuadraturePoints<type>(ghost_type);
  precomputeJacobiansOnQuadraturePoints<type>(nodes, ghost_type);
  checkJacobians<type>(ghost_type);
  constexpr auto polynomial_degree =
      IntegrationOrderFunctor::template getOrder<type>();
  multiplyJacobiansByWeights<type, polynomial_degree>(
      this->jacobians(type, ghost_type));
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
inline void IntegratorGauss<kind, IntegrationOrderFunctor>::initIntegrator(
    const Array<Real> & nodes, ElementType type, GhostType ghost_type) {
  tuple_dispatch<ElementTypes_t<kind>>(
      [&](auto && enum_type) {
        constexpr auto type = std::decay_t<decltype(enum_type)>::value;

        this->template initIntegrator<type>(nodes, ghost_type);
      },
      type);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
void IntegratorGauss<kind, IntegrationOrderFunctor>::
    computeJacobiansOnIntegrationPoints(
        const Array<Real> & nodes, const Matrix<Real> & quad_points,
        Array<Real> & jacobians, ElementType type, GhostType ghost_type,
        const Array<Idx> & filter_elements) const {
  tuple_dispatch<ElementTypes_t<kind>>(
      [&](auto && enum_type) {
        constexpr auto type = std::decay_t<decltype(enum_type)>::value;

        this->template computeJacobiansOnIntegrationPoints<type>(
            nodes, quad_points, jacobians, ghost_type, filter_elements);
      },
      type);
}

} // namespace akantu
