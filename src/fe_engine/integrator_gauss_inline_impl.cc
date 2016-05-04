/**
 * @file   integrator_gauss_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Feb 15 2011
 * @date last modification: Thu Nov 19 2015
 *
 * @brief  inline function of gauss integrator
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

__END_AKANTU__

#include "fe_engine.hh"
#if defined(AKANTU_DEBUG_TOOLS)
#include "aka_debug_tools.hh"
#endif

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
inline void IntegratorGauss<kind>::integrateOnElement(
    const Array<Real> & f, Real * intf, UInt nb_degree_of_freedom,
    const UInt elem, const GhostType & ghost_type) const {
  Array<Real> & jac_loc = jacobians(type, ghost_type);

  UInt nb_quadrature_points = ElementClass<type>::getNbQuadraturePoints();
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == nb_degree_of_freedom,
                      "The vector f do not have the good number of component.");

  Real * f_val = f.storage() + elem * f.getNbComponent();
  Real * jac_val = jac_loc.storage() + elem * nb_quadrature_points;

  integrate(f_val, jac_val, intf, nb_degree_of_freedom, nb_quadrature_points);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
inline Real
IntegratorGauss<kind>::integrate(const Vector<Real> & in_f, UInt index,
                                 const GhostType & ghost_type) const {
  const Array<Real> & jac_loc = jacobians(type, ghost_type);

  UInt nb_quadrature_points =
      GaussIntegrationElement<type>::getNbQuadraturePoints();
  AKANTU_DEBUG_ASSERT(in_f.size() == nb_quadrature_points,
                      "The vector f do not have nb_quadrature_points entries.");

  Real * jac_val = jac_loc.storage() + index * nb_quadrature_points;
  Real intf;

  integrate(in_f.storage(), jac_val, &intf, 1, nb_quadrature_points);

  return intf;
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
inline void IntegratorGauss<kind>::integrate(Real * f, Real * jac, Real * inte,
                                             UInt nb_degree_of_freedom,
                                             UInt nb_quadrature_points) const {
  memset(inte, 0, nb_degree_of_freedom * sizeof(Real));

  Real * cjac = jac;
  for (UInt q = 0; q < nb_quadrature_points; ++q) {
    for (UInt dof = 0; dof < nb_degree_of_freedom; ++dof) {
      inte[dof] += *f * *cjac;
      ++f;
    }
    ++cjac;
  }
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
inline const Matrix<Real> & IntegratorGauss<kind>::getIntegrationPoints(
    const GhostType & ghost_type) const {
  AKANTU_DEBUG_ASSERT(
      quadrature_points.exists(type, ghost_type),
      "Quadrature points for type "
          << quadrature_points.printType(type, ghost_type)
          << " have not been initialized."
          << " Did you use 'computeQuadraturePoints' function ?");
  return quadrature_points(type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
inline UInt IntegratorGauss<kind>::getNbIntegrationPoints(
    const GhostType & ghost_type) const {
  AKANTU_DEBUG_ASSERT(
      quadrature_points.exists(type, ghost_type),
      "Quadrature points for type "
          << quadrature_points.printType(type, ghost_type)
          << " have not been initialized."
          << " Did you use 'computeQuadraturePoints' function ?");
  return quadrature_points(type, ghost_type).cols();
}

template <ElementKind kind>
template <ElementType type, UInt polynomial_degree>
inline Matrix<Real> IntegratorGauss<kind>::getIntegrationPoints() const {
  return GaussIntegrationElement<type,
                                 polynomial_degree>::getQuadraturePoints();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
inline void
IntegratorGauss<kind>::computeQuadraturePoints(const GhostType & ghost_type) {
  Matrix<Real> & quads = quadrature_points(type, ghost_type);
  quads = GaussIntegrationElement<type>::getQuadraturePoints();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
inline void IntegratorGauss<kind>::computeJacobianOnQuadPointsByElement(
    const Matrix<Real> & node_coords, const Matrix<Real> & quad,
    Vector<Real> & jacobians) const {
  // jacobian
  ElementClass<type>::computeJacobian(quad, node_coords, jacobians);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
IntegratorGauss<kind>::IntegratorGauss(const Mesh & mesh, const ID & id,
                                       const MemoryID & memory_id)
    : Integrator(mesh, id, memory_id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void IntegratorGauss<kind>::checkJacobians(const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = this->quadrature_points(type, ghost_type).cols();

  UInt nb_element = mesh.getConnectivity(type, ghost_type).getSize();

  Real * jacobians_val = jacobians(type, ghost_type).storage();

  for (UInt i = 0; i < nb_element * nb_quadrature_points;
       ++i, ++jacobians_val) {
    if (*jacobians_val < 0)
      AKANTU_DEBUG_ERROR(
          "Negative jacobian computed,"
          << " possible problem in the element node ordering (Quadrature Point "
          << i % nb_quadrature_points << ":" << i / nb_quadrature_points << ":"
          << type << ":" << ghost_type << ")");
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void IntegratorGauss<kind>::computeJacobiansOnIntegrationPoints(
    const Array<Real> & nodes, const Matrix<Real> & quad_points,
    Array<Real> & jacobians, const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points = quad_points.cols();

  UInt nb_element = mesh.getNbElement(type, ghost_type);

  jacobians.resize(nb_element * nb_quadrature_points);

  Array<Real>::vector_iterator jacobians_it =
      jacobians.begin_reinterpret(nb_quadrature_points, nb_element);

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type);

  Array<Real>::const_matrix_iterator x_it =
      x_el.begin(spatial_dimension, nb_nodes_per_element);

  //  Matrix<Real> local_coord(spatial_dimension, nb_nodes_per_element);
  for (UInt elem = 0; elem < nb_element; ++elem, ++jacobians_it, ++x_it) {
    const Matrix<Real> & x = *x_it;
    Vector<Real> & J = *jacobians_it;
    computeJacobianOnQuadPointsByElement<type>(x, quad_points, J);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void IntegratorGauss<kind>::precomputeJacobiansOnQuadraturePoints(
    const Array<Real> & nodes, const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  Array<Real> * jacobians_tmp;
  if (!jacobians.exists(type, ghost_type))
    jacobians_tmp = &jacobians.alloc(0, 1, type, ghost_type);
  else {
    jacobians_tmp = &jacobians(type, ghost_type);
  }

  this->computeJacobiansOnIntegrationPoints<type>(
      nodes, quadrature_points(type, ghost_type), *jacobians_tmp, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type, UInt polynomial_degree>
void IntegratorGauss<kind>::multiplyJacobiansByWeights(
    Array<Real> & jacobians) const {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points =
      GaussIntegrationElement<type, polynomial_degree>::getNbQuadraturePoints();
  UInt nb_element = jacobians.getSize() / nb_quadrature_points;

  Vector<Real> weights =
      GaussIntegrationElement<type, polynomial_degree>::getWeights();

  Array<Real>::vector_iterator jacobians_it =
      jacobians.begin_reinterpret(nb_quadrature_points, nb_element);
  Array<Real>::vector_iterator jacobians_end =
      jacobians.end_reinterpret(nb_quadrature_points, nb_element);

  for (; jacobians_it != jacobians_end; ++jacobians_it) {
    Vector<Real> & J = *jacobians_it;
    J *= weights;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
#if defined(AKANTU_COHESIVE_ELEMENT)
template <>
template <ElementType type>
void IntegratorGauss<_ek_cohesive>::precomputeJacobiansOnQuadraturePoints(
    const Array<Real> & nodes, const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points =
      GaussIntegrationElement<type>::getNbQuadraturePoints();

  UInt nb_element = mesh.getNbElement(type, ghost_type);

  Array<Real> * jacobians_tmp;
  if (!jacobians.exists(type, ghost_type))
    jacobians_tmp = &jacobians.alloc(nb_element * nb_quadrature_points, 1, type,
                                     ghost_type);
  else {
    jacobians_tmp = &jacobians(type, ghost_type);
    jacobians_tmp->resize(nb_element * nb_quadrature_points);
  }

  Array<Real>::vector_iterator jacobians_it =
      jacobians_tmp->begin_reinterpret(nb_quadrature_points, nb_element);

  Vector<Real> weights = GaussIntegrationElement<type>::getWeights();

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type);

  Array<Real>::const_matrix_iterator x_it =
      x_el.begin(spatial_dimension, nb_nodes_per_element);

  UInt nb_nodes_per_subelement = nb_nodes_per_element / 2;
  Matrix<Real> x(spatial_dimension, nb_nodes_per_subelement);

  //  Matrix<Real> local_coord(spatial_dimension, nb_nodes_per_element);
  for (UInt elem = 0; elem < nb_element; ++elem, ++jacobians_it, ++x_it) {

    for (UInt s = 0; s < spatial_dimension; ++s)
      for (UInt n = 0; n < nb_nodes_per_subelement; ++n)
        x(s, n) =
            ((*x_it)(s, n) + (*x_it)(s, n + nb_nodes_per_subelement)) * .5;

    Vector<Real> & J = *jacobians_it;

    if (type == _cohesive_1d_2)
      J(0) = 1;
    else
      computeJacobianOnQuadPointsByElement<type>(x, J);
  }

  AKANTU_DEBUG_OUT();
}
#endif

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
void IntegratorGauss<kind>::integrate(const Array<Real> & in_f,
                                      Array<Real> & intf,
                                      UInt nb_degree_of_freedom,
                                      const Array<Real> & jacobians,
                                      UInt nb_element) const {
  AKANTU_DEBUG_IN();

  UInt nb_points = jacobians.getSize() / nb_element;

  Array<Real>::const_matrix_iterator J_it;
  Array<Real>::matrix_iterator inte_it;
  Array<Real>::const_matrix_iterator f_it;

  intf.resize(nb_element);
  f_it = in_f.begin_reinterpret(nb_degree_of_freedom, nb_points, nb_element);
  inte_it = intf.begin_reinterpret(nb_degree_of_freedom, 1, nb_element);

  for (UInt el = 0; el < nb_element; ++el, ++J_it, ++f_it, ++inte_it) {
    const Matrix<Real> & f = *f_it;
    const Matrix<Real> & J = *J_it;
    Matrix<Real> & inte_f = *inte_it;

    inte_f.mul<false, false>(f, J);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void IntegratorGauss<kind>::integrate(
    const Array<Real> & in_f, Array<Real> & intf, UInt nb_degree_of_freedom,
    const GhostType & ghost_type, const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(jacobians.exists(type, ghost_type),
                      "No jacobians for the type "
                          << jacobians.printType(type, ghost_type));

  const Array<Real> & jac_loc = jacobians(type, ghost_type);
  if (filter_elements != empty_filter) {
    UInt nb_element = filter_elements.getSize();
    Array<Real> * filtered_J = new Array<Real>(0, jac_loc.getNbComponent());
    FEEngine::filterElementalData(mesh, jac_loc, *filtered_J, type, ghost_type,
                                  filter_elements);
    this->integrate(in_f, intf, nb_degree_of_freedom, *filtered_J, nb_element);
    delete filtered_J;
  } else {
    UInt nb_element = mesh.getNbElement(type, ghost_type);
    this->integrate(in_f, intf, nb_degree_of_freedom, jac_loc, nb_element);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type, UInt polynomial_degree>
void IntegratorGauss<kind>::integrate(const Array<Real> & in_f,
                                      Array<Real> & intf,
                                      UInt nb_degree_of_freedom,
                                      const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  Matrix<Real> quads = this->getIntegrationPoints<type, polynomial_degree>();

  Array<Real> jacobians;
  this->computeJacobiansOnIntegrationPoints<type>(mesh.getNodes(), quads,
                                                  jacobians, ghost_type);
  this->multiplyJacobiansByWeights<type, polynomial_degree>(jacobians);

  this->integrate(in_f, intf, nb_degree_of_freedom, jacobians,
                  mesh.getNbElement(type, ghost_type));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
Real IntegratorGauss<kind>::integrate(
    const Array<Real> & in_f, const GhostType & ghost_type,
    const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(jacobians.exists(type, ghost_type),
                      "No jacobians for the type "
                          << jacobians.printType(type, ghost_type));

  Array<Real> intfv(0, 1);
  integrate<type>(in_f, intfv, 1, ghost_type, filter_elements);

  UInt nb_values = intfv.getSize();
  if (nb_values == 0)
    return 0.;

  UInt nb_values_to_sum = nb_values >> 1;

  std::sort(intfv.begin(), intfv.end());

  // as long as the half is not empty
  while (nb_values_to_sum) {
    UInt remaining = (nb_values - 2 * nb_values_to_sum);
    if (remaining)
      intfv(nb_values - 2) += intfv(nb_values - 1);

    // sum to consecutive values and store the sum in the first half
    for (UInt i = 0; i < nb_values_to_sum; ++i) {
      intfv(i) = intfv(2 * i) + intfv(2 * i + 1);
    }

    nb_values = nb_values_to_sum;
    nb_values_to_sum >>= 1;
  }

  AKANTU_DEBUG_OUT();
  return intfv(0);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
void IntegratorGauss<kind>::integrateOnIntegrationPoints(
    const Array<Real> & in_f, Array<Real> & intf, UInt nb_degree_of_freedom,
    const Array<Real> & jacobians, UInt nb_element) const {
  AKANTU_DEBUG_IN();

  UInt nb_points = jacobians.getSize() / nb_element;

  Array<Real>::const_scalar_iterator J_it;
  Array<Real>::vector_iterator inte_it;
  Array<Real>::const_vector_iterator f_it;

  intf.resize(nb_element * nb_points);

  J_it = jacobians.begin();
  f_it = in_f.begin(nb_degree_of_freedom);
  inte_it = intf.begin(nb_degree_of_freedom);

  for (UInt el = 0; el < nb_element; ++el, ++J_it, ++f_it, ++inte_it) {
    const Real & J = *J_it;
    const Vector<Real> & f = *f_it;
    Vector<Real> & inte_f = *inte_it;

    inte_f = f;
    inte_f *= J;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void IntegratorGauss<kind>::integrateOnIntegrationPoints(
    const Array<Real> & in_f, Array<Real> & intf, UInt nb_degree_of_freedom,
    const GhostType & ghost_type, const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(jacobians.exists(type, ghost_type),
                      "No jacobians for the type "
                          << jacobians.printType(type, ghost_type));

  const Array<Real> & jac_loc = this->jacobians(type, ghost_type);

  if (filter_elements != empty_filter) {

    UInt nb_element = filter_elements.getSize();
    Array<Real> * filtered_J = new Array<Real>(0, jac_loc.getNbComponent());
    FEEngine::filterElementalData(mesh, jac_loc, *filtered_J, type, ghost_type,
                                  filter_elements);

    this->integrateOnIntegrationPoints(in_f, intf, nb_degree_of_freedom,
                                       *filtered_J, nb_element);
  } else {
    UInt nb_element = mesh.getNbElement(type, ghost_type);
    this->integrateOnIntegrationPoints(in_f, intf, nb_degree_of_freedom,
                                       jac_loc, nb_element);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
#define INIT_INTEGRATOR(type)                                                  \
  computeQuadraturePoints<type>(ghost_type);                                   \
  precomputeJacobiansOnQuadraturePoints<type>(nodes, ghost_type);              \
  checkJacobians<type>(ghost_type);                                            \
  multiplyJacobiansByWeights<type,                                             \
                             ElementClassProperty<type>::polynomial_degree>(   \
      this->jacobians(type, ghost_type));

template <>
inline void
IntegratorGauss<_ek_regular>::initIntegrator(const Array<Real> & nodes,
                                             const ElementType & type,
                                             const GhostType & ghost_type) {
  AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(INIT_INTEGRATOR);
}

#if defined(AKANTU_COHESIVE_ELEMENT)
template <>
inline void
IntegratorGauss<_ek_cohesive>::initIntegrator(const Array<Real> & nodes,
                                              const ElementType & type,
                                              const GhostType & ghost_type) {
  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(INIT_INTEGRATOR);
}
#endif

#if defined(AKANTU_STRUCTURAL_MECHANICS)
template <>
inline void
IntegratorGauss<_ek_structural>::initIntegrator(const Array<Real> & nodes,
                                                const ElementType & type,
                                                const GhostType & ghost_type) {
  AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(INIT_INTEGRATOR);
}
#endif

#undef INIT_INTEGRATOR
