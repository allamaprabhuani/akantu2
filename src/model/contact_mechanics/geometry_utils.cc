/**
 * @file   geometry_utils.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Wed Oct 02 2019
 * @date last modification: Thu Jun 24 2021
 *
 * @brief  Implementation of various utilities needed for contact geometry
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "geometry_utils.hh"
#include "element_class_helper.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void GeometryUtils::normal(const Mesh & mesh, const Array<Real> & positions,
                           const Element & element, Vector<Real> & normal,
                           bool outward) {

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt surface_dimension = spatial_dimension - 1;

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);
  Matrix<Real> coords(spatial_dimension, nb_nodes_per_element);

  UInt * elem_val = mesh.getConnectivity(element.type, _not_ghost).storage();
  mesh.extractNodalValuesFromElement(positions, coords.storage(),
                                     elem_val +
                                         element.element * nb_nodes_per_element,
                                     nb_nodes_per_element, spatial_dimension);

  Matrix<Real> vectors(spatial_dimension, surface_dimension);
  switch (spatial_dimension) {
  case 1: {
    normal[0] = 1;
    break;
  }
  case 2: {
    vectors(0) = Vector<Real>(coords(1)) - Vector<Real>(coords(0));
    Math::normal2(vectors.storage(), normal.storage());
    break;
  }
  case 3: {
    vectors(0) = Vector<Real>(coords(1)) - Vector<Real>(coords(0));
    vectors(1) = Vector<Real>(coords(2)) - Vector<Real>(coords(0));
    Math::normal3(vectors(0).storage(), vectors(1).storage(), normal.storage());
    break;
  }
  default: {
    AKANTU_ERROR("Unknown dimension : " << spatial_dimension);
  }
  }

  // to ensure that normal is always outwards from master surface
  if (outward) {

    const auto & element_to_subelement =
        mesh.getElementToSubelement(element.type)(element.element);

    Vector<Real> outside(spatial_dimension);
    mesh.getBarycenter(element, outside);

    // check if mesh facets exists for cohesive elements contact
    Vector<Real> inside(spatial_dimension);
    if (mesh.isMeshFacets()) {
      mesh.getMeshParent().getBarycenter(element_to_subelement[0], inside);
    } else {
      mesh.getBarycenter(element_to_subelement[0], inside);
    }

    Vector<Real> inside_to_outside = outside - inside;
    auto projection = inside_to_outside.dot(normal);

    if (projection < 0) {
      normal *= -1.0;
    }
  }
}

/* -------------------------------------------------------------------------- */
void GeometryUtils::normal(const Mesh & mesh, const Element & element,
                           Matrix<Real> & tangents, Vector<Real> & normal,
                           bool outward) {
  UInt spatial_dimension = mesh.getSpatialDimension();

  // to ensure that normal is always outwards from master surface we
  // compute a direction vector form inside of element attached to the
  // suurface element

  Vector<Real> inside_to_outside(spatial_dimension);
  if (outward) {

    const auto & element_to_subelement =
        mesh.getElementToSubelement(element.type)(element.element);

    Vector<Real> outside(spatial_dimension);
    mesh.getBarycenter(element, outside);

    // check if mesh facets exists for cohesive elements contact
    Vector<Real> inside(spatial_dimension);
    if (mesh.isMeshFacets()) {
      mesh.getMeshParent().getBarycenter(element_to_subelement[0], inside);
    } else {
      mesh.getBarycenter(element_to_subelement[0], inside);
    }

    inside_to_outside = outside - inside;
    // auto projection = inside_to_outside.dot(normal);

    // if (projection < 0) {
    //  normal *= -1.0;
    //}
  }

  // to ensure that direction of tangents are correct, cross product
  // of tangents should give be in the same direction as of inside to outside
  switch (spatial_dimension) {
  case 2: {
    normal[0] = -tangents(0, 1);
    normal[1] = tangents(0, 0);

    auto ddot = inside_to_outside.dot(normal);
    if (ddot < 0) {
      tangents *= -1.0;
      normal *= -1.0;
    }

    break;
  }
  case 3: {
    auto tang_trans = tangents.transpose();
    auto tang1 = Vector<Real>(tang_trans(0));
    auto tang2 = Vector<Real>(tang_trans(1));

    auto tang1_cross_tang2 = tang1.crossProduct(tang2);
    normal = tang1_cross_tang2 / tang1_cross_tang2.norm();

    auto ddot = inside_to_outside.dot(normal);
    if (ddot < 0) {
      tang_trans(1) *= -1.0;
      normal *= -1.0;
    }

    tangents = tang_trans.transpose();

    break;
  }
  default:
    break;
  }
}

/* -------------------------------------------------------------------------- */
void GeometryUtils::covariantBasis(const Mesh & mesh,
                                   const Array<Real> & positions,
                                   const Element & element,
                                   const Vector<Real> & normal,
                                   Vector<Real> & natural_coord,
                                   Matrix<Real> & tangents) {
  UInt spatial_dimension = mesh.getSpatialDimension();

  const ElementType type = element.type;
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
  UInt * elem_val = mesh.getConnectivity(type, _not_ghost).storage();
  Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);

  mesh.extractNodalValuesFromElement(positions, nodes_coord.storage(),
                                     elem_val +
                                         element.element * nb_nodes_per_element,
                                     nb_nodes_per_element, spatial_dimension);

  auto && dnds = ElementClassHelper<_ek_regular>::getDNDS(natural_coord, type);
  tangents.mul<false, true>(dnds, nodes_coord);

  auto temp_tangents = tangents.transpose();
  for (UInt i = 0; i < spatial_dimension - 1; ++i) {
    auto temp = Vector<Real>(temp_tangents(i));
    temp_tangents(i) = temp.normalize();
  }

  tangents = temp_tangents.transpose();

  // to ensure that direction of tangents are correct, cross product
  // of tangents should give the normal vector computed earlier
  switch (spatial_dimension) {
  case 2: {
    Vector<Real> e_z(3);
    e_z[0] = 0.;
    e_z[1] = 0.;
    e_z[2] = 1.;

    Vector<Real> tangent(3);
    tangent[0] = tangents(0, 0);
    tangent[1] = tangents(0, 1);
    tangent[2] = 0.;

    auto exp_normal = e_z.crossProduct(tangent);

    auto ddot = normal.dot(exp_normal);
    if (ddot < 0) {
      tangents *= -1.0;
    }
    break;
  }
  case 3: {
    auto tang_trans = tangents.transpose();
    auto tang1 = Vector<Real>(tang_trans(0));
    auto tang2 = Vector<Real>(tang_trans(1));

    auto tang1_cross_tang2 = tang1.crossProduct(tang2);
    auto exp_normal = tang1_cross_tang2 / tang1_cross_tang2.norm();

    auto ddot = normal.dot(exp_normal);
    if (ddot < 0) {
      tang_trans(1) *= -1.0;
    }

    tangents = tang_trans.transpose();
    break;
  }
  default:
    break;
  }
}

/* -------------------------------------------------------------------------- */
void GeometryUtils::covariantBasis(const Mesh & mesh,
                                   const Array<Real> & positions,
                                   const Element & element,
                                   Vector<Real> & natural_coord,
                                   Matrix<Real> & tangents) {

  UInt spatial_dimension = mesh.getSpatialDimension();

  const ElementType & type = element.type;
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
  UInt * elem_val = mesh.getConnectivity(type, _not_ghost).storage();
  Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);

  mesh.extractNodalValuesFromElement(positions, nodes_coord.storage(),
                                     elem_val +
                                         element.element * nb_nodes_per_element,
                                     nb_nodes_per_element, spatial_dimension);

  auto && dnds = ElementClassHelper<_ek_regular>::getDNDS(natural_coord, type);
  tangents.mul<false, true>(dnds, nodes_coord);

  auto temp_tangents = tangents.transpose();
  for (UInt i = 0; i < spatial_dimension - 1; ++i) {
    auto temp = Vector<Real>(temp_tangents(i));
    temp_tangents(i) = temp.normalize();
  }

  tangents = temp_tangents.transpose();
}

/* -------------------------------------------------------------------------- */
void GeometryUtils::curvature(const Mesh & mesh, const Array<Real> & positions,
                              const Element & element,
                              const Vector<Real> & natural_coord,
                              Matrix<Real> & curvature) {
  UInt spatial_dimension = mesh.getSpatialDimension();
  const ElementType & type = element.type;

  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
  UInt * elem_val = mesh.getConnectivity(type, _not_ghost).storage();

  auto && d2nds2 =
      ElementClassHelper<_ek_regular>::getD2NDS2(natural_coord, type);

  Matrix<Real> coords(spatial_dimension, nb_nodes_per_element);
  mesh.extractNodalValuesFromElement(positions, coords.storage(),
                                     elem_val +
                                         element.element * nb_nodes_per_element,
                                     nb_nodes_per_element, spatial_dimension);

  curvature.mul<false, true>(coords, d2nds2);
}

/* -------------------------------------------------------------------------- */
UInt GeometryUtils::orthogonalProjection(
    const Mesh & mesh, const Array<Real> & positions,
    const Vector<Real> & slave, const Array<Element> & elements, Real & gap,
    Vector<Real> & natural_projection, Vector<Real> & normal, Real alpha,
    UInt max_iterations, Real projection_tolerance, Real extension_tolerance) {

  UInt index = UInt(-1);
  Real min_gap = std::numeric_limits<Real>::max();

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt surface_dimension = spatial_dimension - 1;

  UInt nb_same_sides{0};
  UInt nb_boundary_elements{0};

  UInt counter{0};

  const auto & contact_group = mesh.getElementGroup("contact_surface");

  for (const auto & element : elements) {
    // filter out elements which are not there in the element group
    // contact surface created by the surface selector and is stored
    // in the mesh or mesh_facet, if a element is not there it
    // returnas UInt(-1)

    const auto & elements_of_type = contact_group.getElements(element.type);
    if (elements_of_type.find(element.element) == UInt(-1)) {
      continue;
    }

    nb_boundary_elements++;

    // find the natural coordinate corresponding to the minimum gap
    // between slave node and master element

    Vector<Real> master(spatial_dimension);

    Vector<Real> xi(natural_projection.size());
    GeometryUtils::naturalProjection(mesh, positions, element, slave, master,
                                     xi, max_iterations, projection_tolerance);

    Matrix<Real> tangent_ele(surface_dimension, spatial_dimension);
    GeometryUtils::covariantBasis(mesh, positions, element, xi, tangent_ele);

    Vector<Real> normal_ele(spatial_dimension);
    GeometryUtils::normal(mesh, element, tangent_ele, normal_ele);

    // if gap between master projection and slave point is zero, then
    // it means that slave point lies on the master element, hence the
    // normal from master to slave cannot be computed in that case

    auto master_to_slave = slave - master;
    Real temp_gap = master_to_slave.norm();

    if (temp_gap != 0) {
      master_to_slave /= temp_gap;
    }

    // also the slave point should lie inside the master surface in
    // case of explicit or outside in case of implicit, one way to
    // check that is by checking the dot product of normal at each
    // master element, if the values of all dot product is same then
    // the slave point lies on the same side of each master element

    // A alpha parameter is introduced which is 1 in case of explicit
    // and -1 in case of implicit, therefor the variation (dot product
    // + alpha) should be close to zero (within tolerance) for both
    // cases

    Real direction_tolerance = 1e-8;
    auto product = master_to_slave.dot(normal_ele);
    auto variation = std::abs(product + alpha);

    if (variation <= direction_tolerance and temp_gap <= min_gap and
        GeometryUtils::isValidProjection(xi, extension_tolerance)) {

      gap = -temp_gap;
      min_gap = temp_gap;
      index = counter;
      natural_projection = xi;
      normal = normal_ele;
    }

    if (temp_gap == 0 or variation <= direction_tolerance) {
      nb_same_sides++;
    }

    counter++;
  }

  // if point is not on the same side of all the boundary elements
  // than it is not consider even if the closet master element is
  // found
  if (nb_same_sides != nb_boundary_elements) {
    index = UInt(-1);
  }

  return index;
}

/* -------------------------------------------------------------------------- */
UInt GeometryUtils::orthogonalProjection(
    const Mesh & mesh, const Array<Real> & positions,
    const Vector<Real> & slave, const Array<Element> & elements, Real & gap,
    Vector<Real> & natural_projection, Vector<Real> & normal,
    Matrix<Real> & tangent, Real /*alpha*/, UInt max_iterations,
    Real projection_tolerance, Real extension_tolerance) {

  UInt index = UInt(-1);
  Real min_gap = std::numeric_limits<Real>::max();

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt surface_dimension = spatial_dimension - 1;

  const auto & contact_group = mesh.getElementGroup("contact_surface");

  for (auto && tuple : enumerate(elements)) {
    auto & counter = std::get<0>(tuple);
    const auto & element = std::get<1>(tuple);
    // filter out elements which are not there in the element group
    // contact surface created by the surface selector and is stored
    // in the mesh or mesh_facet, if a element is not there it
    // returnas UInt(-1)

    const auto & elements_of_type = contact_group.getElements(element.type);
    if (elements_of_type.find(element.element) == UInt(-1)) {
      continue;
    }

    Vector<Real> master(spatial_dimension);

    Vector<Real> xi_ele(natural_projection.size());
    GeometryUtils::naturalProjection(mesh, positions, element, slave, master,
                                     xi_ele, max_iterations,
                                     projection_tolerance);

    Matrix<Real> tangent_ele(surface_dimension, spatial_dimension);
    GeometryUtils::covariantBasis(mesh, positions, element, xi_ele,
                                  tangent_ele);

    Vector<Real> normal_ele(spatial_dimension);
    GeometryUtils::normal(mesh, element, tangent_ele, normal_ele);

    // if gap between master projection and slave point is zero, then
    // it means that slave point lies on the master element, hence the
    // normal from master to slave cannot be computed in that case

    auto master_to_slave = slave - master;
    Real temp_gap = master_to_slave.norm();

    if (temp_gap != 0) {
      master_to_slave /= temp_gap;
    }

    // A alpha parameter is introduced which is 1 in case of explicit
    // and -1 in case of implicit, therefor the variation (dot product
    // + alpha) should be close to zero (within tolerance) for both
    // cases

    auto product = master_to_slave.dot(normal_ele);

    if (product < 0 and temp_gap <= min_gap and
        GeometryUtils::isValidProjection(xi_ele, extension_tolerance)) {

      gap = -temp_gap;
      min_gap = temp_gap;
      index = counter;
      natural_projection = xi_ele;
      normal = normal_ele;
      tangent = tangent_ele;
    }
  }

  return index;
}

/* -------------------------------------------------------------------------- */
void GeometryUtils::realProjection(const Mesh & mesh,
                                   const Array<Real> & positions,
                                   const Vector<Real> & slave,
                                   const Element & element,
                                   const Vector<Real> & normal,
                                   Vector<Real> & projection) {

  UInt spatial_dimension = mesh.getSpatialDimension();

  const ElementType & type = element.type;
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);

  UInt * elem_val = mesh.getConnectivity(type, _not_ghost).storage();
  Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);

  mesh.extractNodalValuesFromElement(positions, nodes_coord.storage(),
                                     elem_val +
                                         element.element * nb_nodes_per_element,
                                     nb_nodes_per_element, spatial_dimension);

  Vector<Real> point(nodes_coord(0));
  Real alpha = (slave - point).dot(normal);

  projection = slave - alpha * normal;
}

/* -------------------------------------------------------------------------- */
void GeometryUtils::realProjection(const Mesh & mesh,
                                   const Array<Real> & positions,
                                   const Element & element,
                                   const Vector<Real> & natural_coord,
                                   Vector<Real> & projection) {

  auto spatial_dimension = mesh.getSpatialDimension();

  const auto & type = element.type;

  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);
  auto shapes =
      ElementClassHelper<_ek_regular>::getN(natural_coord, element.type);

  Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);
  UInt * elem_val = mesh.getConnectivity(type, _not_ghost).storage();
  mesh.extractNodalValuesFromElement(positions, nodes_coord.storage(),
                                     elem_val +
                                         element.element * nb_nodes_per_element,
                                     nb_nodes_per_element, spatial_dimension);

  projection.mul<false>(nodes_coord, shapes);
}

/* -------------------------------------------------------------------------- */
void GeometryUtils::naturalProjection(
    const Mesh & mesh, const Array<Real> & positions, const Element & element,
    const Vector<Real> & slave_coords, Vector<Real> & master_coords,
    Vector<Real> & natural_projection, UInt max_iterations,
    Real projection_tolerance) {

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt surface_dimension = spatial_dimension - 1;

  const ElementType & type = element.type;
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
  UInt * elem_val = mesh.getConnectivity(type, _not_ghost).storage();
  Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);

  mesh.extractNodalValuesFromElement(positions, nodes_coord.storage(),
                                     elem_val +
                                         element.element * nb_nodes_per_element,
                                     nb_nodes_per_element, spatial_dimension);

  // initial guess
  natural_projection.zero();

  // obhjective function  computed on the natural_guess
  Matrix<Real> f(surface_dimension, 1);

  // jacobian matrix computed on the natural_guess
  Matrix<Real> J(surface_dimension, surface_dimension);

  // Jinv = J^{-1}
  Matrix<Real> Jinv(surface_dimension, surface_dimension);

  // dxi = \xi_{k+1} - \xi_{k} in the iterative process
  Matrix<Real> dxi(surface_dimension, 1);

  // gradient at natural projection
  Matrix<Real> gradient(surface_dimension, spatial_dimension);

  // second derivative at natural peojection
  Matrix<Real> double_gradient(surface_dimension, surface_dimension);

  // second derivative of shape function at natural projection
  Matrix<Real> d2nds2(surface_dimension * surface_dimension,
                      nb_nodes_per_element);

  auto compute_double_gradient = [&d2nds2, &nodes_coord, surface_dimension,
                                  spatial_dimension](UInt & alpha,
                                                     UInt & beta) {
    auto index = alpha * surface_dimension + beta;
    Vector<Real> d_alpha_beta(spatial_dimension);

    auto d2nds2_transpose = d2nds2.transpose();
    Vector<Real> d2nds2_alpha_beta(d2nds2_transpose(index));

    d_alpha_beta.mul<false>(nodes_coord, d2nds2_alpha_beta);

    return d_alpha_beta;
  };

  /* --------------------------- */
  /* init before iteration loop  */
  /* --------------------------- */
  // do interpolation
  auto update_f = [&f, &master_coords, &natural_projection, &nodes_coord,
                   &slave_coords, &gradient, surface_dimension,
                   spatial_dimension, type]() {
    // compute real coords on natural projection
    auto && shapes =
        ElementClassHelper<_ek_regular>::getN(natural_projection, type);

    master_coords.mul<false>(nodes_coord, shapes);

    auto distance = slave_coords - master_coords;

    // first derivative of shape function at natural projection
    auto && dnds =
        ElementClassHelper<_ek_regular>::getDNDS(natural_projection, type);
    gradient.mul<false, true>(dnds, nodes_coord);

    // gradient transpose at natural projection
    Matrix<Real> gradient_transpose(surface_dimension, spatial_dimension);

    gradient_transpose = gradient.transpose();

    // loop over surface dimensions
    for (auto alpha : arange(surface_dimension)) {
      Vector<Real> gradient_alpha(gradient_transpose(alpha));
      f(alpha, 0) = -2. * gradient_alpha.dot(distance);
    }

    // compute initial error
    auto error = f.norm<L_2>();
    return error;
  };

  auto projection_error = update_f();

  /* --------------------------- */
  /* iteration loop              */
  /* --------------------------- */
  UInt iterations{0};
  while (projection_tolerance < projection_error and
         iterations < max_iterations) {

    // compute covariant components of metric tensor
    auto a = GeometryUtils::covariantMetricTensor(gradient);

    // computing second derivative at natural projection
    d2nds2 =
        ElementClassHelper<_ek_regular>::getD2NDS2(natural_projection, type);

    // real coord - physical guess
    auto distance = slave_coords - master_coords;

    // computing Jacobian J
    for (auto alpha : arange(surface_dimension)) {
      for (auto beta : arange(surface_dimension)) {
        auto dgrad_alpha_beta = compute_double_gradient(alpha, beta);
        J(alpha, beta) = 2. * (a(alpha, beta) - dgrad_alpha_beta.dot(distance));
      }
    }

    Jinv.inverse(J);

    // compute increment
    dxi.mul<false, false>(Jinv, f, -1.0);

    // update our guess
    natural_projection += Vector<Real>(dxi(0));

    projection_error = update_f();
    iterations++;
  }
}

/* -------------------------------------------------------------------------- */
void GeometryUtils::contravariantBasis(const Matrix<Real> & covariant,
                                       Matrix<Real> & contravariant) {

  auto inv_A = GeometryUtils::contravariantMetricTensor(covariant);
  contravariant.mul<false, false>(inv_A, covariant);
}

/* -------------------------------------------------------------------------- */
Matrix<Real>
GeometryUtils::covariantMetricTensor(const Matrix<Real> & covariant_bases) {
  Matrix<Real> A(covariant_bases.rows(), covariant_bases.rows());
  A.mul<false, true>(covariant_bases, covariant_bases);
  return A;
}

/* -------------------------------------------------------------------------- */
Matrix<Real>
GeometryUtils::contravariantMetricTensor(const Matrix<Real> & covariant_bases) {
  auto A = GeometryUtils::covariantMetricTensor(covariant_bases);
  auto inv_A = A.inverse();
  return inv_A;
}

/* -------------------------------------------------------------------------- */
Matrix<Real> GeometryUtils::covariantCurvatureTensor(
    const Mesh & mesh, const Array<Real> & positions, const Element & element,
    const Vector<Real> & natural_coord, const Vector<Real> & normal) {

  UInt spatial_dimension = mesh.getSpatialDimension();
  auto surface_dimension = spatial_dimension - 1;

  const ElementType & type = element.type;
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt * elem_val = mesh.getConnectivity(type, _not_ghost).storage();

  auto && d2nds2 =
      ElementClassHelper<_ek_regular>::getD2NDS2(natural_coord, type);

  Matrix<Real> coords(spatial_dimension, nb_nodes_per_element);
  mesh.extractNodalValuesFromElement(positions, coords.storage(),
                                     elem_val +
                                         element.element * nb_nodes_per_element,
                                     nb_nodes_per_element, spatial_dimension);

  Matrix<Real> curvature(spatial_dimension,
                         surface_dimension * surface_dimension);
  curvature.mul<false, true>(coords, d2nds2);

  Matrix<Real> H(surface_dimension, surface_dimension);

  UInt i = 0;
  for (auto alpha : arange(surface_dimension)) {
    for (auto beta : arange(surface_dimension)) {
      Vector<Real> temp(curvature(i));
      H(alpha, beta) = temp.dot(normal);
      i++;
    }
  }

  return H;
}

} // namespace akantu
