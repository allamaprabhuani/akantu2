/**
 * Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "element_class_helper.hh"
#include "geometry_utils.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_GEOMETRY_UTILS_INLINE_IMPL_CC__
#define __AKANTU_GEOMETRY_UTILS_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline bool GeometryUtils::isBoundaryElement(const Mesh & mesh,
                                             const Element & subelement) {
  const auto & element_to_subelement =
      mesh.getElementToSubelement(subelement.type)(subelement.element);

  // for regular boundary elements when akantu::SurfaceSelector is set to
  // physical surfaces, the mesh contains only 1 element attached to a
  // boundary sub-element
  if (element_to_subelement.size() == 1 and
      element_to_subelement[0].kind() == _ek_regular) {
    return true;
  }

  // for cohesive interface elements when akantu::SurfaceSelector is set
  // either cohesive surface selector or all surface selector, in this
  // case mesh passed is actually mesh_facet and for boundary or
  // cohesive  interface 2 elements are associated to a sub-element
  // we want only one regular element attached to the sub-element

  Int nb_elements_regular{0};
  // Int nb_elements_cohesive{0};

  for (auto elem : element_to_subelement) {
    if (elem == ElementNull) {
      continue;
    }

    if (elem.kind() == _ek_regular) {
      ++nb_elements_regular;
    }

    // if (elem.kind() == _ek_cohesive) {
    //   ++nb_elements_cohesive;
    // }
  }

  Int nb_elements = element_to_subelement.size();
  return nb_elements_regular < nb_elements;
}

/* -------------------------------------------------------------------------- */
template <class Derived>
inline bool
GeometryUtils::isValidProjection(const Eigen::MatrixBase<Derived> & projection,
                                 Real extension_tolerance) {
  Int nb_xi_inside = 0;
  for (auto xi : projection) {
    if (xi >= -1.0 - extension_tolerance and xi <= 1.0 + extension_tolerance) {
      nb_xi_inside++;
    }
  }

  return nb_xi_inside == projection.size();
}

/* -------------------------------------------------------------------------- */
inline Vector<Real> GeometryUtils::outsideDirection(const Mesh & mesh,
                                                    const Element & element) {
  const auto & element_to_subelement = mesh.getElementToSubelement()(element);

  Vector<Real> outside = mesh.getBarycenter(element);

  // check if mesh facets exists for cohesive elements contact
  Vector<Real> inside;
  if (mesh.isMeshFacets()) {
    inside = mesh.getMeshParent().getBarycenter(element_to_subelement[0]);
  } else {
    inside = mesh.getBarycenter(element_to_subelement[0]);
  }

  return (outside - inside);
}

/* -------------------------------------------------------------------------- */
template <class Derived>
Vector<Real> GeometryUtils::normal(const Mesh & mesh,
                                   const Eigen::MatrixBase<Derived> & coords,
                                   const Element & element, bool outward) {
  Int spatial_dimension = coords.rows();
  Vector<Real> normal(spatial_dimension);

  switch (spatial_dimension) {
  case 1: {
    normal[0] = 1;
    break;
  }
  case 2: {
    normal = Math::normal(coords(1) - coords(0));
    break;
  }
  case 3: {
    normal = Math::normal(coords(1) - coords(0), coords(2) - coords(0));
    break;
  }
  default: {
    AKANTU_ERROR("Unknown dimension : " << spatial_dimension);
  }
  }

  // to ensure that normal is always outwards from master surface
  if (outward) {
    auto projection = outsideDirection(mesh, element).dot(normal);
    if (projection < 0) {
      normal *= -1.0;
    }
  }
  return normal;
}

/* -------------------------------------------------------------------------- */
template <class Derived>
Vector<Real> GeometryUtils::normal(const Mesh & mesh, const Element & element,
                                   Eigen::MatrixBase<Derived> & tangents,
                                   bool outward) {
  auto spatial_dimension = mesh.getSpatialDimension();
  // to ensure that normal is always outwards from master surface we
  // compute a direction vector form inside of element attached to the
  // suurface element

  Vector<Real> normal(spatial_dimension);
  // to ensure that direction of tangents are correct, cross product
  // of tangents should give be in the same direction as of inside to outside
  switch (spatial_dimension) {
  case 2: {
    normal(0) = -tangents(1, 0);
    normal(1) = tangents(0, 0);
    break;
  }
  case 3: {
    VectorProxy<Real, 3> tangent1(tangents(0).data());
    VectorProxy<Real, 3> tangent2(tangents(1).data());
    normal = (tangent1.cross(tangent2)).normalized();
    break;
  }
  default:
    break;
  }

  if (outward) {
    auto ddot = outsideDirection(mesh, element).dot(normal);
    if (ddot < 0) {
      tangents *= -1.0;
      normal *= -1.0;
    }
  }
  return normal;
}

/* -------------------------------------------------------------------------- */
template <class Derived1, class Derived2>
inline Matrix<Real>
GeometryUtils::covariantBasis(const Eigen::MatrixBase<Derived1> & coords,
                              const Element & element,
                              Eigen::MatrixBase<Derived2> & natural_coord) {
  auto && dnds =
      ElementClassHelper<_ek_regular>::getDNDS(natural_coord, element.type);

  Matrix<Real> tangents_transpose = coords * dnds.transpose();
  for (auto && vect : tangents_transpose) {
    vect = vect.normalized();
  }

  return tangents_transpose;
}

/* -------------------------------------------------------------------------- */
template <class Derived1, class Derived2, class Derived3>
inline Matrix<Real>
GeometryUtils::covariantBasis(const Eigen::MatrixBase<Derived1> & coords,
                              const Element & element,
                              const Eigen::MatrixBase<Derived2> & normal,
                              Eigen::MatrixBase<Derived3> & natural_coord) {
  auto tangents = covariantBasis(coords, element, natural_coord);

  // to ensure that direction of tangents are correct, cross product
  // of tangents should give the normal vector computed earlier
  Int spatial_dimension = coords.rows();
  Vector<Real, 3> exp_normal;

  switch (spatial_dimension) {
  case 2: {
    Vector<Real, 3> e_z{0., 0., 1.};

    Vector<Real, 3> tangent;
    tangent[0] = tangents(0, 0);
    tangent[1] = tangents(1, 0);
    tangent[2] = 0.;

    exp_normal = e_z.cross(tangent);
    break;
  }
  case 3: {
    VectorProxy<Real, 3> tangent1(tangents(0).data());
    VectorProxy<Real, 3> tangent2(tangents(1).data());
    exp_normal = (tangent1.cross(tangent2)).normalized();
    break;
  }
  default:
    AKANTU_TO_IMPLEMENT();
  }

  auto ddot = normal.dot(exp_normal);
  if (ddot < 0) {
    tangents(1) *= -1.0;
  }

  return tangents;
}

/* -------------------------------------------------------------------------- */
template <class Derived1, class Derived2>
inline Matrix<Real>
GeometryUtils::curvature(const Eigen::MatrixBase<Derived1> & coords,
                         const Element & element,
                         const Eigen::MatrixBase<Derived2> & natural_coord) {
  auto && d2nds2 =
      ElementClassHelper<_ek_regular>::getD2NDS2(natural_coord, element.type);
  return coords * d2nds2.transpose();
}

/* -------------------------------------------------------------------------- */
template <class Derived1, class Derived2, class Derived3, class Derived4,
          class ElementList>
Element GeometryUtils::orthogonalProjection(
    const Mesh & mesh, const Array<Real> & positions,
    const Eigen::MatrixBase<Derived1> & slave, const ElementList & elements,
    Real & gap, Eigen::MatrixBase<Derived2> & natural_projection,
    Eigen::MatrixBase<Derived3> & normal, Eigen::MatrixBase<Derived4> & tangent,
    Real /*alpha*/, Int max_iterations, Real projection_tolerance,
    Real extension_tolerance) {

  auto found_element = ElementNull;
  auto min_gap = std::numeric_limits<Real>::max();

  const auto & contact_group = mesh.getElementGroup("contact_surface");

  for (auto && element : elements) {
    // filter out elements which are not there in the element group
    // contact surface created by the surface selector and is stored
    // in the mesh or mesh_facet, if a element is not there it
    // returnas UInt(-1)

    const auto & elements_of_type = contact_group.getElements(element.type);
    if (elements_of_type.find(element.element) == -1) {
      continue;
    }

    auto coords = mesh.extractNodalValuesFromElement(positions, element);

    auto && [xi_ele, master] = GeometryUtils::naturalProjection(
        coords, element, slave, max_iterations, projection_tolerance);

    auto && tangent_ele =
        GeometryUtils::covariantBasis(coords, element, xi_ele);

    auto && normal_ele = GeometryUtils::normal(mesh, element, tangent_ele);

    // if gap between master projection and slave point is zero, then
    // it means that slave point lies on the master element, hence the
    // normal from master to slave cannot be computed in that case
    auto master_to_slave = (slave - master).eval();
    auto temp_gap = master_to_slave.norm();

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
      found_element = element;
      natural_projection = xi_ele;
      normal = normal_ele;
      tangent = tangent_ele;
    }
  }

  return found_element;
}

/* -------------------------------------------------------------------------- */
template <class Derived1, class Derived2, class Derived3>
Vector<Real>
GeometryUtils::realProjection(const Eigen::MatrixBase<Derived1> & coords,
                              const Eigen::MatrixBase<Derived2> & slave,
                              const Eigen::MatrixBase<Derived3> & normal) {
  auto alpha = (slave - coords(0)).dot(normal);
  return slave - alpha * normal;
}

/* -------------------------------------------------------------------------- */
template <class Derived1, class Derived2>
Vector<Real> GeometryUtils::realProjection(
    const Eigen::MatrixBase<Derived1> & coords, const Element & element,
    const Eigen::MatrixBase<Derived2> & natural_coord) {
  auto shapes =
      ElementClassHelper<_ek_regular>::getN(natural_coord, element.type);
  return coords * shapes;
}

/* -------------------------------------------------------------------------- */
template <class Derived1, class Derived2>
std::pair<Vector<Real>, Vector<Real>> GeometryUtils::naturalProjection(
    const Eigen::MatrixBase<Derived1> & coords, const Element & element,
    const Eigen::MatrixBase<Derived2> & slave_coords, Int max_iterations,
    Real projection_tolerance) {

  auto spatial_dimension = coords.rows();
  auto surface_dimension = spatial_dimension - 1;

  auto type = element.type;

  Vector<Real> master_coords(spatial_dimension);
  Vector<Real> natural_projection(surface_dimension);

  // initial guess
  natural_projection.zero();

  // obhjective function  computed on the natural_guess
  Vector<Real> f(surface_dimension);

  // jacobian matrix computed on the natural_guess
  Matrix<Real> J(surface_dimension, surface_dimension);

  // dxi = \xi_{k+1} - \xi_{k} in the iterative process
  Vector<Real> dxi(surface_dimension);

  // gradient at natural projection
  Matrix<Real> gradient(surface_dimension, spatial_dimension);

  // second derivative at natural peojection
  Matrix<Real> double_gradient(surface_dimension, surface_dimension);

  // second derivative of shape function at natural projection
  Matrix<Real> d2nds2(surface_dimension * surface_dimension, coords.cols());

  auto compute_double_gradient = [&d2nds2, &coords, surface_dimension,
                                  spatial_dimension](Int & alpha, Int & beta) {
    auto index = alpha * surface_dimension + beta;
    Vector<Real> d_alpha_beta(spatial_dimension);

    d_alpha_beta = coords * d2nds2.transpose()(index);

    return d_alpha_beta;
  };

  /* --------------------------- */
  /* init before iteration loop  */
  /* --------------------------- */
  // do interpolation
  auto update_f = [&f, &master_coords, &natural_projection, &coords,
                   &slave_coords, &gradient, surface_dimension, type]() {
    // compute real coords on natural projection
    auto && shapes =
        ElementClassHelper<_ek_regular>::getN(natural_projection, type);

    master_coords = coords * shapes;
    auto distance = slave_coords - master_coords;

    // first derivative of shape function at natural projection
    auto && dnds =
        ElementClassHelper<_ek_regular>::getDNDS(natural_projection, type);
    gradient = dnds * coords.transpose();

    // loop over surface dimensions
    for (auto alpha : arange(surface_dimension)) {
      f(alpha) = -2. * gradient.transpose()(alpha).dot(distance);
    }

    // compute initial error
    return f.norm();
  };

  auto projection_error = update_f();

  /* --------------------------- */
  /* iteration loop              */
  /* --------------------------- */
  Int iterations{0};
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

    // compute increment
    dxi = -1 * J.inverse() * f;

    // update our guess
    natural_projection += dxi;

    projection_error = update_f();
    iterations++;
  }

  return std::make_pair(natural_projection, master_coords);
}

/* -------------------------------------------------------------------------- */
template <class Derived>
Matrix<Real> GeometryUtils::contravariantBasis(
    const Eigen::MatrixBase<Derived> & covariant) {
  auto && inv_A = GeometryUtils::contravariantMetricTensor(covariant);
  return inv_A * covariant;
}

/* -------------------------------------------------------------------------- */
template <class Derived>
Matrix<Real> GeometryUtils::covariantMetricTensor(
    const Eigen::MatrixBase<Derived> & covariant_bases) {
  auto A = covariant_bases.transpose() * covariant_bases;
  return A;
}

/* -------------------------------------------------------------------------- */
template <class Derived>
Matrix<Real> GeometryUtils::contravariantMetricTensor(
    const Eigen::MatrixBase<Derived> & covariant_bases) {
  Matrix<Real> A_inv = GeometryUtils::covariantMetricTensor(covariant_bases);
  return A_inv.inverse();
}

/* -------------------------------------------------------------------------- */
template <class Derived1, class Derived2, class Derived3>
Matrix<Real> GeometryUtils::covariantCurvatureTensor(
    const Eigen::MatrixBase<Derived1> & coords, const Element & element,
    const Eigen::MatrixBase<Derived2> & natural_coord,
    const Eigen::MatrixBase<Derived3> & normal) {

  auto spatial_dimension = coords.rows();
  auto surface_dimension = spatial_dimension - 1;

  auto type = element.type;

  auto && d2nds2 =
      ElementClassHelper<_ek_regular>::getD2NDS2(natural_coord, type);

  Matrix<Real> curvature = coords * d2nds2.transpose();
  Matrix<Real> H(surface_dimension, surface_dimension);

  Int i = 0;
  for (auto alpha : arange(surface_dimension)) {
    for (auto beta : arange(surface_dimension)) {
      H(alpha, beta) = curvature(i).dot(normal);
      i++;
    }
  }

  return H;
}

} // namespace akantu

#endif
