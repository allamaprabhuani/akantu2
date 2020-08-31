/**
 * @file   geometry_utils.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Mmay 20 2019
 * @date last modification: Mon May 20 2019
 *
 * @brief  Implementation of various utilities needed for contact geometry 
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "geometry_utils.hh"
/* -------------------------------------------------------------------------- */


namespace akantu {
  
/* -------------------------------------------------------------------------- */
void GeometryUtils::normal(const Mesh & mesh, const Array<Real> & positions,
			   const Element & element, Vector<Real> & normal, bool outward) {

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt surface_dimension = spatial_dimension - 1;
  
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);
  Matrix<Real> coords(spatial_dimension, nb_nodes_per_element);

  UInt * elem_val = mesh.getConnectivity(element.type, _not_ghost).storage();
  mesh.extractNodalValuesFromElement(positions, coords.storage(),
                                     elem_val + element.element * nb_nodes_per_element,
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
  default: { AKANTU_ERROR("Unknown dimension : " << spatial_dimension); }
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
void GeometryUtils::covariantBasis(const Mesh & mesh, const Array<Real> & positions,
				   const Element & element, const Vector<Real> & normal,
				   Vector<Real> & natural_coord,
				   Matrix<Real> & tangents) {

  UInt spatial_dimension = mesh.getSpatialDimension();
  
  const ElementType & type = element.type;
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
  UInt * elem_val = mesh.getConnectivity(type, _not_ghost).storage();
  Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);

  mesh.extractNodalValuesFromElement(positions, nodes_coord.storage(),
                                     elem_val + element.element * nb_nodes_per_element,
                                     nb_nodes_per_element, spatial_dimension);

  UInt surface_dimension = spatial_dimension - 1;
  Matrix<Real> dnds(surface_dimension, nb_nodes_per_element);
  
#define GET_SHAPE_DERIVATIVES_NATURAL(type)				\
  ElementClass<type>::computeDNDS(natural_coord, dnds)
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_DERIVATIVES_NATURAL);
#undef GET_SHAPE_DERIVATIVES_NATURAL

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
void GeometryUtils::curvature(const Mesh & mesh, const Array<Real> & positions,
			      const Element & element, const Vector<Real> & natural_coord,
			      Matrix<Real> & curvature) {
  UInt spatial_dimension = mesh.getSpatialDimension();
  auto surface_dimension = spatial_dimension - 1;

  const ElementType & type = element.type;
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
  UInt * elem_val = mesh.getConnectivity(type, _not_ghost).storage();
  
  Matrix<Real> dn2ds2(surface_dimension * surface_dimension,
		      Mesh::getNbNodesPerElement(type));
  
#define GET_SHAPE_SECOND_DERIVATIVES_NATURAL(type)			\
    ElementClass<type>::computeDN2DS2(natural_coord, dn2ds2)
    AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_SECOND_DERIVATIVES_NATURAL);
#undef GET_SHAPE_SECOND_DERIVATIVES_NATURAL

    Matrix<Real> coords(spatial_dimension, nb_nodes_per_element);
    mesh.extractNodalValuesFromElement(positions, coords.storage(),
				       elem_val + element.element * nb_nodes_per_element,
				       nb_nodes_per_element, spatial_dimension);

    curvature.mul<false, true>(coords, dn2ds2);
}

  
/* -------------------------------------------------------------------------- */
UInt GeometryUtils::orthogonalProjection(const Mesh & mesh, const Array<Real> & positions,
					 const Vector<Real> & slave,
					 const Array<Element> & elements,
					 Real & gap, Vector<Real> & natural_projection,
					 Vector<Real> & normal, Real alpha,
					 Real projection_tolerance, 
					 Real extension_tolerance) {

  UInt index = UInt(-1);
  Real min_gap = std::numeric_limits<Real>::max();

  UInt spatial_dimension = mesh.getSpatialDimension();

  UInt nb_same_sides{0};
  UInt nb_boundary_elements{0};

  UInt counter{0};
  for (auto & element : elements) {
    if (!GeometryUtils::isBoundaryElement(mesh, element))
      continue;

    nb_boundary_elements++;

    Vector<Real> normal_ele(spatial_dimension);
    GeometryUtils::normal(mesh, positions, element, normal_ele);
    
    Vector<Real> master(spatial_dimension);
    GeometryUtils::realProjection(mesh, positions, slave, element, normal_ele, master);

    Vector<Real> xi(natural_projection.size());
    GeometryUtils::naturalProjection(mesh, positions, element, master, xi,
				     projection_tolerance);

    // if gap between master projection and slave point is zero, then
    // it means that slave point lies on the master element, hence the
    // normal from master to slave cannot be computed in that case
    
    auto master_to_slave = slave - master;
    Real temp_gap = master_to_slave.norm();
    
    if (temp_gap != 0) 
      master_to_slave /= temp_gap;
    
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
         
    if (variation <= direction_tolerance and
	temp_gap <= min_gap and
	GeometryUtils::isValidProjection(xi, extension_tolerance)) {

      gap     = -temp_gap;
      min_gap = temp_gap;
      index   = counter;
      natural_projection = xi;
      normal = normal_ele;
      
    }

    if(temp_gap == 0 or variation <= direction_tolerance)
      nb_same_sides++;
    
    counter++;
  }

  // if point is not on the same side of all the boundary elements
  // than it is not consider even if the closet master element is
  // found
  if(nb_same_sides != nb_boundary_elements)
    index = UInt(-1);
  
  return index;
}

/* -------------------------------------------------------------------------- */
void GeometryUtils::realProjection(const Mesh & mesh, const Array<Real> & positions,
				   const Vector<Real> & slave, const Element & element,
				   const Vector<Real> & normal, 
				   Vector<Real> & projection) {

  UInt spatial_dimension = mesh.getSpatialDimension();

  const ElementType & type = element.type;
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);

  UInt * elem_val = mesh.getConnectivity(type, _not_ghost).storage();
  Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);

  mesh.extractNodalValuesFromElement(positions, nodes_coord.storage(),
                                     elem_val + element.element * nb_nodes_per_element,
                                     nb_nodes_per_element, spatial_dimension);
  
  Vector<Real> point(nodes_coord(0));
  Real alpha = (slave - point).dot(normal);

  projection = slave - alpha * normal;
}

/* -------------------------------------------------------------------------- */
void GeometryUtils::realProjection(const Mesh & mesh, const Array<Real> & positions,
				   const Element & element, const Vector<Real> & natural_coord,
				   Vector<Real> & projection) {

  UInt spatial_dimension = mesh.getSpatialDimension();

  const ElementType & type = element.type;

  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(element.type);
  
  Vector<Real> shapes(nb_nodes_per_element);  
#define GET_SHAPE_NATURAL(type)		\
  ElementClass<type>::computeShapes(natural_coord, shapes)
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_NATURAL);
#undef GET_SHAPE_NATURAL

  Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);
  UInt * elem_val = mesh.getConnectivity(type, _not_ghost).storage();
  mesh.extractNodalValuesFromElement(positions, nodes_coord.storage(),
                                     elem_val + element.element * nb_nodes_per_element,
                                     nb_nodes_per_element, spatial_dimension);

  projection.mul<false>(nodes_coord, shapes);
}

/* -------------------------------------------------------------------------- */
void GeometryUtils::naturalProjection(const Mesh & mesh, const Array<Real> & positions,
				      const Element & element,
				      Vector<Real> & real_projection,
				      Vector<Real> & natural_projection,
				      Real projection_tolerance) {

  UInt spatial_dimension = mesh.getSpatialDimension();
  
  const ElementType & type = element.type;
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
  UInt * elem_val = mesh.getConnectivity(type, _not_ghost).storage();
  Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);

  mesh.extractNodalValuesFromElement(positions, nodes_coord.storage(),
                                     elem_val + element.element * nb_nodes_per_element,
                                     nb_nodes_per_element, spatial_dimension);

#define GET_NATURAL_COORDINATE(type)                                           \
  ElementClass<type>::inverseMap(real_projection, nodes_coord,                 \
                                 natural_projection, projection_tolerance)
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_NATURAL_COORDINATE);
#undef GET_NATURAL_COORDINATE
}


/* -------------------------------------------------------------------------- */
void GeometryUtils::contravariantBasis(const Matrix<Real> & covariant,
				       Matrix<Real> & contravariant) {
  
  auto inv_A = GeometryUtils::contravariantMetricTensor(covariant);
  contravariant.mul<false, false>(inv_A, covariant); 
}


/* -------------------------------------------------------------------------- */
Matrix<Real> GeometryUtils::covariantMetricTensor(const Matrix<Real> & covariant_bases) {
  
  Matrix<Real> A(covariant_bases.rows(), covariant_bases.rows());
  A.mul<false, true>(covariant_bases, covariant_bases);
  return A;
}
  

/* -------------------------------------------------------------------------- */
Matrix<Real> GeometryUtils::contravariantMetricTensor(const Matrix<Real> & covariant_bases)  {

  auto A = GeometryUtils::covariantMetricTensor(covariant_bases);
  auto inv_A = A.inverse();
  return inv_A;
}

/* -------------------------------------------------------------------------- */
Matrix<Real> GeometryUtils::covariantCurvatureTensor(const Mesh & mesh,
						     const Array<Real> & positions,
						     const Element & element,
						     const Vector<Real> & natural_coord,
						     const Vector<Real> & normal) {

  UInt spatial_dimension = mesh.getSpatialDimension();
  auto surface_dimension = spatial_dimension - 1;

  const ElementType & type = element.type;
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt * elem_val = mesh.getConnectivity(type, _not_ghost).storage();
  
  Matrix<Real> dn2ds2(surface_dimension * surface_dimension,
		      nb_nodes_per_element);
  
#define GET_SHAPE_SECOND_DERIVATIVES_NATURAL(type)			\
  ElementClass<type>::computeDN2DS2(natural_coord, dn2ds2)
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_SECOND_DERIVATIVES_NATURAL);
#undef GET_SHAPE_SECOND_DERIVATIVES_NATURAL

  Matrix<Real> coords(spatial_dimension, nb_nodes_per_element);
  mesh.extractNodalValuesFromElement(positions, coords.storage(),
				     elem_val + element.element * nb_nodes_per_element,
				     nb_nodes_per_element, spatial_dimension);

  Matrix<Real> curvature(spatial_dimension, surface_dimension*surface_dimension);
  curvature.mul<false, true>(coords, dn2ds2);

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

}
