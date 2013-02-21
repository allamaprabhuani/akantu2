/**
 * @file   element_class_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Nov 14 16:40:51 2012
 *
 * @brief Implementation of the inline templated function of the element class
 * descriptions
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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


/* -------------------------------------------------------------------------- */
/* GaussIntegrationElement                                                    */
/* -------------------------------------------------------------------------- */
template<ElementType element_type>
const types::Matrix<Real> GaussIntegrationElement<element_type>::getQuadraturePoints() {
  types::Matrix<Real> quads(quad,
			    ElementClass<element_type>::getNaturalSpaceDimension(),
			    nb_quadrature_points);
  return quads;
}

/* -------------------------------------------------------------------------- */
template<ElementType element_type>
const types::Vector<Real> GaussIntegrationElement<element_type>::getWeights() {
  return types::Vector<Real>(weights, nb_quadrature_points);
}

/* -------------------------------------------------------------------------- */
/* GeometricalElement                                                         */
/* -------------------------------------------------------------------------- */
template<GeometricalType geometrical_type, GeometricalShapeType shape>
inline const types::Matrix<UInt>
GeometricalElement<geometrical_type, shape>::getFacetLocalConnectivityPerElement() {
  return types::Matrix<UInt>(facet_connectivity, nb_facets, nb_nodes_per_facet);
}

/* -------------------------------------------------------------------------- */
template<GeometricalType geometrical_type, GeometricalShapeType shape>
inline bool
GeometricalElement<geometrical_type, shape>::contains(const types::Vector<Real> & coords) {
  return GeometricalShapeContains<shape>::contains(coords);
}

/* -------------------------------------------------------------------------- */
template<>
inline bool
GeometricalShapeContains<_gst_point>::contains(const types::Vector<Real> & coords) {
  return (coords(0) < std::numeric_limits<Real>::epsilon());
}

/* -------------------------------------------------------------------------- */
template<>
inline bool
GeometricalShapeContains<_gst_square>::contains(const types::Vector<Real> & coords) {
  bool in = true;
  for (UInt i = 0; i < coords.size() && in; ++i)
    in &= ((coords(i) >= -1.) && (coords(i) <= 1.));
  return in;
}

/* -------------------------------------------------------------------------- */
template<>
inline bool
GeometricalShapeContains<_gst_triangle>::contains(const types::Vector<Real> & coords) {
  bool in = true;
  Real sum = 0;
    for (UInt i = 0; (i < coords.size()) && in; ++i) {
    in &= ((coords(i) >= 0) && (coords(i) <= 1.));
    sum += coords(i);
  }
  if(in) return (in && (sum <= 1));
  return in;
}

/* -------------------------------------------------------------------------- */
/* InterpolationElement                                                       */
/* -------------------------------------------------------------------------- */
template<InterpolationType interpolation_type, InterpolationKind kind>
inline void
InterpolationElement<interpolation_type, kind>::computeShapes(const types::Matrix<Real> & natural_coord,
							      types::Matrix<Real> & N) {
  UInt nb_points = natural_coord.cols();
  for (UInt p = 0; p < nb_points; ++p) {
    types::Vector<Real> Np = N(p);
    computeShapes(natural_coord(p), Np);
  }
}

/* -------------------------------------------------------------------------- */
template<InterpolationType interpolation_type, InterpolationKind kind>
inline void
InterpolationElement<interpolation_type, kind>::computeDNDS(const types::Matrix<Real> & natural_coord,
							    types::Tensor3<Real> & dnds) {
  UInt nb_points = natural_coord.cols();
  for (UInt p = 0; p < nb_points; ++p) {
    types::Matrix<Real> dnds_p = dnds(p);
    computeDNDS(natural_coord(p), dnds_p);
  }
}

/* -------------------------------------------------------------------------- */
/**
 * interpolate the field on a point given in natural coordinates the field which
 * values are given on the node of the element
 *
 * @param natural_coords natural coordinates of point where to interpolate \xi
 * @param nodal_values values of the function per node @f$ f_{ij} = f_{n_i j} @f$ so it should be a matrix of size nb_nodes_per_element @f$\times@f$ nb_degree_of_freedom
 * @param interpolated interpolated value of f @f$ f_j(\xi) = \sum_i f_{n_i j} N_i @f$
 */
template<InterpolationType interpolation_type, InterpolationKind kind>
inline void
InterpolationElement<interpolation_type, kind>::interpolateOnNaturalCoordinates(const types::Vector<Real> & natural_coords,
										const types::Matrix<Real> & nodal_values,
										types::Vector<Real> & interpolated) {
  types::Vector<Real> shapes(nb_nodes_per_element);
  computeShapes(natural_coords, shapes);

  types::Matrix<Real> interpm(interpolated.storage(), nodal_values.rows(), 1);
  types::Matrix<Real> shapesm(shapes.storage(), nb_nodes_per_element, 1);
  interpm.mul<false, false>(nodal_values, shapesm);
}


/* -------------------------------------------------------------------------- */
/// @f$ gradient_{ij} = \frac{\partial f_j}{\partial s_i} = \sum_k \frac{\partial N_k}{\partial s_i}f_{j n_k} @f$
template<InterpolationType interpolation_type, InterpolationKind kind>
inline void
InterpolationElement<interpolation_type, kind>::gradientOnNaturalCoordinates(const types::Vector<Real> & natural_coords,
									     const types::Matrix<Real> & f,
									     types::Matrix<Real> & gradient) {
  types::Matrix<Real> dnds(natural_space_dimension, nb_nodes_per_element);
  computeDNDS(natural_coords, dnds);
  gradient.mul<false, true>(f, dnds);
}

/* -------------------------------------------------------------------------- */
/* ElementClass                                                               */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
inline void ElementClass<type, kind>::computeJMat(const types::Tensor3<Real> & dnds,
						  const types::Matrix<Real> & node_coords,
						  types::Tensor3<Real> & J) {
  UInt nb_points = dnds.size(2);
  for (UInt p = 0; p < nb_points; ++p) {
    types::Matrix<Real> J_p = J(p);
    computeJMat(dnds(p), node_coords, J_p);
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
inline void ElementClass<type, kind>::computeJMat(const types::Matrix<Real> & dnds,
						  const types::Matrix<Real> & node_coords,
						  types::Matrix<Real> & J) {
  /// @f$ J = dxds = dnds * x @f$
  J.mul<false, true>(dnds, node_coords);
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
inline void ElementClass<type, kind>::computeJacobian(const types::Matrix<Real> & natural_coords,
						      const types::Matrix<Real> & node_coords,
						      types::Vector<Real> & jacobians) {
  UInt nb_points = natural_coords.cols();
  types::Matrix<Real> dnds(interpolation_element::natural_space_dimension,
			   interpolation_element::nb_nodes_per_element);
  types::Matrix<Real> J(natural_coords.rows(),
			node_coords.rows());

  for (UInt p = 0; p < nb_points; ++p) {
    interpolation_element::computeDNDS(natural_coords(p), dnds);
    computeJMat(dnds, node_coords, J);
    computeJacobian(J, jacobians(p));
  }
}



/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
inline void ElementClass<type, kind>::computeJacobian(const types::Tensor3<Real> & J,
						      types::Vector<Real> & jacobians) {
  UInt nb_points = J.size(2);
  for (UInt p = 0; p < nb_points; ++p) {
    computeJacobian(J(p), jacobians(p));
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
inline void ElementClass<type, kind>::computeJacobian(const types::Matrix<Real> & J,
						      Real & jacobians) {
  if(J.rows() == J.cols()) {
    jacobians = Math::det<element_property::spatial_dimension>(J.storage());
  } else {
    interpolation_element::computeSpecialJacobian(J, jacobians);
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
inline void
ElementClass<type, kind>::computeShapeDerivatives(const types::Tensor3<Real> & J,
						  const types::Tensor3<Real> & dnds,
						  types::Tensor3<Real> & shape_deriv) {
  UInt nb_points = J.size(2);
  for (UInt p = 0; p < nb_points; ++p) {
    types::Matrix<Real> shape_deriv_p = shape_deriv(p);
    computeShapeDerivatives(J(p), dnds(p), shape_deriv_p);
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
inline void
ElementClass<type, kind>::computeShapeDerivatives(const types::Matrix<Real> & J,
						  const types::Matrix<Real> & dnds,
						  types::Matrix<Real> & shape_deriv) {
  types::Matrix<Real> inv_J(J.rows(), J.cols());
  Math::inv<element_property::spatial_dimension>(J.storage(), inv_J.storage());

  shape_deriv.mul<false, false>(inv_J, dnds);
}


/* -------------------------------------------------------------------------- */
template<ElementType type, ElementKind kind>
inline void
ElementClass<type, kind>::computeNormalsOnNaturalCoordinates(const types::Matrix<Real> & coord,
							     types::Matrix<Real> & f,
							     types::Matrix<Real> & normals) {
  UInt dimension = normals.rows();
  UInt nb_points = coord.cols();

  AKANTU_DEBUG_ASSERT((dimension - 1) == interpolation_element::natural_space_dimension,
		      "cannot extract a normal because of dimension mismatch "
		      << dimension << " " << interpolation_element::natural_space_dimension);

  types::Matrix<Real> J(dimension, interpolation_element::natural_space_dimension);
  for (UInt p = 0; p < nb_points; ++p) {
    interpolation_element::gradientOnNaturalCoordinates(coord(p), f, J);
    if (dimension == 2) {
      Math::normal2(J.storage(), normals(p).storage());
    }
    if (dimension == 3){
      Math::normal3(J(0).storage(), J(1).storage(), normals(p).storage());
    }
  }
}

/* ------------------------------------------------------------------------- */
/**
 * In the non linear cases we need to iterate to find the natural coordinates @f$\xi@f$
 * provided real coordinates @f$x@f$.
 *
 * We want to solve: @f$ x- \phi(\xi) = 0@f$ with @f$\phi(\xi) = \sum_I N_I(\xi) x_I@f$
 * the mapping function which uses the nodal coordinates @f$x_I@f$.
 *
 * To that end we use the Newton method and the following series:
 *
 * @f$ \frac{\partial \phi(x_k)}{\partial \xi} \left( \xi_{k+1} - \xi_k \right) = x - \phi(x_k)@f$
 *
 * When we consider elements embedded in a dimension higher than them (2D triangle in a 3D space for example)
 * @f$ J = \frac{\partial \phi(\xi_k)}{\partial \xi}@f$ is of dimension @f$dim_{space} \times dim_{elem}@f$ which
 * is not invertible in most cases. Rather we can solve the problem:
 *
 * @f$ J^T J \left( \xi_{k+1} - \xi_k \right) = J^T \left( x - \phi(\xi_k) \right) @f$
 *
 * So that
 *
 * @f$ d\xi = \xi_{k+1} - \xi_k = (J^T J)^{-1} J^T \left( x - \phi(\xi_k) \right) @f$
 *
 * So that if the series converges we have:
 *
 * @f$ 0 = J^T \left( \phi(\xi_\infty) - x \right) @f$
 *
 * And we see that this is ill-posed only if @f$ J^T x = 0@f$ which means that the vector provided
 * is normal to any tangent which means it is outside of the element itself.
 *
 * @param real_coords: the real coordinates the natural coordinates are sought for
 * @param node_coords: the coordinates of the nodes forming the element
 * @param natural_coords: output->the sought natural coordinates
 * @param spatial_dimension: spatial dimension of the problem
 *
 **/
template <ElementType type, ElementKind kind>
inline void ElementClass<type, kind>::inverseMap(const types::Vector<Real> & real_coords,
						 const types::Matrix<Real> & node_coords,
						 types::RVector & natural_coords,
						 Real tolerance) {
  UInt spatial_dimension = real_coords.size();
  UInt dimension = natural_coords.size();

  //matrix copy of the real_coords
  types::RMatrix mreal_coords(real_coords.storage(), spatial_dimension, 1);

  //initial guess
  //  types::RMatrix natural_guess(natural_coords.storage(), dimension, 1);
  natural_coords.clear();

  // real space coordinates provided by initial guess
  types::RMatrix physical_guess(dimension, 1);

  // objective function f = real_coords - physical_guess
  types::RMatrix f(dimension, 1);

  // dnds computed on the natural_guess
  //  types::RMatrix dnds(interpolation_element::nb_nodes_per_element, spatial_dimension);

  // J Jacobian matrix computed on the natural_guess
  types::RMatrix J(spatial_dimension, dimension);

  // G = J^t * J
  types::RMatrix G(spatial_dimension, spatial_dimension);

  // Ginv = G^{-1}
  types::RMatrix Ginv(spatial_dimension, spatial_dimension);

  // J = Ginv * J^t
  types::RMatrix F(spatial_dimension, dimension);

  // dxi = \xi_{k+1} - \xi in the iterative process
  types::RMatrix dxi(spatial_dimension, 1);

  /* --------------------------- */
  /* init before iteration loop  */
  /* --------------------------- */
  // do interpolation
  types::Vector<Real> physical_guess_v(physical_guess.storage(), dimension);
  interpolation_element::interpolateOnNaturalCoordinates(natural_coords,
							 node_coords,
							 physical_guess_v);

  // compute initial objective function value f = real_coords - physical_guess
  f  = mreal_coords;
  f -= physical_guess;

  // compute initial error
  Real inverse_map_error = f.norm();

  /* --------------------------- */
  /* iteration loop              */
  /* --------------------------- */
  while(tolerance < inverse_map_error) {
    //compute J^t
    interpolation_element::gradientOnNaturalCoordinates(natural_coords, node_coords, J);

    //compute G
    G.mul<true, false>(J, J);

    // inverse G
    Ginv.inverse(G);

    //compute F
    F.mul<false, true>(Ginv, J);

    //compute increment
    dxi.mul<false,false>(F, f);

    //update our guess
    natural_coords += dxi(0);

    //interpolate
    interpolation_element::interpolateOnNaturalCoordinates(natural_coords,
							   node_coords,
							   physical_guess_v);

    // compute error
    f  = mreal_coords;
    f -= physical_guess;
    inverse_map_error = f.norm();
  }
  //  memcpy(natural_coords.storage(), natural_guess.storage(), sizeof(Real) * natural_coords.size());
}
