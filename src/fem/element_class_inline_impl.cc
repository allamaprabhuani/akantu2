/**
 * @file   element_class_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @date   Thu Jul 15 10:28:28 2010
 *
 * @brief  Implementation of the inline functions of the class element_class
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
template<ElementType type> inline Real * ElementClass<type>::getQuadraturePoints() {
  return quad;
}

/* -------------------------------------------------------------------------- */
template<ElementType type> inline UInt ElementClass<type>::getShapeSize() {
  return nb_nodes_per_element;
}

/* -------------------------------------------------------------------------- */
template<ElementType type> inline UInt ElementClass<type>::getShapeDerivativesSize() {
  return nb_nodes_per_element * spatial_dimension;
}

/* -------------------------------------------------------------------------- */
template<ElementType type>
void ElementClass<type>::preComputeStandards(const Real * coord,
                                             const UInt dimension,
                                             Real * shape,
                                             Real * dshape,
                                             Real * jacobians) {
  // ask for computation of shapes
  computeShapes(quad, nb_quadrature_points, shape);

  // compute dnds
  Real dnds[nb_nodes_per_element * spatial_dimension * nb_quadrature_points];
  computeDNDS(quad, nb_quadrature_points, dnds);

  // compute dxds
  Real dxds[dimension * spatial_dimension * nb_quadrature_points];
  computeDXDS(dnds, nb_quadrature_points, coord, dimension, dxds);

  // jacobian
  computeJacobian(dxds, nb_quadrature_points, dimension, jacobians);
  // if dimension == spatial_dimension compute shape derivatives
  if (dimension == spatial_dimension) {
    computeShapeDerivatives(dxds, dnds, nb_quadrature_points, dimension, dshape);
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ElementClass<type>::computeShapes(const Real * natural_coords,
                                              const UInt nb_points,
                                              Real * shapes) {
  Real * cpoint = const_cast<Real *>(natural_coords);
  for (UInt p = 0; p < nb_points; ++p) {
    computeShapes(cpoint, shapes);
    shapes += nb_nodes_per_element;
    cpoint += spatial_dimension;
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ElementClass<type>::computeShapes(const Real * natural_coords,
                                              const UInt nb_points,
                                              Real * shapes,
                                              const Real * local_coord,
                                              UInt id) {
  Real * cpoint = const_cast<Real *>(natural_coords);
  for (UInt p = 0; p < nb_points; ++p) {
    computeShapes(cpoint, shapes, local_coord, id);
    shapes += nb_nodes_per_element;
    cpoint += spatial_dimension;
  }
}

/* -------------------------------------------------------------------------- */
template<ElementType type>
inline void ElementClass<type>::computeDNDS(const Real * natural_coords,
                                            const UInt nb_points,
                                            Real * dnds) {
  Real * cpoint = const_cast<Real *>(natural_coords);
  Real * cdnds = dnds;
  for (UInt p = 0; p < nb_points; ++p) {
    computeDNDS(cpoint, cdnds);
    cpoint += spatial_dimension;
    cdnds += nb_nodes_per_element * spatial_dimension;
  }
}

/* -------------------------------------------------------------------------- */
template<ElementType type>
inline void ElementClass<type>::computeDXDS(const Real * dnds,
                                            const UInt nb_points,
                                            const Real * node_coords,
                                            const UInt dimension,
                                            Real * dxds) {
  Real * cdnds = const_cast<Real *>(dnds);
  Real * cdxds = dxds;
  for (UInt p = 0; p < nb_points; ++p) {
    computeDXDS(cdnds, node_coords, dimension, cdxds);
    cdnds += nb_nodes_per_element * spatial_dimension;
    cdxds += spatial_dimension * dimension;
  }
}
/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ElementClass<type>::computeDXDS(const Real * dnds,
                                            const Real * node_coords,
                                            const UInt dimension,
                                            Real * dxds) {
  /// @f$ J = dxds = dnds * x @f$
  Math::matrix_matrix(spatial_dimension, dimension, nb_nodes_per_element,
                      dnds, node_coords, dxds);
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ElementClass<type>::computeJacobian(const Real * dxds,
                                                const UInt nb_points,
                                                const UInt dimension,
                                                Real * jac) {
  Real * cdxds = const_cast<Real *>(dxds);
  Real * cjac = jac;
  for (UInt p = 0; p < nb_points; ++p) {
    computeJacobian(cdxds, dimension, *cjac);
    // AKANTU_DEBUG_ASSERT((cjac[0] > 0),
    //                  "Negative jacobian computed, possible problem in the element node order.");
    cdxds += spatial_dimension * dimension;
    cjac++;
  }
}
/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ElementClass<type>::computeShapeDerivatives(const Real * dxds,
                                                        const Real * dnds,
                                                        const UInt nb_points,
                                                        __attribute__ ((unused)) const UInt dimension,
                                                        Real * shape_deriv) {
  AKANTU_DEBUG_ASSERT(dimension == spatial_dimension,"gradient in space "
                      << dimension
                      << " cannot be evaluated for element of dimension "
                      << spatial_dimension);

  Real * cdxds = const_cast<Real *>(dxds);
  Real * cdnds = const_cast<Real *>(dnds);
  for (UInt p = 0; p < nb_points; ++p) {
    computeShapeDerivatives(cdxds, cdnds, shape_deriv);
    cdnds += spatial_dimension * nb_nodes_per_element;
    cdxds += spatial_dimension * spatial_dimension;
    shape_deriv += nb_nodes_per_element * spatial_dimension;
  }
}


/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ElementClass<type>::computeShapeDerivatives(const Real * natural_coords,
                                                        const UInt  nb_points,
                                                        const UInt dimension,
                                                        Real * shape_deriv,
                                                        const Real * local_coord,
                                                        UInt id) {
  AKANTU_DEBUG_ASSERT(dimension == spatial_dimension,"Gradient in space "
                      << dimension
                      << " cannot be evaluated for element of dimension "
                      << spatial_dimension);

  Real * cpoint = const_cast<Real *>(natural_coords);
  for (UInt p = 0; p < nb_points; ++p) {
    computeShapeDerivatives(cpoint, shape_deriv, local_coord, id);
    shape_deriv += nb_nodes_per_element * dimension;
    cpoint += dimension;
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ElementClass<type>::computeShapeDerivatives(const Real * dxds,
                                                        const Real * dnds,
                                                        Real * shape_deriv) {
  /// @f$ dxds = J^{-1} @f$
  Real inv_dxds[spatial_dimension * spatial_dimension];
  if (spatial_dimension == 1) inv_dxds[0] = 1./dxds[0];
  if (spatial_dimension == 2) Math::inv2(dxds, inv_dxds);
  if (spatial_dimension == 3) Math::inv3(dxds, inv_dxds);

  Math::matrixt_matrixt(nb_nodes_per_element, spatial_dimension, spatial_dimension,
                        dnds, inv_dxds, shape_deriv);
}
/* -------------------------------------------------------------------------- */
template<ElementType type>
inline Real ElementClass<type>::getInradius(__attribute__ ((unused)) const Real * coord) {
  AKANTU_DEBUG_ERROR("Function not implemented for type : " << type);
  return 0;
}
/* -------------------------------------------------------------------------- */
template<ElementType type>
inline void ElementClass<type>::computeNormalsOnQuadPoint(const Real * coord,
                                                          const UInt dimension,
                                                          Real * normals) {
  AKANTU_DEBUG_ASSERT((dimension - 1) == spatial_dimension,
                      "cannot extract a normal because of dimension mismatch "
                      << dimension << " " << spatial_dimension);

  Real * cpoint = const_cast<Real *>(quad);
  Real * cnormals = normals;
  Real dnds[spatial_dimension*nb_nodes_per_element];
  Real dxds[spatial_dimension*dimension];

  for (UInt p = 0; p < nb_quadrature_points; ++p) {
    computeDNDS(cpoint,dnds);
    computeDXDS(dnds,coord,dimension,dxds);
    if (dimension == 2) {
      Math::normal2(dxds,cnormals);
    }
    if (dimension == 3){
      Math::normal3(dxds,dxds+dimension,cnormals);
    }
    cpoint += spatial_dimension;
    cnormals += dimension;
  }
}
/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ElementClass<type>::interpolateOnNaturalCoordinates(const Real * natural_coords,
                                                                const Real * nodal_values,
                                                                UInt dimension,
                                                                Real * interpolated){


  Real shapes[nb_nodes_per_element];
  computeShapes(natural_coords,shapes);
  Math::matrix_matrix(1, dimension, nb_nodes_per_element,
                      shapes, nodal_values, interpolated);
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ElementClass<type>::inverseMap(const types::RVector & real_coords,
                                           const types::RMatrix & node_coords,
                                           UInt dimension,
                                           types::RVector & natural_coords,
                                           Real tolerance){

  //matric copy of the real_coords
  types::RMatrix mreal_coords(real_coords.storage(),1,spatial_dimension);
  //initial guess
  types::RMatrix natural_guess(1,dimension,0.);
  // realspace coordinates provided by initial guess
  types::RMatrix physical_guess(1,dimension);
  // objective function f = real_coords - physical_guess
  types::RMatrix f(1,dimension);
  // dnds computed on the natural_guess
  types::RMatrix dnds(nb_nodes_per_element,spatial_dimension);
  // dxds computed on the natural_guess
  types::RMatrix dxds(spatial_dimension,dimension);
  // transposed dxds computed on the natural_guess
  types::RMatrix dxds_t(dimension,spatial_dimension);
  // G = dxds * dxds_t
  types::RMatrix G(spatial_dimension,spatial_dimension);
  // Ginv = G^{-1}
  types::RMatrix Ginv(spatial_dimension,spatial_dimension);
  // J = Ginv * dxds
  types::RMatrix J(spatial_dimension,dimension);
  // dxi = \xi_{k+1} - \xi in the iterative process
  types::RMatrix dxi(1,spatial_dimension);

  /* --------------------------- */
  // init before iteration loop
  /* --------------------------- */

  // do interpolation
  interpolateOnNaturalCoordinates(natural_guess.storage(),node_coords.storage(),dimension,physical_guess.storage());
  // compute initial objective function value f = real_coords - physical_guess
  f = physical_guess;
  f*= -1.;
  f+= mreal_coords;
  // compute initial error
  Real inverse_map_error = f.norm();


  /* --------------------------- */
  // iteration loop
  /* --------------------------- */

  while(tolerance < inverse_map_error)
    {
      //compute dxds
      computeDNDS(natural_guess.storage(), dnds.storage());
      computeDXDS(dnds.storage(),node_coords.storage(),dimension,dxds.storage());
      //compute G
      dxds_t = dxds;
      dxds_t.transpose();
      G.mul<false,false>(dxds,dxds_t);
      // inverse G
      if      (spatial_dimension == 1) Ginv[0] = 1./G[0];
      else if (spatial_dimension == 2) Math::inv2(G.storage(),Ginv.storage());
      else if (spatial_dimension == 3) Math::inv3(G.storage(),Ginv.storage());

      //compute J
      J.mul<false,false>(Ginv,dxds);

      //compute increment
      dxi.mul<false,false>(f,J);


      //update our guess
      natural_guess += dxi;
      //interpolate
      interpolateOnNaturalCoordinates(natural_guess.storage(),node_coords.storage(),dimension,physical_guess.storage());
      // compute error
      f = physical_guess;
      f*= -1.;
      f+= mreal_coords;
      inverse_map_error = f.norm();
    }
    memcpy(natural_coords.storage(),natural_guess.storage(),sizeof(Real)*natural_coords.size());
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ElementClass<type>::computeShapes(const Real * natural_coords,
                                              Real * shapes) {
  computeShapes(natural_coords, shapes, NULL,0);
}
/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ElementClass<type>::computeShapes(__attribute__ ((unused)) const Real * natural_coords,
                                              __attribute__ ((unused)) Real * shapes,
                                              __attribute__ ((unused)) const Real * local_coord,
                                              __attribute__ ((unused)) UInt id) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}
/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ElementClass<type>::computeDNDS(__attribute__((unused)) const Real * natural_coords,
                                            __attribute__((unused)) Real * dnds){
  //  computeDNDS(natural_coords, dnds);
  AKANTU_DEBUG_TO_IMPLEMENT();
}
/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ElementClass<type>::computeShapeDerivatives(__attribute__ ((unused)) const Real * natural_coords,
                                                        __attribute__ ((unused)) Real * shape_deriv,
                                                        __attribute__ ((unused)) const Real * local_coord,
                                                        __attribute__ ((unused)) UInt id) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}


/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ElementClass<type>::computeJacobian(__attribute__ ((unused)) const Real * dxds,
                                                __attribute__ ((unused)) const UInt dimension,
                                                __attribute__ ((unused)) Real & jac) {
  //AKANTU_DEBUG_ERROR("Function not implemented for type : " << type);
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
inline Real * ElementClass<type>::getGaussIntegrationWeights() {
  return gauss_integration_weights;
}

/* -------------------------------------------------------------------------- */


template <ElementType type>
inline bool ElementClass<type>::contains(__attribute__ ((unused)) const types::RVector & natural_coords) {
  AKANTU_DEBUG_TO_IMPLEMENT();
  return false;
}

/* -------------------------------------------------------------------------- */
#include "element_classes/element_class_segment_2_inline_impl.cc"
#include "element_classes/element_class_segment_3_inline_impl.cc"
#include "element_classes/element_class_triangle_3_inline_impl.cc"
#include "element_classes/element_class_triangle_6_inline_impl.cc"
#include "element_classes/element_class_tetrahedron_4_inline_impl.cc"
#include "element_classes/element_class_tetrahedron_10_inline_impl.cc"
#include "element_classes/element_class_quadrangle_4_inline_impl.cc"
#include "element_classes/element_class_quadrangle_8_inline_impl.cc"
#include "element_classes/element_class_hexahedron_8_inline_impl.cc"
#include "element_classes/element_class_bernoulli_beam_2_inline_impl.cc"
