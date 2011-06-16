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
template<ElementType type> inline UInt ElementClass<type>::getNbQuadraturePoint() {
  return nb_quadrature_points;
}

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
    // 			"Negative jacobian computed, possible problem in the element node order.");
    cdxds += spatial_dimension * dimension;
    cjac++;
  }
}
/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ElementClass<type>::computeShapeDerivatives(const Real * dxds,
							const Real * dnds,
							const UInt nb_points,
							const UInt dimension,
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
inline void ElementClass<type>::computeShapes(__attribute__ ((unused)) const Real * natural_coords,
					      __attribute__ ((unused)) Real * shapes) {
  AKANTU_DEBUG_ERROR("Function not implemented for type : " << type);
}
/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ElementClass<type>::computeDNDS(__attribute__ ((unused)) const Real * natural_coords,
					    __attribute__ ((unused)) Real * dnds) {
  AKANTU_DEBUG_ERROR("Function not implemented for type : " << type);
}


/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ElementClass<type>::computeJacobian(__attribute__ ((unused)) const Real * dxds,
						__attribute__ ((unused)) const UInt dimension,
						__attribute__ ((unused)) Real & jac) {
  AKANTU_DEBUG_ERROR("Function not implemented for type : " << type);
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
inline Real * ElementClass<type>::getGaussIntegrationWeights() {
  return gauss_integration_weights;
}


#include "element_classes/element_class_segment_2.cc"
#include "element_classes/element_class_segment_3.cc"
#include "element_classes/element_class_triangle_3.cc"
#include "element_classes/element_class_triangle_6.cc"
#include "element_classes/element_class_tetrahedron_4.cc"
#include "element_classes/element_class_tetrahedron_10.cc"
#include "element_classes/element_class_quadrangle_4.cc"
#include "element_classes/element_class_hexahedron_8.cc"
