/**
 * @file   integrator_gauss_inline_impl.cc
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Feb 10 20:43:52 2011
 *
 * @brief  inline function of gauss integrator
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
template <ElementType type>
inline void IntegratorGauss::integrateOnElement(const Vector<Real> & f,
						Real * intf,
						UInt nb_degree_of_freedom,
						const UInt elem,
						const GhostType & ghost_type) const {

  Vector<Real> * jac_loc = jacobians(type, ghost_type);

  UInt nb_quadrature_points = ElementClass<type>::getNbQuadraturePoints();
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == nb_degree_of_freedom ,
		      "The vector f do not have the good number of component.");

  Real * f_val    = f.values + elem * f.getNbComponent();
  Real * jac_val  = jac_loc->values + elem * nb_quadrature_points;

  integrate(f_val, jac_val, intf, nb_degree_of_freedom, nb_quadrature_points);
}


/* -------------------------------------------------------------------------- */
template <ElementType type>
inline Real IntegratorGauss::integrate(const types::RVector & in_f,
				       UInt index,
				       const GhostType & ghost_type) const{

  
  const Vector<Real> & jac_loc = jacobians(type, ghost_type);

  UInt nb_quadrature_points = ElementClass<type>::getNbQuadraturePoints();
  AKANTU_DEBUG_ASSERT(in_f.size() == nb_quadrature_points ,
		      "The vector f do not have nb_quadrature_points entries.");

  Real * jac_val  = jac_loc.storage() + index * nb_quadrature_points;
  Real intf;

  integrate(in_f.storage(), jac_val, &intf, 1, nb_quadrature_points);

  return intf;
}


/* -------------------------------------------------------------------------- */

inline void IntegratorGauss::integrate(Real *f, Real *jac, Real * inte,
			   UInt nb_degree_of_freedom,
			   UInt nb_quadrature_points) const {
  memset(inte, 0, nb_degree_of_freedom * sizeof(Real));

  Real *cjac = jac;
  for (UInt q = 0; q < nb_quadrature_points; ++q) {
    for (UInt dof = 0; dof < nb_degree_of_freedom; ++dof) {
      inte[dof] += *f * *cjac;
      ++f;
    }
    ++cjac;
  }
}




/* -------------------------------------------------------------------------- */
template <ElementType type>
inline const Vector<Real> & IntegratorGauss::getQuadraturePoints(const GhostType & ghost_type) const {
  AKANTU_DEBUG_ASSERT(quadrature_points.exists(type, ghost_type),
		      "Quadrature points for type "
		      << quadrature_points.printType(type, ghost_type)
 		      << " have not been initialized."
		      << " Did you use 'computeQuadraturePoints' function ?");
  return quadrature_points(type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void IntegratorGauss::computeQuadraturePoints(const GhostType & ghost_type) {
  UInt n_coords = ElementClass<type>::getNbQuadraturePoints();
  UInt dim = ElementClass<type>::getSpatialDimension();

  quadrature_points.alloc(0, dim, type, ghost_type);
  Vector<Real> & quad_points = quadrature_points(type, ghost_type);

  Real * coord_val = ElementClass<type>::getQuadraturePoints();
  for (UInt i = 0; i < n_coords; ++i) {
    quad_points.push_back(coord_val);
    coord_val += dim;
  }
}
/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void IntegratorGauss::
computeJacobianOnQuadPointsByElement(UInt spatial_dimension,
				     Real * node_coords,
				     UInt nb_nodes_per_element,
				     Real * jacobians) {

  Real * quad = ElementClass<type>::getQuadraturePoints();
  const UInt nb_quad_points = ElementClass<type>::getNbQuadraturePoints();
  // compute dnds
  Real * dnds = new Real[nb_nodes_per_element * spatial_dimension * nb_quad_points];
  ElementClass<type>::computeDNDS(quad,
				  nb_quad_points,
				  dnds);
  // compute dxds
  const UInt element_dimension = ElementClass<type>::getSpatialDimension();
  Real * dxds = new Real[element_dimension * spatial_dimension * nb_quad_points];
  ElementClass<type>::computeDXDS(dnds, nb_quad_points,
				  node_coords, spatial_dimension, dxds);

  delete [] dnds;

  // jacobian
  ElementClass<type>::computeJacobian(dxds, nb_quad_points,
				      spatial_dimension, jacobians);

  delete [] dxds;
}
/* -------------------------------------------------------------------------- */
