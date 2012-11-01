/**
 * @file   integrator_cohesive_inline_impl.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Thu Feb 23 17:58:45 2012
 *
 * @brief  IntegratorCohesive inline implementation
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
template <class Inte>
IntegratorCohesive<Inte>::IntegratorCohesive(const Mesh & mesh,
					     const ID & id,
					     const MemoryID & memory_id) :
  Inte(mesh, id, memory_id) {
  AKANTU_DEBUG_IN();

  std::stringstream sstr;
  sstr << id << "sub_integrator";
  //sub_type_integrator = new Inte(mesh, sstr.str(), memory_id);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<class Inte>
template <ElementType type>
void IntegratorCohesive<Inte>::precomputeJacobiansOnQuadraturePoints(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();
  const ElementType sub_type = ElementType(CohesiveElementSubElementType<type>::value);
  //  sub_type_integrator->precomputeJacobiansOnQuadraturePoints<sub_type>(ghost_type);

  /* ---------------------------------*/

  UInt spatial_dimension = Inte::mesh->getSpatialDimension();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points = CohesiveElement<type>::getNbQuadraturePoints();
  Real * weights = CohesiveElement<type>::getGaussIntegrationWeights();

  UInt * elem_val = Inte::mesh->getConnectivity(type,ghost_type).storage();;
  UInt nb_element = Inte::mesh->getConnectivity(type,ghost_type).getSize();

  Vector<Real> & jacobians_tmp = Inte::jacobians.alloc(nb_element*nb_quadrature_points,
								      1,
								      type,
								      ghost_type);

  Real * jacobians_val = jacobians_tmp.storage();

  Real local_coord[spatial_dimension * nb_nodes_per_element/2];
  for (UInt elem = 0; elem < nb_element; ++elem) {

    // extract the coordinates of the first line nodes first
    UInt * connectivity = elem_val+elem*nb_nodes_per_element;
    for (UInt n = 0; n < nb_nodes_per_element/2; ++n) {
      memcpy(local_coord + n * spatial_dimension,
  	     Inte::mesh->getNodes().storage() + connectivity[n] * spatial_dimension,
  	     spatial_dimension * sizeof(Real));
    }

    static_cast<Inte*>(this)->computeJacobianOnQuadPointsByElement<sub_type>(spatial_dimension,
									     local_coord,
									     nb_nodes_per_element/2,
									     jacobians_val);

    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      *jacobians_val++ *= weights[q];
    }
    //    jacobians_val += nb_quadrature_points;
  }


  /* ---------------------------------*/

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<class Inte>
template <ElementType type>
inline void IntegratorCohesive<Inte>::integrateOnElement(const Vector<Real> & f,
							 Real * intf,
							 UInt nb_degree_of_freedom,
							 const UInt elem,
							 const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();
  const ElementType sub_type = ElementType(CohesiveElementSubElementType<type>::value);
  static_cast<Inte*>(this)->integrateOnElement<sub_type>(f, intf, nb_degree_of_freedom, elem, ghost_type);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<class Inte>
template <ElementType type>
const Vector<Real> &
IntegratorCohesive<Inte>::getQuadraturePoints(const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();
  const ElementType sub_type = ElementType(CohesiveElementSubElementType<type>::value);
  AKANTU_DEBUG_OUT();
  return Inte::template getQuadraturePoints<sub_type>(ghost_type);
}

/* -------------------------------------------------------------------------- */
template<class Inte>
template <ElementType type>
void IntegratorCohesive<Inte>::computeQuadraturePoints(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();
  const ElementType sub_type = ElementType(CohesiveElementSubElementType<type>::value);
  static_cast<Inte*>(this)->template computeQuadraturePoints<sub_type>(ghost_type);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<class Inte>
template <ElementType type>
inline void IntegratorCohesive<Inte>::
computeJacobianOnQuadPointsByElement(UInt spatial_dimension,
				     Real * node_coords,
				     UInt nb_nodes_per_element,
				     Real * jacobians) {
  AKANTU_DEBUG_IN();
  const ElementType sub_type = ElementType(CohesiveElementSubElementType<type>::value);
  static_cast<Inte*>(this)->computeJacobianOnQuadPointsByElement<sub_type>(spatial_dimension,
									   node_coords,
									   nb_nodes_per_element,
									   jacobians);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<class Inte>
template <ElementType type>
void IntegratorCohesive<Inte>::integrate(const Vector<Real> & in_f,
					 Vector<Real> &intf,
					 UInt nb_degree_of_freedom,
					 const GhostType & ghost_type,
					 const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();
  //  const ElementType sub_type = ElementType(CohesiveElementSubElementType<type>::value);
  //  Inte::template integrate<sub_type>(in_f, intf, nb_degree_of_freedom, ghost_type, filter_elements);

  /* ------------------------------------------------------------------------ */

  AKANTU_DEBUG_ASSERT(Inte::jacobians.exists(type, ghost_type),
		      "No jacobians for the type "
		      << Inte::jacobians.printType(type, ghost_type));

  UInt nb_element = Inte::mesh->getNbElement(type,ghost_type);
  const Vector<Real> & jac_loc = Inte::jacobians(type, ghost_type);

  UInt nb_quadrature_points = CohesiveElement<type>::getNbQuadraturePoints();

  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    nb_element      = filter_elements->getSize();
    filter_elem_val = filter_elements->values;
  }

  Real * in_f_val = in_f.storage();
  Real * intf_val = intf.storage();
  Real * jac_val  = jac_loc.storage();

  UInt offset_in_f = in_f.getNbComponent()*nb_quadrature_points;
  UInt offset_intf = intf.getNbComponent();

  Real * jac      = jac_val;

  for (UInt el = 0; el < nb_element; ++el) {
    if(filter_elements != NULL) {
      jac      = jac_val  + filter_elem_val[el] * nb_quadrature_points;
    }

    Inte::integrate(in_f_val, jac, intf_val, nb_degree_of_freedom, nb_quadrature_points);

    in_f_val += offset_in_f;
    intf_val += offset_intf;
    if(filter_elements == NULL) {
      jac      += nb_quadrature_points;
    }
  }

  /* ------------------------------------------------------------------------ */

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template<class Inte>
template <ElementType type>
Real IntegratorCohesive<Inte>::integrate(const Vector<Real> & in_f,
					 const GhostType & ghost_type,
					 const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(Inte::jacobians.exists(type, ghost_type),
		      "No jacobians for the type "
		      << Inte::jacobians.printType(type, ghost_type));

  UInt nb_element = Inte::mesh->getNbElement(type, ghost_type);
  const Vector<Real> & jac_loc = Inte::jacobians(type, ghost_type);

  UInt nb_quadrature_points = CohesiveElement<type>::getNbQuadraturePoints();

  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    nb_element      = filter_elements->getSize();
    filter_elem_val = filter_elements->values;
  }

  Real intf = 0.;
  Real * in_f_val  = in_f.storage();
  Real * jac_val   = jac_loc.storage();
  UInt offset_in_f = in_f.getNbComponent() * nb_quadrature_points;
  Real * jac       = jac_val;

  for (UInt el = 0; el < nb_element; ++el) {
    if(filter_elements != NULL) {
      jac = jac_val  + filter_elem_val[el] * nb_quadrature_points;
    }

    Real el_intf = 0;
    Inte::integrate(in_f_val, jac, &el_intf, 1, nb_quadrature_points);
    intf += el_intf;

    in_f_val += offset_in_f;
    if(filter_elements == NULL) {
      jac += nb_quadrature_points;
    }
  }

  AKANTU_DEBUG_OUT();
  return intf;
}


/* -------------------------------------------------------------------------- */
// template<class Inte>
// template <ElementType type>
// Real IntegratorCohesive<Inte>::integrate(const Vector<Real> & in_f,
// 					 const GhostType & ghost_type,
// 					 const Vector<UInt> * filter_elements) const {
//   AKANTU_DEBUG_IN();
//   const ElementType sub_type = ElementType(CohesiveElementSubElementType<type>::value);
//   AKANTU_DEBUG_OUT();
//   return  Inte::template integrate<sub_type>(in_f, ghost_type, filter_elements);
// }


/* -------------------------------------------------------------------------- */
template<class Inte>
template <ElementType type>
void IntegratorCohesive<Inte>::integrateOnQuadraturePoints(const Vector<Real> & in_f,
							   Vector<Real> &intf,
							   UInt nb_degree_of_freedom,
							   const GhostType & ghost_type,
							   const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();
  const ElementType sub_type = ElementType(CohesiveElementSubElementType<type>::value);
  Inte::template integrateOnQuadraturePoints<sub_type>(in_f,
						       intf,
						       nb_degree_of_freedom,
						       ghost_type,
						       filter_elements);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<class Inte>
template<ElementType type>
void IntegratorCohesive<Inte>::checkJacobians(const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();
  //  const ElementType sub_type = ElementType(CohesiveElementSubElementType<type>::value);
  //  sub_type_integrator->checkJacobians<sub_type>(ghost_type);

  UInt nb_quadrature_points = CohesiveElement<type>::getNbQuadraturePoints();

  UInt nb_element;

  nb_element = Inte::mesh->getConnectivity(type,ghost_type).getSize();

  Real * jacobians_val = Inte::jacobians(type, ghost_type).storage();

  for (UInt i = 0; i < nb_element*nb_quadrature_points; ++i,++jacobians_val){
    AKANTU_DEBUG_ASSERT(*jacobians_val >0,
			"Negative jacobian computed,"
			<< " possible problem in the element node order");
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/* specialization                                                             */
/* -------------------------------------------------------------------------- */

template<>
template<>
inline void IntegratorCohesive<IntegratorGauss>
::precomputeJacobiansOnQuadraturePoints<_not_defined>(__attribute__((unused)) const GhostType & ghost_type) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template<>
template<>
inline void IntegratorCohesive<IntegratorGauss>
::integrateOnElement<_not_defined>(__attribute__((unused)) const Vector<Real> & f,
				   __attribute__((unused)) Real * intf,
				   __attribute__((unused)) UInt nb_degree_of_freedom,
				   __attribute__((unused)) const UInt elem,
				   __attribute__((unused)) const GhostType & ghost_type) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template<>
template<>
inline const Vector<Real> &
IntegratorCohesive<IntegratorGauss>
::getQuadraturePoints<_not_defined>(__attribute__((unused)) const GhostType & ghost_type) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template<>
template<>
inline void IntegratorCohesive<IntegratorGauss>
::computeQuadraturePoints<_not_defined>(__attribute__((unused)) const GhostType & ghost_type) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template<>
template<>
inline void IntegratorCohesive<IntegratorGauss>::
computeJacobianOnQuadPointsByElement<_not_defined>(__attribute__((unused)) UInt spatial_dimension,
						   __attribute__((unused)) Real * node_coords,
						   __attribute__((unused)) UInt nb_nodes_per_element,
						   __attribute__((unused)) Real * jacobians) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template<>
template<>
inline void IntegratorCohesive<IntegratorGauss>
::integrate<_not_defined>(__attribute__((unused)) const Vector<Real> & in_f,
			  __attribute__((unused)) Vector<Real> &intf,
			  __attribute__((unused)) UInt nb_degree_of_freedom,
			  __attribute__((unused)) const GhostType & ghost_type,
			  __attribute__((unused)) const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}


/* -------------------------------------------------------------------------- */
template<>
template<>
inline Real IntegratorCohesive<IntegratorGauss>
::integrate<_not_defined>(__attribute__((unused)) const Vector<Real> & in_f,
			  __attribute__((unused)) const GhostType & ghost_type,
			  __attribute__((unused)) const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}


/* -------------------------------------------------------------------------- */
template<>
template<>
inline void IntegratorCohesive<IntegratorGauss>
::integrateOnQuadraturePoints<_not_defined>(__attribute__((unused)) const Vector<Real> & in_f,
					    __attribute__((unused)) Vector<Real> &intf,
					    __attribute__((unused)) UInt nb_degree_of_freedom,
					    __attribute__((unused)) const GhostType & ghost_type,
					    __attribute__((unused)) const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template<>
template<>
inline void IntegratorCohesive<IntegratorGauss>
::checkJacobians<_not_defined>(__attribute__((unused)) const GhostType & ghost_type) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

