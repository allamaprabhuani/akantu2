/**
 * @file   integrator_cohesive_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Feb 22 17:22:21 2012
 *
 * @brief  
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
  Integrator(mesh, id, memory_id) {
  AKANTU_DEBUG_IN();

  std::stringstream sstr;
  sstr << id << "sub_integrator";
  sub_type_integrator = new Inte(mesh, sstr.str(), memory_id);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<class Inte>
template <ElementType type>
void IntegratorCohesive<Inte>::precomputeJacobiansOnQuadraturePoints(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();
  const ElementType sub_type = ElementType(CohesiveElementSubElementType<type>::value);
  sub_type_integrator->precomputeJacobiansOnQuadraturePoints<sub_type>(ghost_type);
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
  sub_type_integrator->integrateOnElement<sub_type>(f, intf, nb_degree_of_freedom, elem, ghost_type);
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
  return sub_type_integrator->getQuadraturePoints<sub_type>(ghost_type);
}

/* -------------------------------------------------------------------------- */
template<class Inte>
template <ElementType type>
void IntegratorCohesive<Inte>::computeQuadraturePoints(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();
  const ElementType sub_type = ElementType(CohesiveElementSubElementType<type>::value);
  sub_type_integrator->template computeQuadraturePoints<sub_type>(ghost_type);
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
  sub_type_integrator->computeJacobianOnQuadPointsByElement<sub_type>(spatial_dimension,
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
  const ElementType sub_type = ElementType(CohesiveElementSubElementType<type>::value);
  sub_type_integrator->integrate<sub_type>(in_f, intf, nb_degree_of_freedom, ghost_type, filter_elements);
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template<class Inte>
template <ElementType type>
Real IntegratorCohesive<Inte>::integrate(const Vector<Real> & in_f,
					       const GhostType & ghost_type,
					       const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();
  const ElementType sub_type = ElementType(CohesiveElementSubElementType<type>::value);
  AKANTU_DEBUG_OUT();
  return  sub_type_integrator->integrate<sub_type>(in_f, ghost_type, filter_elements);
}


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
  sub_type_integrator->integrateOnQuadraturePoints<sub_type>(in_f,
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
  const ElementType sub_type = ElementType(CohesiveElementSubElementType<type>::value);
  sub_type_integrator->checkJacobians<sub_type>(ghost_type);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/* specialization                                                             */
/* -------------------------------------------------------------------------- */

template<>
template<>
inline void IntegratorCohesive<IntegratorGauss>::precomputeJacobiansOnQuadraturePoints<_not_defined>(__attribute__((unused)) const GhostType & ghost_type) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template<>
template<>
inline void IntegratorCohesive<IntegratorGauss>::integrateOnElement<_not_defined>(__attribute__((unused)) const Vector<Real> & f,
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
IntegratorCohesive<IntegratorGauss>::getQuadraturePoints<_not_defined>(__attribute__((unused)) const GhostType & ghost_type) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template<>
template<>
inline void IntegratorCohesive<IntegratorGauss>::computeQuadraturePoints<_not_defined>(__attribute__((unused)) const GhostType & ghost_type) {
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
inline void IntegratorCohesive<IntegratorGauss>::integrate<_not_defined>(__attribute__((unused)) const Vector<Real> & in_f,
									 __attribute__((unused)) Vector<Real> &intf,
									 __attribute__((unused)) UInt nb_degree_of_freedom,
									 __attribute__((unused)) const GhostType & ghost_type,
									 __attribute__((unused)) const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}


/* -------------------------------------------------------------------------- */
template<>
template<>
inline Real IntegratorCohesive<IntegratorGauss>::integrate<_not_defined>(__attribute__((unused)) const Vector<Real> & in_f,
									 __attribute__((unused)) const GhostType & ghost_type,
									 __attribute__((unused)) const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}


/* -------------------------------------------------------------------------- */
template<>
template<>
inline void IntegratorCohesive<IntegratorGauss>::integrateOnQuadraturePoints<_not_defined>(__attribute__((unused)) const Vector<Real> & in_f,
											   __attribute__((unused)) Vector<Real> &intf,
											   __attribute__((unused)) UInt nb_degree_of_freedom,
											   __attribute__((unused)) const GhostType & ghost_type,
											   __attribute__((unused)) const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template<>
template<>
inline void IntegratorCohesive<IntegratorGauss>::checkJacobians<_not_defined>(__attribute__((unused)) const GhostType & ghost_type) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

