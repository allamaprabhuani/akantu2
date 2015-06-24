/**
 * @file   mesh_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jul 15 2010
 * @date last modification: Fri Sep 05 2014
 *
 * @brief  Implementation of the inline functions of the mesh class
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#if defined(AKANTU_COHESIVE_ELEMENT)
#  include "cohesive_element.hh"
#endif

#ifndef __AKANTU_MESH_INLINE_IMPL_CC__
#define __AKANTU_MESH_INLINE_IMPL_CC__


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline RemovedNodesEvent::RemovedNodesEvent(const Mesh & mesh) :
  new_numbering(mesh.getNbNodes(), 1, "new_numbering") {
}

/* -------------------------------------------------------------------------- */
inline RemovedElementsEvent::RemovedElementsEvent(const Mesh & mesh) :
  new_numbering("new_numbering", mesh.getID()) {
}

/* -------------------------------------------------------------------------- */
template <>
inline void Mesh::sendEvent<RemovedElementsEvent>(RemovedElementsEvent & event) {
  connectivities.onElementsRemoved(event.getNewNumbering());
  EventHandlerManager<MeshEventHandler>::sendEvent(event);
}

/* -------------------------------------------------------------------------- */
template <>
inline void Mesh::sendEvent<RemovedNodesEvent>(RemovedNodesEvent & event) {
  if(created_nodes)    removeNodesFromArray(*nodes           , event.getNewNumbering());
  if(nodes_global_ids) removeNodesFromArray(*nodes_global_ids, event.getNewNumbering());
  if(nodes_type.getSize() != 0)  removeNodesFromArray(nodes_type      , event.getNewNumbering());

  EventHandlerManager<MeshEventHandler>::sendEvent(event);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void Mesh::removeNodesFromArray(Array<T> & vect, const Array<UInt> & new_numbering) {
  Array<T> tmp(vect.getSize(), vect.getNbComponent());
  UInt nb_component = vect.getNbComponent();
  UInt new_nb_nodes = 0;
  for (UInt i = 0; i < new_numbering.getSize(); ++i) {
    UInt new_i = new_numbering(i);
    if(new_i != UInt(-1)) {
      memcpy(tmp.storage() + new_i * nb_component,
	     vect.storage() + i * nb_component,
	     nb_component * sizeof(T));
      ++new_nb_nodes;
    }
  }

  tmp.resize(new_nb_nodes);
  vect.copy(tmp);
}


/* -------------------------------------------------------------------------- */
#ifdef AKANTU_CORE_CXX11
template <typename... Args>
inline void Mesh::translate(Args... params) {
  // check that the number of parameters corresponds to the dimension
  AKANTU_DEBUG_ASSERT(sizeof...(Args) <= spatial_dimension , "Number of arguments greater than dimension.");

  // unpack parameters
  Real s[] = { params... };

  Array<Real>& nodes = getNodes();
  for (UInt i = 0; i < nodes.getSize(); ++i)
    for (UInt k = 0; k < sizeof...(Args); ++k)
      nodes(i, k) += s[k];
}
#endif

/* -------------------------------------------------------------------------- */
inline UInt Mesh::elementToLinearized(const Element & elem) const {
  AKANTU_DEBUG_ASSERT(elem.type < _max_element_type &&
		      elem.element < types_offsets.storage()[elem.type+1],
		      "The element " << elem
		      << "does not exists in the mesh " << getID());

  return types_offsets.storage()[elem.type] + elem.element;
}

/* -------------------------------------------------------------------------- */
inline Element Mesh::linearizedToElement (UInt linearized_element) const {

  UInt t;

  for (t = _not_defined;
       t != _max_element_type && linearized_element >= types_offsets(t);
       ++t);

  AKANTU_DEBUG_ASSERT(linearized_element < types_offsets(t),
   		      "The linearized element " << linearized_element
   		      << "does not exists in the mesh " << getID());

  --t;
  ElementType type = ElementType(t);
  return Element(type,
		 linearized_element - types_offsets.storage()[t],
		 _not_ghost,
		 getKind(type));
}

/* -------------------------------------------------------------------------- */
inline void Mesh::updateTypesOffsets(const GhostType & ghost_type) {
  Array<UInt> * types_offsets_ptr = &this->types_offsets;
  if(ghost_type == _ghost) types_offsets_ptr = &this->ghost_types_offsets;
  Array<UInt> & types_offsets = *types_offsets_ptr;

  types_offsets.clear();
  type_iterator it   = firstType(_all_dimensions, ghost_type, _ek_not_defined);
  type_iterator last = lastType(_all_dimensions, ghost_type, _ek_not_defined);

  for (; it != last; ++it)
    types_offsets(*it) = connectivities(*it, ghost_type).getSize();

  for (UInt t = _not_defined + 1; t < _max_element_type; ++t)
    types_offsets(t) += types_offsets(t - 1);
  for (UInt t = _max_element_type; t > _not_defined; --t)
    types_offsets(t) = types_offsets(t - 1);
  types_offsets(0) = 0;
}

/* -------------------------------------------------------------------------- */
inline const Mesh::ConnectivityTypeList & Mesh::getConnectivityTypeList(const GhostType & ghost_type) const {
  if (ghost_type == _not_ghost)
    return type_set;
  else
    return ghost_type_set;
}

/* -------------------------------------------------------------------------- */
inline Array<UInt> * Mesh::getNodesGlobalIdsPointer() {
  AKANTU_DEBUG_IN();
  if(nodes_global_ids == NULL) {
    std::stringstream sstr; sstr << getID() << ":nodes_global_ids";
    nodes_global_ids = &(alloc<UInt>(sstr.str(), nodes->getSize(), 1));
  }
  AKANTU_DEBUG_OUT();
  return nodes_global_ids;
}

/* -------------------------------------------------------------------------- */
inline Array<Int> * Mesh::getNodesTypePointer() {
  AKANTU_DEBUG_IN();
  if(nodes_type.getSize() == 0) {
    nodes_type.resize(nodes->getSize());
    nodes_type.set(-1);
  }
  AKANTU_DEBUG_OUT();
  return &nodes_type;
}

/* -------------------------------------------------------------------------- */
inline Array<UInt> * Mesh::getConnectivityPointer(const ElementType & type,
						  const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  Array<UInt> * tmp;
  if(!connectivities.exists(type, ghost_type)) {
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    tmp = &(connectivities.alloc(0, nb_nodes_per_element,
				 type, ghost_type));

    AKANTU_DEBUG_INFO("The connectivity vector for the type "
		      << type << " created");

    if (ghost_type == _not_ghost) type_set.insert(type);
    else ghost_type_set.insert(type);

    updateTypesOffsets(ghost_type);
  } else {
    tmp = &connectivities(type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
  return tmp;
}

/* -------------------------------------------------------------------------- */
inline Array< std::vector<Element> > * Mesh::getElementToSubelementPointer(const ElementType & type,
                                                                           const GhostType & ghost_type) {
  Array< std::vector<Element> > * tmp =
    getDataPointer< std::vector<Element> >("element_to_subelement", type, ghost_type, 1, true);
  return tmp;
}

/* -------------------------------------------------------------------------- */
inline Array<Element > * Mesh::getSubelementToElementPointer(const ElementType & type,
                                                             const GhostType & ghost_type) {
  Array<Element> * tmp =
    getDataPointer<Element>("subelement_to_element", type, ghost_type,
			    getNbFacetsPerElement(type), true, is_mesh_facets);
  return tmp;
}


/* -------------------------------------------------------------------------- */
inline const Array< std::vector<Element> > & Mesh::getElementToSubelement(const ElementType & type,
                                                                          const GhostType & ghost_type) const {
  return getData< std::vector<Element> >("element_to_subelement", type, ghost_type);
}

/* -------------------------------------------------------------------------- */
inline Array< std::vector<Element> > & Mesh::getElementToSubelement(const ElementType & type,
                                                                    const GhostType & ghost_type) {
  return getData< std::vector<Element> >("element_to_subelement", type, ghost_type);
}

/* -------------------------------------------------------------------------- */
inline const Array<Element> & Mesh::getSubelementToElement(const ElementType & type,
                                                           const GhostType & ghost_type) const {
  return getData<Element>("subelement_to_element", type, ghost_type);
}

/* -------------------------------------------------------------------------- */
inline Array<Element> & Mesh::getSubelementToElement(const ElementType & type,
                                                     const GhostType & ghost_type) {
  return getData<Element>("subelement_to_element", type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline Array<T> * Mesh::getDataPointer(const std::string & data_name,
				       const ElementType & el_type,
                                       const GhostType & ghost_type,
                                       UInt nb_component,
				       bool size_to_nb_element,
				       bool resize_with_parent) {
  Array<T> & tmp = mesh_data.getElementalDataArrayAlloc<T>(data_name,
                                                           el_type, ghost_type,
                                                           nb_component);

  if (size_to_nb_element) {
    if (resize_with_parent)
      tmp.resize(mesh_parent->getNbElement(el_type, ghost_type));
    else
      tmp.resize(this->getNbElement(el_type, ghost_type));
  } else {
    tmp.resize(0);
  }

  return &tmp;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline const Array<T> & Mesh::getData(const std::string & data_name,
				      const ElementType & el_type,
				      const GhostType & ghost_type) const {
  return mesh_data.getElementalDataArray<T>(data_name, el_type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline Array<T> & Mesh::getData(const std::string & data_name,
				const ElementType & el_type,
                                const GhostType & ghost_type) {
  return mesh_data.getElementalDataArray<T>(data_name, el_type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline const ElementTypeMapArray<T> & Mesh::getData(const std::string & data_name) const {
  return mesh_data.getElementalData<T>(data_name);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline ElementTypeMapArray<T> & Mesh::getData(const std::string & data_name) {
  return mesh_data.getElementalData<T>(data_name);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline ElementTypeMapArray<T> & Mesh::registerData(const std::string & data_name) {
  this->mesh_data.registerElementalData<T>(data_name);
  return this->getData<T>(data_name);
}



/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbElement(const ElementType & type,
			       const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  try {
    const Array<UInt> & conn = connectivities(type, ghost_type);
    AKANTU_DEBUG_OUT();
    return conn.getSize();
  } catch (...) {
    AKANTU_DEBUG_OUT();
    return 0;
  }
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbElement(const UInt spatial_dimension,
			       const GhostType & ghost_type,
			       const ElementKind & kind) const {
  AKANTU_DEBUG_IN();
  UInt nb_element = 0;

  type_iterator it   = firstType(spatial_dimension, ghost_type, kind);
  type_iterator last = lastType(spatial_dimension, ghost_type, kind);
  for (; it != last; ++it) nb_element += getNbElement(*it, ghost_type);

  AKANTU_DEBUG_OUT();
  return nb_element;
}


/* -------------------------------------------------------------------------- */
inline void Mesh::getBarycenter(UInt element, const ElementType & type,
				Real * barycenter,
				GhostType ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt * conn_val = getConnectivity(type, ghost_type).storage();
  UInt nb_nodes_per_element = getNbNodesPerElement(type);

  Real local_coord[spatial_dimension * nb_nodes_per_element];

  UInt offset = element * nb_nodes_per_element;
  for (UInt n = 0; n < nb_nodes_per_element; ++n) {
    memcpy(local_coord + n * spatial_dimension,
	   nodes->storage() + conn_val[offset + n] * spatial_dimension,
	   spatial_dimension*sizeof(Real));
  }

  Math::barycenter(local_coord, nb_nodes_per_element, spatial_dimension, barycenter);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void Mesh::getBarycenter(const Element & element, Vector<Real> & barycenter) const {
  getBarycenter(element.element, element.type, barycenter.storage(), element.ghost_type);
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbNodesPerElement(const ElementType & type) {
  UInt nb_nodes_per_element = 0;
#define GET_NB_NODES_PER_ELEMENT(type)					\
  nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement()
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_NB_NODES_PER_ELEMENT);
#undef GET_NB_NODES_PER_ELEMENT
  return nb_nodes_per_element;
}

/* -------------------------------------------------------------------------- */
inline ElementType Mesh::getP1ElementType(const ElementType & type) {
  ElementType p1_type = _not_defined;
#define GET_P1_TYPE(type)				\
  p1_type = ElementClass<type>::getP1ElementType()

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_P1_TYPE);
#undef GET_P1_TYPE
  return p1_type;
}

/* -------------------------------------------------------------------------- */
inline ElementKind Mesh::getKind(const ElementType & type) {
  ElementKind kind = _ek_not_defined;
#define GET_KIND(type)				\
  kind = ElementClass<type>::getKind()
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_KIND);
#undef GET_KIND
  return kind;
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getSpatialDimension(const ElementType & type) {
  UInt spatial_dimension = 0;
#define GET_SPATIAL_DIMENSION(type)				\
  spatial_dimension = ElementClass<type>::getSpatialDimension()
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SPATIAL_DIMENSION);
#undef GET_SPATIAL_DIMENSION

  return spatial_dimension;
}

/* -------------------------------------------------------------------------- */
inline ElementType Mesh::getFacetType(const ElementType & type) {
  ElementType surface_type = _not_defined;
#define GET_FACET_TYPE(type)					\
  surface_type = ElementClass<type>::getFacetType()
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_FACET_TYPE);
#undef GET_FACET_TYPE

  return surface_type;
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbFacetsPerElement(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt n_facet = 0;
#define GET_NB_FACET(type)				\
  n_facet = ElementClass<type>::getNbFacetsPerElement()

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_NB_FACET);
#undef GET_NB_FACET

  AKANTU_DEBUG_OUT();
  return n_facet;
}

/* -------------------------------------------------------------------------- */
inline MatrixProxy<UInt> Mesh::getFacetLocalConnectivity(const ElementType & type) {
  AKANTU_DEBUG_IN();

#define GET_FACET_CON(type)						\
  AKANTU_DEBUG_OUT();							\
  return ElementClass<type>::getFacetLocalConnectivityPerElement()

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_FACET_CON);
#undef GET_FACET_CON

  AKANTU_DEBUG_OUT();
  return Matrix<UInt>(); // This avoid a compilation warning but will certainly
			 // also cause a segfault if reached
}

/* -------------------------------------------------------------------------- */
inline Matrix<UInt> Mesh::getFacetConnectivity(UInt element,
					       const ElementType & type,
					       const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  Matrix<UInt> local_facets(getFacetLocalConnectivity(type), false);
  Matrix<UInt> facets(local_facets.rows(), local_facets.cols());

  const Array<UInt> & conn = connectivities(type, ghost_type);

  for (UInt f = 0; f < facets.rows(); ++f) {
    for (UInt n = 0; n < facets.cols(); ++n) {
      facets(f, n) = conn(element, local_facets(f, n));
    }
  }

  AKANTU_DEBUG_OUT();
  return facets;
}


/* -------------------------------------------------------------------------- */
template<typename T>
inline void Mesh::extractNodalValuesFromElement(const Array<T> & nodal_values,
						T * local_coord,
						UInt * connectivity,
						UInt n_nodes,
						UInt nb_degree_of_freedom) const {
  for (UInt n = 0; n < n_nodes; ++n) {
    memcpy(local_coord + n * nb_degree_of_freedom,
	   nodal_values.storage() + connectivity[n] * nb_degree_of_freedom,
	   nb_degree_of_freedom * sizeof(T));
  }
}

/* -------------------------------------------------------------------------- */
inline void Mesh::addConnectivityType(const ElementType & type,
				      const GhostType & ghost_type){
  getConnectivityPointer(type, ghost_type);
}

/* -------------------------------------------------------------------------- */
inline bool Mesh::isPureGhostNode(UInt n) const {
  return nodes_type.getSize() ? (nodes_type(n) == -3) : false;
}

/* -------------------------------------------------------------------------- */
inline bool Mesh::isLocalOrMasterNode(UInt n) const {
  return nodes_type.getSize() ? (nodes_type(n) == -2) || (nodes_type(n) == -1) : true;
}


/* -------------------------------------------------------------------------- */
inline bool Mesh::isLocalNode(UInt n) const {
  return nodes_type.getSize() ? nodes_type(n) == -1 : true;
}

/* -------------------------------------------------------------------------- */
inline bool Mesh::isMasterNode(UInt n) const {
  return nodes_type.getSize() ? nodes_type(n) == -2 : false;
}

/* -------------------------------------------------------------------------- */
inline bool Mesh::isSlaveNode(UInt n) const {
  return nodes_type.getSize() ? nodes_type(n) >= 0 : false;
}

/* -------------------------------------------------------------------------- */
inline Int Mesh::getNodeType(UInt local_id) const {
  return nodes_type.getSize() ? nodes_type(local_id) : -1;
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNodeGlobalId(UInt local_id) const {
  return nodes_global_ids ? (*nodes_global_ids)(local_id) : local_id;
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbGlobalNodes() const {
  return nodes_global_ids ? nb_global_nodes : nodes->getSize();
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__


#endif /* __AKANTU_MESH_INLINE_IMPL_CC__ */
