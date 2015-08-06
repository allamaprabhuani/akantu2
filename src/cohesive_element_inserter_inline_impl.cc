/**
 * @file   cohesive_element_inserter_inline_impl.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 *
 * @brief  Cohesive element inserter inline functions
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
inline UInt CohesiveElementInserter::getNbDataForElements(const Array<Element> & elements,
							  SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;

  if (tag == _gst_ce_inserter) {
    UInt nb_nodes = 0;

    Array<Element>::const_iterator<Element> it  = elements.begin();
    Array<Element>::const_iterator<Element> end = elements.end();
    for (; it != end; ++it) {
      const Element & el = *it;
      nb_nodes += Mesh::getNbNodesPerElement(el.type);
    }

    size += nb_nodes * sizeof(UInt);
  }

  if (tag == _gst_ce_groups) {

    size = elements.getSize() * (sizeof(bool) + sizeof(unsigned int));

  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void CohesiveElementInserter::packElementData(CommunicationBuffer & buffer,
						     const Array<Element> & elements,
						     SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  if (tag == _gst_ce_inserter)
    packUnpackGlobalConnectivity<true>(buffer, elements);

  if (tag == _gst_ce_groups)
    packUnpackGroupedInsertionData<true>(buffer,elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void CohesiveElementInserter::unpackElementData(CommunicationBuffer & buffer,
						       const Array<Element> & elements,
						       SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  if (tag == _gst_ce_inserter)
    packUnpackGlobalConnectivity<false>(buffer, elements);

  if (tag == _gst_ce_groups)
    packUnpackGroupedInsertionData<false>(buffer,elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<bool pack_mode>
inline void CohesiveElementInserter::packUnpackGlobalConnectivity(CommunicationBuffer & buffer,
								  const Array<Element> & elements) const {
  AKANTU_DEBUG_IN();

  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;

  Array<UInt>::iterator<Vector<UInt> > conn_begin;
  UInt nb_nodes_per_elem = 0;
  UInt index;

  Array<Element>::const_iterator<Element> it  = elements.begin();
  Array<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;

    if (el.type != current_element_type || el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type   = el.ghost_type;

      nb_nodes_per_elem = Mesh::getNbNodesPerElement(current_element_type);

      conn_begin = mesh.connectivities(current_element_type,
				       current_ghost_type).begin(nb_nodes_per_elem);
    }

    /// get element connectivity
    Vector<UInt> current_conn = conn_begin[el.element];

    /// loop on all connectivity nodes
    for (UInt n = 0; n < nb_nodes_per_elem; ++n) {
      UInt node = current_conn(n);

      if (pack_mode) {
	/// if node is local or master pack its global id, otherwise
	/// dummy data
	if (mesh.isLocalOrMasterNode(node))
	  index = mesh.getNodeGlobalId(node);
	else
	  index = UInt(-1);

	buffer << index;
      }
      else {
	buffer >> index;

	/// update slave nodes' index
	if (index != UInt(-1) && mesh.isSlaveNode(node))
	  (*mesh.nodes_global_ids)(node) = index;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<bool pack_mode>
inline void CohesiveElementInserter::packUnpackGroupedInsertionData(CommunicationBuffer & buffer,
								    const Array<Element> & elements) const {

  AKANTU_DEBUG_IN();
  
  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;
  ElementTypeMapArray<UInt> & physical_names = mesh_facets.registerData<UInt>("physical_names");

  Array<bool> *vect = NULL;
  Array<unsigned int> *vect2 = NULL;

  Array<Element>::const_iterator<Element> it  = elements.begin();
  Array<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    if(el.type != current_element_type || el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type   = el.ghost_type;
      vect = &const_cast<Array<bool> &>(insertion_facets(el.type, el.ghost_type));
      vect2 = &(physical_names(el.type, el.ghost_type));
    }

    Vector<bool> data(vect->storage() + el.element, 1);
    Vector<unsigned int> data2(vect2->storage() + el.element, 1);

    if(pack_mode) {
      buffer << data;
      buffer << data2;
    }
    else {
      buffer >> data;
      buffer >> data2;
    }
  }

  AKANTU_DEBUG_OUT();
}
