/**
 * @file   model_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date   Wed Aug 25 08:50:54 2010
 *
 * @brief  inline implementation of the model class
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
inline SynchronizerRegistry & Model::getSynchronizerRegistry(){
  AKANTU_DEBUG_ASSERT(synch_registry,"synchronizer registry not initialized:"
		      << " did you call createSynchronizerRegistry ?");
  return  *synch_registry;
}

/* -------------------------------------------------------------------------- */
template <typename FEMClass>
inline FEMClass & Model::getFEMClassBoundary(std::string name) {
  AKANTU_DEBUG_IN();

  if (name == "") name = default_fem;

  FEMMap::const_iterator it_boun = fems_boundary.find(name);

  FEMClass * tmp_fem_boundary;

  if (it_boun == fems_boundary.end()){
    AKANTU_DEBUG_INFO("Creating FEM boundary " << name);

    FEMMap::const_iterator it = fems.find(name);
    AKANTU_DEBUG_ASSERT(it != fems.end(), "The FEM " << name << " is not registered");

    UInt spatial_dimension = it->second->getElementDimension();
    std::stringstream sstr; sstr << id << ":fem_boundary:" << name;

    tmp_fem_boundary = new FEMClass(it->second->getMesh(),
				    spatial_dimension-1,
				    sstr.str(),
				    memory_id);
    fems_boundary[name] = tmp_fem_boundary;
  } else {
    tmp_fem_boundary = dynamic_cast<FEMClass *>(it_boun->second);
  }

  AKANTU_DEBUG_OUT();
  return *tmp_fem_boundary;
}


/* -------------------------------------------------------------------------- */
template <typename FEMClass>
inline FEMClass & Model::getFEMClass(std::string name) const{
  AKANTU_DEBUG_IN();

  if (name == "") name = default_fem;

  FEMMap::const_iterator it = fems.find(name);
  AKANTU_DEBUG_ASSERT(it != fems.end(), "The FEM " << name << " is not registered");

  AKANTU_DEBUG_OUT();
  return dynamic_cast<FEMClass &>(*(it->second));
}

/* -------------------------------------------------------------------------- */

inline void Model::unRegisterFEMObject(const std::string & name){

  FEMMap::iterator it = fems.find(name);
  AKANTU_DEBUG_ASSERT(it != fems.end(), "FEM object with name "
		      << name << " was not found");

  delete((*it).second);
  fems.erase(it);
  if (!fems.empty())
    default_fem = (*fems.begin()).first;
}

/* -------------------------------------------------------------------------- */

template <typename FEMClass>
inline void Model::registerFEMObject(const std::string & name,
				     Mesh & mesh,
				     UInt spatial_dimension){
  if (fems.size() == 0) default_fem = name;

#ifndef AKANTU_NDEBUG
  FEMMap::iterator it = fems.find(name);
  AKANTU_DEBUG_ASSERT(it == fems.end(), "FEM object with name "
		      << name << " was already created");
#endif

  std::stringstream sstr; sstr << id << ":fem:" << name;
  fems[name] = new FEMClass(mesh, spatial_dimension, sstr.str(), memory_id);

  // MeshUtils::buildFacets(fems[name]->getMesh());

  // std::stringstream sstr2; sstr2 << id << ":fem_boundary:" << name;
  // fems_boundary[name] = new FEMClass(mesh, spatial_dimension-1, sstr2.str(), memory_id);
}

/* -------------------------------------------------------------------------- */
inline FEM & Model::getFEM(std::string name) const{
  AKANTU_DEBUG_IN();

  if (name == "") name = default_fem;

  FEMMap::const_iterator it = fems.find(name);
  AKANTU_DEBUG_ASSERT(it != fems.end(),"The FEM " << name << " is not registered");

  AKANTU_DEBUG_OUT();
  return *(it->second);
}


/* -------------------------------------------------------------------------- */
inline FEM & Model::getFEMBoundary(std::string name){
  AKANTU_DEBUG_IN();

  if (name == "") name = default_fem;

  FEMMap::const_iterator it = fems_boundary.find(name);
  AKANTU_DEBUG_ASSERT(it != fems_boundary.end(),
		      "The FEM boundary  " << name << " is not registered");
  AKANTU_DEBUG_ASSERT(it->second != NULL,
		      "The FEM boundary " << name << " was not created");

  AKANTU_DEBUG_OUT();
  return *(it->second);
}
/* -------------------------------------------------------------------------- */
/// @todo : should merge with a single function which handles local and global
inline void Model::changeLocalEquationNumberforPBC(std::map<UInt,UInt> & pbc_pair,
					    UInt dimension){
  for (std::map<UInt,UInt>::iterator it = pbc_pair.begin();
       it != pbc_pair.end();++it) {
    Int node_master = (*it).second;
    Int node_slave = (*it).first;
    for (UInt i = 0; i < dimension; ++i) {
      (*dof_synchronizer->getLocalDOFEquationNumbersPointer())
	(node_slave*dimension+i) = dimension*node_master+i;
      (*dof_synchronizer->getGlobalDOFEquationNumbersPointer())
	(node_slave*dimension+i) = dimension*node_master+i;
    }
  }
}
/* -------------------------------------------------------------------------- */
inline bool Model::getIsPBCSlaveNode(const UInt node) {
  // if no pbc is defined, is_pbc_slave_node is of size zero
  if (is_pbc_slave_node.getSize() == 0)
    return false;
  else
    return is_pbc_slave_node(node);
}

/* -------------------------------------------------------------------------- */
inline UInt Model::getNbQuadraturePoints(const Vector<Element> & elements) const {
  UInt nb_quad = 0;
  Vector<Element>::const_iterator<Element> it  = elements.begin();
  Vector<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    nb_quad += getFEM().getNbQuadraturePoints(el.type, el.ghost_type);
  }
  return nb_quad;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void Model::packElementalDataHelper(const ByElementTypeVector<T> & data_to_pack,
                                           CommunicationBuffer & buffer,
                                           const Vector<Element> & elements,
                                           bool per_quadrature_point_data) const {
  packUnpackElementalDataHelper<T, true>(const_cast<ByElementTypeVector<T> &>(data_to_pack),
                                         buffer,
                                         elements,
                                         per_quadrature_point_data);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void Model::unpackElementalDataHelper(ByElementTypeVector<T> & data_to_unpack,
                                             CommunicationBuffer & buffer,
                                             const Vector<Element> & elements,
                                             bool per_quadrature_point_data) const {
  packUnpackElementalDataHelper<T, false>(data_to_unpack, buffer, elements, per_quadrature_point_data);
}

/* -------------------------------------------------------------------------- */
template<typename T, bool pack_helper>
inline void Model::packUnpackElementalDataHelper(ByElementTypeVector<T> & data_to_pack,
                                                 CommunicationBuffer & buffer,
                                                 const Vector<Element> & element,
                                                 bool per_quadrature_point_data) const {
  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;
  UInt nb_quad_per_elem = 0;
  UInt nb_component = 0;

  Vector<T> * vect = NULL;

  Vector<Element>::const_iterator<Element> it  = element.begin();
  Vector<Element>::const_iterator<Element> end = element.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    if(el.type != current_element_type || el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type   = el.ghost_type;
      vect = &data_to_pack(el.type, el.ghost_type);
      if(per_quadrature_point_data)
        nb_quad_per_elem = this->getFEM().getNbQuadraturePoints(el.type, el.ghost_type);
      else nb_quad_per_elem = 1;
      nb_component = vect->getNbComponent();
    }

    types::Vector<T> data(vect->storage() + el.element * nb_component * nb_quad_per_elem,
			  nb_component * nb_quad_per_elem);
    if(pack_helper)
      buffer << data;
    else
      buffer >> data;
  }
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void Model::packNodalDataHelper(Vector<T> & data_to_pack,
                                       CommunicationBuffer & buffer,
                                       const Vector<Element> & element) const {
  packUnpackNodalDataHelper<T, true>(data_to_pack, buffer, element);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void Model::unpackNodalDataHelper(Vector<T> & data_to_unpack,
                                         CommunicationBuffer & buffer,
                                         const Vector<Element> & element) const {
  packUnpackNodalDataHelper<T, false>(data_to_unpack, buffer, element);
}

/* -------------------------------------------------------------------------- */
template<typename T, bool pack_helper>
inline void Model::packUnpackNodalDataHelper(Vector<T> & data,
                                             CommunicationBuffer & buffer,
                                             const Vector<Element> & elements) const {
  UInt nb_component = data.getNbComponent();
  UInt nb_nodes_per_element = 0;

  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;
  UInt * conn = NULL;

  Vector<Element>::const_iterator<Element> it  = elements.begin();
  Vector<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    if(el.type != current_element_type || el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type   = el.ghost_type;
      conn = mesh.getConnectivity(el.type, el.ghost_type).storage();
      nb_nodes_per_element = Mesh::getNbNodesPerElement(el.type);
    }

    UInt el_offset  = el.element * nb_nodes_per_element;
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n];
      types::Vector<T> data_vect(data.storage() + offset_conn * nb_component,
				 nb_component);

      if(pack_helper)
	buffer << data_vect;
      else
	buffer >> data_vect;
    }
  }
}
