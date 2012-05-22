/**
 * @file   model_inline_impl.cc
 * @author Guillaume ANCIAUX <guillaume.anciaux@epfl.ch>
 * @date   Fri Aug 20 17:18:08 2010
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

  FEMMap::const_iterator it = fems.find(name);
  AKANTU_DEBUG_ASSERT(it == fems.end(), "FEM object with name "
		      << name << " was already created");

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
