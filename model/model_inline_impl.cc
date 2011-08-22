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
inline Model::Model(const ID & id,
		    const MemoryID & memory_id) :
  Memory(memory_id), id(id),synch_registry(NULL) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
inline SynchronizerRegistry & Model::getSynchronizerRegistry(){
  AKANTU_DEBUG_ASSERT(synch_registry,"synchronizer registry not initialized:"
		      << " did you call createSynchronizerRegistry ?");
  return  *synch_registry;
}
/* -------------------------------------------------------------------------- */
inline void Model::createSynchronizerRegistry(DataAccessor * data_accessor){
  synch_registry = new SynchronizerRegistry(*data_accessor);
}

/* -------------------------------------------------------------------------- */
inline void Model::initPBC(UInt x, UInt y, UInt z){
  Mesh & mesh = getFEM().getMesh();
  mesh.computeBoundingBox();
  if (x) MeshUtils::computePBCMap(mesh,0,pbc_pair);
  if (y) MeshUtils::computePBCMap(mesh,1,pbc_pair);
  if (z) MeshUtils::computePBCMap(mesh,2,pbc_pair);

  std::map<UInt,UInt>::iterator it = pbc_pair.begin();
  std::map<UInt,UInt>::iterator end = pbc_pair.end();
  
  Real * coords = mesh.getNodes().values;
  UInt dim = mesh.getSpatialDimension();
  while(it != end){
    UInt i1 = (*it).first;
    UInt i2 = (*it).second;
    
    AKANTU_DEBUG_INFO("pairing " << i1 << "(" 
		      << coords[dim*i1] << "," << coords[dim*i1+1] << "," 
		      << coords[dim*i1+2]
		      << ") with"
		      << i2 << "(" 
		      << coords[dim*i2] << "," << coords[dim*i2+1] << "," 
		      << coords[dim*i2+2]
		      << ")");	
    ++it;
  }

}
/* -------------------------------------------------------------------------- */
inline Synchronizer & Model::createParallelSynch(MeshPartition * partition,
						 __attribute__((unused)) DataAccessor * data_accessor){
  AKANTU_DEBUG_IN();
  /* ------------------------------------------------------------------------ */
  /* Parallel initialization                                                  */
  /* ------------------------------------------------------------------------ */
  StaticCommunicator * comm = 
    StaticCommunicator::getStaticCommunicator();
  Int prank = comm->whoAmI();
  
  DistributedSynchronizer * synch = NULL;
  if(prank == 0) 
    synch = 
      DistributedSynchronizer::createDistributedSynchronizerMesh(getFEM().getMesh(),
								 partition);
  else 
    synch = 
      DistributedSynchronizer::createDistributedSynchronizerMesh(getFEM().getMesh(),
								 NULL);

  AKANTU_DEBUG_OUT();
  return *synch;
}


/* -------------------------------------------------------------------------- */
inline Model::~Model() {
  AKANTU_DEBUG_IN();

  FEMMap::iterator it;
  for (it = fems.begin(); it != fems.end(); ++it) {
    if(it->second) delete it->second;
  }

  for (it = fems_boundary.begin(); it != fems_boundary.end(); ++it) {
    if(it->second) delete it->second;
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <typename FEMClass>
inline FEMClass & Model::getFEMClassBoundary(std::string name){
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


    MeshUtils::buildFacets(it->second->getMesh());

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
inline void Model::changeLocalEquationNumberforPBC(std::map<UInt,UInt> & pbc_pair,
					    UInt dimension){
  for (std::map<UInt,UInt>::iterator it = pbc_pair.begin(); 
       it != pbc_pair.end();++it) {
    Int node_master = (*it).second;
    Int node_slave = (*it).first;
    for (UInt i = 0; i < dimension; ++i) {
      (*dof_synchronizer->getLocalDOFEquationNumbersPointer())
	(node_slave*dimension+i) = dimension*node_master+i;
    }
  }
}
/* -------------------------------------------------------------------------- */

