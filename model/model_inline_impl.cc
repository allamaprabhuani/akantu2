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
inline Model::Model(const ModelID & id,
		    const MemoryID & memory_id) :
  Memory(memory_id), id(id) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline Model::~Model() {
  AKANTU_DEBUG_IN();
  
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <typename FEMClass>
inline FEMClass & Model::getFEMBoundary(std::string name){
  AKANTU_DEBUG_IN();
  if (name == "") name = default_fem;
  if (fems_boundary[name] == NULL){
    UInt spatial_dimension = fems[name]->getElementDimension();
    std::stringstream sstr; sstr << id << ":femboundary:" << name;
    MeshUtils::buildFacets(fems[name]->getMesh());
    fems_boundary[name] = new FEMClass(fems[name]->getMesh(), 
				       spatial_dimension-1, 
				       sstr.str(), 
				       memory_id); 
  }
  AKANTU_DEBUG_OUT();
  return dynamic_cast<FEMClass>(*fems_boundary[name]);
}


/* -------------------------------------------------------------------------- */
template <typename FEMClass>
inline FEMClass & Model::getFEM(std::string name) const{
  AKANTU_DEBUG_IN();
  if (name == "") name = default_fem;
  AKANTU_DEBUG_OUT();
  std::map<std::string,FEM*>::const_iterator it = fems.find(name);
  AKANTU_DEBUG_ASSERT(it != fems.end(),"fem  " << name << " not registered");
  return dynamic_cast<FEMClass &>(*(it->second));
}

/* -------------------------------------------------------------------------- */

template <typename FEMClass>
inline void Model::registerFEMObject(const std::string & name,
				     Mesh & mesh, 
				     UInt spatial_dimension){
  if (fems.size() == 0) default_fem = name;
  AKANTU_DEBUG_ASSERT(fems.count(name) == 0,"fem object with name " 
		      << name 
		      << " was already created: abort");

  std::stringstream sstr; sstr << id << ":fem:" << name;
  fems[name] = new FEMClass(mesh, spatial_dimension, sstr.str(), memory_id); 
  MeshUtils::buildFacets(fems[name]->getMesh());
  std::stringstream sstr2; sstr2 << id << ":femboundary" << name;
  fems_boundary[name] = new FEMClass(mesh, spatial_dimension-1, sstr2.str(), memory_id); 
}

/* -------------------------------------------------------------------------- */
inline FEM & Model::getFEM(std::string name) const{
  AKANTU_DEBUG_IN();
  if (name == "") name = default_fem;
  AKANTU_DEBUG_OUT();
  std::map<std::string,FEM*>::const_iterator it = fems.find(name);
  AKANTU_DEBUG_ASSERT(it != fems.end(),"fem  " << name << " not registered");
  return *(it->second);
}


/* -------------------------------------------------------------------------- */
inline FEM & Model::getFEMBoundary(std::string name){
  AKANTU_DEBUG_IN();
  if (name == "") name = default_fem;
  std::map<std::string,FEM*>::const_iterator it = fems_boundary.find(name);
  AKANTU_DEBUG_ASSERT(it != fems_boundary.end(),"fem_boundary  " << name << " not registered");
  AKANTU_DEBUG_ASSERT(fems_boundary[name] != NULL,"fem boundary was not created");
  AKANTU_DEBUG_OUT();
  return *fems_boundary[name];
}
