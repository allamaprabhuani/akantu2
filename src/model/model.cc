/**
 * @file   model.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Oct 03 12:24:07 2011
 *
 * @brief  implementation of model common parts
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

#include "sub_boundary.hh"
#include "model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
Model::Model(Mesh& m, UInt dim, const ID & id,
	     const MemoryID & memory_id) :
  Memory(memory_id), mesh(m), id(id),
  spatial_dimension(dim == _all_dimensions ? m.getSpatialDimension() : dim),
  synch_registry(NULL),is_pbc_slave_node(0,1,"is_pbc_slave_node") {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Model::createSynchronizerRegistry(DataAccessor * data_accessor){
  synch_registry = new SynchronizerRegistry(*data_accessor);
}

/* -------------------------------------------------------------------------- */
void Model::setPBC(UInt x, UInt y, UInt z){
  Mesh & mesh = getFEM().getMesh();
  mesh.computeBoundingBox();
  if (x) MeshUtils::computePBCMap(mesh, 0, pbc_pair);
  if (y) MeshUtils::computePBCMap(mesh, 1, pbc_pair);
  if (z) MeshUtils::computePBCMap(mesh, 2, pbc_pair);
}

/* -------------------------------------------------------------------------- */
void Model::setPBC(SurfacePairList & surface_pairs,
		   ElementType surface_e_type){
  Mesh & mesh = getFEM().getMesh();

  SurfacePairList::iterator s_it;
  for(s_it = surface_pairs.begin(); s_it != surface_pairs.end(); ++s_it) {
    MeshUtils::computePBCMap(mesh, *s_it, surface_e_type, pbc_pair);
  }
}

/* -------------------------------------------------------------------------- */
void Model::initPBC() {
  Mesh & mesh = getFEM().getMesh();

  std::map<UInt,UInt>::iterator it = pbc_pair.begin();
  std::map<UInt,UInt>::iterator end = pbc_pair.end();

  is_pbc_slave_node.resize(mesh.getNbNodes());
#ifndef AKANTU_NDEBUG
  Real * coords = mesh.getNodes().values;
  UInt dim = mesh.getSpatialDimension();
#endif
  while(it != end){
    UInt i1 = (*it).first;

    is_pbc_slave_node(i1) = true;

#ifndef AKANTU_NDEBUG
    UInt i2 = (*it).second;
    AKANTU_DEBUG_INFO("pairing " << i1 << " ("
		      << coords[dim*i1] << "," << coords[dim*i1+1] << ","
		      << coords[dim*i1+2]
		      << ") with "
		      << i2 << " ("
		      << coords[dim*i2] << "," << coords[dim*i2+1] << ","
		      << coords[dim*i2+2]
		      << ")");
#endif
    ++it;
  }
}

/* -------------------------------------------------------------------------- */
Synchronizer & Model::createParallelSynch(MeshPartition * partition,
					  __attribute__((unused)) DataAccessor * data_accessor){
  AKANTU_DEBUG_IN();
  /* ------------------------------------------------------------------------ */
  /* Parallel initialization                                                  */
  /* ------------------------------------------------------------------------ */
  StaticCommunicator & comm =
    StaticCommunicator::getStaticCommunicator();
  Int prank = comm.whoAmI();

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
Model::~Model() {
  AKANTU_DEBUG_IN();

  FEMMap::iterator it;
  for (it = fems.begin(); it != fems.end(); ++it) {
    if(it->second) delete it->second;
  }

  for (it = fems_boundary.begin(); it != fems_boundary.end(); ++it) {
    if(it->second) delete it->second;
  }

  delete synch_registry;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Model::dumpBoundary(const std::string & boundary_name) {
  Boundary & boundary =  mesh.getBoundary();
  Boundary::iterator bit = boundary.find(boundary_name);
  Boundary::iterator bend = boundary.end();

  if(bit == bend) {
    AKANTU_EXCEPTION("No boundary named " << boundary_name << "!");
  }
  bit->dump();

}

/* -------------------------------------------------------------------------- */
void Model::dumpBoundary() {
  Boundary & boundary =  mesh.getBoundary();
  Boundary::iterator bit = boundary.begin();
  Boundary::iterator bend = boundary.end();

  for(; bit != bend; ++bit) {
    bit->dump();
  }
}

/* -------------------------------------------------------------------------- */
void Model::setBoundaryDirectory(const std::string & directory) {
  Boundary & boundary     = mesh.getBoundary();
  Boundary::iterator bit  = boundary.begin();
  Boundary::iterator bend = boundary.end();

  for (; bit != bend; ++bit) {
    bit->setDirectory(directory);
  }
}

/* -------------------------------------------------------------------------- */
void Model::setBoundaryDirectory(const std::string & directory,
				 const std::string & boundary_name) {
  Boundary & boundary =  mesh.getBoundary();
  Boundary::iterator bit = boundary.find(boundary_name);
  Boundary::iterator bend = boundary.end();
  
  if(bit == bend) {
    AKANTU_EXCEPTION("No boundary named " << boundary_name << "!");
  }
  bit->setDirectory(directory);
}

/* -------------------------------------------------------------------------- */
void Model::setBoundaryBaseName(const std::string & basename,
				const std::string & boundary_name) {
  Boundary & boundary =  mesh.getBoundary();
  Boundary::iterator bit = boundary.find(boundary_name);
  Boundary::iterator bend = boundary.end();
  
  if(bit == bend) {
    AKANTU_EXCEPTION("No boundary named " << boundary_name << "!");
  }
  bit->setBaseName(basename);
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
