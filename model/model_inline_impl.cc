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
inline Model::Model(Mesh & mesh,
	     UInt spatial_dimension,
	     const ModelID & id,
	     const MemoryID & memory_id) :
  Memory(memory_id), id(id) {
  AKANTU_DEBUG_IN();
  this->spatial_dimension = (spatial_dimension == 0) ? mesh.getSpatialDimension() : spatial_dimension;
  std::stringstream sstr; sstr << id << ":fem";
  this->fem = new FEM(mesh, spatial_dimension, sstr.str(), memory_id);
  this->fem_boundary = NULL;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline Model::~Model() {
  AKANTU_DEBUG_IN();
  delete fem;

  if(fem_boundary) delete fem_boundary;

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
inline FEM & Model::getFEMBoundary(){
  AKANTU_DEBUG_IN();

  if (!fem_boundary){
    MeshUtils::buildFacets(fem->getMesh());
    std::stringstream sstr; sstr << id << ":femboundary";
    this->fem_boundary = new FEM(fem->getMesh(), spatial_dimension-1, sstr.str(), memory_id);
  }

  AKANTU_DEBUG_OUT();
  return *this->fem_boundary;
}
