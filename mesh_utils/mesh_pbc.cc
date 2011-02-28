/**
 * @file   mesh_pbc.cc
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @date   Tue Feb  8 10:48:01 2011
 *
 * @brief  periodic boundary condition connectivity tweak
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

#include "mesh_utils.hh"
#include <map>

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/// class that sorts a set of nodes of same coordinates in 'dir' direction
class CoordinatesComparison {
public:
  CoordinatesComparison (const UInt dimension, 
			 const UInt dirx, const UInt diry, 
			 Real * coords):
    dim(dimension),dir_x(dirx),dir_y(diry),coordinates(coords){}
  
  bool operator() (UInt n1, UInt n2){
    Real p1_x = coordinates[dim*n1+dir_x];
    Real p2_x = coordinates[dim*n2+dir_x];
    Real p1_y = coordinates[dim*n1+dir_y];
    Real p2_y = coordinates[dim*n2+dir_y];
    Real diff_x = p1_x - p2_x;
    if (fabs(diff_x) > Math::tolerance)
      return diff_x > 0.0 ? false : true;
    else {
      Real diff_y = p1_y - p2_y;
      return diff_y > 0 ? false : true;
    }
  }
private:
  UInt dim;
  UInt dir_x;
  UInt dir_y;
  Real * coordinates;
};

/* -------------------------------------------------------------------------- */
void MeshUtils::tweakConnectivityForPBC(Mesh & mesh,
					UInt flag_x,UInt flag_y,UInt flag_z){
  std::map<UInt,UInt> pbc_pair;
  mesh.computeBoundingBox();

  if (flag_x) computePBCMap(mesh,0,pbc_pair);
  if (flag_y) computePBCMap(mesh,1,pbc_pair);
  if (flag_z) computePBCMap(mesh,2,pbc_pair);

  {
    std::map<UInt,UInt>::iterator it = pbc_pair.begin();
    std::map<UInt,UInt>::iterator end = pbc_pair.end();
    
    Real * coords = mesh.nodes->values;
    UInt dim = mesh.getSpatialDimension();
    while(it != end){
      UInt i1 = (*it).first;
      UInt i2 = (*it).second;
      
      AKANTU_DEBUG_INFO("pairing " << i1 << "(" 
			<< coords[dim*i1] << "," << coords[dim*i1+1] << "," << coords[dim*i1+2]
			<< ") with"
			<< i2 << "(" 
			<< coords[dim*i2] << "," << coords[dim*i2+1] << "," << coords[dim*i2+2]
			<< ")");	
      ++it;
    }
  }

  //allocate and initialize list of reversed elements
  mesh.initByElementTypeUIntVector(mesh.reversed_elements_pbc,1,0,mesh.id,"reversed");
  // now loop over the elements to change the connectivity of some elements
  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    UInt nb_elem = mesh.getNbElement(type);
    UInt nb_nodes_per_elem = mesh.getNbNodesPerElement(type);
    UInt * conn = mesh.getConnectivityPointer(*it)->values;
    UInt index = 0;
    Vector<UInt> & list = *(mesh.reversed_elements_pbc[type]);
    for (UInt el = 0; el < nb_elem; el++) {
      for (UInt k = 0; k < nb_nodes_per_elem; ++k,++index){
	if (pbc_pair.count(conn[index])){
	  conn[index] = pbc_pair[conn[index]];
	}	
      }
      list.push_back(el);
    }
  }
}

/* -------------------------------------------------------------------------- */

void MeshUtils::computePBCMap(const Mesh & mymesh,const UInt dir,
			      std::map<UInt,UInt> & pbc_pair){
  std::vector<UInt> selected_left;
  std::vector<UInt> selected_right;

  Real * coords = mymesh.nodes->values;
  const UInt nb_nodes = mymesh.nodes->getSize();
  const UInt dim = mymesh.getSpatialDimension();
  AKANTU_DEBUG_INFO("min " << mymesh.xmin[dir]);
  AKANTU_DEBUG_INFO("max " << mymesh.xmax[dir]);

  for (UInt i = 0; i < nb_nodes; ++i) {
    AKANTU_DEBUG_TRACE("treating " << coords[dim*i+dir]);
    if(Math::are_float_equal(coords[dim*i+dir],mymesh.xmin[dir])){
      AKANTU_DEBUG_TRACE("pushing node " << i << " on the left side");
      selected_left.push_back(i);
    }
    else if(Math::are_float_equal(coords[dim*i+dir],mymesh.xmax[dir])){
      selected_right.push_back(i);
      AKANTU_DEBUG_TRACE("pushing node " << i << " on the right side");
    }
  }

  AKANTU_DEBUG_INFO("found " << selected_left.size() << " and " << selected_right.size() 
	       << " nodes at each boundary for direction " << dir);

  UInt dir_x,dir_y;

  if (dir == 0){
    dir_x = 1;dir_y = 2;
  }
  else if (dir == 1){
    dir_x = 0;dir_y = 2;
  }
  else if (dir == 2){
    dir_x = 0;dir_y = 1;
  }

  CoordinatesComparison compare_nodes(dim,dir_x,dir_y,coords);

  std::sort(selected_left.begin(),selected_left.end(),compare_nodes);
  std::sort(selected_right.begin(),selected_right.end(),compare_nodes);
  
  std::vector<UInt>::iterator it_left = selected_left.begin();
  std::vector<UInt>::iterator end_left = selected_left.end();

  std::vector<UInt>::iterator it_right = selected_right.begin();
  std::vector<UInt>::iterator end_right = selected_right.end();

  while ((it_left != end_left) && (it_right != end_right) ){
    UInt i1 = *it_left;
    UInt i2 = *it_right;

    AKANTU_DEBUG_TRACE("do I pair? " << i1 << "(" 
		      << coords[dim*i1] << "," << coords[dim*i1+1] << "," << coords[dim*i1+2]
		      << ") with"
		      << i2 << "(" 
		      << coords[dim*i2] << "," << coords[dim*i2+1] << "," << coords[dim*i2+2]
		      << ") in direction " << dir);	

    
    Real dx = coords[dim*i1+dir_x] - coords[dim*i2+dir_x];
    Real dy = coords[dim*i1+dir_y] - coords[dim*i2+dir_y];

    if (fabs(dx*dx+dy*dy) < Math::tolerance)
      {
  	//then i match these pairs
	if (pbc_pair.count(i2)){
	  i2 = pbc_pair[i2];
	}
  	pbc_pair[i1] = i2;

	AKANTU_DEBUG_TRACE("pairing " << i1 << "(" 
			  << coords[dim*i1] << "," << coords[dim*i1+1] << "," << coords[dim*i1+2]
			  << ") with"
			  << i2 << "(" 
			  << coords[dim*i2] << "," << coords[dim*i2+1] << "," << coords[dim*i2+2]
			  << ") in direction " << dir);	
  	++it_left;
  	++it_right;
      }
    else if (fabs(dy) < Math::tolerance && dx > 0) ++it_right;
    else if (fabs(dy) < Math::tolerance && dx < 0) ++it_left;
    else if (dy > 0) ++it_right;
    else if (dy < 0) ++it_left;
    else {
      AKANTU_DEBUG_ERROR("this should not append");	
    }
  }
  AKANTU_DEBUG_INFO("found " <<  pbc_pair.size() << " pairs for direction " << dir);

}

__END_AKANTU__
