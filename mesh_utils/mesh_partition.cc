/**
 * @file   mesh_partition.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Aug 16 17:16:59 2010
 *
 * @brief  implementation of common part of all partitioner
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "mesh_partition.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MeshPartition::MeshPartition(const Mesh & mesh, UInt spatial_dimension,
			     const MemoryID & memory_id) :
  Memory(memory_id), mesh(mesh), spatial_dimension(spatial_dimension),
  dxadj(NULL), dadjncy(NULL) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MeshPartition::~MeshPartition() {
  if(dxadj) delete dxadj;
  if(dadjncy) delete dadjncy;
}

/* -------------------------------------------------------------------------- */
void MeshPartition::buildDualGraph() {
  AKANTU_DEBUG_IN();

  // UInt nb_nodes = mesh.getNbNodes();

  // const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  // Mesh::ConnectivityTypeList::const_iterator it;

  // UInt nb_types = type_list.size();
  // UInt nb_good_types = 0;

  // UInt nb_nodes_per_element[nb_types];
  // UInt nb_nodes_per_element_p1[nb_types];

  // UInt magic_number[nb_types];

  // UInt * conn_val[nb_types];
  // UInt nb_element[nb_types];

  // for(it = type_list.begin(); it != type_list.end(); ++it) {
  //   ElementType type = *it;
  //   if(mesh.getSpatialDimension(type) != spatial_dimension) continue;

  //   nb_nodes_per_element[nb_good_types]    = mesh.getNbNodesPerElement(type);
  //   nb_nodes_per_element_p1[nb_good_types] = mesh.getNbNodesPerElementP1(type);

  //   magic_number[nb_good_types] = mesh.getNbNodesPerElementP1(mesh.getSurfaceElementType(type));

  //   conn_val[nb_good_types] = mesh.getConnectivity(type).values;
  //   nb_element[nb_good_types] = mesh.getConnectivity(type).getSize();
  //   nb_good_types++;
  // }

  // /// array for the node-element list
  // UInt * node_offset = new UInt[nb_nodes + 1];
  // UInt * node_index;

  // /// count number of occurrence of each node
  // for (UInt t = 0; t < nb_good_types; ++t) {
  //   memset(node_offset, 0, (nb_nodes + 1)*sizeof(UInt));
  //   for (UInt el = 0; el < nb_element[t]; ++el) {
  //     UInt el_offset = el*nb_nodes_per_element[t];
  //     for (UInt n = 0; n < nb_nodes_per_element_p1[t]; ++n) {
  // 	node_offset[conn_val[t][el_offset + n]]++;
  //     }
  //   }
  // }

  // /// convert the occurrence array in a csr one
  // for (UInt i = 1; i < nb_nodes; ++i) node_offset[i] += node_offset[i-1];
  // for (UInt i = nb_nodes; i > 0; --i) node_offset[i]  = node_offset[i-1];
  // node_offset[0] = 0;

  // /// rearrange element to get the node-element list
  // node_index  = new UInt[node_offset[nb_nodes]];
  // for (UInt t = 0, linerized_el = 0; t < nb_good_types; ++t)
  //   for (UInt el = 0; el < nb_element[t]; ++el, ++linerized_el) {
  //     UInt el_offset = el*nb_nodes_per_element[t];
  //     for (UInt n = 0; n < nb_nodes_per_element_p1[t]; ++n)
  // 	node_index[node_offset[conn_val[t][el_offset + n]]++] = linerized_el;
  //   }

  // for (UInt i = nb_nodes; i > 0; --i) node_offset[i]  = node_offset[i-1];
  // node_offset[0] = 0;

  /* ****************************************** */

  // UInt nb_total_element = 0;
  // UInt nb_total_node_element = 0;
  // for (UInt t = 0; t < nb_good_types; ++t) {
  //   nb_total_element += nb_element[t];
  //   nb_total_node_element += nb_element[t]*nb_nodes_per_element_p1[t];
  // }

  // std::stringstream sstr_dxadj; sstr_dxadj << mesh.getID() << ":dxadj";
  // dxadj = new Vector<UInt>(nb_total_element, 1, sstr_dxadj.str());

  // std::stringstream sstr_dadjncy; sstr_dadjncy << mesh.getID() << ":dadjncy";
  // dadjncy = new Vector<UInt>(nb_total_node_element, 1, sstr_dxadj.str());

  // UInt * dxadj_val = dxadj->values;
  // UInt * dadjncy_val = dadjncy->values;

  // /// initialize the dxadj array
  // for (UInt t = 0, linerized_el = 0; t < nb_good_types; ++t)
  //   for (UInt el = 0; el < nb_element[t]; ++el, ++linerized_el)
  //     dxadj_val[linerized_el] = nb_nodes_per_element_p1[t];



  // /// convert the dxadj_val array in a csr one
  // for (UInt i = 1; i < nb_total_element; ++i) dxadj_val[i] += dxadj_val[i-1];
  // for (UInt i = nb_total_element; i > 0; --i) dxadj_val[i]  = dxadj_val[i-1];
  // dxadj_val[0] = 0;

  // /// weight map to determine adjacency
  // UInt index[200], weight[200];   /// key, value
  // UInt mask = (1 << 11) - 1;      /// hash function
  // Int * mark = new Int[mask + 1]; /// collision detector
  // for (UInt i = 0; i < mask + 1; ++i) mark[i] = -1;


  // for (UInt t = 0, linerized_el = 0; t < nb_good_types; ++t) {
  //   for (UInt el = 0; el < nb_element[t]; ++el, ++linerized_el) {
  //     UInt el_offset = el*nb_nodes_per_element[t];

  //     /// fill the weight map
  //     UInt m = 0;
  //     for (UInt n = 0; n < nb_nodes_per_element_p1[t]; ++n) {
  // 	UInt node = conn_val[t][el_offset + n];

  // 	for (UInt k = node_offset[node + 1] - 1;
  // 	     k >= node_offset[node];
  // 	     --k) {
  // 	  UInt current_el = node_index[k];
  // 	  if(current_el <= linerized_el) break;

  // 	  UInt mark_offset  = current_el & mask;
  // 	  Int current_mark = mark[mark_offset];
  // 	  if(current_mark == -1) { /// if element not in map
  // 	    index[m] = current_el;
  // 	    weight[m] = 1;
  // 	    mark[mark_offset] = m++;
  // 	  } else if (index[current_mark] == current_el) { /// if element already in map
  // 	    weight[current_mark]++;
  // 	  } else { /// if element already in map but collision in the keys
  // 	    UInt i;
  // 	    for (i = 0; i < m; ++i) {
  // 	      if(index[i] == current_el) {
  // 		weight[i]++;
  // 		break;
  // 	      }
  // 	    }
  // 	    if(i == m) {
  // 	      index[m] = current_el;
  // 	      weight[m++] = 1;
  // 	    }
  // 	  }
  // 	}
  //     }

  //     /// each element with a weight of the size of a facet are adjacent
  //     for (UInt n = 0; n < m; ++n) {
  // 	if(weight[n] == magic_number[t]) {
  // 	  UInt adjacent_el = index[n];
  // 	  dadjncy_val[dxadj_val[linerized_el]++] = adjacent_el;
  // 	  dadjncy_val[dxadj_val[adjacent_el ]++] = linerized_el;
  // 	}
  // 	mark[index[n] & mask] = -1;
  //     }
  //   }
  // }

  // UInt k_start = 0;
  // for (UInt t = 0, linerized_el = 0, j = 0; t < nb_good_types; ++t)
  //   for (UInt el = 0; el < nb_element[t]; ++el, ++linerized_el) {
  //     for (UInt k = k_start; k < dxadj_val[linerized_el]; ++k, ++j)
  // 	dadjncy_val[j] = dadjncy_val[k];
  //     dxadj_val[linerized_el] = j;
  //     k_start += nb_nodes_per_element_p1[t];
  //   }

  // for (UInt i = nb_total_element; i > 0; ++i) dxadj_val[i] = dxadj_val[i-1];
  // dxadj_val[0] = 0;

  // delete [] node_index;
  // delete [] node_offset;
  // delete [] mark;


  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
