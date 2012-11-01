/**
 * @file   mesh_partition.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Aug 17 16:19:43 2010
 *
 * @brief  implementation of common part of all partitioner
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
#include "mesh_partition.hh"
#include "mesh_utils.hh"
#include "aka_types.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MeshPartition::MeshPartition(const Mesh & mesh, UInt spatial_dimension,
			     const MemoryID & memory_id) :
  Memory(memory_id), id("MeshPartitioner"),
  mesh(mesh), spatial_dimension(spatial_dimension),
  partitions             ("partition"             , id),
  ghost_partitions       ("ghost_partition"       , id),
  ghost_partitions_offset("ghost_partition_offset", id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MeshPartition::~MeshPartition() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * conversion in c++ of a the  GENDUALMETIS (mesh.c) function wrote by George in
 * Metis (University of Minnesota)
 */
void MeshPartition::buildDualGraph(Vector<Int> & dxadj, Vector<Int> & dadjncy) {
  AKANTU_DEBUG_IN();

  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  UInt nb_types = type_list.size();
  UInt nb_good_types = 0;

  UInt nb_nodes_per_element[nb_types];
  UInt nb_nodes_per_element_p1[nb_types];

  UInt magic_number[nb_types];

  UInt * conn_val[nb_types];
  UInt nb_element[nb_types];

  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    if(Mesh::getSpatialDimension(type) != mesh.getSpatialDimension()) continue;

    ElementType type_p1 = Mesh::getP1ElementType(type);

    nb_nodes_per_element[nb_good_types]    = Mesh::getNbNodesPerElement(type);
    nb_nodes_per_element_p1[nb_good_types] = Mesh::getNbNodesPerElement(type_p1);
    conn_val[nb_good_types]                = mesh.getConnectivity(type, _not_ghost).values;
    nb_element[nb_good_types]              = mesh.getConnectivity(type, _not_ghost).getSize();
    magic_number[nb_good_types]            =
      Mesh::getNbNodesPerElement(Mesh::getFacetElementType(type_p1));
    nb_good_types++;
  }


  // Vector<UInt> node_offset;
  // Vector<UInt> node_index;

  CSR<UInt> node_to_elem;

  //  MeshUtils::buildNode2Elements(mesh, node_offset, node_index);
  MeshUtils::buildNode2Elements(mesh, node_to_elem);

  // UInt * node_offset_val = node_offset.values;
  // UInt * node_index_val = node_index.values;


  UInt nb_total_element = 0;
  UInt nb_total_node_element = 0;
  for (UInt t = 0; t < nb_good_types; ++t) {
    nb_total_element += nb_element[t];
    nb_total_node_element += nb_element[t]*nb_nodes_per_element_p1[t];
  }

  dxadj.resize(nb_total_element + 1);

  Int * dxadj_val = dxadj.values;
  /// initialize the dxadj array
  for (UInt t = 0, linerized_el = 0; t < nb_good_types; ++t)
    for (UInt el = 0; el < nb_element[t]; ++el, ++linerized_el)
      dxadj_val[linerized_el] = nb_nodes_per_element_p1[t];


  /// convert the dxadj_val array in a csr one
  for (UInt i = 1; i < nb_total_element; ++i) dxadj_val[i] += dxadj_val[i-1];
  for (UInt i = nb_total_element; i > 0; --i) dxadj_val[i]  = dxadj_val[i-1];
  dxadj_val[0] = 0;

  dadjncy.resize(2*dxadj_val[nb_total_element]);
  Int * dadjncy_val = dadjncy.values;

  /// weight map to determine adjacency
  unordered_map<UInt, UInt>::type weight_map;

  // UInt index[200], weight[200];   /// key, value
  // UInt mask = (1 << 11) - 1;      /// hash function
  // Int * mark = new Int[mask + 1]; /// collision detector
  // for (UInt i = 0; i < mask + 1; ++i) mark[i] = -1;


  for (UInt t = 0, linerized_el = 0; t < nb_good_types; ++t) {
    for (UInt el = 0; el < nb_element[t]; ++el, ++linerized_el) {
      UInt el_offset = el*nb_nodes_per_element[t];

      /// fill the weight map
      //      UInt m = 0;
      for (UInt n = 0; n < nb_nodes_per_element_p1[t]; ++n) {
  	UInt node = conn_val[t][el_offset + n];
	CSR<UInt>::iterator k;
	for (k = node_to_elem.rbegin(node); k != node_to_elem.rend(node); --k) {
  	// for (UInt k = node_offset_val[node + 1] - 1;
  	//      k >= node_offset_val[node];
  	//      --k) {
  	  UInt current_el = *k;//node_index_val[k];
  	  if(current_el <= linerized_el) break;

	  unordered_map<UInt, UInt>::type::iterator it_w;
	  it_w = weight_map.find(current_el);

	  if(it_w == weight_map.end()) {
	    weight_map[current_el] = 1;
	  } else {
	    it_w->second++;
	  }

  	  // UInt mark_offset  = current_el & mask;
  	  // Int current_mark = mark[mark_offset];
  	  // if(current_mark == -1) { /// if element not in map
  	  //   index[m] = current_el;
  	  //   weight[m] = 1;
  	  //   mark[mark_offset] = m++;
  	  // } else if (index[current_mark] == current_el) { /// if element already in map
  	  //   weight[current_mark]++;
  	  // } else { /// if element already in map but collision in the keys
  	  //   UInt i;
  	  //   for (i = 0; i < m; ++i) {
  	  //     if(index[i] == current_el) {
  	  // 	weight[i]++;
  	  // 	break;
  	  //     }
  	  //   }
  	  //   if(i == m) {
  	  //     index[m] = current_el;
  	  //     weight[m++] = 1;
  	  //   }
  	  // }
  	}
      }
      /// each element with a weight of the size of a facet are adjacent
      unordered_map<UInt, UInt>::type::iterator it_w;
      for(it_w = weight_map.begin(); it_w != weight_map.end(); ++it_w) {
	if(it_w->second == magic_number[t]) {
  	  UInt adjacent_el = it_w->first;

	  UInt index_adj = dxadj_val[adjacent_el ]++;
	  UInt index_lin = dxadj_val[linerized_el]++;

  	  dadjncy_val[index_lin] = adjacent_el;
  	  dadjncy_val[index_adj] = linerized_el;
	}
      }

      // for (UInt n = 0; n < m; ++n) {
      // 	if(weight[n] == magic_number[t]) {
      // 	  UInt adjacent_el = index[n];
      // 	  dadjncy_val[dxadj_val[linerized_el]++] = adjacent_el;
      // 	  dadjncy_val[dxadj_val[adjacent_el ]++] = linerized_el;
      // 	}
      // 	mark[index[n] & mask] = -1;
      // }

      weight_map.clear();
    }
  }

  Int k_start = 0;
  for (UInt t = 0, linerized_el = 0, j = 0; t < nb_good_types; ++t)
    for (UInt el = 0; el < nb_element[t]; ++el, ++linerized_el) {
      for (Int k = k_start; k < dxadj_val[linerized_el]; ++k, ++j)
	dadjncy_val[j] = dadjncy_val[k];
      dxadj_val[linerized_el] = j;
      k_start += nb_nodes_per_element_p1[t];
    }

  for (UInt i = nb_total_element; i > 0; --i) dxadj_val[i] = dxadj_val[i-1];
  dxadj_val[0] = 0;

  //  delete [] mark;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshPartition::fillPartitionInformation(const Mesh & mesh,
					     const Int * linearized_partitions) {
  AKANTU_DEBUG_IN();

  // Vector<UInt> node_offset;
  // Vector<UInt> node_index;
  CSR<UInt> node_to_elem;

  // MeshUtils::buildNode2Elements(mesh, node_offset, node_index);
  MeshUtils::buildNode2Elements(mesh, node_to_elem);
  // UInt * node_offset_val = node_offset.values;
  // UInt * node_index_val = node_index.values;

  Mesh::type_iterator it  = mesh.firstType(spatial_dimension);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension);

  UInt linearized_el = 0;
  for(; it != end; ++it) {
    ElementType type = *it;

    UInt nb_element = mesh.getNbElement(*it);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    partitions             .alloc(nb_element,     1, type, _not_ghost);
    ghost_partitions_offset.alloc(nb_element + 1, 1, type, _ghost);
    ghost_partitions       .alloc(0,              1, type, _ghost);

    const Vector<UInt> & connectivity = mesh.getConnectivity(type, _not_ghost);

    for (UInt el = 0; el < nb_element; ++el, ++linearized_el) {
      UInt part = linearized_partitions[linearized_el];

      partitions(type, _not_ghost)(el) = part;
      std::list<UInt> list_adj_part;
      for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	UInt node = connectivity.values[el * nb_nodes_per_element + n];
	CSR<UInt>::iterator ne;
	//	for (UInt ne = node_offset_val[node]; ne < node_offset_val[node + 1]; ++ne) {
	for (ne = node_to_elem.begin(node); ne != node_to_elem.end(node); ++ne) {
	  UInt adj_el = *ne;//node_index_val[ne];
	  UInt adj_part = linearized_partitions[adj_el];
	  if(part != adj_part) {
	    list_adj_part.push_back(adj_part);
	  }
	}
      }

      list_adj_part.sort();
      list_adj_part.unique();

      for(std::list<UInt>::iterator adj_it = list_adj_part.begin();
	  adj_it != list_adj_part.end();
	  ++adj_it) {
	ghost_partitions(type, _ghost).push_back(*adj_it);
	ghost_partitions_offset(type, _ghost)(el)++;
      }
    }

    /// convert the ghost_partitions_offset array in an offset array
    Vector<UInt> & ghost_partitions_offset_ptr = ghost_partitions_offset(type, _ghost);
    for (UInt i = 1; i < nb_element; ++i)
      ghost_partitions_offset_ptr(i) += ghost_partitions_offset_ptr(i-1);
    for (UInt i = nb_element; i > 0; --i)
      ghost_partitions_offset_ptr(i)  = ghost_partitions_offset_ptr(i-1);
    ghost_partitions_offset_ptr(0) = 0;
  }

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
