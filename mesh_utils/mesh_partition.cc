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
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MeshPartition::MeshPartition(const Mesh & mesh, UInt spatial_dimension,
			     const MemoryID & memory_id) :
  Memory(memory_id), id("MeshPartitioner"),
  mesh(mesh), spatial_dimension(spatial_dimension) {
  AKANTU_DEBUG_IN();

  for(UInt type = _not_defined; type < _max_element_type; ++type) {
    partitions[type]              = NULL;
    ghost_partitions[type]        = NULL;
    ghost_partitions_offset[type] = NULL;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MeshPartition::~MeshPartition() {
  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    if(Mesh::getSpatialDimension(type) != mesh.getSpatialDimension()) continue;

    if(partitions[type]) {
      dealloc(partitions[type]->getID());
      dealloc(ghost_partitions[type]->getID());
      dealloc(ghost_partitions_offset[type]->getID());
    }
  }
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
    conn_val[nb_good_types]                = mesh.getConnectivity(type).values;
    nb_element[nb_good_types]              = mesh.getConnectivity(type).getSize();
    magic_number[nb_good_types]            =
      Mesh::getNbNodesPerElement(Mesh::getFacetElementType(type_p1));
    nb_good_types++;
  }


  Vector<UInt> node_offset;
  Vector<UInt> node_index;

  MeshUtils::buildNode2Elements(mesh, node_offset, node_index);

  UInt * node_offset_val = node_offset.values;
  UInt * node_index_val = node_index.values;


  UInt nb_total_element = 0;
  UInt nb_total_node_element = 0;
  for (UInt t = 0; t < nb_good_types; ++t) {
    nb_total_element += nb_element[t];
    nb_total_node_element += nb_element[t]*nb_nodes_per_element_p1[t];
  }

  dxadj.resize(nb_total_element + 1);
  dadjncy.resize(nb_total_node_element);

  Int * dxadj_val = dxadj.values;
  Int * dadjncy_val = dadjncy.values;

  /// initialize the dxadj array
  for (UInt t = 0, linerized_el = 0; t < nb_good_types; ++t)
    for (UInt el = 0; el < nb_element[t]; ++el, ++linerized_el)
      dxadj_val[linerized_el] = nb_nodes_per_element_p1[t];


  /// convert the dxadj_val array in a csr one
  for (UInt i = 1; i < nb_total_element; ++i) dxadj_val[i] += dxadj_val[i-1];
  for (UInt i = nb_total_element; i > 0; --i) dxadj_val[i]  = dxadj_val[i-1];
  dxadj_val[0] = 0;

  /// weight map to determine adjacency
  UInt index[200], weight[200];   /// key, value
  UInt mask = (1 << 11) - 1;      /// hash function
  Int * mark = new Int[mask + 1]; /// collision detector
  for (UInt i = 0; i < mask + 1; ++i) mark[i] = -1;


  for (UInt t = 0, linerized_el = 0; t < nb_good_types; ++t) {
    for (UInt el = 0; el < nb_element[t]; ++el, ++linerized_el) {
      UInt el_offset = el*nb_nodes_per_element[t];

      /// fill the weight map
      UInt m = 0;
      for (UInt n = 0; n < nb_nodes_per_element_p1[t]; ++n) {
  	UInt node = conn_val[t][el_offset + n];

  	for (UInt k = node_offset_val[node + 1] - 1;
  	     k >= node_offset_val[node];
  	     --k) {
  	  UInt current_el = node_index_val[k];
  	  if(current_el <= linerized_el) break;

  	  UInt mark_offset  = current_el & mask;
  	  Int current_mark = mark[mark_offset];
  	  if(current_mark == -1) { /// if element not in map
  	    index[m] = current_el;
  	    weight[m] = 1;
  	    mark[mark_offset] = m++;
  	  } else if (index[current_mark] == current_el) { /// if element already in map
  	    weight[current_mark]++;
  	  } else { /// if element already in map but collision in the keys
  	    UInt i;
  	    for (i = 0; i < m; ++i) {
  	      if(index[i] == current_el) {
  		weight[i]++;
  		break;
  	      }
  	    }
  	    if(i == m) {
  	      index[m] = current_el;
  	      weight[m++] = 1;
  	    }
  	  }
  	}
      }

      /// each element with a weight of the size of a facet are adjacent
      for (UInt n = 0; n < m; ++n) {
  	if(weight[n] == magic_number[t]) {
  	  UInt adjacent_el = index[n];
  	  dadjncy_val[dxadj_val[linerized_el]++] = adjacent_el;
  	  dadjncy_val[dxadj_val[adjacent_el ]++] = linerized_el;
  	}
  	mark[index[n] & mask] = -1;
      }
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

  delete [] mark;

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
