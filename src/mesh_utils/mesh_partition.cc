/**
 * @file   mesh_partition.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Aug 16 17:16:59 2010
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
  partitions             ("partition"             , id, memory_id),
  ghost_partitions       ("ghost_partition"       , id, memory_id),
  ghost_partitions_offset("ghost_partition_offset", id, memory_id) {
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
 * conversion in c++ of the GENDUALMETIS (mesh.c) function wrote by George in
 * Metis (University of Minnesota)
 */
void MeshPartition::buildDualGraph(Vector<Int> & dxadj, Vector<Int> & dadjncy,
				   Vector<Int> & edge_loads,
				   const EdgeLoadFunctor & edge_load_func,
				   const Vector<UInt> & pairs) {
  AKANTU_DEBUG_IN();

  // tweak mesh;
  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  UInt nb_types = type_list.size();
  UInt nb_good_types = 0;

  UInt nb_nodes_per_element[nb_types];
  UInt nb_nodes_per_element_p1[nb_types];

  UInt magic_number[nb_types];

  //  UInt * conn_val[nb_types];
  UInt nb_element[nb_types];

  Vector<UInt> * conn[nb_types];
  Vector<UInt> * conn_tmp[nb_types];

  Vector<Element> lin_to_element;

  Element el;
  el.ghost_type = _not_ghost;

  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    if(Mesh::getSpatialDimension(type) != mesh.getSpatialDimension()) continue;
    el.type = type;

    ElementType type_p1 = Mesh::getP1ElementType(type);

    nb_nodes_per_element[nb_good_types]    = Mesh::getNbNodesPerElement(type);
    nb_nodes_per_element_p1[nb_good_types] = Mesh::getNbNodesPerElement(type_p1);
    nb_element[nb_good_types]              = mesh.getConnectivity(type, _not_ghost).getSize();
    magic_number[nb_good_types]            =
      Mesh::getNbNodesPerElement(Mesh::getFacetElementType(type_p1));

    conn[nb_good_types] = &const_cast<Vector<UInt> &>(mesh.getConnectivity(type, _not_ghost));

    for (UInt i = 0; i < nb_element[nb_good_types]; ++i) {
      el.element = i;
      lin_to_element.push_back(el);
    }


    if(pairs.getSize() != 0) {
      conn_tmp[nb_good_types] = new Vector<UInt>(mesh.getConnectivity(type, _not_ghost));
      for (UInt i = 0; i < pairs.getSize(); ++i) {
	for (UInt el = 0; el < nb_element[nb_good_types]; ++el) {
	  for (UInt n = 0; n < nb_nodes_per_element[nb_good_types]; ++n) {
	    if(pairs(i, 1) == (*conn[nb_good_types])(el, n))
	      (*conn[nb_good_types])(el, n) = pairs(i, 0);
	  }
	}
      }
    }

    nb_good_types++;
  }

  CSR<UInt> node_to_elem;

  MeshUtils::buildNode2Elements(mesh, node_to_elem);

  UInt nb_total_element = 0;
  UInt nb_total_node_element = 0;
  for (UInt t = 0; t < nb_good_types; ++t) {
    nb_total_element += nb_element[t];
    nb_total_node_element += nb_element[t]*nb_nodes_per_element_p1[t];
  }

  dxadj.resize(nb_total_element + 1);

  /// initialize the dxadj array
  for (UInt t = 0, linerized_el = 0; t < nb_good_types; ++t)
    for (UInt el = 0; el < nb_element[t]; ++el, ++linerized_el)
      dxadj(linerized_el) = nb_nodes_per_element_p1[t];

  /// convert the dxadj_val array in a csr one
  for (UInt i = 1; i < nb_total_element; ++i) dxadj(i) += dxadj(i-1);
  for (UInt i = nb_total_element; i > 0; --i) dxadj(i)  = dxadj(i-1);
  dxadj(0) = 0;

  dadjncy.resize(2*dxadj(nb_total_element));
  edge_loads.resize(2*dxadj(nb_total_element));

  /// weight map to determine adjacency
  unordered_map<UInt, UInt>::type weight_map;

  for (UInt t = 0, linerized_el = 0; t < nb_good_types; ++t) {
    for (UInt el = 0; el < nb_element[t]; ++el, ++linerized_el) {
      /// fill the weight map
      for (UInt n = 0; n < nb_nodes_per_element_p1[t]; ++n) {
  	UInt node = (*conn[t])(el, n);
	CSR<UInt>::iterator k;
	for (k = node_to_elem.rbegin(node); k != node_to_elem.rend(node); --k) {
  	  UInt current_el = *k;
  	  if(current_el <= linerized_el) break;

	  unordered_map<UInt, UInt>::type::iterator it_w;
	  it_w = weight_map.find(current_el);

	  if(it_w == weight_map.end()) {
	    weight_map[current_el] = 1;
	  } else {
	    it_w->second++;
	  }
  	}
      }
      /// each element with a weight of the size of a facet are adjacent
      unordered_map<UInt, UInt>::type::iterator it_w;
      for(it_w = weight_map.begin(); it_w != weight_map.end(); ++it_w) {
	if(it_w->second == magic_number[t]) {
  	  UInt adjacent_el = it_w->first;

	  UInt index_adj = dxadj(adjacent_el )++;
	  UInt index_lin = dxadj(linerized_el)++;

  	  dadjncy(index_lin) = adjacent_el;
  	  dadjncy(index_adj) = linerized_el;
	}
      }

      weight_map.clear();
    }
  }

  if(pairs.getSize() != 0) {
    for (UInt i = 0; i < nb_good_types; ++i) {
      conn[i]->copy(*conn_tmp[i]);
      delete conn_tmp[i];
    }
  }

  Int k_start = 0;
  for (UInt t = 0, linerized_el = 0, j = 0; t < nb_good_types; ++t)
    for (UInt el = 0; el < nb_element[t]; ++el, ++linerized_el) {
      for (Int k = k_start; k < dxadj(linerized_el); ++k, ++j)
	dadjncy(j) = dadjncy(k);
      dxadj(linerized_el) = j;
      k_start += nb_nodes_per_element_p1[t];
    }

  for (UInt i = nb_total_element; i > 0; --i) dxadj(i) = dxadj(i - 1);
  dxadj(0) = 0;

  UInt adj = 0;
  for (UInt i = 0; i < nb_total_element; ++i) {
    UInt nb_adj = dxadj(i + 1) - dxadj(i);
    for (UInt j = 0; j < nb_adj; ++j, ++adj) {
      Int el_adj_id = dadjncy(dxadj(i) + j);
      Element el     = lin_to_element(i);
      Element el_adj = lin_to_element(el_adj_id);

      Int load = edge_load_func(el, el_adj);
      edge_loads(adj) = load;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshPartition::fillPartitionInformation(const Mesh & mesh,
					     const Int * linearized_partitions) {
  AKANTU_DEBUG_IN();

  CSR<UInt> node_to_elem;

  MeshUtils::buildNode2Elements(mesh, node_to_elem);

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
	for (ne = node_to_elem.begin(node); ne != node_to_elem.end(node); ++ne) {
	  UInt adj_el = *ne;
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
