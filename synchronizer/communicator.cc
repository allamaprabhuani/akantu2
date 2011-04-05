/**
 * @file   communicator.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 26 16:36:21 2010
 *
 * @brief implementation of a  communicator using a static_communicator for real
 * send/receive
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
#include <map>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "ghost_synchronizer.hh"
#include "communicator.hh"
#include "static_communicator.hh"
#include "mesh_utils.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
Communicator::Communicator(SynchronizerID id,
			   MemoryID memory_id) :
  Synchronizer(id, memory_id),
  static_communicator(StaticCommunicator::getStaticCommunicator()) {
  AKANTU_DEBUG_IN();

  for (UInt t = _not_defined; t < _max_element_type; ++t) {
    element_to_send_offset   [t] = NULL;
    element_to_send	     [t] = NULL;
    element_to_receive_offset[t] = NULL;
    element_to_receive       [t] = NULL;
  }

  nb_proc = static_communicator->getNbProc();
  rank    = static_communicator->whoAmI();

  send_buffer = new Vector<Real>[nb_proc];
  recv_buffer = new Vector<Real>[nb_proc];

  send_element = new std::vector<Element>[nb_proc];
  recv_element = new std::vector<Element>[nb_proc];

  size_to_send.clear();
  size_to_receive.clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Communicator::~Communicator() {
  AKANTU_DEBUG_IN();

  for (UInt p = 0; p < nb_proc; ++p) {
    send_buffer[p].resize(0);
    recv_buffer[p].resize(0);

    send_element[p].clear();
    recv_element[p].clear();
  }

  delete [] send_buffer;
  delete [] recv_buffer;

  delete [] send_element;
  delete [] recv_element;


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Communicator * Communicator::createCommunicatorDistributeMesh(Mesh & mesh,
							      const MeshPartition * partition,
							      UInt root,
							      SynchronizerID id,
							      MemoryID memory_id) {
  AKANTU_DEBUG_IN();

  const UInt TAG_SIZES        = 0;
  const UInt TAG_CONNECTIVITY = 1;
  const UInt TAG_DATA         = 2;
  const UInt TAG_PARTITIONS   = 3;
  const UInt TAG_NB_NODES     = 4;
  const UInt TAG_NODES        = 5;
  const UInt TAG_COORDINATES  = 6;
  const UInt TAG_NODES_TYPE   = 7;

  StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator();
  UInt nb_proc = comm->getNbProc();
  UInt my_rank = comm->whoAmI();

  Communicator * communicator = new Communicator(id, memory_id);

  if(nb_proc == 1) return communicator;

  UInt * local_connectivity = NULL;
  UInt * local_partitions = NULL;
  UInt * local_data = NULL;
  Vector<UInt> * old_nodes = mesh.getNodesGlobalIdsPointer();
  Vector<Real> * nodes = mesh.getNodesPointer();

  UInt spatial_dimension = nodes->getNbComponent();

  /* ------------------------------------------------------------------------ */
  /*  Local (rank == root)                                                    */
  /* ------------------------------------------------------------------------ */
  if(my_rank == root) {
    AKANTU_DEBUG_ASSERT(partition->getNbPartition() == nb_proc,
			"The number of partition does not match the number of processors");

    /**
     * connectivity and communications scheme construction
     */
    const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
    Mesh::ConnectivityTypeList::const_iterator it;
    for(it = type_list.begin(); it != type_list.end(); ++it) {
      ElementType type = *it;

      if(Mesh::getSpatialDimension(type) != mesh.getSpatialDimension()) continue;

      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
      UInt nb_element = mesh.getNbElement(*it);
      UInt nb_local_element[nb_proc];
      UInt nb_ghost_element[nb_proc];
      UInt nb_element_to_send[nb_proc];

      memset(nb_local_element, 0, nb_proc*sizeof(UInt));
      memset(nb_ghost_element, 0, nb_proc*sizeof(UInt));
      memset(nb_element_to_send, 0, nb_proc*sizeof(UInt));

      UInt * partition_num = partition->getPartition(type).values;

      UInt * ghost_partition = partition->getGhostPartition(type).values;
      UInt * ghost_partition_offset = partition->getGhostPartitionOffset(type).values;

      Mesh::UIntDataMap & uint_data_map = mesh.getUIntDataMap(type);
      UInt nb_tags = uint_data_map.size();

      /* -------------------------------------------------------------------- */
      /// constructing the reordering structures
      for (UInt el = 0; el < nb_element; ++el) {
	nb_local_element[partition_num[el]]++;
	for (UInt part = ghost_partition_offset[el];
	     part < ghost_partition_offset[el + 1];
	     ++part) {
	  nb_ghost_element[ghost_partition[part]]++;
	}
	nb_element_to_send[partition_num[el]] +=
	  ghost_partition_offset[el + 1] - ghost_partition_offset[el] + 1;
      }

      /// allocating buffers
      UInt * buffers[nb_proc];
      UInt * buffers_tmp[nb_proc];
      for (UInt p = 0; p < nb_proc; ++p) {
	UInt size = nb_nodes_per_element * (nb_local_element[p] +
					    nb_ghost_element[p]);
	buffers[p] = new UInt[size];
	buffers_tmp[p] = buffers[p];
      }

      /// copying the local connectivity
      UInt * conn_val = mesh.getConnectivity(type).values;
      for (UInt el = 0; el < nb_element; ++el) {
	memcpy(buffers_tmp[partition_num[el]],
	       conn_val + el * nb_nodes_per_element,
	       nb_nodes_per_element * sizeof(UInt));
	buffers_tmp[partition_num[el]] += nb_nodes_per_element;
      }

      /// copying the connectivity of ghost element
      for (UInt el = 0; el < nb_element; ++el) {
	for (UInt part = ghost_partition_offset[el];
	     part < ghost_partition_offset[el + 1];
	     ++part) {
	  UInt proc = ghost_partition[part];
	  memcpy(buffers_tmp[proc],
		 conn_val + el * nb_nodes_per_element,
		 nb_nodes_per_element * sizeof(UInt));
	  buffers_tmp[proc] += nb_nodes_per_element;
	}
      }

      UInt names_size = 0;
      Mesh::UIntDataMap::iterator it_data;
      for(it_data = uint_data_map.begin(); it_data != uint_data_map.end(); ++it_data) {
	names_size += it_data->first.size() + 1;
      }

      /* -------->>>>-SIZE + CONNECTIVITY------------------------------------ */
       /// send all connectivity and ghost information to all processors
      std::vector<CommunicationRequest *> requests;
      for (UInt p = 0; p < nb_proc; ++p) {
	if(p != root) {
	  UInt size[6];
	  size[0] = (UInt) type;
	  size[1] = nb_local_element[p];
	  size[2] = nb_ghost_element[p];
	  size[3] = nb_element_to_send[p];
	  size[4] = nb_tags;
	  size[5] = names_size;
	  comm->send(size, 6, p, TAG_SIZES);

	  AKANTU_DEBUG_INFO("Sending connectivities to proc " << p);
	  requests.push_back(comm->asyncSend(buffers[p],
					    nb_nodes_per_element * (nb_local_element[p] +
								    nb_ghost_element[p]),
					    p, TAG_CONNECTIVITY));
	} else {
	  local_connectivity = buffers[p];
	}
      }

      /// create the renumbered connectivity
      AKANTU_DEBUG_INFO("Renumbering local connectivities");
      MeshUtils::renumberMeshNodes(mesh,
				   local_connectivity,
				   nb_local_element[root],
				   nb_ghost_element[root],
				   type,
				   *old_nodes);

      comm->waitAll(requests);
      comm->freeCommunicationRequest(requests);
      requests.clear();

      for (UInt p = 0; p < nb_proc; ++p) {
	delete [] buffers[p];
      }

      /* -------------------------------------------------------------------- */
      for (UInt p = 0; p < nb_proc; ++p) {
	buffers[p] = new UInt[nb_ghost_element[p] + nb_element_to_send[p]];
	buffers_tmp[p] = buffers[p];
      }

      /// splitting the partition information to send them to processors
      UInt count_by_proc[nb_proc];
      memset(count_by_proc, 0, nb_proc*sizeof(UInt));
      for (UInt el = 0; el < nb_element; ++el) {
	*(buffers_tmp[partition_num[el]]++) = ghost_partition_offset[el + 1] - ghost_partition_offset[el];
	for (UInt part = ghost_partition_offset[el], i = 0;
	     part < ghost_partition_offset[el + 1];
	     ++part, ++i) {
	  *(buffers_tmp[partition_num[el]]++) = ghost_partition[part];
	}
      }

      for (UInt el = 0; el < nb_element; ++el) {
	for (UInt part = ghost_partition_offset[el], i = 0;
	     part < ghost_partition_offset[el + 1];
	     ++part, ++i) {
	  *(buffers_tmp[ghost_partition[part]]++) = partition_num[el];
	}
      }

      /* -------->>>>-PARTITIONS--------------------------------------------- */
      /// last data to compute the communication scheme
      for (UInt p = 0; p < nb_proc; ++p) {
	if(p != root) {
	  AKANTU_DEBUG_INFO("Sending partition informations to proc " << p);
	  requests.push_back(comm->asyncSend(buffers[p],
					     nb_element_to_send[p] + nb_ghost_element[p],
					     p, TAG_PARTITIONS));
	} else {
	  local_partitions = buffers[p];
	}
      }

      AKANTU_DEBUG_INFO("Creating communications scheme");
      communicator->fillCommunicationScheme(local_partitions,
					    nb_local_element[root],
					    nb_ghost_element[root],
					    type);

      comm->waitAll(requests);
      comm->freeCommunicationRequest(requests);
      requests.clear();

      for (UInt p = 0; p < nb_proc; ++p) {
	delete [] buffers[p];
      }

      /* -------------------------------------------------------------------- */
      /// send int data assossiated to the mesh
      if(nb_tags) {
	UInt uint_names_size = names_size / sizeof(UInt) + (names_size % sizeof(UInt) ? 1 : 0);
	for (UInt p = 0; p < nb_proc; ++p) {
	  UInt size = nb_tags * (nb_local_element[p] + nb_ghost_element[p])
	    + uint_names_size;
	  buffers[p] = new UInt[size];
	  std::fill_n(buffers[p], size, 0);
	  buffers_tmp[p] = buffers[p];
	}

	char * names = new char[names_size];
	char * names_tmp = names;
	memset(names, 0, names_size);
	for(it_data = uint_data_map.begin(); it_data != uint_data_map.end(); ++it_data) {
	  UInt * data = it_data->second->values;

	  memcpy(names_tmp, it_data->first.data(), it_data->first.size());
	  names_tmp += it_data->first.size() + 1;

	  /// copying data for the local element
	  for (UInt el = 0; el < nb_element; ++el) {
	    UInt proc = partition_num[el];
	    *(buffers_tmp[proc]) = data[el];
	    buffers_tmp[proc]++;
	  }

	  /// copying the data for the ghost element
	  for (UInt el = 0; el < nb_element; ++el) {
	    for (UInt part = ghost_partition_offset[el];
		 part < ghost_partition_offset[el + 1];
		 ++part) {
	      UInt proc = ghost_partition[part];
	      *(buffers_tmp[proc]) = data[el];
	      buffers_tmp[proc]++;
	    }
	  }
	}

	/* -------->>>>-TAGS------------------------------------------------- */
	for (UInt p = 0; p < nb_proc; ++p) {
	  memcpy((char *)buffers_tmp[p], names, names_size);
	  if(p != root) {
	    UInt size = nb_tags * (nb_local_element[p] + nb_ghost_element[p])
	      + uint_names_size;
	    AKANTU_DEBUG_INFO("Sending associated data to proc " << p);
	    requests.push_back(comm->asyncSend(buffers[p], size, p, TAG_DATA));
	  } else {
	    local_data = buffers[p];
	  }
	}
	MeshUtils::setUIntData(mesh, local_data, nb_tags, type);

	comm->waitAll(requests);
	comm->freeCommunicationRequest(requests);
	requests.clear();

	for (UInt p = 0; p < nb_proc; ++p) delete [] buffers[p];
      }
    }

    /* -------->>>>-SIZE----------------------------------------------------- */
    for (UInt p = 0; p < nb_proc; ++p) {
      if(p != root) {
	UInt size[6];
	size[0] = (UInt) _not_defined;
	size[1] = 0;
	size[2] = 0;
	size[3] = 0;
	size[4] = 0;
	size[5] = 0;
	comm->send(size, 6, p, TAG_SIZES);
      }
    }

    /* ---------------------------------------------------------------------- */
    /* ---------------------------------------------------------------------- */
    /**
     * Nodes coordinate construction and synchronization
     */

    std::multimap< UInt, std::pair<UInt, UInt> > nodes_to_proc;
    /// get the list of nodes to send and send them
    Real * local_nodes = NULL;
    UInt nb_nodes_per_proc[nb_proc];

    /* --------<<<<-NB_NODES + NODES----------------------------------------- */
    for (UInt p = 0; p < nb_proc; ++p) {
      UInt nb_nodes;
      UInt * buffer;
      if(p != root) {
	AKANTU_DEBUG_INFO("Receiving list of nodes from proc " << p);
	comm->receive(&nb_nodes, 1, p, TAG_NB_NODES);
	buffer = new UInt[nb_nodes];
	nb_nodes_per_proc[p] = nb_nodes;
	comm->receive(buffer, nb_nodes, p, TAG_NODES);
      } else {
	nb_nodes = old_nodes->getSize();
	nb_nodes_per_proc[p] = nb_nodes;
	buffer = old_nodes->values;
      }

      /// get the coordinates for the selected nodes
      Real * nodes_to_send = new Real[nb_nodes * spatial_dimension];
      Real * nodes_to_send_tmp = nodes_to_send;
      for (UInt n = 0; n < nb_nodes; ++n) {
	memcpy(nodes_to_send_tmp,
	       nodes->values + spatial_dimension * buffer[n],
	       spatial_dimension * sizeof(Real));
	nodes_to_proc.insert(std::make_pair(buffer[n], std::make_pair(p, n)));
	nodes_to_send_tmp += spatial_dimension;
      }

      /* -------->>>>-COORDINATES-------------------------------------------- */
      if(p != root) { /// send them for distant processors
	delete [] buffer;
	AKANTU_DEBUG_INFO("Sending coordinates to proc " << p);
	comm->send(nodes_to_send, nb_nodes * spatial_dimension, p, TAG_COORDINATES);
	delete [] nodes_to_send;
      } else { /// save them for local processor
	local_nodes = nodes_to_send;
      }
    }


    /// construct the local nodes coordinates
    UInt nb_nodes = old_nodes->getSize();
    nodes->resize(nb_nodes);
    memcpy(nodes->values, local_nodes, nb_nodes * spatial_dimension * sizeof(Real));
    delete [] local_nodes;

//     Vector<Int> * nodes_type_per_proc[nb_proc];
//     for (UInt p = 0; p < nb_proc; ++p) {
//       nodes_type_per_proc[p] = new Vector<Int>(nb_nodes);
//     }

//     std::multimap< UInt, std::pair<UInt, UInt> >::iterator it_node;
//     std::pair< std::multimap< UInt, std::pair<UInt, UInt> >::iterator,
//       std::multimap< UInt, std::pair<UInt, UInt> >::iterator > it_range;
//     for (UInt i = 0; i < mesh.nb_global_nodes; ++i) {
//       it_range = nodes_to_proc.equal_range(i);
//       Int node_type = (it_range.first == it_range.second) ?
// 	-1 : (it_range.first)->second.first;
//       for (it_node = it_range.first; it_node != it_range.second; ++it_node) {
// 	nodes_type_per_proc[it_node->second.first]->values[it_node->second.second] = node_type;
//       }
//     }
    /* -------->>>>-NODES_TYPE----------------------------------------------- */
//     std::vector<CommunicationRequest *> requests;
//     for (UInt p = 0; p < nb_proc; ++p) {
//       if(p != root) {
// 	AKANTU_DEBUG_INFO("Sending nodes types to proc " << p);
// 	requests.push_back(comm->asyncSend(nodes_type_per_proc[p]->values,
// 					   nb_nodes_per_proc[p], p, TAG_NODES_TYPE));
//       }
//     }

//     communicator->fillNodesType(nodes_type_per_proc[root]->values, mesh);

//     comm->waitAll(requests);
//     comm->freeCommunicationRequest(requests);
//     requests.clear();

//     for (UInt p = 0; p < nb_proc; ++p) {
//       delete nodes_type_per_proc[p];
//     }

    /* ---------------------------------------------------------------------- */
    /*  Distant (rank != root)                                                */
    /* ---------------------------------------------------------------------- */
  } else {
    /**
     * connectivity and communications scheme construction on distant processors
     */
    ElementType type = _not_defined;
    do {
      /* --------<<<<-SIZE--------------------------------------------------- */
      UInt size[6];
      comm->receive(size, 6, root, TAG_SIZES);

      type          = (ElementType) size[0];
      UInt nb_local_element     = size[1];
      UInt nb_ghost_element     = size[2];
      UInt nb_element_to_send   = size[3];
      UInt nb_tags              = size[4];
      UInt names_size           = size[5];
      UInt uint_names_size = names_size / sizeof(UInt) + (names_size % sizeof(UInt) ? 1 : 0);

      if(type != _not_defined) {
	UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
	/* --------<<<<-CONNECTIVITY----------------------------------------- */
	local_connectivity = new UInt[(nb_local_element + nb_ghost_element) *
				      nb_nodes_per_element];
	AKANTU_DEBUG_INFO("Receiving connectivities from proc " << root);
	comm->receive(local_connectivity, nb_nodes_per_element * (nb_local_element +
								  nb_ghost_element),
			   root, TAG_CONNECTIVITY);

	AKANTU_DEBUG_INFO("Renumbering local connectivities");
	MeshUtils::renumberMeshNodes(mesh,
				     local_connectivity,
				     nb_local_element,
				     nb_ghost_element,
				     type,
				     *old_nodes);

	delete [] local_connectivity;

	/* --------<<<<-PARTITIONS--------------------------------------------- */
	local_partitions = new UInt[nb_element_to_send + nb_ghost_element * 2];
	AKANTU_DEBUG_INFO("Receiving partition informations from proc " << root);
	comm->receive(local_partitions,
			   nb_element_to_send + nb_ghost_element * 2,
			   root, TAG_PARTITIONS);

	AKANTU_DEBUG_INFO("Creating communications scheme");
	communicator->fillCommunicationScheme(local_partitions,
					     nb_local_element,
					     nb_ghost_element,
					     type);

	delete [] local_partitions;

	/* --------<<<<-TAGS------------------------------------------------- */
	if(nb_tags) {
	  AKANTU_DEBUG_INFO("Receiving associated data from proc " << root);
	  UInt size = (nb_local_element + nb_ghost_element) * nb_tags
	    + uint_names_size;
	  local_data = new UInt[size];
	  comm->receive(local_data, size, root, TAG_DATA);

	  MeshUtils::setUIntData(mesh, local_data, nb_tags, type);

	  delete [] local_data;
	}
      }
    } while(type != _not_defined);

    /**
     * Nodes coordinate construction and synchronization on distant processors
     */
    /* -------->>>>-NB_NODES + NODES----------------------------------------- */
    AKANTU_DEBUG_INFO("Sending list of nodes to proc " << root);
    UInt nb_nodes = old_nodes->getSize();
    comm->send(&nb_nodes, 1, root, TAG_NB_NODES);
    comm->send(old_nodes->values, nb_nodes, root, TAG_NODES);

    /* --------<<<<-COORDINATES---------------------------------------------- */
    nodes->resize(nb_nodes);
    AKANTU_DEBUG_INFO("Receiving coordinates from proc " << root);
    comm->receive(nodes->values, nb_nodes * spatial_dimension, root, TAG_COORDINATES);

    /* --------<<<<-NODES_TYPE----------------------------------------------- */
//     Int * nodes_types = new Int[nb_nodes];
//     AKANTU_DEBUG_INFO("Receiving nodes types from proc " << root);
//     comm->receive(nodes_types, nb_nodes,
// 		  root, TAG_NODES_TYPE);

//     communicator->fillNodesType(nodes_types, mesh);
  }

  comm->broadcast(&(mesh.nb_global_nodes), 1, root);

  AKANTU_DEBUG_OUT();
  return communicator;
}

/* -------------------------------------------------------------------------- */
void Communicator::fillNodesType(Int * nodes_type_tmp, Mesh & mesh) {
  AKANTU_DEBUG_IN();

  StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator();
  UInt my_rank = comm->whoAmI();


  UInt nb_nodes = mesh.getNbNodes();
  std::cout << nb_nodes << std::endl;
  Int * nodes_type = mesh.getNodesTypePointer()->values;


  memcpy(nodes_type, nodes_type_tmp, nb_nodes * sizeof(Int));

  UInt * nodes_set = new UInt[nb_nodes];
  std::fill_n(nodes_set, nb_nodes, 0);

  Mesh::ConnectivityTypeList::const_iterator it;

  const UInt NORMAL_SET = 1;
  const UInt GHOST_SET  = 2;

  bool * already_seen = new bool[nb_nodes];
  std::fill_n(already_seen, nb_nodes, false);

  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList(_not_ghost);
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    UInt nb_element = mesh.getNbElement(type);
    UInt * conn_val = mesh.getConnectivity(type).values;

    for (UInt e = 0; e < nb_element; ++e) {
      for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	if(!already_seen[*conn_val]) {
	  nodes_set[*conn_val] += NORMAL_SET;
	  already_seen[*conn_val] = true;
	}
	conn_val++;
      }
    }
  }

  std::fill_n(already_seen, nb_nodes, false);
  const Mesh::ConnectivityTypeList & ghost_type_list = mesh.getConnectivityTypeList(_ghost);
  for(it = ghost_type_list.begin(); it != ghost_type_list.end(); ++it) {
    ElementType type = *it;
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    UInt nb_element = mesh.getNbGhostElement(type);
    UInt * conn_val = mesh.getGhostConnectivity(type).values;

    for (UInt e = 0; e < nb_element; ++e) {
      for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	if(!already_seen[*conn_val]) {
	  nodes_set[*conn_val] += GHOST_SET;
	  already_seen[*conn_val] = true;
	}
	conn_val++;
      }
    }
  }

  delete [] already_seen;

  for (UInt i = 0; i < nb_nodes; ++i) {
    if(nodes_set[i] == NORMAL_SET) nodes_type[i] = -1;
    if(nodes_set[i] == GHOST_SET) nodes_type[i] = -3;
    if(nodes_set[i] == (GHOST_SET + NORMAL_SET)) {
      if(nodes_type[i] == (Int) my_rank) nodes_type[i] = -2;
    }
  }

  delete [] nodes_set;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Communicator::fillCommunicationScheme(UInt * partition,
					   UInt nb_local_element,
					   UInt nb_ghost_element,
					   ElementType type) {
  AKANTU_DEBUG_IN();

  Element element;
  element.type = type;

  UInt * part = partition;

  part = partition;
  for (UInt lel = 0; lel < nb_local_element; ++lel) {
    UInt nb_send = *part; part++;
    element.element = lel;
    for (UInt p = 0; p < nb_send; ++p) {
      UInt proc = *part; part++;

      AKANTU_DEBUG(dblAccessory, "Must send : " << element << " to proc " << proc);
      (send_element[proc]).push_back(element);
    }
  }

  for (UInt gel = 0; gel < nb_ghost_element; ++gel) {
    UInt proc = *part; part++;
    element.element = gel;

    AKANTU_DEBUG(dblAccessory, "Must recv : " << element << " from proc " << proc);
    recv_element[proc].push_back(element);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Communicator::synchronize(GhostSynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  asynchronousSynchronize(tag);

  waitEndSynchronize(tag);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Communicator::asynchronousSynchronize(GhostSynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(send_requests.size() == 0,
		      "There must be some pending sending communications");

  for (UInt p = 0; p < nb_proc; ++p) {
    UInt ssize = (size_to_send[tag]).values[p];
    if(p == rank || ssize == 0) continue;

    send_buffer[p].resize(ssize);
    Real * buffer = send_buffer[p].values;

    Element * elements = &(send_element[p].at(0));
    UInt nb_elements   =  send_element[p].size();
    AKANTU_DEBUG_INFO("Packing data for proc " << p
		      << " (" << ssize << "/" << nb_elements
		      <<" data to send/elements)");
    for (UInt el = 0; el < nb_elements; ++el) {
      ghost_synchronizer->packData(&buffer, *elements, tag);
      elements++;
    }
    std::cerr << std::dec;
    AKANTU_DEBUG_INFO("Posting send to proc " << p);
    send_requests.push_back(static_communicator->asyncSend(send_buffer[p].values,
						       ssize,
						       p,
						       (Int) tag));
  }

  AKANTU_DEBUG_ASSERT(recv_requests.size() == 0,
		      "There must be some pending receive communications");

  for (UInt p = 0; p < nb_proc; ++p) {
    UInt rsize = (size_to_receive[tag]).values[p];
    if(p == rank || rsize == 0) continue;
    recv_buffer[p].resize(rsize);

    AKANTU_DEBUG_INFO("Posting receive from proc " << p << " (" << rsize << " data to receive)");
    recv_requests.push_back(static_communicator->asyncReceive(recv_buffer[p].values,
							      rsize,
							      p,
							      (Int) tag));
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Communicator::waitEndSynchronize(GhostSynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  std::vector<CommunicationRequest *> req_not_finished;
  std::vector<CommunicationRequest *> * req_not_finished_tmp = &req_not_finished;
  std::vector<CommunicationRequest *> * recv_requests_tmp = &recv_requests;
  while(!recv_requests_tmp->empty()) {

    for (std::vector<CommunicationRequest *>::iterator req_it = recv_requests_tmp->begin();
	 req_it != recv_requests_tmp->end() ; ++req_it) {
      CommunicationRequest * req = *req_it;

      if(static_communicator->testRequest(req)) {
	UInt proc = req->getSource();
	AKANTU_DEBUG_INFO("Unpacking data coming from proc " << proc);
	Real * buffer = recv_buffer[proc].values;

	Element * elements = &recv_element[proc].at(0);
	UInt nb_elements   =  recv_element[proc].size();
	for (UInt el = 0; el < nb_elements; ++el) {
	  ghost_synchronizer->unpackData(&buffer, *elements, tag);
	  elements++;
	}

	static_communicator->freeCommunicationRequest(req);
      } else {
	req_not_finished_tmp->push_back(req);
      }
    }

    std::vector<CommunicationRequest *> * swap = req_not_finished_tmp;
    req_not_finished_tmp = recv_requests_tmp;
    recv_requests_tmp = swap;

    req_not_finished_tmp->clear();
  }


  static_communicator->waitAll(send_requests);
  static_communicator->freeCommunicationRequest(send_requests);
  send_requests.clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Communicator::allReduce(Real * values, UInt nb_values, const SynchronizerOperation & op) {
  AKANTU_DEBUG_IN();
  static_communicator->allReduce(values, nb_values, op);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Communicator::registerTag(GhostSynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(size_to_send.find(tag) == size_to_send.end(),
		      "The GhostSynchronizationTag " << tag
		      << "is already registered in " << id);

  size_to_send   [tag].resize(nb_proc);
  size_to_receive[tag].resize(nb_proc);

  for (UInt p = 0; p < nb_proc; ++p) {
    UInt ssend    = 0;
    UInt sreceive = 0;
    if(p != rank) {
      for (std::vector<Element>::const_iterator sit = send_element[p].begin();
	   sit != send_element[p].end();
	   ++sit) {
	ssend += ghost_synchronizer->getNbDataToPack(*sit, tag);
      }

      for (std::vector<Element>::const_iterator rit = recv_element[p].begin();
	   rit != recv_element[p].end();
	   ++rit) {
	sreceive += ghost_synchronizer->getNbDataToUnpack(*rit, tag);
      }
      AKANTU_DEBUG_INFO("I have " << ssend << " data to send to " << p
			<< " and " << sreceive << " data to receive for tag " << tag);
    }

    size_to_send   [tag].values[p] = ssend;
    size_to_receive[tag].values[p] = sreceive;
  }

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
