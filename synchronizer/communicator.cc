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
 * <insert license here>
 *
 */

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

  UInt nb_proc = static_communicator->getNbProc();
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

  UInt psize = static_communicator->getNbProc();
  for (UInt p = 0; p < psize; ++p) {
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

  StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator();

  UInt nb_proc = comm->getNbProc();
  UInt my_rank = comm->whoAmI();

  UInt * local_connectivity;
  UInt * local_partitions;
  Vector<UInt> old_nodes;
  Vector<Real> * nodes = mesh.getNodesPointer();

  UInt spatial_dimension = nodes->getNbComponent();

  Communicator * communicator = new Communicator();

  if(nb_proc == 1) return communicator;

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

      /// constricting the reordering structures
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

      /// send all connectivity and ghost information to all processors
      std::vector<CommunicationRequest *> requests;
      for (UInt p = 0; p < nb_proc; ++p) {
	if(p != root) {
	  UInt size[4];
	  size[0] = (UInt) type;
	  size[1] = nb_local_element[p];
	  size[2] = nb_ghost_element[p];
	  size[3] = nb_element_to_send[p];
	  comm->send(size, 4, p, 0);
	  AKANTU_DEBUG_INFO("Sending connectivities to proc " << p);
	  requests.push_back(comm->asyncSend(buffers[p],
					    nb_nodes_per_element * (nb_local_element[p] +
								    nb_ghost_element[p]),
					    p, 1));
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
				   old_nodes);

      comm->waitAll(requests);
      comm->freeCommunicationRequest(requests);
      requests.clear();
      for (UInt p = 0; p < nb_proc; ++p) {
	delete [] buffers[p];
      }


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

      /// last data to compute the communication scheme
      for (UInt p = 0; p < nb_proc; ++p) {
	if(p != root) {
	  AKANTU_DEBUG_INFO("Sending partition informations to proc " << p);
	  requests.push_back(comm->asyncSend(buffers[p],
					     nb_element_to_send[p] + nb_ghost_element[p],
					     p, 2));
	} else {
	  local_partitions = buffers[p];
	}
      }

      AKANTU_DEBUG_INFO("Creating communications scheme");
      communicator->fillCommunicationScheme(local_partitions,
					    nb_local_element[root],
					    nb_ghost_element[root],
					    nb_element_to_send[root],
					    type);

      comm->waitAll(requests);
      comm->freeCommunicationRequest(requests);
      requests.clear();

      for (UInt p = 0; p < nb_proc; ++p) {
	delete [] buffers[p];
      }
    }

    for (UInt p = 0; p < nb_proc; ++p) {
      if(p != root) {
	UInt size[4];
	size[0] = (UInt) _not_defined;
	size[1] = 0;
	size[2] = 0;
	size[3] = 0;
	comm->send(size, 4, p, 0);
      }
    }

    /**
     * Nodes coordinate construction and synchronization
     */

    /// get the list of nodes to send and send them
    Real * local_nodes;
    for (UInt p = 0; p < nb_proc; ++p) {
      UInt nb_nodes;
      UInt * buffer;
      if(p != root) {
	AKANTU_DEBUG_INFO("Receiving list of nodes from proc " << p);
	comm->receive(&nb_nodes, 1, p, 0);
	buffer = new UInt[nb_nodes];
	comm->receive(buffer, nb_nodes, p, 3);
      } else {
	nb_nodes = old_nodes.getSize();
	buffer = old_nodes.values;
      }

      /// get the coordinates for the selected nodes
      Real * nodes_to_send = new Real[nb_nodes * spatial_dimension];
      Real * nodes_to_send_tmp = nodes_to_send;
      for (UInt n = 0; n < nb_nodes; ++n) {
	memcpy(nodes_to_send_tmp,
	       nodes->values + spatial_dimension * buffer[n],
	       spatial_dimension * sizeof(Real));
	nodes_to_send_tmp += spatial_dimension;
      }

      if(p != root) { /// send them for distant processors
	delete [] buffer;
	AKANTU_DEBUG_INFO("Sending coordinates to proc " << p);
	comm->send(nodes_to_send, nb_nodes * spatial_dimension, p, 4);
	delete [] nodes_to_send;
      } else { /// save them for local processor
	local_nodes = nodes_to_send;
      }
    }

    /// construct the local nodes coordinates
    UInt nb_nodes = old_nodes.getSize();
    nodes->resize(nb_nodes);
    memcpy(nodes->values, local_nodes, nb_nodes * spatial_dimension * sizeof(Real));
    delete [] local_nodes;

    /* ---------------------------------------------------------------------- */
    /*  Distant (rank != root)                                                */
    /* ---------------------------------------------------------------------- */
  } else {
    /**
     * connectivity and communications scheme construction on distant processors
     */
    ElementType type = _not_defined;
    do {
      UInt size[4];
      comm->receive(size, 4, root, 0);

      type          = (ElementType) size[0];
      UInt nb_local_element     = size[1];
      UInt nb_ghost_element     = size[2];
      UInt nb_element_to_send   = size[3];

      if(type != _not_defined) {
	UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

	local_connectivity = new UInt[(nb_local_element + nb_ghost_element) *
				      nb_nodes_per_element];
	AKANTU_DEBUG_INFO("Receiving connectivities from proc " << root);
	comm->receive(local_connectivity, nb_nodes_per_element * (nb_local_element +
								  nb_ghost_element),
			   root, 1);

	AKANTU_DEBUG_INFO("Renumbering local connectivities");
	MeshUtils::renumberMeshNodes(mesh,
				     local_connectivity,
				     nb_local_element,
				     nb_ghost_element,
				     type,
				     old_nodes);

	delete [] local_connectivity;

	local_partitions = new UInt[nb_element_to_send + nb_ghost_element * 2];
	AKANTU_DEBUG_INFO("Receiving partition informations from proc " << root);
	comm->receive(local_partitions,
			   nb_element_to_send + nb_ghost_element * 2,
			   root, 2);

	AKANTU_DEBUG_INFO("Creating communications scheme");
	communicator->fillCommunicationScheme(local_partitions,
					     nb_local_element,
					     nb_ghost_element,
					     nb_element_to_send,
					     type);

	delete [] local_partitions;
      }
    } while(type != _not_defined);

    /**
     * Nodes coordinate construction and synchronization on distant processors
     */
    AKANTU_DEBUG_INFO("Sending list of nodes to proc " << root);
    UInt nb_nodes = old_nodes.getSize();
    comm->send(&nb_nodes, 1, root, 0);
    comm->send(old_nodes.values, nb_nodes, root, 3);

    nodes->resize(nb_nodes);
    AKANTU_DEBUG_INFO("Receiving coordinates from proc " << root);
    comm->receive(nodes->values, nb_nodes * spatial_dimension, root, 4);
  }

  AKANTU_DEBUG_OUT();
  return communicator;
}

/* -------------------------------------------------------------------------- */
void Communicator::fillCommunicationScheme(UInt * partition,
					   UInt nb_local_element,
					   UInt nb_ghost_element,
					   UInt nb_element_to_send,
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

  UInt psize = static_communicator->getNbProc();
  UInt prank = static_communicator->whoAmI();

  AKANTU_DEBUG_ASSERT(send_requests.size() == 0,
		      "There must be some pending sending communications");

  for (UInt p = 0; p < psize; ++p) {
    UInt ssize = (size_to_send[tag]).values[p];
    if(p == prank || ssize == 0) continue;

    AKANTU_DEBUG_INFO("Packing data for proc" << p);
    send_buffer[p].resize(ssize);
    Real * buffer = send_buffer[p].values;


    Element * elements = &(send_element[p].at(0));
    UInt nb_elements   =  send_element[p].size();
    for (UInt el = 0; el < nb_elements; ++el) {
      ghost_synchronizer->packData(&buffer, *elements, tag);
      elements++;
    }

    AKANTU_DEBUG_INFO("Posting send to proc " << p);
    send_requests.push_back(static_communicator->asyncSend(send_buffer[p].values,
						       ssize,
						       p,
						       (Int) tag));
  }

  AKANTU_DEBUG_ASSERT(recv_requests.size() == 0,
		      "There must be some pending receive communications");

  for (UInt p = 0; p < psize; ++p) {
    UInt rsize = size_to_receive[tag].values[p];
    if(p == prank || rsize == 0) continue;
    recv_buffer[p].resize(rsize);

    AKANTU_DEBUG_INFO("Posting receive from proc " << p);
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
void Communicator::registerTag(GhostSynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  UInt psize = static_communicator->getNbProc();
  UInt prank = static_communicator->whoAmI();

  AKANTU_DEBUG_ASSERT(size_to_send.find(tag) == size_to_send.end(),
		      "The GhostSynchronizationTag " << tag
		      << "is already registered in " << id);

  size_to_send   [tag].resize(psize);
  size_to_receive[tag].resize(psize);

  for (UInt p = 0; p < psize; ++p) {
    UInt ssend    = 0;
    UInt sreceive = 0;
    if(p != prank) {
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
			<< " and " << sreceive << " data to receive");
    }

    size_to_send   [tag].values[p] = ssend;
    size_to_receive[tag].values[p] = sreceive;
  }

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
