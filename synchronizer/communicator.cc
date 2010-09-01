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
static const Communicator & Communicator::createCommunicator(UInt root,
							     Mesh & mesh,
							     const MeshPartition & partition) {
  AKANTU_DEBUG_IN();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();

  UInt nb_proc = comm.getNbProc();
  UInt my_rank = comm.getWhoAmI();

  UInt * local_connectivity;
  Vector<UInt> old_nodes;
  Vector<Real> * nodes = mesh.getNodesPointer();

  UInt spatial_dimension = nodes.getNbComponent();

  Communicator * communicator = new Communicator(nb_proc);

  if(my_rank == root) {
    AKANTU_DEBUG_ASSERT(partition.getNbPartition() == nb_proc,
			"The number of partition does not match the number of processors");

    /**
     * connectivity and communications scheme construction
     */
    const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
    Mesh::ConnectivityTypeList::const_iterator it;
    for(it = type_list.begin(); it != type_list.end(); ++it) {
      ElementType type = *it;

      if(Mesh::getSpatialDimension(type) != mesh.getSpatialDimension()) continue;

      // Vector<UInt> ** proc_element_to_send    = element_to_send   [type];
      // Vector<UInt> ** proc_element_to_receive = element_to_receive[type];

      UInt nb_element = mesh.getNbElement(*it);
      UInt nb_local_element[nb_proc];
      UInt nb_ghost_element[nb_proc];
      UInt nb_element_to_send[nb_proc];

      UInt * partition_num = partition.getPartition(type).values;
      UInt * reordering = new UInt[nb_element];

      memset(nb_local_element, 0, nb_proc*sizeof(UInt));
      memset(nb_ghost_element, 0, nb_proc*sizeof(UInt));
      memset(nb_element_to_send, 0, nb_proc*sizeof(UInt));

      UInt * ghost_partition = partition.getGhostPartition(type).values;
      UInt size_ghost_partition = partition.getGhostPartition(type).getSize();
      UInt * ghost_partition_offset = partition.getGhostPartitionOffset(type).values;
      UInt * ghost_ordering = new UInt[size_ghost_partition];

      /// constricting the reordering structures
      for (UInt el = 0; el < nb_element; ++el) {
	reordering[el] = nb_local_element[partition_num[el]]++;

	for (UInt part = ghost_partition_offset[el];
	     part < ghost_partition_offset[el + 1];
	     ++part) {
	  ghost_ordering[part] = nb_ghost_element[ghost_partition[part]]++;
	}
	nb_element_to_send[partition_num[el]] +=
	  ghost_partition_offset[el + 1] - ghost_partition_offset[el];
      }

      /// allocating buffers
      UInt * buffers[nb_proc];
      UInt * buffers_tmp[nb_proc];
      for (UInt p = 0; p < nb_proc; ++p) {
	buffers = new UInt[nb_nodes_per_element * (nb_local_element[p] +
						   nb_ghost_element[p])];
	buffers_tmp = buffers;
      }

      /// copying the local connectivity
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
      std::vector<CommunicationRequest> requests;
      for (UInt p = 0; p < nb_proc; ++p) {
	if(p != root) {
	  UInt size[4];
	  size[0] = (UInt) type;
	  size[1] = nb_local_element[p];
	  size[2] = nb_ghost_element[p];
	  size[3] = nb_element_to_send[p];
	  comm.send<UInt>(size, 4, p, 0);
	  requests.push_back(comm.asyncSend<UInt>(buffers[p],
						  nb_nodes_per_element * (nb_local_element[p] +
									  nb_ghost_element[p]),
						  p, 1));
	} else {
	  local_connectivity = buffers[p];
	}
      }

      /// create the renumbered connectivity
      renumberMeshNodes(local_connectivity,
			mesh.getConnectivityPointer(type),
p			old_nodes);

      comm.waitAll(requests);
      for (UInt p = 0; p < nb_proc; ++p) {
	delete [] buffers[p];
      }

      /// splitting the partition information to send them to processors
      for (UInt el = 0; el < nb_element; ++el) {
	reordering[el] = nb_local_element[partition_num[el]]++;

	for (UInt part = ghost_partition_offset[el];
	     part < ghost_partition_offset[el + 1];
	     ++part) {
	  ghost_ordering[part] = nb_ghost_element[ghost_partition[part]]++;
	}
	nb_element_to_send[partition_num[el]] +=
	  ghost_partition_offset[el + 1] - ghost_partition_offset[el];
      }


      /// last data to compute the communication scheme
      for (UInt p = 0; p < nb_proc; ++p) {
	if(p != root) {
	  requests.push_back(comm.asyncSend<UInt>(buffers[p],
						  nb_nodes_per_element * (nb_local_element[p] +
									  nb_ghost_element[p]),
						  p, 1));
	} else {
	  local_connectivity = buffers[p];
	}
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
	comm.receive<UInt>(&nb_nodes, 1, p, 0);
	buffer = new UInt[nb_nodes];
	comm.receive<UInt>(buffer, nb_nodes, p, 1);
      } else {
	nb_nodes = old_nodes.getSize();
	buffer = old_nodes.values;
      }

      /// get the coordinates for the selected nodes
      Real * nodes_to_send = new Real[nb_nodes * spatial_dimension];
      Real * nodes_to_send_tmp = nodes_to_send;
      for (UInt n = 0; n < nb_nodes; ++n) {
	memcpy(nodes_to_send_tmp,
	       nodes.values + spatial_dimension * buffer[n],
	       spatial_dimension * sizeof(Real));
	nodes_to_send_tmp += spatial_dimension;
      }

      if(p != root) { /// send them for distant processors
	delete [] buffer;
	comm.send<Real>(nodes_to_send, nb_nodes * spatial_dimension, p, 2);
	delete [] nodes_to_send;
      } else { /// save them for local processor
	local_nodes = nodes_to_send;
      }
    }

    /// construct the local nodes coordinates
    nb_nodes = old_nodes.getSize();
    nodes.resize(nb_nodes);
    memcpy(nodes, local_nodes, nb_nodes * spatial_dimension * sizeof(Real));
    delete [] local_nodes;

  } else {
    /**
     * connectivity and communications scheme construction on distant processors
     */
    UInt size[4];
    do {
      comm.receive<UInt>(size, 4, root, 0);

      ElementType type        = (ElementType) size[0];
      UInt nb_local_element   = size[1];
      UInt nb_ghost_element   = size[2];
      UInt nb_element_to_send = size[3];

      if(type != _not_defined) {
	local_connectivity = new UInt[nb_element * nb_nodes_per_element];
	comm.receive<UInt>(buffer, nb_nodes_per_element (nb_local_element +
							 nb_ghost_element),
			   root, 1);

	renumberMeshNodes(local_connectivity,
			  mesh.getConnectivityPointer(type),
			  old_nodes);

	delete [] local_connectivity;
      }
    } while(type != _not_defined);


    /**
     * Nodes coordinate construction and synchronization on distant processors
     */
    UInt nb_nodes = old_nodes.getSize();
    comm.send<UInt>(&nb_nodes, 1, root, 0);
    comm.send<UInt>(old_nodes.values, nb_nodes, root, 1);

    nodes.resize(nb_nodes);
    comm.receive<UInt>(nodes->values, nb_nodes * spatial_dimension, root, 2);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void communicationScheme(local_connectivity, nb_local_element, nb_ghost_element)
