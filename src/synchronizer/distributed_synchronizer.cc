/**
 * @file   distributed_synchronizer.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Thu Jun 16 16:36:52 2011
 *
 * @brief  implementation of a  communicator using a static_communicator for real
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
#include "aka_common.hh"
#include "distributed_synchronizer.hh"
#include "static_communicator.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include <map>
#include <iostream>

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
DistributedSynchronizer::DistributedSynchronizer(Mesh & mesh,
						 SynchronizerID id,
						 MemoryID memory_id) :
  Synchronizer(id, memory_id),
  mesh(mesh), static_communicator(&StaticCommunicator::getStaticCommunicator())
{
  AKANTU_DEBUG_IN();

  nb_proc = static_communicator->getNbProc();
  rank    = static_communicator->whoAmI();

  send_element = new Vector<Element>[nb_proc];
  recv_element = new Vector<Element>[nb_proc];

  mesh.registerEventHandler(*this);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
DistributedSynchronizer::~DistributedSynchronizer() {
  AKANTU_DEBUG_IN();

  for (UInt p = 0; p < nb_proc; ++p) {
  //   send_buffer[p].resize(0);
  //   recv_buffer[p].resize(0);

    send_element[p].clear();
    recv_element[p].clear();
  }

  // delete [] send_buffer;
  // delete [] recv_buffer;

  delete [] send_element;
  delete [] recv_element;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
DistributedSynchronizer * DistributedSynchronizer::
createDistributedSynchronizerMesh(Mesh & mesh,
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

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  UInt nb_proc = comm.getNbProc();
  UInt my_rank = comm.whoAmI();

  DistributedSynchronizer & communicator = *(new DistributedSynchronizer(mesh, id, memory_id));

  if(nb_proc == 1) return &communicator;

  UInt * local_connectivity = NULL;
  UInt * local_partitions = NULL;
  UInt * local_data = NULL;
  Vector<UInt> * old_nodes = mesh.getNodesGlobalIdsPointer();
  old_nodes->resize(0);
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
    Mesh::type_iterator it  = mesh.firstType(mesh.getSpatialDimension(),
					     _not_ghost,
					     _ek_not_defined);
    Mesh::type_iterator end = mesh.lastType(mesh.getSpatialDimension(),
					     _not_ghost,
					     _ek_not_defined);
    UInt count = 0;
    for(; it != end; ++it) {
      ElementType type = *it;

      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
      UInt nb_element = mesh.getNbElement(*it);
      UInt nb_local_element[nb_proc];
      UInt nb_ghost_element[nb_proc];
      UInt nb_element_to_send[nb_proc];

      memset(nb_local_element, 0, nb_proc*sizeof(UInt));
      memset(nb_ghost_element, 0, nb_proc*sizeof(UInt));
      memset(nb_element_to_send, 0, nb_proc*sizeof(UInt));

      UInt * partition_num = partition->getPartition(type, _not_ghost).values;

      UInt * ghost_partition = partition->getGhostPartition(type, _ghost).values;
      UInt * ghost_partition_offset = partition->getGhostPartitionOffset(type, _ghost).values;

      UIntDataMap & uint_data_map = mesh.getUIntDataMap(type, _not_ghost);
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
      UInt * conn_val = mesh.getConnectivity(type, _not_ghost).values;
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
      UIntDataMap::iterator it_data;
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
          AKANTU_DEBUG_INFO("Sending connectivities informations to proc " << p << " TAG("<< Tag::genTag(my_rank, count, TAG_SIZES) <<")");
          comm.send(size, 6, p, Tag::genTag(my_rank, count, TAG_SIZES));

          AKANTU_DEBUG_INFO("Sending connectivities to proc " << p << " TAG("<< Tag::genTag(my_rank, count, TAG_CONNECTIVITY) <<")");
          requests.push_back(comm.asyncSend(buffers[p],
                                            nb_nodes_per_element * (nb_local_element[p] +
                                                                    nb_ghost_element[p]),
                                            p, Tag::genTag(my_rank, count, TAG_CONNECTIVITY)));
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

      comm.waitAll(requests);
      comm.freeCommunicationRequest(requests);
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
          AKANTU_DEBUG_INFO("Sending partition informations to proc " << p << " TAG("<< Tag::genTag(my_rank, count, TAG_PARTITIONS) <<")");
          requests.push_back(comm.asyncSend(buffers[p],
                                             nb_element_to_send[p] + nb_ghost_element[p],
					    p, Tag::genTag(my_rank, count, TAG_PARTITIONS)));
        } else {
          local_partitions = buffers[p];
        }
      }

      AKANTU_DEBUG_INFO("Creating communications scheme");
      communicator.fillCommunicationScheme(local_partitions,
					   nb_local_element[root],
					   nb_ghost_element[root],
					   type);

      comm.waitAll(requests);
      comm.freeCommunicationRequest(requests);
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
            AKANTU_DEBUG_INFO("Sending associated data to proc " << p << " TAG("<< Tag::genTag(my_rank, count, TAG_DATA) <<")");
            requests.push_back(comm.asyncSend(buffers[p], size, p, Tag::genTag(my_rank, count, TAG_DATA)));
          } else {
            local_data = buffers[p];
          }
        }
        MeshUtils::setUIntData(mesh, local_data, nb_tags, type);

        comm.waitAll(requests);
        comm.freeCommunicationRequest(requests);
        requests.clear();

        for (UInt p = 0; p < nb_proc; ++p) delete [] buffers[p];
        delete [] names;
      }
      ++count;
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
	AKANTU_DEBUG_INFO("Sending empty connectivities informations to proc " << p << " TAG("<< Tag::genTag(my_rank, count, TAG_SIZES) <<")");
        comm.send(size, 6, p, Tag::genTag(my_rank, count, TAG_SIZES));
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
    UInt * nodes_per_proc[nb_proc];
    /* --------<<<<-NB_NODES + NODES----------------------------------------- */
    for (UInt p = 0; p < nb_proc; ++p) {
      UInt nb_nodes = 0;
      //      UInt * buffer;
      if(p != root) {
        AKANTU_DEBUG_INFO("Receiving number of nodes from proc " << p << " TAG("<< Tag::genTag(p, 0, TAG_NB_NODES) <<")");
        comm.receive(&nb_nodes, 1, p, Tag::genTag(p, 0, TAG_NB_NODES));
        nodes_per_proc[p] = new UInt[nb_nodes];
        nb_nodes_per_proc[p] = nb_nodes;
        AKANTU_DEBUG_INFO("Receiving list of nodes from proc " << p << " TAG("<< Tag::genTag(p, 0, TAG_NODES) <<")");
        comm.receive(nodes_per_proc[p], nb_nodes, p, Tag::genTag(p, 0, TAG_NODES));
      } else {
        nb_nodes = old_nodes->getSize();
        nb_nodes_per_proc[p] = nb_nodes;
        nodes_per_proc[p] = old_nodes->values;
      }

      /// get the coordinates for the selected nodes
      Real * nodes_to_send = new Real[nb_nodes * spatial_dimension];
      Real * nodes_to_send_tmp = nodes_to_send;
      for (UInt n = 0; n < nb_nodes; ++n) {
        memcpy(nodes_to_send_tmp,
               nodes->values + spatial_dimension * nodes_per_proc[p][n],
               spatial_dimension * sizeof(Real));
        // nodes_to_proc.insert(std::make_pair(buffer[n], std::make_pair(p, n)));
        nodes_to_send_tmp += spatial_dimension;
      }

      /* -------->>>>-COORDINATES-------------------------------------------- */
      if(p != root) { /// send them for distant processors
        AKANTU_DEBUG_INFO("Sending coordinates to proc " << p << " TAG("<< Tag::genTag(my_rank, 0, TAG_COORDINATES) <<")");
        comm.send(nodes_to_send, nb_nodes * spatial_dimension, p, Tag::genTag(my_rank, 0, TAG_COORDINATES));
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

    Vector<Int> * nodes_type_per_proc[nb_proc];
    for (UInt p = 0; p < nb_proc; ++p) {
      nodes_type_per_proc[p] = new Vector<Int>(nb_nodes_per_proc[p]);
    }

    communicator.fillNodesType(mesh);

    /* --------<<<<-NODES_TYPE-1--------------------------------------------- */
    for (UInt p = 0; p < nb_proc; ++p) {
      if(p != root) {
        AKANTU_DEBUG_INFO("Receiving first nodes types from proc " << p << " TAG("<< Tag::genTag(my_rank, count, TAG_NODES_TYPE) <<")");
        comm.receive(nodes_type_per_proc[p]->values,
		     nb_nodes_per_proc[p], p, Tag::genTag(p, 0, TAG_NODES_TYPE));
      } else {
        nodes_type_per_proc[p]->copy(mesh.getNodesType());
      }
      for (UInt n = 0; n < nb_nodes_per_proc[p]; ++n) {
        if((*nodes_type_per_proc[p])(n) == -2)
          nodes_to_proc.insert(std::make_pair(nodes_per_proc[p][n], std::make_pair(p, n)));
      }
    }

    std::multimap< UInt, std::pair<UInt, UInt> >::iterator it_node;
    std::pair< std::multimap< UInt, std::pair<UInt, UInt> >::iterator,
      std::multimap< UInt, std::pair<UInt, UInt> >::iterator > it_range;
    for (UInt i = 0; i < mesh.nb_global_nodes; ++i) {
      it_range = nodes_to_proc.equal_range(i);
      if(it_range.first == nodes_to_proc.end() || it_range.first->first != i) continue;

      UInt node_type = (it_range.first)->second.first;
      for (it_node = it_range.first; it_node != it_range.second; ++it_node) {
        UInt proc = it_node->second.first;
        UInt node = it_node->second.second;
        if(proc != node_type)
          nodes_type_per_proc[proc]->values[node] = node_type;
      }
    }

    /* -------->>>>-NODES_TYPE-2--------------------------------------------- */
    std::vector<CommunicationRequest *> requests;
    for (UInt p = 0; p < nb_proc; ++p) {
      if(p != root) {
        AKANTU_DEBUG_INFO("Sending nodes types to proc " << p << " TAG("<< Tag::genTag(my_rank, 0, TAG_NODES_TYPE) <<")");
        requests.push_back(comm.asyncSend(nodes_type_per_proc[p]->values,
					  nb_nodes_per_proc[p], p, Tag::genTag(my_rank, 0, TAG_NODES_TYPE)));
      } else {
        mesh.getNodesTypePointer()->copy(*nodes_type_per_proc[p]);
      }
    }

    comm.waitAll(requests);
    comm.freeCommunicationRequest(requests);
    requests.clear();

    for (UInt p = 0; p < nb_proc; ++p) {
      if(p != root) delete [] nodes_per_proc[p];
      delete nodes_type_per_proc[p];
    }

    /* ---------------------------------------------------------------------- */
    /*  Distant (rank != root)                                                */
    /* ---------------------------------------------------------------------- */
  } else {
    /**
     * connectivity and communications scheme construction on distant processors
     */
    ElementType type = _not_defined;
    UInt count = 0;
    do {
      /* --------<<<<-SIZE--------------------------------------------------- */
      UInt size[6] = { 0 };
      comm.receive(size, 6, root, Tag::genTag(root, count, TAG_SIZES));

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
        comm.receive(local_connectivity, nb_nodes_per_element * (nb_local_element +
                                                                  nb_ghost_element),
		     root, Tag::genTag(root, count, TAG_CONNECTIVITY));

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
        comm.receive(local_partitions,
		     nb_element_to_send + nb_ghost_element * 2,
		     root, Tag::genTag(root, count, TAG_PARTITIONS));

        AKANTU_DEBUG_INFO("Creating communications scheme");
        communicator.fillCommunicationScheme(local_partitions,
                                             nb_local_element,
                                             nb_ghost_element,
                                             type);

        delete [] local_partitions;

        /* --------<<<<-TAGS------------------------------------------------- */
        if(nb_tags) {
          AKANTU_DEBUG_INFO("Receiving associated data from proc " << root);
          UInt size_data = (nb_local_element + nb_ghost_element) * nb_tags
            + uint_names_size;
          local_data = new UInt[size_data];
          comm.receive(local_data, size_data, root, Tag::genTag(root, count, TAG_DATA));

          MeshUtils::setUIntData(mesh, local_data, nb_tags, type);

          delete [] local_data;
        }
      }
      ++count;
    } while(type != _not_defined);

    /**
     * Nodes coordinate construction and synchronization on distant processors
     */
    /* -------->>>>-NB_NODES + NODES----------------------------------------- */
    AKANTU_DEBUG_INFO("Sending list of nodes to proc " << root);
    UInt nb_nodes = old_nodes->getSize();
    comm.send(&nb_nodes, 1, root, Tag::genTag(my_rank, 0, TAG_NB_NODES));
    comm.send(old_nodes->values, nb_nodes, root, Tag::genTag(my_rank, 0, TAG_NODES));

    /* --------<<<<-COORDINATES---------------------------------------------- */
    nodes->resize(nb_nodes);
    AKANTU_DEBUG_INFO("Receiving coordinates from proc " << root);
    comm.receive(nodes->values, nb_nodes * spatial_dimension, root, Tag::genTag(root, 0, TAG_COORDINATES));

    communicator.fillNodesType(mesh);
    /* --------<<<<-NODES_TYPE-2--------------------------------------------- */
    Int * nodes_types = mesh.getNodesTypePointer()->values;
    AKANTU_DEBUG_INFO("Sending first nodes types to proc " << root);
    comm.send(nodes_types, nb_nodes,
	      root, Tag::genTag(my_rank, 0, TAG_NODES_TYPE));

    /* --------<<<<-NODES_TYPE-2--------------------------------------------- */
    AKANTU_DEBUG_INFO("Receiving nodes types from proc " << root);
    comm.receive(nodes_types, nb_nodes,
		 root, Tag::genTag(root, 0, TAG_NODES_TYPE));
  }

  comm.broadcast(&(mesh.nb_global_nodes), 1, root);

  AKANTU_DEBUG_OUT();
  return &communicator;
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::fillNodesType(Mesh & mesh) {
  AKANTU_DEBUG_IN();

  //  StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator();
  //  UInt my_rank = comm.whoAmI();

  UInt nb_nodes = mesh.getNbNodes();
  //  std::cout << nb_nodes << std::endl;
  Int * nodes_type = mesh.getNodesTypePointer()->values;

  //memcpy(nodes_type, nodes_type_tmp, nb_nodes * sizeof(Int));

  UInt * nodes_set = new UInt[nb_nodes];
  std::fill_n(nodes_set, nb_nodes, 0);

  //  Mesh::ConnectivityTypeList::const_iterator it;

  const UInt NORMAL_SET = 1;
  const UInt GHOST_SET  = 2;

  bool * already_seen = new bool[nb_nodes];

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;
    UInt set = NORMAL_SET;
    if (gt == _ghost) set = GHOST_SET;

    std::fill_n(already_seen, nb_nodes, false);
    Mesh::type_iterator it  = mesh.firstType(0, gt, _ek_not_defined);
    Mesh::type_iterator end = mesh.lastType(0, gt, _ek_not_defined);
    for(; it != end; ++it) {
      ElementType type = *it;

      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
      UInt nb_element = mesh.getNbElement(type, gt);
      UInt * conn_val = mesh.getConnectivity(type, gt).values;

      for (UInt e = 0; e < nb_element; ++e) {
        for (UInt n = 0; n < nb_nodes_per_element; ++n) {
          AKANTU_DEBUG_ASSERT(*conn_val < nb_nodes, "Node " << *conn_val
                              << " bigger than number of nodes " << nb_nodes);
          if(!already_seen[*conn_val]) {
            nodes_set[*conn_val] += set;
            already_seen[*conn_val] = true;
          }
          conn_val++;
        }
      }
    }
  }

  delete [] already_seen;

  for (UInt i = 0; i < nb_nodes; ++i) {
    if(nodes_set[i] == NORMAL_SET) nodes_type[i] = -1;
    else if(nodes_set[i] == GHOST_SET) nodes_type[i] = -3;
    else if(nodes_set[i] == (GHOST_SET + NORMAL_SET)) nodes_type[i] = -2;
  }

  delete [] nodes_set;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::fillCommunicationScheme(UInt * partition,
                                                      UInt nb_local_element,
                                                      UInt nb_ghost_element,
                                                      ElementType type) {
  AKANTU_DEBUG_IN();

  Element element;
  element.type = type;
  element.kind = Mesh::getKind(type);

  UInt * part = partition;

  part = partition;
  for (UInt lel = 0; lel < nb_local_element; ++lel) {
    UInt nb_send = *part; part++;
    element.element = lel;
    element.ghost_type = _not_ghost;
    for (UInt p = 0; p < nb_send; ++p) {
      UInt proc = *part; part++;

      AKANTU_DEBUG(dblAccessory, "Must send : " << element << " to proc " << proc);
      (send_element[proc]).push_back(element);
    }
  }

  for (UInt gel = 0; gel < nb_ghost_element; ++gel) {
    UInt proc = *part; part++;
    element.element = gel;
    element.ghost_type = _ghost;
    AKANTU_DEBUG(dblAccessory, "Must recv : " << element << " from proc " << proc);
    recv_element[proc].push_back(element);
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::asynchronousSynchronize(DataAccessor & data_accessor,
                                		      SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  if (communications.find(tag) == communications.end()) {
    communications[tag].resize(nb_proc);
    computeBufferSize(data_accessor, tag);
  }

  Communication & communication = communications[tag];

  AKANTU_DEBUG_ASSERT(communication.send_requests.size() == 0,
		      "There must be some pending sending communications. Tag is " << tag);


  for (UInt p = 0; p < nb_proc; ++p) {
    UInt ssize = communication.size_to_send[p];
    if(p == rank || ssize == 0) continue;

    CommunicationBuffer & buffer = communication.send_buffer[p];
    buffer.resize(ssize);
#ifndef AKANTU_NDEBUG
    UInt nb_elements   =  send_element[p].getSize();
    AKANTU_DEBUG_INFO("Packing data for proc " << p
		      << " (" << ssize << "/" << nb_elements
		      <<" data to send/elements)");
#endif

    data_accessor.packElementData(buffer, send_element[p], tag);

    AKANTU_DEBUG_ASSERT(buffer.getPackedSize() == ssize,
			"a problem have been introduced with "
			<< "false sent sizes declaration "
			<< buffer.getPackedSize() << " != " << ssize);
    std::cerr << std::dec;
    AKANTU_DEBUG_INFO("Posting send to proc " << p);
    communication.send_requests.push_back(static_communicator->asyncSend(buffer.storage(),
									 ssize,
									 p,
									 (Int) tag));
  }

  AKANTU_DEBUG_ASSERT(communication.recv_requests.size() == 0,
		      "There must be some pending receive communications");

  for (UInt p = 0; p < nb_proc; ++p) {
    UInt rsize = communication.size_to_receive[p];
    if(p == rank || rsize == 0) continue;
    CommunicationBuffer & buffer = communication.recv_buffer[p];
    buffer.resize(rsize);

    AKANTU_DEBUG_INFO("Posting receive from proc " << p << " (" << rsize << " data to receive)");
    communication.recv_requests.push_back(static_communicator->asyncReceive(buffer.storage(),
									    rsize,
									    p,
									    (Int) tag));
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::waitEndSynchronize(DataAccessor & data_accessor,
						 SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(communications.find(tag) != communications.end(), "No communication with the tag \""
		      << tag <<"\" started");

  Communication & communication = communications[tag];

  std::vector<CommunicationRequest *> req_not_finished;
  std::vector<CommunicationRequest *> * req_not_finished_tmp = &req_not_finished;
  std::vector<CommunicationRequest *> * recv_requests_tmp = &(communication.recv_requests);

  //  static_communicator->waitAll(recv_requests);

  while(!recv_requests_tmp->empty()) {

    for (std::vector<CommunicationRequest *>::iterator req_it = recv_requests_tmp->begin();
	 req_it != recv_requests_tmp->end() ; ++req_it) {
      CommunicationRequest * req = *req_it;

      if(static_communicator->testRequest(req)) {
	UInt proc = req->getSource();
	AKANTU_DEBUG_INFO("Unpacking data coming from proc " << proc);
	CommunicationBuffer & buffer = communication.recv_buffer[proc];

	data_accessor.unpackElementData(buffer, recv_element[proc], tag);
	buffer.resize(0);

	AKANTU_DEBUG_ASSERT(buffer.getLeftToUnpack() == 0,
			    "all data have not been unpacked: "
			    << buffer.getLeftToUnpack() << " bytes left");
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

  static_communicator->waitAll(communication.send_requests);
  for (std::vector<CommunicationRequest *>::iterator req_it = communication.send_requests.begin();
       req_it != communication.send_requests.end() ; ++req_it) {
    CommunicationRequest & req = *(*req_it);

    if(static_communicator->testRequest(&req)) {
      UInt proc = req.getDestination();
      CommunicationBuffer & buffer = communication.send_buffer[proc];
      buffer.resize(0);
      static_communicator->freeCommunicationRequest(&req);
    }
  }
  communication.send_requests.clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::computeBufferSize(DataAccessor & data_accessor,
						SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  for (UInt p = 0; p < nb_proc; ++p) {
    UInt ssend    = 0;
    UInt sreceive = 0;
    if(p != rank) {
      ssend    = data_accessor.getNbDataForElements(send_element[p], tag);
      sreceive = data_accessor.getNbDataForElements(recv_element[p], tag);

      AKANTU_DEBUG_INFO("I have " << ssend << "(" << ssend / 1024.
			<< "kB - "<< send_element[p].getSize() <<" element(s)) data to send to " << p
			<< " and " << sreceive << "(" << sreceive / 1024.
			<< "kB - "<< recv_element[p].getSize() <<" element(s)) data to receive for tag "
			<< tag);
    }

    communications[tag].size_to_send   [p] = ssend;
    communications[tag].size_to_receive[p] = sreceive;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  Int prank = StaticCommunicator::getStaticCommunicator().whoAmI();
  Int psize = StaticCommunicator::getStaticCommunicator().getNbProc();
  stream << "[" << prank << "/" << psize << "]" << space << "DistributedSynchronizer [" << std::endl;
  for (UInt p = 0; p < nb_proc; ++p) {
    if (p == UInt(prank)) continue;
    stream << "[" << prank << "/" << psize << "]" << space << " + Communication to proc " << p << " [" << std::endl;
    if(AKANTU_DEBUG_TEST(dblDump)) {
      stream << "[" << prank << "/" << psize << "]" << space << "    - Element to send to proc " << p << " [" << std::endl;

      Vector<Element>::iterator<Element> it_el  = send_element[p].begin();
      Vector<Element>::iterator<Element> end_el = send_element[p].end();
      for(;it_el != end_el; ++it_el)
	stream << "[" << prank << "/" << psize << "]" << space << "       " << *it_el << std::endl;
      stream << "[" << prank << "/" << psize << "]" << space << "   ]" << std::endl;

      stream << "[" << prank << "/" << psize << "]" << space << "    - Element to recv from proc " << p << " [" << std::endl;

      it_el  = recv_element[p].begin();
      end_el = recv_element[p].end();
      for(;it_el != end_el; ++it_el)
	stream << "[" << prank << "/" << psize << "]" << space << "       " << *it_el << std::endl;stream << "[" << prank << "/" << psize << "]" << space << "   ]" << std::endl;
    }

    std::map< SynchronizationTag, Communication>::const_iterator it = communications.begin();
    std::map< SynchronizationTag, Communication>::const_iterator end = communications.end();
    for (; it != end; ++it) {
      const SynchronizationTag & tag = it->first;
      const Communication & communication = it->second;
      UInt ssend    = communication.size_to_send[p];
      UInt sreceive = communication.size_to_receive[p];
      stream << "[" << prank << "/" << psize << "]" << space << "     - Tag " << tag << " -> " << ssend << "byte(s) -- <- " << sreceive << "byte(s)" << std::endl;
    }

    std::cout << "[" << prank << "/" << psize << "]" << space << " ]" << std::endl;
  }
  std::cout << "[" << prank << "/" << psize << "]" << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::onElementsRemoved(const Vector<Element> & element_to_remove,
						const ByElementTypeUInt & new_numbering,
						__attribute__((unused)) const RemovedElementsEvent & event) {
  AKANTU_DEBUG_IN();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  UInt psize = comm.getNbProc();
  UInt prank = comm.whoAmI();

  std::vector<CommunicationRequest *> isend_requests;
  Vector<UInt> * list_of_el = new Vector<UInt>[nb_proc];
  // Handling ghost elements
  for (UInt p = 0; p < psize; ++p) {
    if (p == prank) continue;

    Vector<Element> & recv = recv_element[p];
    if(recv.getSize() == 0) continue;

    Vector<Element>::iterator<Element> recv_begin = recv.begin();
    Vector<Element>::iterator<Element> recv_end   = recv.end();

    Vector<Element>::const_iterator<Element> er_it  = element_to_remove.begin();
    Vector<Element>::const_iterator<Element> er_end = element_to_remove.end();
    // for (;it != end; ++it) {
    //   const Element & el = *it;
    //   if(el.ghost_type == _ghost) {
    // 	Vector<Element>::iterator<Element> pos = std::find(recv_begin, recv_end, el);
    // 	if(pos != recv_end) {
    // 	  UInt i = pos - recv_begin;
    // 	  AKANTU_DEBUG_ASSERT(i < recv.getSize(), "The element is out of the list");
    // 	  list_of_el_to_remove.push_back(i);
    // 	}
    //   } else {
    // 	/// \todo handle the removal of element to send from the ghost of other procs
    // 	AKANTU_EXCEPTION("DistributedSynchronizer do not know how to remove _not_ghost element, ask a developer to finish is job");
    //   }
    // }

    Vector<UInt> & list = list_of_el[p];
    for (UInt i = 0; recv_begin != recv_end; ++i, ++recv_begin) {
      const Element & el = *recv_begin;
      Vector<Element>::const_iterator<Element> pos = std::find(er_it, er_end, el);
      if(pos == er_end) {
	list.push_back(i);
      }
    }

    // if (list_of_el_to_remove.getSize() != 0) {
    //   std::sort(list_of_el_to_remove.begin(), list_of_el_to_remove.end());
    //   UInt lr = 0;
    //   for (UInt l = 0; l < recv.getSize(); ++l) {
    // 	if(lr < list_of_el_to_remove.getSize() && list_of_el_to_remove(lr) == l) ++lr;
    // 	else list.push_back(l);
    //   }

    //   AKANTU_DEBUG_INFO("Sending list of elements (" << list.getSize() << " elements) of proc " << p << " TAG("<< Tag::genTag(p, 0, 0) <<")");
    //   list.push_back(UInt(-1));
    // } else {
    //   list.push_back(UInt(0));
    // }

    if(list.getSize() == recv.getSize())
      list.push_back(UInt(0));
    else list.push_back(UInt(-1));

    isend_requests.push_back(comm.asyncSend(list.storage(), list.getSize(),
					    p, Tag::genTag(prank, 0, 0)));

    list.erase(list.getSize() - 1);
    if(list.getSize() == recv.getSize()) continue;

    Vector<Element> new_recv;
    for (UInt nr = 0; nr < list.getSize(); ++nr) {
      Element & el = recv(list(nr));
      el.element = new_numbering(el.type, el.ghost_type)(el.element);
      new_recv.push_back(el);
    }

    AKANTU_DEBUG_INFO("I had " << recv.getSize() << " elements to recv from proc " << p << " and "
		      << list.getSize() << " elements to keep. I have "
		      << new_recv.getSize() << " elements left.");
    recv.copy(new_recv);
  }

  for (UInt p = 0; p < psize; ++p) {
    if (p == prank) continue;
    Vector<Element> & send = send_element[p];

    if(send.getSize() == 0) continue;

    CommunicationStatus status;
    AKANTU_DEBUG_INFO("Getting number of elements of proc " << p << " not needed anymore TAG("<< Tag::genTag(p, 0, 0) <<")");
    comm.probe<UInt>(p, Tag::genTag(p, 0, 0), status);
    Vector<UInt> list(status.getSize());

    AKANTU_DEBUG_INFO("Receiving list of elements (" << status.getSize() - 1 << " elements) no longer needed by proc " << p << " TAG("<< Tag::genTag(p, 0, 0) <<")");
    comm.receive(list.storage(), list.getSize(),
		 p, Tag::genTag(p, 0, 0));

    if(list.getSize() == 1 && list(0) == 0) continue;

    list.erase(list.getSize() - 1);

    Vector<Element> new_send;
    for (UInt ns = 0; ns < list.getSize(); ++ns) {
      new_send.push_back(send(list(ns)));
    }

    AKANTU_DEBUG_INFO("I had " << send.getSize() << " elements to send to proc " << p << " and "
		      << list.getSize() << " elements to keep. I have "
		      << new_send.getSize() << " elements left.");
    send.copy(new_send);
  }

  comm.waitAll(isend_requests);
  comm.freeCommunicationRequest(isend_requests);

  delete [] list_of_el;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
