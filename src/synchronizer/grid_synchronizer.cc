/**
 * @file   grid_synchronizer.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Sep 20 10:30:37 2011
 *
 * @brief  implementation of the grid synchronizer
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
#include "grid_synchronizer.hh"
#include "aka_grid.hh"
#include "mesh.hh"
#include "fem.hh"
#include "static_communicator.hh"

#include "mesh_io.hh"


/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
GridSynchronizer::GridSynchronizer(const ID & id,
                                   MemoryID memory_id) :
  DistributedSynchronizer(id, memory_id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class E>
GridSynchronizer * GridSynchronizer::createGridSynchronizer(Mesh & mesh,
                                                            const RegularGrid<E> & grid,
                                                            SynchronizerID id,
                                                            MemoryID memory_id) {
  AKANTU_DEBUG_IN();

  StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator();
  UInt nb_proc = comm->getNbProc();
  UInt my_rank = comm->whoAmI();

  GridSynchronizer * communicator = new GridSynchronizer(id, memory_id);
  if(nb_proc == 1) return communicator;

  UInt spatial_dimension = mesh.getSpatialDimension();

  Real * bounding_boxes = new Real[2 * spatial_dimension * nb_proc];
  Real * my_bounding_box = bounding_boxes + 2 * spatial_dimension * my_rank;

  mesh.getLocalLowerBounds(my_bounding_box);
  mesh.getLocalUpperBounds(my_bounding_box + spatial_dimension);

  AKANTU_DEBUG_INFO("Exchange of bounding box to detect the overlapping regions.");

  comm->allGather(bounding_boxes, spatial_dimension * 2);

  bool * intersects_proc = new bool[nb_proc];
  std::fill_n(intersects_proc, nb_proc, true);

  UInt * first_cells = new UInt[3 * nb_proc];
  UInt * last_cells = new UInt[3 * nb_proc];
  std::fill_n(first_cells, 3 * nb_proc, 0);
  std::fill_n(first_cells, 3 * nb_proc, 0);

  ByElementTypeUInt ** element_per_proc = new ByElementTypeUInt* [nb_proc];
  for (UInt p = 0; p < nb_proc; ++p) element_per_proc[p] = NULL;

  for (UInt p = 0; p < nb_proc; ++p) {
    if(p == my_rank) continue;

    Real * proc_bounding_box = bounding_boxes + 2 * spatial_dimension * p;

    bool intersects = false;
    UInt * first_cell_p = first_cells + p * spatial_dimension;
    UInt * last_cell_p =  last_cells  + p * spatial_dimension;
    for (UInt s = 0; s < spatial_dimension; ++s) {
      intersects = Math::intersects(my_bounding_box[s],
                                    my_bounding_box[spatial_dimension + s],
                                    proc_bounding_box[s],
                                    proc_bounding_box[spatial_dimension + s]);

      intersects_proc[p] &= intersects;

      if(intersects) {
        AKANTU_DEBUG_INFO("I intersects with processor " << p << " in direction " << s);

        // is point 1 of proc p in the dimension s in the range ?
        bool point1 = Math::is_in_range(proc_bounding_box[s],
                                        my_bounding_box[s],
                                        my_bounding_box[s+spatial_dimension]);

        // is point 2 of proc p in the dimension s in the range ?
        bool point2 = Math::is_in_range(proc_bounding_box[s+spatial_dimension],
                                        my_bounding_box[s],
                                        my_bounding_box[s+spatial_dimension]);

        Real start = 0.;
        Real end = 0.;

        if(point1 && !point2) {
          /* |-----------|         my_bounding_box(i)
           *       |-----------|   proc_bounding_box(i)
           *       1           2
           */
          start = proc_bounding_box[s];
          end   = my_bounding_box[s+spatial_dimension];
        } else if(point1 && point2) {
          /* |-----------------|   my_bounding_box(i)
           *   |-----------|       proc_bounding_box(i)
           *   1           2
           */
          start = proc_bounding_box[s];
          end   = proc_bounding_box[s+spatial_dimension];
        } else if(!point1 && point2) {
          /*       |-----------|  my_bounding_box(i)
           * |-----------|        proc_bounding_box(i)
           * 1           2
           */
          start = my_bounding_box[s];
          end   = proc_bounding_box[s+spatial_dimension];
        }


        first_cell_p[s] = grid.getCell(start, s);
        last_cell_p [s] = grid.getCell(end, s);
      }
    }

    typedef typename RegularGrid<E>::Cell Cell;

    std::vector<Cell> * cells = new std::vector<Cell>;

    if(intersects_proc[p]) {
      AKANTU_DEBUG_INFO("I intersects with processor " << p);

      Cell cell(grid);

      for (UInt i = 0; i < spatial_dimension; ++i) {
        if(first_cell_p[i] != 0) --first_cell_p[i];
        if(last_cell_p[i] != 0) ++last_cell_p[i];
      }

      for (UInt fd = first_cell_p[0]; fd < last_cell_p[0]; ++fd) {
        cell.position[0] = fd;
        if(spatial_dimension == 1) {
	  cell.updateID();
          cells->push_back(cell);
	}
        else {
          for (UInt sd = first_cell_p[1]; sd < last_cell_p[1] ; ++sd) {
            cell.position[1] = sd;
            if(spatial_dimension == 2) {
	      cell.updateID();
              cells->push_back(cell);
            } else {
              for (UInt ld = first_cell_p[2]; fd < last_cell_p[2] ; ++ld) {
                cell.position[2] = ld;
		cell.updateID();
                cells->push_back(cell);
              }
            }
          }
        }
      }

      typename std::vector<Cell>::iterator cur_cell = cells->begin();
      typename std::vector<Cell>::iterator last_cell = cells->end();
      std::set<E, CompElementLess> * to_send = new std::set<E, CompElementLess>();

      for (; cur_cell != last_cell; ++cur_cell) {
        typename RegularGrid<E>::const_iterator cur_elem = grid.beginCell(*cur_cell);
        typename RegularGrid<E>::const_iterator last_elem = grid.endCell(*cur_cell);

        for (; cur_elem != last_elem; ++cur_elem) {
          to_send->insert(*cur_elem);
        }
      }

      AKANTU_DEBUG_INFO("I have prepared " << to_send->size() << " elements to send to processor " << p);

      delete cells;

      std::stringstream sstr; sstr << "element_per_proc_" << p;
      element_per_proc[p] = new ByElementTypeUInt(sstr.str(), id);
      ByElementTypeUInt & elempproc = *(element_per_proc[p]);

      typename std::set<E>::iterator elem = to_send->begin();
      typename std::set<E>::iterator last_elem = to_send->end();
      for (; elem != last_elem; ++elem) {
        ElementType type = elem->type;
        UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);

        // /!\ this part must be slow due to the access in the ByElementTypeUInt
        if(!elempproc.exists(type))
          elempproc.alloc(0, nb_nodes_per_element, type, _not_ghost);

        UInt global_connect[nb_nodes_per_element];
        UInt * local_connect = mesh.getConnectivity(type).storage() + elem->element * nb_nodes_per_element;
        for (UInt i = 0; i < nb_nodes_per_element; ++i) {
          global_connect[i] = mesh.getNodeGlobalId(local_connect[i]);
        }

        elempproc(type).push_back(global_connect);
        communicator->send_element[p].push_back(*elem);
      }

      delete to_send;
    }
  }

  AKANTU_DEBUG_INFO("I have finished to compute intersection,"
                    << " no it's time to communicate with my neighbors");



#define SIZE_TAG       0
#define DATA_TAG       1
#define ASK_NODES_TAG  2
#define SEND_NODES_TAG 3
  /* tag = |__________21_________|___8____|_3_|
   *       |          proc       | num mes| ct|
   */
#define GEN_TAG(proc, msg_count, tag) (((proc & 0x1FFFFF) << 11) + ((msg_count & 0xFF)<< 3) + (tag & 0x7))

  /**
   * Sending loop, sends the connectivity asynchronously to all concerned proc
   */
  std::vector<CommunicationRequest *> isend_requests;
  for (UInt p = 0; p < nb_proc; ++p) {
    if(p == my_rank) continue;
    if(intersects_proc[p]) {
      ByElementTypeUInt & elempproc = *(element_per_proc[p]);

      ByElementTypeUInt::type_iterator it_type   = elempproc.firstType(0, _not_ghost);
      ByElementTypeUInt::type_iterator last_type = elempproc.lastType(0, _not_ghost);

      UInt count = 0;
      for (; it_type != last_type; ++it_type) {
        Vector<UInt> & conn = elempproc(*it_type, _not_ghost);
        UInt info[2];
        info[0] = (UInt) *it_type;
        info[1] = conn.getSize() * conn.getNbComponent();

        AKANTU_DEBUG_INFO("I have " << conn.getSize() << " elements of type " << *it_type
                          << " to send to processor " << p
                          << " (communication tag : " << GEN_TAG(my_rank, count, DATA_TAG) << ")");

        isend_requests.push_back(comm->asyncSend(info, 2, p, GEN_TAG(my_rank, count, SIZE_TAG)));
        isend_requests.push_back(comm->asyncSend(conn.storage(),
                                                conn.getSize() * conn.getNbComponent(),
                                                p, GEN_TAG(my_rank, count, DATA_TAG)));

        ++count;
      }

      UInt info[2];
      info[0] = (UInt) _not_defined;
      info[1] = 0;
      isend_requests.push_back(comm->asyncSend(info, 2, p, GEN_TAG(my_rank, count, SIZE_TAG)));
    }
  }

  /**
   * Receives the connectivity and store them in the ghosts elements
   */
  Vector<UInt> global_nodes_ids = mesh.getGlobalNodesIds();
  std::vector<CommunicationRequest *> isend_nodes_requests;
  UInt nb_nodes_to_recv[nb_proc];
  UInt nb_total_nodes_to_recv = 0;
  UInt nb_current_nodes = global_nodes_ids.getSize();

  for (UInt p = 0; p < nb_proc; ++p) {
    nb_nodes_to_recv[p] = 0;
    if(p == my_rank) continue;

    Vector<UInt> ask_nodes;
    UInt count = 0;
    if(intersects_proc[p]) {
      ElementType type = _not_defined;
      do {
        UInt info[2];
        comm->receive(info, 2, p, GEN_TAG(p, count, SIZE_TAG));

        type = (ElementType) info[0];
        if(type != _not_defined) {
          UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);;
          UInt nb_element = info[1] / nb_nodes_per_element;

          Vector<UInt> tmp_conn(nb_element, nb_nodes_per_element);

          comm->receive(tmp_conn.storage(), info[1], p, GEN_TAG(p, count, DATA_TAG));

          AKANTU_DEBUG_INFO("I will receive " << nb_element << " elements of type " << ElementType(info[0])
                          << " from processor " << p
                          << " (communication tag : " << GEN_TAG(p, count, DATA_TAG) << ")");


          Vector<UInt> & ghost_connectivity = const_cast<Vector<UInt> &>(mesh.getConnectivity(type, _ghost));
          UInt nb_ghost_element = ghost_connectivity.getSize();
          Element element(type, 0, _ghost);

          UInt conn[nb_nodes_per_element];
          for (UInt el = 0; el < nb_element; ++el) {
            UInt nb_node_to_ask_for_elem = 0;

            for (UInt n = 0; n < nb_nodes_per_element; ++n) {
              UInt gn = tmp_conn(el, n);
              UInt ln = global_nodes_ids.find(gn);
              if(ln == UInt(-1)) {
                global_nodes_ids.push_back(gn);
                ln = nb_current_nodes;
                ++nb_current_nodes;
                ask_nodes.push_back(gn);
                ++nb_node_to_ask_for_elem;
              }

              conn[n] = ln;
            }

            // all the nodes are already known locally, the element should already exists
            UInt c = UInt(-1);
            if(nb_node_to_ask_for_elem != 0) {
              c = ghost_connectivity.find(conn);
            }

            if(c == UInt(-1)) {
              element.element = nb_ghost_element;
              ++nb_ghost_element;
              ghost_connectivity.push_back(conn);
            }
            communicator->recv_element[p].push_back(element);
          }
        }
        count++;
      } while(type != _not_defined);

      AKANTU_DEBUG_INFO("I have " << ask_nodes.getSize()
                        << " missing nodes for elements coming from processor " << p
                        << " (communication tag : " << GEN_TAG(my_rank, 0, ASK_NODES_TAG) << ")");


      isend_nodes_requests.push_back(comm->asyncSend(ask_nodes.storage(), ask_nodes.getSize(),
                                                     p, GEN_TAG(my_rank, 0, ASK_NODES_TAG)));
      nb_nodes_to_recv[p] = ask_nodes.getSize();
      nb_total_nodes_to_recv += ask_nodes.getSize();
    }
  }


  Vector<Real> & nodes = const_cast<Vector<Real> &>(mesh.getNodes());
  UInt nb_nodes = nodes.getSize();
  std::vector<CommunicationRequest *> irecv_nodes_requests;
  nodes.resize(nb_total_nodes_to_recv + nb_nodes);
  for (UInt p = 0; p < nb_proc; ++p) {
    if((p != my_rank) && (nb_nodes_to_recv[p] > 0)) {
      irecv_nodes_requests.push_back(comm->asyncReceive(nodes.storage() + nb_nodes * spatial_dimension,
                                                        nb_nodes_to_recv[p] * spatial_dimension,
                                                        p, GEN_TAG(p, 0, SEND_NODES_TAG)));
      nb_nodes += nb_nodes_to_recv[p];
    }
  }


  comm->waitAll(isend_requests);
  comm->freeCommunicationRequest(isend_requests);

  for (UInt p = 0; p < nb_proc; ++p) {
    if(element_per_proc[p]) delete element_per_proc[p];
  }
  delete [] element_per_proc;

  /**
   * Sends requested nodes to proc
   */
  for (UInt p = 0; p < nb_proc; ++p) {
    if(p == my_rank) continue;

    Vector<UInt> asked_nodes;
    CommunicationStatus status;
    comm->probe<UInt>(p, GEN_TAG(p, 0, ASK_NODES_TAG), status);
    UInt nb_nodes_to_send = status.getSize();
    asked_nodes.resize(nb_nodes_to_send);

    AKANTU_DEBUG_INFO("I have " << nb_nodes_to_send
                      << " to send to processor " << p
                      << "(communication tag : " << GEN_TAG(p, 0, ASK_NODES_TAG) << ")");

    comm->receive(asked_nodes.storage(), nb_nodes_to_send, p, GEN_TAG(p, 0, ASK_NODES_TAG));
    Vector<Real> nodes_to_send(0, spatial_dimension);
    for (UInt n = 0; n < nb_nodes_to_send; ++n) {
      UInt ln = global_nodes_ids.find(asked_nodes(n));
      nodes_to_send.push_back(nodes.storage() + ln * spatial_dimension);
    }
    comm->send(nodes_to_send.storage(), nb_nodes_to_send * spatial_dimension, p, GEN_TAG(my_rank, 0, SEND_NODES_TAG));
  }

  comm->waitAll(isend_nodes_requests);
  comm->freeCommunicationRequest(isend_nodes_requests);

  comm->waitAll(irecv_nodes_requests);
  comm->freeCommunicationRequest(irecv_nodes_requests);


  AKANTU_DEBUG_OUT();
  return communicator;
}
/* -------------------------------------------------------------------------- */


template GridSynchronizer *
GridSynchronizer::createGridSynchronizer<QuadraturePoint>(Mesh & mesh,
                                                          const RegularGrid<QuadraturePoint> & grid,
                                                          SynchronizerID id,
                                                          MemoryID memory_id);
template GridSynchronizer *
GridSynchronizer::createGridSynchronizer<Element>(Mesh & mesh,
                                                  const RegularGrid<Element> & grid,
                                                  SynchronizerID id,
                                                  MemoryID memory_id);

__END_AKANTU__
