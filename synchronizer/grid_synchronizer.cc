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
#include "static_communicator.hh"

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

  Real * bounding_boxes = new Real[spatial_dimension * nb_proc];
  Real * my_bounding_box = bounding_boxes + spatial_dimension * my_rank;

  mesh.getLowerBounds(my_bounding_box);
  mesh.getUpperBounds(my_bounding_box + spatial_dimension);

  comm->allGather(bounding_boxes, 3);


  bool * intersects_proc = new bool[nb_proc];
  std::fill_n(intersects_proc, nb_proc, false);

  UInt * first_cells = new UInt[3 * nb_proc];
  UInt * last_cells = new UInt[3 * nb_proc];
  std::fill_n(first_cells, 3 * nb_proc, 0);
  std::fill_n(first_cells, 3 * nb_proc, 0);

  ByElementTypeUInt ** element_per_proc = new ByElementTypeUInt* [nb_proc];


  for (UInt p = 0; p < nb_proc; ++p) {
    if(p == my_rank) continue;

    Real * proc_bounding_box = bounding_boxes + spatial_dimension * p;

    bool intersects = false;
    UInt * first_cell = first_cells + p * spatial_dimension;
    UInt * last_cell =  last_cells  + p * spatial_dimension;
    for (UInt s = 0; s < spatial_dimension; ++s) {
      intersects = Math::intersects(my_bounding_box[s],
                                    my_bounding_box[spatial_dimension + s],
                                    proc_bounding_box[s],
                                    proc_bounding_box[spatial_dimension + s]);

      intersects_proc[p] &= intersects;

      if(intersects) {
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

        first_cell[s] = grid.getCell(start, s);
        last_cell [s] = grid.getCell(end, s);
      }
    }

    std::vector<UInt> * cells = new std::vector<UInt>();

    if(intersects_proc[p]) {
      UInt cell[3] = { 0 };
      if(first_cell[0] != 0) --first_cell[0];
      if(last_cell[0] != 0) ++last_cell[0];
      for (UInt fd = first_cell[0]; fd <= last_cell[0]; ++fd) {
        cell[0] = fd;

        if(first_cell[1] != 0) --first_cell[1];
        if(last_cell[1] != 0) ++last_cell[1];
        for (UInt sd = first_cell[1]; sd <= last_cell[1] ; ++sd) {
          cell[1] = sd;

          if(first_cell[2] != 0) --first_cell[2];
          if(last_cell[2] != 0) ++last_cell[2];
          for (UInt ld = first_cell[2]; fd <= last_cell[2] ; ++ld) {
            cell[2] = ld;
            cells->push_back(grid.getCell(cell));
          }
        }
      }

      std::vector<UInt>::iterator cur_cell = cells->begin();
      std::vector<UInt>::iterator last_cell = cells->end();
      std::set<E, CompElementLess> * to_send = new std::set<E, CompElementLess>();

      for (; cur_cell != last_cell; ++cur_cell) {
        typename RegularGrid<E>::const_iterator cur_elem = grid.beginCell(*cur_cell);
        typename RegularGrid<E>::const_iterator last_elem = grid.endCell(*cur_cell);

        for (; cur_elem != last_elem; ++cur_elem) {
          to_send->insert(*cur_elem);
        }
      }

      delete cells;

      std::stringstream sstr; sstr << "element_per_proc_" << p;
      element_per_proc[p] = new ByElementTypeUInt(sstr.str(), id);
      ByElementTypeUInt & elempproc = *(element_per_proc[p]);

      typename std::set<E>::iterator elem = to_send->begin();
      typename std::set<E>::iterator last_elem = to_send->end();
      for (; elem != last_elem; ++elem) {
        ElementType type = elem->type;
        UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);

        // /!\ this part is slow due to the access in the ByElementTypeUInt
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


  std::vector<CommunicationRequest *> isend_requests;
  for (UInt p = 0; p < nb_proc; ++p) {
    if(intersects_proc[p]) {
      for (UInt i = 0; i < N; ++i) {
        
      }

      UInt info[2];
      info[0] = (UInt) 
      isend_requests.push_back(static_communicator->asyncSend(info,
                                                              ssize,
                                                              p,
                                                              (Int) tag));

      isend_requests.push_back(static_communicator->asyncSend(buffer.storage(),
                                                              ssize,
                                                              p,
                                                              (Int) tag));

    }
  }



  AKANTU_DEBUG_OUT();
  return communicator;
}
/* -------------------------------------------------------------------------- */



__END_AKANTU__
