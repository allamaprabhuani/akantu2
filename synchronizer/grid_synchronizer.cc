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
#include "mesh.hh"
#include "aka_grid.hh"
#include "static_communicator.hh"

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
GridSynchronizer::GridSynchronizer(ID & id = "grid_synchronizer",
				   MemoryID memory_id = 0) :
  Synchronizer(id, memory_id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
GridSynchronizer * GridSynchronizer::createGridSynchronizer(Mesh & mesh,
							    const RegularGrid & grid,
							    SynchronizerID id = "grid_synchronizer",
							    MemoryID memory_id = 0) {
  AKANTU_DEBUG_IN();

  StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator();
  UInt nb_proc = comm->getNbProc();
  UInt my_rank = comm->whoAmI();

  GridSynchronizer * synchronizer = new GridSynchronizer(id, memory_id);

  if(nb_proc == 1) return synchronizer;
  UInt spatial_dimension = mesh.getSpatialDimension();

  Real * bounding_boxes = new Real[spatial_dimension * nb_proc];
  Real * my_bounding_box = bounding_boxes + spatial_dimension * my_rank;

  mesh.getLowerBounds(my_bounding_box);
  mesh.getUpperBounds(my_bounding_box + spatial_dimension);

  comm->allGather(bounding_boxes, 3);


  bool intersects_proc = new bool[nb_proc];
  std::fill_n(intersects_proc, nb_proc, false);

  UInt first_cells = new UInt[spatial_dimension * nb_proc];
  UInt last_cells = new UInt[spatial_dimension * nb_proc];

  for (UInt p = 0; p < nb_proc; ++p) {
    if(p == my_rank) continue;

    Real * proc_bounding_box = bounding_boxes + spatial_dimension * p;

    bool intersects = false;
    UInt first_cell = first_cells + p * spatial_dimension;
    UInt last_cell =  last_cells  + p * spatial_dimension;
    for (UInt s = 0; s < spatial_dimension; ++s) {
      intersects = Math::instersects(my_bounding_box[s],
				     my_bounding_box[spatial_dimension + s],
				     proc_bounding_box[s],
				     proc_bounding_box[spatial_dimension + s]);

      intersects_proc[p] *= intersects;

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

	first_cell[s] = grid.getCell(start);
	last_cell [s] = grid.getCell(end);
      }
    }

    

  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */



__END_AKANTU__
