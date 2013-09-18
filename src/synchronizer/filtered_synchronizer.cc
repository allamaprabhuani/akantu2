/**
 * @file   filtered_synchronizer.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Sep 17 11:17:30 2013
 *
 * @brief  filtered synchronizer
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
#include "filtered_synchronizer.hh"

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
FilteredSynchronizer::FilteredSynchronizer(Mesh & mesh,
					   SynchronizerID id,
					   MemoryID memory_id) :
  DistributedSynchronizer(mesh, id, memory_id) {
  AKANTU_DEBUG_IN();
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
FilteredSynchronizer * FilteredSynchronizer::
createFilteredSynchronizer(const DistributedSynchronizer & d_synchronizer,
			   SynchElementFilter & filter) {
  AKANTU_DEBUG_IN();
  
  FilteredSynchronizer & f_synchronizer = 
    *(new FilteredSynchronizer(d_synchronizer.mesh,
			       d_synchronizer.id + ":filtered",
			       d_synchronizer.memory_id));
  
  f_synchronizer.setupSynchronizer(d_synchronizer,
				   filter);

  AKANTU_DEBUG_OUT();
  return &f_synchronizer;
}

/* -------------------------------------------------------------------------- */
void FilteredSynchronizer::setupSynchronizer(const DistributedSynchronizer & d_synchronizer,
					     SynchElementFilter & filter) {
  AKANTU_DEBUG_IN();
  
  Array<Element> * d_send_element = d_synchronizer.send_element;
  Array<Element> * d_recv_element = d_synchronizer.recv_element;

  this->updateElementList(d_send_element,
			  this->send_element,
			  filter);
  this->updateElementList(d_recv_element,
			  this->recv_element,
			  filter);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FilteredSynchronizer::updateElementList(Array<Element> * source_elements,
					     Array<Element> * destination_elements,
					     SynchElementFilter & filter) {
  AKANTU_DEBUG_IN();
  
  // loop over procs
  for (UInt p = 0; p < this->nb_proc; ++p) {
    if (p == this->rank) continue;

    // access the element for this proc
    const Array<Element> & unfiltered_elements = source_elements[p];
    Array<Element> & filtered_elements = destination_elements[p];

    // iterator to loop over all source elements
    Array<Element>::const_iterator<Element> it  = unfiltered_elements.begin();
    Array<Element>::const_iterator<Element> end = unfiltered_elements.end();
    
    // if filter accepts this element, push it into the destination elements
    for (; it != end; ++it) {
      const Element & element = *it;
      if (filter(element)) {
	filtered_elements.push_back(element);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
