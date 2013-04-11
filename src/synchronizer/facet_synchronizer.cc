/**
 * @file   facet_synchronizer.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Tue Mar 26 09:55:38 2013
 *
 * @brief  Facet synchronizer for parallel simulations with cohesive elments
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

#include "facet_synchronizer.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
FacetSynchronizer::FacetSynchronizer(DistributedSynchronizer & distributed_synchronizer,
				     Mesh & mesh_facets,
				     SynchronizerID id,
				     MemoryID memory_id) :
  DistributedSynchronizer(mesh_facets, id, memory_id),
  distributed_synchronizer(distributed_synchronizer),
  mesh_facets(mesh_facets) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
FacetSynchronizer * FacetSynchronizer::
createFacetSynchronizer(DistributedSynchronizer & distributed_synchronizer,
			Mesh & mesh_facets,
			SynchronizerID id,
			MemoryID memory_id) {
  AKANTU_DEBUG_IN();

  FacetSynchronizer & f_synchronizer = *(new FacetSynchronizer(distributed_synchronizer,
							       mesh_facets,
							       id, memory_id));

  f_synchronizer.setupFacetSynchronization();

  AKANTU_DEBUG_OUT();
  return &f_synchronizer;
}

/* -------------------------------------------------------------------------- */
void FacetSynchronizer::updateDistributedSynchronizer(DataAccessor & data_accessor,
						      const ByElementTypeUInt & cohesive_el_to_facet) {
  AKANTU_DEBUG_IN();

  Array<Element> * distrib_send_element = distributed_synchronizer.send_element;
  Array<Element> * distrib_recv_element = distributed_synchronizer.recv_element;

  updateElementList(distrib_send_element, send_element, cohesive_el_to_facet);
  updateElementList(distrib_recv_element, recv_element, cohesive_el_to_facet);

  std::map<SynchronizationTag, Communication>::iterator it
    = distributed_synchronizer.communications.begin();

  std::map<SynchronizationTag, Communication>::iterator end
    = distributed_synchronizer.communications.end();

  for (; it != end; ++it) {
    SynchronizationTag tag = it->first;
    distributed_synchronizer.computeBufferSize(data_accessor, tag);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FacetSynchronizer::updateElementList(Array<Element> * elements,
					  const Array<Element> * facets,
					  const ByElementTypeUInt & cohesive_el_to_facet) {
  AKANTU_DEBUG_IN();

  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;

  ElementType current_coh_element_type = _not_defined;

  UInt max_uint = std::numeric_limits<UInt>::max();
  const Array<UInt> * cohesive_el_to_f = NULL;

  Element cohesive_element(_not_defined, 0, _not_ghost, _ek_cohesive);

  for (UInt p = 0; p < nb_proc; ++p) {
    const Array<Element> & fa = facets[p];
    Array<Element> & el = elements[p];

    Array<Element>::const_iterator<Element> it = fa.begin();
    Array<Element>::const_iterator<Element> end = fa.end();

    for (; it != end; ++it) {
      const Element & facet = *it;

      if(facet.type != current_element_type || facet.ghost_type != current_ghost_type) {
	current_element_type = facet.type;
	current_coh_element_type = FEM::getCohesiveElementType(current_element_type);
	current_ghost_type   = facet.ghost_type;
	cohesive_el_to_f = &(cohesive_el_to_facet(current_element_type,
						  current_ghost_type));
	cohesive_element.type = current_coh_element_type;
	cohesive_element.ghost_type = current_ghost_type;
      }

      cohesive_element.element = (*cohesive_el_to_f)(facet.element);

      if (cohesive_element.element != max_uint)
	el.push_back(cohesive_element);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FacetSynchronizer::setupFacetSynchronization() {
  AKANTU_DEBUG_IN();

  Array<Element> * distrib_send_element = distributed_synchronizer.send_element;
  Array<Element> * distrib_recv_element = distributed_synchronizer.recv_element;

  /// build rank to facet correspondance
  ByElementTypeUInt rank_to_facet("rank_to_facet", id);
  initRankToFacet(rank_to_facet);
  buildRankToFacet(rank_to_facet, distrib_recv_element);

  /// allocate send/recv barycenters arrays
  std::list<ElementBarycenter> * send_elbary =
    new std::list<ElementBarycenter>[nb_proc];

  std::list<ElementBarycenter> * recv_elbary =
    new std::list<ElementBarycenter>[nb_proc];

  /// compute barycenters
  getFacetBarycentersPerElement<_not_ghost>(rank_to_facet,
					    distrib_send_element,
					    send_elbary);

  getFacetBarycentersPerElement<_ghost>(rank_to_facet,
					distrib_recv_element,
					recv_elbary);

  /// build send/recv facet arrays
  buildFacetList(recv_elbary, recv_element);
  buildFacetList(send_elbary, send_element);

  /// delete temporary data
  delete [] recv_elbary;
  delete [] send_elbary;


#ifndef AKANTU_NDEBUG
  /// count recv facets for each processor
  Array<UInt> nb_facets_recv(nb_proc);
  nb_facets_recv.clear();

  UInt spatial_dimension = mesh_facets.getSpatialDimension();

  Mesh::type_iterator first = mesh_facets.firstType(spatial_dimension - 1, _ghost);
  Mesh::type_iterator last  = mesh_facets.lastType(spatial_dimension - 1, _ghost);

  for (; first != last; ++first) {
    const Array<UInt> & r_to_f = rank_to_facet(*first, _ghost);
    UInt nb_facet = r_to_f.getSize();

    for (UInt f = 0; f < nb_facet; ++f) {
      UInt proc = r_to_f(f);
      if (proc != rank)
	++nb_facets_recv(proc);
    }
  }

  for (UInt p = 0; p < nb_proc; ++p) {
    AKANTU_DEBUG_ASSERT(nb_facets_recv(p) == recv_element[p].getSize(),
			"Wrong number of recv facets");
  }

#endif

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FacetSynchronizer::buildFacetList(std::list<ElementBarycenter> * element_barycenter,
				       Array<Element> * final_facets) {
  AKANTU_DEBUG_IN();

  for (UInt p = 0; p < nb_proc; ++p) {
    if (p == rank) continue;

    std::list<ElementBarycenter> & elbary = element_barycenter[p];
    if (elbary.size() == 0) continue;

    Array<Element> & final = final_facets[p];

    elbary.sort();

    /// introduce lists' elements in the array
    std::list<ElementBarycenter>::iterator elbary_pre = elbary.begin();
    std::list<ElementBarycenter>::iterator elbary_it = elbary.begin();
    std::list<ElementBarycenter>::iterator elbary_end = elbary.end();

    final.push_back((*elbary_it).elem);

    for (++elbary_it; elbary_it != elbary_end; ++elbary_it, ++elbary_pre) {
      /// skip double facets
      if ( (*elbary_pre).bary == (*elbary_it).bary ) continue;

      final.push_back((*elbary_it).elem);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FacetSynchronizer::initRankToFacet(ByElementTypeUInt & rank_to_facet) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh_facets.getSpatialDimension();

  mesh_facets.initByElementTypeArray(rank_to_facet, 1, spatial_dimension - 1);

  GhostType gt = _ghost;

  Mesh::type_iterator first = mesh_facets.firstType(spatial_dimension - 1, gt);
  Mesh::type_iterator last  = mesh_facets.lastType(spatial_dimension - 1, gt);

  for (; first != last; ++first) {
    Array<UInt> & rank_to_f = rank_to_facet(*first, gt);
    UInt nb_facet = mesh_facets.getNbElement(*first, gt);
    rank_to_f.resize(nb_facet);

    for (UInt f = 0; f < nb_facet; ++f)
      rank_to_f(f) = rank;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FacetSynchronizer::buildRankToFacet(ByElementTypeUInt & rank_to_facet,
					 const Array<Element> * elements) {
  AKANTU_DEBUG_IN();

  for (UInt p = 0; p < nb_proc; ++p) {
    if (p == rank) continue;
    const Array<Element> & elem = elements[p];
    UInt nb_element = elem.getSize();

    for (UInt el = 0; el < nb_element; ++el) {
      ElementType type = elem(el).type;
      GhostType gt = elem(el).ghost_type;
      UInt el_index = elem(el).element;

      const Array<Element> & facet_to_element =
	mesh_facets.getSubelementToElement(type, gt);
      UInt nb_facets_per_element = Mesh::getNbFacetsPerElement(type);
      ElementType facet_type = Mesh::getFacetType(type);

      for (UInt f = 0; f < nb_facets_per_element; ++f) {
	UInt facet_index = facet_to_element(el_index, f).element;
	GhostType facet_gt = facet_to_element(el_index, f).ghost_type;

	if (facet_gt == _not_ghost) continue;

	Array<UInt> & t_to_f = rank_to_facet(facet_type, facet_gt);
	if ((p < t_to_f(facet_index)) || (t_to_f(facet_index) == rank))
	  t_to_f(facet_index) = p;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
