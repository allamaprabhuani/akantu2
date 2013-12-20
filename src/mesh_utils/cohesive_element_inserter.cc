/**
 * @file   cohesive_element_inserter.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Fri Nov 29 17:16:03 2013
 *
 * @brief  Cohesive element inserter functions
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
#include <algorithm>
#include <limits>
#include "cohesive_element_inserter.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

CohesiveElementInserter::CohesiveElementInserter(Mesh & mesh,
						 Mesh & mesh_facets,
						 const ID & id) :
  id(id),
  mesh(mesh),
  mesh_facets(&mesh_facets),
  insertion_facets("insertion_facets", id),
  is_extrinsic(true),
  insertion_limits(mesh.getSpatialDimension(), 2),
  check_facets("check_facets", id) { }

/* -------------------------------------------------------------------------- */
CohesiveElementInserter::CohesiveElementInserter(Mesh & mesh,
						 const ID & id) :
  id(id),
  mesh(mesh),
  insertion_facets("insertion_facets", id),
  is_extrinsic(false),
  insertion_limits(mesh.getSpatialDimension(), 2),
  check_facets("check_facets", id) {

  mesh_facets = new Mesh(mesh.initMeshFacets());
  MeshUtils::buildAllFacets(mesh, *mesh_facets);

  init();
}

/* -------------------------------------------------------------------------- */
CohesiveElementInserter::~CohesiveElementInserter() {
  if (!is_extrinsic) delete mesh_facets;

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  delete distributed_synchronizer;
#endif
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::init() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  MeshUtils::resetFacetToDouble(*mesh_facets);

  /// initialize facet insertion array
  mesh_facets->initByElementTypeArray(insertion_facets, 1,
				      spatial_dimension - 1,
				      false,
				      _ek_regular,
				      true);

  /// init insertion limits
  for (UInt dim = 0; dim < spatial_dimension; ++dim) {
    insertion_limits(dim, 0) = std::numeric_limits<Real>::max() * (-1.);
    insertion_limits(dim, 1) = std::numeric_limits<Real>::max();
  }

  if (is_extrinsic) {
    mesh_facets->initByElementTypeArray(check_facets, 1, spatial_dimension - 1);
    initFacetsCheck();
  }

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  facet_synchronizer = NULL;
  distributed_synchronizer = NULL;
#endif

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::initFacetsCheck() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType facet_gt = *gt;
    Mesh::type_iterator it   = mesh_facets->firstType(spatial_dimension - 1, facet_gt);
    Mesh::type_iterator last = mesh_facets->lastType(spatial_dimension - 1, facet_gt);

    for (; it != last; ++it) {
      ElementType facet_type = *it;

      Array<bool> & f_check = check_facets(facet_type, facet_gt);

      const Array< std::vector<Element> > & element_to_facet
	= mesh_facets->getElementToSubelement(facet_type, facet_gt);

      UInt nb_facet = element_to_facet.getSize();
      f_check.resize(nb_facet);

      for (UInt f = 0; f < nb_facet; ++f) {
	if (element_to_facet(f)[1] == ElementNull ||
	    (element_to_facet(f)[0].ghost_type == _ghost &&
	     element_to_facet(f)[1].ghost_type == _ghost)) {
	  f_check(f) = false;
	}
	else f_check(f) = true;
      }
    }
  }

 AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::limitCheckFacets() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  Vector<Real> bary_facet(spatial_dimension);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end();
       ++gt) {
    GhostType ghost_type = *gt;

    Mesh::type_iterator it  = mesh_facets->firstType(spatial_dimension - 1, ghost_type);
    Mesh::type_iterator end = mesh_facets->lastType(spatial_dimension - 1, ghost_type);
    for(; it != end; ++it) {
      ElementType type = *it;
      Array<bool> & f_check = check_facets(type, ghost_type);
      UInt nb_facet = mesh_facets->getNbElement(type, ghost_type);

      for (UInt f = 0; f < nb_facet; ++f) {
	if (f_check(f)) {

	  mesh_facets->getBarycenter(f, type, bary_facet.storage(), ghost_type);

	  UInt coord_in_limit = 0;

	  while (coord_in_limit < spatial_dimension &&
		 bary_facet(coord_in_limit) > insertion_limits(coord_in_limit, 0) &&
		 bary_facet(coord_in_limit) < insertion_limits(coord_in_limit, 1))
	    ++coord_in_limit;

	  if (coord_in_limit != spatial_dimension)
	    f_check(f) = false;

	}
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::setLimit(char direction,
				       Real first_limit,
				       Real second_limit) {
  UInt direction_index;

  switch(direction) {
  case 'x':
    direction_index = 0; break;
  case 'y':
    direction_index = 1; break;
  case 'z':
    direction_index = 2; break;
  default:
    /// assign dummy value
    direction_index = 10;
  }

  AKANTU_DEBUG_ASSERT(direction_index < mesh.getSpatialDimension(),
		      "Specify the direction as 'x', 'y' or 'z' according to the mesh spatial dimension");

  insertion_limits(direction_index, 0) = std::min(first_limit, second_limit);
  insertion_limits(direction_index, 1) = std::max(first_limit, second_limit);
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::insertIntrinsicElements() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  Vector<Real> bary_facet(spatial_dimension);

  Mesh::type_iterator it  = mesh_facets->firstType(spatial_dimension - 1);
  Mesh::type_iterator end = mesh_facets->lastType(spatial_dimension - 1);

  for(; it != end; ++it) {
    const ElementType type_facet = *it;
    Array<bool> & f_insertion = insertion_facets(type_facet);
    Array<std::vector<Element> > & element_to_facet
      = mesh_facets->getElementToSubelement(type_facet);

    UInt nb_facet = mesh_facets->getNbElement(type_facet);

    for (UInt f = 0; f < nb_facet; ++f) {

      if (element_to_facet(f)[1] == ElementNull) continue;

      mesh_facets->getBarycenter(f, type_facet, bary_facet.storage());

      UInt coord_in_limit = 0;

      while (coord_in_limit < spatial_dimension &&
	     bary_facet(coord_in_limit) > insertion_limits(coord_in_limit, 0) &&
	     bary_facet(coord_in_limit) < insertion_limits(coord_in_limit, 1))
	++coord_in_limit;

      if (coord_in_limit == spatial_dimension)
	f_insertion(f) = true;
    }
  }

  NewNodesEvent node_event;
  node_event.getList().extendComponentsInterlaced(2, 1);
  NewElementsEvent element_event;

  MeshUtils::insertCohesiveElements(mesh,
				    *mesh_facets,
				    insertion_facets,
				    node_event.getList(),
				    element_event.getList());

  UInt nb_new_nodes = node_event.getList().getSize();

  mesh.nb_global_nodes += nb_new_nodes;
  mesh_facets->nb_global_nodes += nb_new_nodes;

  mesh.sendEvent(node_event);
  mesh.updateTypesOffsets(_not_ghost);
  mesh.sendEvent(element_event);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::insertExtrinsicElements() {
  AKANTU_DEBUG_IN();

  NewNodesEvent node_event;
  node_event.getList().extendComponentsInterlaced(2, 1);
  NewElementsEvent element_event;

  MeshUtils::insertCohesiveElements(mesh,
				    *mesh_facets,
				    insertion_facets,
				    node_event.getList(),
				    element_event.getList());

  UInt nb_new_nodes = node_event.getList().getSize();
  UInt nb_new_elements = element_event.getList().getSize();

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  if (mesh.nodes_type) {
    /// update global ids
    nb_new_nodes = updateGlobalIDs(node_event);

    /// compute total number of new elements
    StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
    comm.allReduce(&nb_new_elements, 1, _so_sum);
  }
#endif

  if (nb_new_nodes > 0) {
    mesh.nb_global_nodes += nb_new_nodes;
    mesh_facets->nb_global_nodes += nb_new_nodes;
    mesh.sendEvent(node_event);
  }

  if (nb_new_elements > 0) {
    updateInsertionFacets();
    mesh.updateTypesOffsets(_not_ghost);
    mesh.sendEvent(element_event);
    MeshUtils::resetFacetToDouble(*mesh_facets);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::updateInsertionFacets() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType facet_gt = *gt;
    Mesh::type_iterator it   = mesh_facets->firstType(spatial_dimension - 1, facet_gt);
    Mesh::type_iterator last = mesh_facets->lastType(spatial_dimension - 1, facet_gt);

    for (; it != last; ++it) {
      ElementType facet_type = *it;

      Array<bool> & f_check = check_facets(facet_type, facet_gt);
      Array<bool> & ins_facets = insertion_facets(facet_type, facet_gt);

      UInt nb_facets = f_check.getSize();

      for (UInt f = 0; f < ins_facets.getSize(); ++f) {
	if (ins_facets(f)) {
	  ++nb_facets;
	  ins_facets(f) = false;
	  f_check(f) = false;
	}
      }

      f_check.resize(nb_facets);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "CohesiveElementInserter [" << std::endl;

  stream << space << " + mesh [" << std::endl;
  mesh.printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << " + mesh_facets [" << std::endl;
  mesh_facets->printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << " + is_extrinsic : " << is_extrinsic << std::endl;

  stream << space << "]" << std::endl;
}



__END_AKANTU__
