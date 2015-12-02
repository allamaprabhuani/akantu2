/**
 * @file   neighborhood_max_criterion.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Wed Oct 14 21:31:07 2015
 *
 * @brief  Implementation of class NeighborhoodMaxCriterion
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
#include "neighborhood_max_criterion.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
NeighborhoodMaxCriterion::NeighborhoodMaxCriterion(const SolidMechanicsModel & model, 
						   const ElementTypeMapReal & quad_coordinates,
						   const ID & criterion_id,
						   const ID & id,
						   const MemoryID & memory_id)  :
  NeighborhoodBase(model, quad_coordinates, id, memory_id),
  Parsable(_st_non_local, id),
  is_highest("is_highest", id),
  criterion(criterion_id, id)
 { 
   
   AKANTU_DEBUG_IN();
   
   this->registerParam("radius"       , neighborhood_radius             , 100.,
		      _pat_parsable | _pat_readable  , "Non local radius");
  
   Mesh & mesh = this->model.getMesh();
   /// allocate the element type map arrays for _not_ghosts: One entry per quad
   GhostType ghost_type = _not_ghost;
   Mesh::type_iterator it        = mesh.firstType(spatial_dimension, ghost_type);
   Mesh::type_iterator last_type = mesh.lastType (spatial_dimension, ghost_type);
   
   for(; it != last_type; ++it) {
     UInt new_size = this->quad_coordinates(*it, ghost_type).getSize();
     this->is_highest.alloc(new_size, 1, *it, ghost_type, true);
     this->criterion.alloc(new_size, 1, *it, ghost_type, true);
   }
   
   /// criterion needs allocation also for ghost
   ghost_type = _ghost;
   it = mesh.firstType(spatial_dimension, ghost_type);
   last_type = mesh.lastType(spatial_dimension, ghost_type);
   for(; it != last_type; ++it) {
     UInt new_size = this->quad_coordinates(*it, ghost_type).getSize();
     this->criterion.alloc(new_size, 1, *it, ghost_type, true);
   }

     AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
NeighborhoodMaxCriterion::~NeighborhoodMaxCriterion() {
  AKANTU_DEBUG_IN();


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NeighborhoodMaxCriterion::initNeighborhood() {
  AKANTU_DEBUG_IN();

  /// parse the input parameter
  const Parser & parser = getStaticParser();
  const ParserSection & section_neighborhood = *(parser.getSubSections(_st_neighborhood).first);
  this->parseSection(section_neighborhood);

  AKANTU_DEBUG_INFO("Creating the grid");
  this->createGrid();

  /// insert the non-ghost quads into the grid
  this->insertAllQuads(_not_ghost);

  /// store the number of current ghost elements for each type in the mesh
  ElementTypeMap<UInt> nb_ghost_protected;
  Mesh & mesh = this->model.getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, _ghost);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _ghost);
  for(; it != last_type; ++it)
    nb_ghost_protected(mesh.getNbElement(*it, _ghost), *it, _ghost);

  /// create the grid synchronizer
  this->createGridSynchronizer();

  /// insert the ghost quads into the grid
  this->insertAllQuads(_ghost);

  /// create the pair lists
  this->updatePairList();
  
  /// remove the unneccessary ghosts
  this->cleanupExtraGhostElements(nb_ghost_protected);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NeighborhoodMaxCriterion::createGridSynchronizer() {
  this->is_creating_grid = true;
  std::set<SynchronizationTag> tags;
  tags.insert(_gst_nh_criterion);

  std::stringstream sstr; sstr << getID() << ":grid_synchronizer";
  this->grid_synchronizer = GridSynchronizer::createGridSynchronizer(this->model.getMesh(),
								     *spatial_grid,
								     sstr.str(),
								     synch_registry,
								     tags, 0, false);
  this->is_creating_grid = false;
}

/* -------------------------------------------------------------------------- */
void NeighborhoodMaxCriterion::insertAllQuads(const GhostType & ghost_type) {
  IntegrationPoint q;
  q.ghost_type = ghost_type;
  Mesh & mesh = this->model.getMesh();

  Mesh::type_iterator it        = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType (spatial_dimension, ghost_type);
  for(; it != last_type; ++it) {
    UInt nb_element = mesh.getNbElement(*it, ghost_type);
    UInt nb_quad    = this->model.getFEEngine().getNbIntegrationPoints(*it, ghost_type);

    const Array<Real> & quads = this->quad_coordinates(*it, ghost_type);
    q.type = *it;

    Array<Real>::const_vector_iterator quad = quads.begin(spatial_dimension);

    for (UInt e = 0; e < nb_element; ++e) {
      q.element = e;
      for (UInt nq = 0; nq < nb_quad; ++nq) {
	q.num_point = nq;
	q.global_num = q.element * nb_quad + nq;
	spatial_grid->insert(q, *quad);
	++quad;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void NeighborhoodMaxCriterion::findMaxQuads(std::vector<IntegrationPoint> & max_quads) {
  AKANTU_DEBUG_IN();

  /// clear the element type maps
  this->is_highest.clear();
  this->criterion.clear();
  /// update the values of the criterion
  const ID field_name = criterion.getName();
  for (UInt m = 0; m < this->model.getNbMaterials(); ++m) {
    const Material & material = this->model.getMaterial(m);
    if (material.isInternal<Real>(field_name, _ek_regular))
      for(UInt g = _not_ghost; g <= _ghost; ++g) {
	GhostType ghost_type = (GhostType) g;
	material.flattenInternal(field_name, criterion, ghost_type, _ek_regular);
      }
  }

  /// start the exchange the value of the criterion on the ghost elements
  SynchronizerRegistry & synch_registry = this->model.getSynchronizerRegistry();
  synch_registry.asynchronousSynchronize(_gst_nh_criterion);

  /// compare to not-ghost neighbors
  checkNeighbors(_not_ghost);

  /// finish the exchange
  synch_registry.waitEndSynchronize(_gst_nh_criterion);

  /// compare to ghost neighbors
  checkNeighbors(_ghost);


  /// extract the quads with highest criterion in their neighborhood
  IntegrationPoint quad; 
  quad.ghost_type = _not_ghost;
  Mesh & mesh = this->model.getMesh();
  Mesh::type_iterator it        = mesh.firstType(spatial_dimension, _not_ghost);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _not_ghost);
  for(; it != last_type; ++it) {
    quad.type = *it;
    Array<bool>::const_scalar_iterator is_highest_it  = is_highest(*it, _not_ghost).begin();
    Array<bool>::const_scalar_iterator is_highest_end = is_highest(*it, _not_ghost).end();

    UInt nb_quadrature_points = this->model.getFEEngine().getNbIntegrationPoints(*it, _not_ghost);
    UInt q = 0;

    /// loop over is_highest for the current element type
    for (;is_highest_it != is_highest_end; ++is_highest_it, ++q) {
      if (*is_highest_it) {
	/// gauss point has the highest stress in his neighbourhood
	quad.element = q / nb_quadrature_points;
	quad.global_num = q;
	quad.num_point = q % nb_quadrature_points;
	max_quads.push_back(quad);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NeighborhoodMaxCriterion::checkNeighbors(const GhostType & ghost_type2) {                               
  AKANTU_DEBUG_IN();

  PairList::const_iterator first_pair = pair_list[ghost_type2].begin();
  PairList::const_iterator last_pair  = pair_list[ghost_type2].end();

  // Compute the weights
  for(;first_pair != last_pair; ++first_pair) {
    const IntegrationPoint & lq1 = first_pair->first;
    const IntegrationPoint & lq2 = first_pair->second;

    Array<bool> & has_highest_eq_stress_1 = is_highest(lq1.type, lq1.ghost_type);

    const Array<Real> & criterion_1 = this->criterion(lq1.type, lq1.ghost_type);
    const Array<Real> & criterion_2 = this->criterion(lq2.type, lq2.ghost_type);

    if(criterion_1(lq1.global_num) < criterion_2(lq2.global_num))
      has_highest_eq_stress_1(lq1.global_num) = false;
    else if (ghost_type2 != _ghost) {
      Array<bool> & has_highest_eq_stress_2 = is_highest(lq2.type, lq2.ghost_type);
      has_highest_eq_stress_2(lq2.global_num) = false;
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NeighborhoodMaxCriterion::cleanupExtraGhostElements(const ElementTypeMap<UInt> & nb_ghost_protected) {
                             
  Mesh & mesh = this->model.getMesh();
  /// create remove elements event
  RemovedElementsEvent remove_elem(mesh);
  /// create set of ghosts to keep
  std::set<Element> relevant_ghost_elements;
  PairList::const_iterator first_pair = pair_list[_ghost].begin();
  PairList::const_iterator last_pair  = pair_list[_ghost].end();
  for(;first_pair != last_pair; ++first_pair) {
    const IntegrationPoint & q2 = first_pair->second;
    relevant_ghost_elements.insert(q2);
  }

  Array<Element> ghosts_to_erase(0);
  Mesh::type_iterator it        = mesh.firstType(spatial_dimension, _ghost);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _ghost);

  Element element;  
  element.ghost_type = _ghost;

  std::set<Element>::const_iterator end = relevant_ghost_elements.end();
  for(; it != last_type; ++it) {
    element.type = *it;
    UInt nb_ghost_elem = mesh.getNbElement(*it, _ghost);
    UInt nb_ghost_elem_protected = 0;
    try {
      nb_ghost_elem_protected = nb_ghost_protected(*it, _ghost);
    } catch (...) {}

    if(!remove_elem.getNewNumbering().exists(*it, _ghost))
      remove_elem.getNewNumbering().alloc(nb_ghost_elem, 1, *it, _ghost);
    else remove_elem.getNewNumbering(*it, _ghost).resize(nb_ghost_elem);
    Array<UInt> & new_numbering = remove_elem.getNewNumbering(*it, _ghost);
    for (UInt g = 0; g < nb_ghost_elem; ++g) {
      element.element = g;
      if (element.element >= nb_ghost_elem_protected && relevant_ghost_elements.find(element) == end) {
	ghosts_to_erase.push_back(element);
	new_numbering(element.element) = UInt(-1);
      }
    }
    /// renumber remaining ghosts
    UInt ng = 0;
    for (UInt g = 0; g < nb_ghost_elem; ++g) {
      if (new_numbering(g) != UInt(-1)) {
	new_numbering(g) = ng;
	++ng;
      }
    }
  }

  mesh.sendEvent(remove_elem);
  this->onElementsRemoved(ghosts_to_erase, remove_elem.getNewNumbering(), remove_elem);

}


__END_AKANTU__
