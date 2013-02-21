/**
 * @file   material_non_local_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Aug 31 11:09:48 2011
 *
 * @brief  Non-local inline implementation
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
__END_AKANTU__

/* -------------------------------------------------------------------------- */
#include "aka_types.hh"
#include "grid_synchronizer.hh"
#include "synchronizer_registry.hh"
#include "integrator.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt DIM, template <UInt> class WeightFunction>
MaterialNonLocal<DIM, WeightFunction>::MaterialNonLocal(SolidMechanicsModel & model,
							const ID & id)  :
  Material(model, id), radius(100.), weight_func(NULL), cell_list(NULL),
  update_weights(0), compute_stress_calls(0), is_creating_grid(false), grid_synchronizer(NULL) {
  AKANTU_DEBUG_IN();

  this->registerParam("radius"         , radius        , 100., _pat_readable  , "Non local radius");
  this->registerParam("UpdateWeights"  , update_weights, 0U  , _pat_modifiable, "Update weights frequency");
  this->registerParam("Weight function", weight_func   ,       _pat_internal);

  this->is_non_local = true;

  this->weight_func = new WeightFunction<DIM>(*this);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
MaterialNonLocal<spatial_dimension, WeightFunction>::~MaterialNonLocal() {
  AKANTU_DEBUG_IN();

  delete cell_list;
  delete weight_func;
  delete grid_synchronizer;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::initMaterial() {
  AKANTU_DEBUG_IN();
  //  Material::initMaterial();
  Mesh & mesh = this->model->getFEM().getMesh();

  ByElementTypeReal quadrature_points_coordinates("quadrature_points_coordinates_tmp_nl", id);
  this->initInternalVector(quadrature_points_coordinates, spatial_dimension, true);
  computeQuadraturePointsCoordinates(quadrature_points_coordinates, _not_ghost);

  ByElementType<UInt> nb_ghost_protected;
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, _ghost);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _ghost);
  for(; it != last_type; ++it)
    nb_ghost_protected(mesh.getNbElement(*it, _ghost), *it, _ghost);


  createCellList(quadrature_points_coordinates);
  updatePairList(quadrature_points_coordinates);

  cleanupExtraGhostElement(nb_ghost_protected);

  weight_func->setRadius(radius);
  weight_func->init();

  computeWeights(quadrature_points_coordinates);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::cleanupExtraGhostElement(const ByElementType<UInt> & nb_ghost_protected) {
  AKANTU_DEBUG_IN();

  // Create list of element to keep
  std::set<Element> relevant_ghost_element;

  pair_type::const_iterator first_pair_types = existing_pairs[1].begin();
  pair_type::const_iterator last_pair_types = existing_pairs[1].end();
  for (; first_pair_types != last_pair_types; ++first_pair_types) {
    ElementType type2 = first_pair_types->second;
    GhostType ghost_type2 = _ghost;
    UInt nb_quad2 = this->model->getFEM().getNbQuadraturePoints(type2);
    Vector<UInt> & elem_filter = element_filter(type2, ghost_type2);

    const Vector<UInt> & pairs =
      pair_list(first_pair_types->first, _not_ghost)(first_pair_types->second, ghost_type2);
    Vector<UInt>::const_iterator< types::Vector<UInt> > first_pair = pairs.begin(2);
    Vector<UInt>::const_iterator< types::Vector<UInt> > last_pair  = pairs.end(2);
    for(;first_pair != last_pair; ++first_pair) {
      UInt _q2 = (*first_pair)(1);
      QuadraturePoint q2(type2, elem_filter(_q2 / nb_quad2), _q2 % nb_quad2, ghost_type2);
      relevant_ghost_element.insert(q2);
    }
  }

  // Create list of element to remove and new numbering for element to keep
  Mesh & mesh = this->model->getFEM().getMesh();
  std::set<Element> ghost_to_erase;

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, _ghost);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _ghost);

  RemovedElementsEvent remove_elem(mesh);

  Element element;
  element.ghost_type = _ghost;
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

    Vector<UInt> & elem_filter = element_filter(*it, _ghost);
    Vector<UInt> & new_numbering = remove_elem.getNewNumbering(*it, _ghost);
    UInt ng = 0;
    for (UInt g = 0; g < nb_ghost_elem; ++g) {
      element.element = elem_filter(g);
      if(element.element >= nb_ghost_elem_protected &&
	 (std::find(relevant_ghost_element.begin(),
		    relevant_ghost_element.end(),
		    element) == relevant_ghost_element.end())) {
	ghost_to_erase.insert(element);
	remove_elem.getList().push_back(element);

	new_numbering(g) = UInt(-1);
      } else {
	new_numbering(g) = ng;
	++ng;
      }
    }
  }

  mesh.sendEvent(remove_elem);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::createCellList(ByElementTypeReal & quadrature_points_coordinates) {
  AKANTU_DEBUG_IN();

  const Real safety_factor = 1.2; // for the cell grid spacing
  Mesh & mesh = this->model->getFEM().getMesh();
  mesh.computeBoundingBox();

  Real lower_bounds[spatial_dimension];
  Real upper_bounds[spatial_dimension];
  mesh.getLocalLowerBounds(lower_bounds);
  mesh.getLocalUpperBounds(upper_bounds);

  Real spacing[spatial_dimension];
  for (UInt i = 0; i < spatial_dimension; ++i) {
    spacing[i] = radius * safety_factor;
  }

  cell_list = new RegularGrid<QuadraturePoint>(spatial_dimension,
					       lower_bounds,
					       upper_bounds,
					       spacing);


  fillCellList(quadrature_points_coordinates, _not_ghost);

  is_creating_grid = true;
  SynchronizerRegistry & synch_registry = this->model->getSynchronizerRegistry();
  std::stringstream sstr; sstr << id << ":grid_synchronizer";
  grid_synchronizer = GridSynchronizer::createGridSynchronizer(mesh,
							       *cell_list,
							       sstr.str());
  synch_registry.registerSynchronizer(*grid_synchronizer, _gst_mnl_for_average);
  synch_registry.registerSynchronizer(*grid_synchronizer, _gst_mnl_weight);
  is_creating_grid = false;

  this->computeQuadraturePointsCoordinates(quadrature_points_coordinates, _ghost);
  fillCellList(quadrature_points_coordinates, _ghost);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::fillCellList(const ByElementTypeReal & quadrature_points_coordinates,
								       const GhostType & ghost_type) {
  Mesh & mesh = this->model->getFEM().getMesh();

  QuadraturePoint q;
  q.ghost_type = ghost_type;

  Mesh::type_iterator it        = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType (spatial_dimension, ghost_type);
  for(; it != last_type; ++it) {
    Vector<UInt> & elem_filter = element_filter(*it, ghost_type);
    UInt nb_element = elem_filter.getSize();
    UInt nb_quad    = this->model->getFEM().getNbQuadraturePoints(*it, ghost_type);

    const Vector<Real> & quads = quadrature_points_coordinates(*it, ghost_type);
    q.type = *it;

    Vector<Real>::const_iterator<types::RVector> quad = quads.begin(spatial_dimension);
    UInt * elem = elem_filter.storage();

    for (UInt e = 0; e < nb_element; ++e) {
      q.element = *elem;
      for (UInt nq = 0; nq < nb_quad; ++nq) {
	q.num_point = nq;
	q.setPosition(*quad);
	cell_list->insert(q, *quad);
	++quad;
      }
      ++elem;
    }
  }
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::updatePairList(const ByElementTypeReal & quadrature_points_coordinates) {
  AKANTU_DEBUG_IN();

  Mesh & mesh = this->model->getFEM().getMesh();

  GhostType ghost_type = _not_ghost;

  // generate the pair of neighbor depending of the cell_list
  Mesh::type_iterator it        = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);

  for(; it != last_type; ++it) {
    // Preparing datas
    const Vector<Real> & quads = quadrature_points_coordinates(*it, ghost_type);
    Vector<Real>::const_iterator<types::RVector> first_quad = quads.begin(spatial_dimension);
    Vector<Real>::const_iterator<types::RVector> last_quad  = quads.end(spatial_dimension);

    ByElementTypeUInt & pairs = pair_list(ByElementTypeUInt("pairs", id, memory_id),
					  *it,
					  ghost_type);

    ElementType current_element_type = _not_defined;
    GhostType current_ghost_type = _casper;
    UInt existing_pairs_num = 0;

    Vector<UInt> * neighbors = NULL;
    Vector<UInt> * element_index_material2 = NULL;

    UInt my_num_quad = 0;

    // loop over quad points
    for(;first_quad != last_quad; ++first_quad, ++my_num_quad) {
      RegularGrid<QuadraturePoint>::Cell cell = cell_list->getCell(*first_quad);

      RegularGrid<QuadraturePoint>::neighbor_cells_iterator first_neigh_cell =
      cell_list->beginNeighborCells(cell);
      RegularGrid<QuadraturePoint>::neighbor_cells_iterator last_neigh_cell =
      cell_list->endNeighborCells(cell);

      // loop over neighbors cells of the one containing the current quadrature
      // point
      for (; first_neigh_cell != last_neigh_cell; ++first_neigh_cell) {
	RegularGrid<QuadraturePoint>::iterator first_neigh_quad =
	cell_list->beginCell(*first_neigh_cell);
	RegularGrid<QuadraturePoint>::iterator last_neigh_quad =
	cell_list->endCell(*first_neigh_cell);

	// loop over the quadrature point in the current cell of the cell list
	for (;first_neigh_quad != last_neigh_quad; ++first_neigh_quad){
	  QuadraturePoint & quad = *first_neigh_quad;
	  UInt nb_quad_per_elem =
	    this->model->getFEM().getNbQuadraturePoints(quad.type,
							quad.ghost_type);

	  // little optimization to not search in the map at each quad points
	  if(quad.type != current_element_type ||
	     quad.ghost_type != current_ghost_type) {

	    current_element_type = quad.type;
	    current_ghost_type   = quad.ghost_type;
	    existing_pairs_num = quad.ghost_type == _not_ghost ? 0 : 1;
	    if(!pairs.exists(current_element_type, current_ghost_type)) {
	      neighbors = &(pairs.alloc(0, 2,
					current_element_type,
					current_ghost_type));
	    } else {
	      neighbors = &(pairs(current_element_type,
				  current_ghost_type));
	    }
	    existing_pairs[existing_pairs_num].insert(std::pair<ElementType,
								ElementType>(*it,
											  current_element_type));
	    element_index_material2 =
	      &(this->model->getElementIndexByMaterial(current_element_type,
						       current_ghost_type));
	  }

	  UInt neigh_num_quad =
	    (*element_index_material2)(quad.element) * nb_quad_per_elem +
	    quad.num_point;

	  const types::RVector & neigh_quad = quad.getPosition();

	  Real distance = first_quad->distance(neigh_quad);
	  if(distance <= radius &&
	     (current_ghost_type == _ghost ||
	      (current_ghost_type == _not_ghost && my_num_quad <= neigh_num_quad))) { // sotring only half lists
	    UInt pair[2];
	    pair[0] = my_num_quad;
	    pair[1] = neigh_num_quad;

	    neighbors->push_back(pair);
	  }

	}
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::computeWeights(const ByElementTypeReal & quadrature_points_coordinates) {
  AKANTU_DEBUG_IN();

  GhostType ghost_type1;
  ghost_type1 = _not_ghost;

  ByElementTypeReal quadrature_points_volumes("quadrature_points_volumes", id, memory_id);
  this->initInternalVector(quadrature_points_volumes, 1, true);
  resizeInternalVector(quadrature_points_volumes);

  const FEM & fem = this->model->getFEM();

  weight_func->updateInternals(quadrature_points_volumes);

  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type2 = (GhostType) gt;
    UInt existing_pairs_num = gt - _not_ghost;
    pair_type::iterator first_pair_types = existing_pairs[existing_pairs_num].begin();
    pair_type::iterator last_pair_types = existing_pairs[existing_pairs_num].end();

    // Compute the weights
    for (; first_pair_types != last_pair_types; ++first_pair_types) {
      ElementType type1 = first_pair_types->first;
      ElementType type2 = first_pair_types->second;

      const Vector<UInt> & pairs = pair_list(type1, ghost_type1)(type2, ghost_type2);

      std::string ghost_id = "";
      if (ghost_type1 == _ghost) ghost_id = ":ghost";

      ByElementTypeReal & weights_type_1 = pair_weight(type1, ghost_type1);
      std::stringstream sstr; sstr << id << ":pair_weight:" << type1 << ghost_id;
      weights_type_1.setID(sstr.str());

      Vector<Real> * tmp_weight = NULL;
      if(!weights_type_1.exists(type2, ghost_type2)) {
	tmp_weight = &(weights_type_1.alloc(0, 2, type2, ghost_type2));
      } else {
	tmp_weight = &(weights_type_1(type2, ghost_type2));
      }
      Vector<Real> & weights = *tmp_weight;
      weights.resize(pairs.getSize());
      weights.clear();

      const Vector<Real> & jacobians_1 = fem.getIntegratorInterface().getJacobians(type1, ghost_type1);
      const Vector<Real> & jacobians_2 = fem.getIntegratorInterface().getJacobians(type2, ghost_type2);

      UInt nb_quad1 = fem.getNbQuadraturePoints(type1);
      UInt nb_quad2 = fem.getNbQuadraturePoints(type2);

      Vector<Real> & quads_volumes1 = quadrature_points_volumes(type1, ghost_type1);
      Vector<Real> & quads_volumes2 = quadrature_points_volumes(type2, ghost_type2);

      Vector<Real>::const_iterator<types::RVector> iquads1;
      Vector<Real>::const_iterator<types::RVector> iquads2;
      iquads1 = quadrature_points_coordinates(type1, ghost_type1).begin(spatial_dimension);
      iquads2 = quadrature_points_coordinates(type2, ghost_type2).begin(spatial_dimension);

      Vector<UInt>::const_iterator< types::Vector<UInt> > first_pair = pairs.begin(2);
      Vector<UInt>::const_iterator< types::Vector<UInt> > last_pair  = pairs.end(2);
      Vector<Real>::iterator<types::RVector> weight  = weights.begin(2);

      this->weight_func->selectType(type1, ghost_type1, type2, ghost_type2);

      // Weight function
      for(;first_pair != last_pair; ++first_pair, ++weight) {
	UInt _q1 = (*first_pair)(0);
	UInt _q2 = (*first_pair)(1);
	const types::RVector & pos1 = iquads1[_q1];
	const types::RVector & pos2 = iquads2[_q2];
	QuadraturePoint q1(_q1 / nb_quad1, _q1 % nb_quad1, _q1, pos1, type1, ghost_type1);
	QuadraturePoint q2(_q2 / nb_quad2, _q2 % nb_quad2, _q2, pos2, type2, ghost_type2);

	Real r = pos1.distance(pos2);

	Real w2J2 = jacobians_2(_q2);
	(*weight)(0) = w2J2 * this->weight_func->operator()(r, q1, q2);
	if(ghost_type2 != _ghost && _q1 != _q2) {
	  Real w1J1 = jacobians_1(_q1);
	  (*weight)(1) = w1J1 * this->weight_func->operator()(r, q2, q1);
	} else
	  (*weight)(1) = 0;

	quads_volumes1(_q1) += (*weight)(0);
	if(ghost_type2 != _ghost) quads_volumes2(_q2) += (*weight)(1);
      }
    }
  }

  //normalize the weights
  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type2 = (GhostType) gt;
    UInt existing_pairs_num = gt - _not_ghost;
    pair_type::iterator first_pair_types = existing_pairs[existing_pairs_num].begin();
    pair_type::iterator last_pair_types = existing_pairs[existing_pairs_num].end();

    first_pair_types = existing_pairs[existing_pairs_num].begin();
    for (; first_pair_types != last_pair_types; ++first_pair_types) {
      ElementType type1 = first_pair_types->first;
      ElementType type2 = first_pair_types->second;

      const Vector<UInt> & pairs = pair_list(type1, ghost_type1)(type2, ghost_type2);
      Vector<Real> & weights = pair_weight(type1, ghost_type1)(type2, ghost_type2);

      Vector<Real> & quads_volumes1 = quadrature_points_volumes(type1, ghost_type1);
      Vector<Real> & quads_volumes2 = quadrature_points_volumes(type2, ghost_type2);

      Vector<UInt>::const_iterator< types::Vector<UInt> > first_pair = pairs.begin(2);
      Vector<UInt>::const_iterator< types::Vector<UInt> > last_pair  = pairs.end(2);
      Vector<Real>::iterator<types::RVector> weight  = weights.begin(2);

      for(;first_pair != last_pair; ++first_pair, ++weight) {
	UInt q1 = (*first_pair)(0);
	UInt q2 = (*first_pair)(1);

	(*weight)(0) *= 1. / quads_volumes1(q1);
	if(ghost_type2 != _ghost) (*weight)(1) *= 1. / quads_volumes2(q2);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
template<typename T>
void MaterialNonLocal<spatial_dimension, WeightFunction>::weightedAvergageOnNeighbours(const ByElementTypeVector<T> & to_accumulate,
										       ByElementTypeVector<T> & accumulated,
										       UInt nb_degree_of_freedom,
										       GhostType ghost_type2) const {
  AKANTU_DEBUG_IN();

  UInt existing_pairs_num = 0;
  if (ghost_type2 == _ghost) existing_pairs_num = 1;

  pair_type::const_iterator first_pair_types = existing_pairs[existing_pairs_num].begin();
  pair_type::const_iterator last_pair_types = existing_pairs[existing_pairs_num].end();

  GhostType ghost_type1 = _not_ghost; // does not make sens the ghost vs ghost so this should always by _not_ghost

  for (; first_pair_types != last_pair_types; ++first_pair_types) {
    const Vector<UInt> & pairs =
      pair_list(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);
    const Vector<Real> & weights =
      pair_weight(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);


    const Vector<T> & to_acc = to_accumulate(first_pair_types->second, ghost_type2);
    Vector<T> & acc = accumulated(first_pair_types->first, ghost_type1);

    if(ghost_type2 == _not_ghost) acc.clear();

    Vector<UInt>::const_iterator< types::Vector<UInt> > first_pair = pairs.begin(2);
    Vector<UInt>::const_iterator< types::Vector<UInt> > last_pair  = pairs.end(2);
    Vector<Real>::const_iterator< types::Vector<Real> > pair_w = weights.begin(2);

    for(;first_pair != last_pair; ++first_pair, ++pair_w) {
      UInt q1 = (*first_pair)(0);
      UInt q2 = (*first_pair)(1);
      for(UInt d = 0; d < nb_degree_of_freedom; ++d){
	acc(q1, d) += (*pair_w)(0) * to_acc(q2, d);
	if(ghost_type2 != _ghost) acc(q2, d) += (*pair_w)(1) * to_acc(q1, d);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::updateResidual(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  // Update the weights for the non local variable averaging
  if(ghost_type == _not_ghost &&
     this->update_weights &&
     (this->compute_stress_calls % this->update_weights == 0)) {
    ByElementTypeReal quadrature_points_coordinates("quadrature_points_coordinates", id);
    Mesh & mesh = this->model->getFEM().getMesh();
    mesh.initByElementTypeVector(quadrature_points_coordinates, spatial_dimension, 0);
    computeQuadraturePointsCoordinates(quadrature_points_coordinates, _not_ghost);
    computeQuadraturePointsCoordinates(quadrature_points_coordinates, _ghost);
    computeWeights(quadrature_points_coordinates);
  }
  if(ghost_type == _not_ghost) ++this->compute_stress_calls;

  computeAllStresses(ghost_type);

  computeNonLocalStresses(ghost_type);
  assembleResidual(ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::computeAllNonLocalStresses(GhostType ghost_type) {
  // Update the weights for the non local variable averaging
  if(ghost_type == _not_ghost) {
    if(this->update_weights && (this->compute_stress_calls % this->update_weights == 0)) {
      this->model->getSynchronizerRegistry().asynchronousSynchronize(_gst_mnl_weight);

      ByElementTypeReal quadrature_points_coordinates("quadrature_points_coordinates", id);
      Mesh & mesh = this->model->getFEM().getMesh();
      mesh.initByElementTypeVector(quadrature_points_coordinates, spatial_dimension, 0);
      computeQuadraturePointsCoordinates(quadrature_points_coordinates, _not_ghost);
      computeQuadraturePointsCoordinates(quadrature_points_coordinates, _ghost);

      this->model->getSynchronizerRegistry().waitEndSynchronize(_gst_mnl_weight);

      computeWeights(quadrature_points_coordinates);
    }

    typename std::map<ID, NonLocalVariable>::iterator it = non_local_variables.begin();
    typename std::map<ID, NonLocalVariable>::iterator end = non_local_variables.end();
    for(;it != end; ++it) {
      NonLocalVariable & non_local_variable = it->second;

      resizeInternalVector(*non_local_variable.non_local_variable);
      this->weightedAvergageOnNeighbours(*non_local_variable.local_variable, *non_local_variable.non_local_variable,
					 non_local_variable.non_local_variable_nb_component, _not_ghost);
    }

    ++this->compute_stress_calls;
  } else {

    typename std::map<ID, NonLocalVariable>::iterator it = non_local_variables.begin();
    typename std::map<ID, NonLocalVariable>::iterator end = non_local_variables.end();
    for(;it != end; ++it) {
      NonLocalVariable & non_local_variable = it->second;
      this->weightedAvergageOnNeighbours(*non_local_variable.local_variable, *non_local_variable.non_local_variable,
					 non_local_variable.non_local_variable_nb_component, _ghost);
    }

    computeNonLocalStresses(_not_ghost);
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
bool MaterialNonLocal<spatial_dimension, WeightFunction>::parseParam(const std::string & key,
								   const std::string & value,
								   __attribute__((unused)) const ID & id) {
  std::stringstream sstr(value);
  if(key == "radius") { sstr >> radius; }
  else if(key == "UpdateWeights") { sstr >> update_weights; }
  else if(!weight_func->parseParam(key, value)) return Material::parseParam(key, value, id);
  return true;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::savePairs(const std::string & filename) const {
  std::ofstream pout;

  std::stringstream sstr;

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int prank = comm.whoAmI();
  sstr << filename << "." << prank;

  pout.open(sstr.str().c_str());

  GhostType ghost_type1;
  ghost_type1 = _not_ghost;

  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type2 = (GhostType) gt;
    UInt existing_pairs_num = gt - _not_ghost;

    pair_type::const_iterator first_pair_types = existing_pairs[existing_pairs_num].begin();
    pair_type::const_iterator last_pair_types = existing_pairs[existing_pairs_num].end();

    for (; first_pair_types != last_pair_types; ++first_pair_types) {
      const Vector<UInt> & pairs =
	pair_list(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);
      const Vector<Real> & weights =
	pair_weight(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);

      pout << "Types : " << first_pair_types->first << " (" << ghost_type1 << ") - " << first_pair_types->second << " (" << ghost_type2 << ")" << std::endl;

      Vector<UInt>::const_iterator< types::Vector<UInt> > first_pair = pairs.begin(2);
      Vector<UInt>::const_iterator< types::Vector<UInt> > last_pair  = pairs.end(2);
      Vector<Real>::const_iterator<types::RVector> pair_w = weights.begin(2);

      for(;first_pair != last_pair; ++first_pair, ++pair_w) {
	UInt q1 = (*first_pair)(0);
	UInt q2 = (*first_pair)(1);
	pout << q1 << " " << q2 << " "<< (*pair_w)(0) << " " << (*pair_w)(1) << std::endl;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::neighbourhoodStatistics(const std::string & filename) const {
  std::ofstream pout;
  pout.open(filename.c_str());

  const Mesh & mesh = this->model->getFEM().getMesh();

  GhostType ghost_type1;
  ghost_type1 = _not_ghost;

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int prank = comm.whoAmI();

  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type2 = (GhostType) gt;
    UInt existing_pairs_num = gt - _not_ghost;

    ByElementTypeUInt nb_neighbors("nb_neighbours", id, memory_id);
    mesh.initByElementTypeVector(nb_neighbors, 1, spatial_dimension);
    resizeInternalVector(nb_neighbors);

    pair_type::const_iterator first_pair_types = existing_pairs[existing_pairs_num].begin();
    pair_type::const_iterator last_pair_types = existing_pairs[existing_pairs_num].end();

    for (; first_pair_types != last_pair_types; ++first_pair_types) {
      const Vector<UInt> & pairs =
	pair_list(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);
      if(prank == 0) {
	pout << ghost_type2 << " ";
	pout << "Types : " << first_pair_types->first << " " << first_pair_types->second << std::endl;
      }
      Vector<UInt>::const_iterator< types::Vector<UInt> > first_pair = pairs.begin(2);
      Vector<UInt>::const_iterator< types::Vector<UInt> > last_pair  = pairs.end(2);
      Vector<UInt> & nb_neigh_1 = nb_neighbors(first_pair_types->first, ghost_type1);
      Vector<UInt> & nb_neigh_2 = nb_neighbors(first_pair_types->second, ghost_type2);
      for(;first_pair != last_pair; ++first_pair) {
	UInt q1 = (*first_pair)(0);
	UInt q2 = (*first_pair)(1);
	++(nb_neigh_1(q1));
	if(q1 != q2) ++(nb_neigh_2(q2));
      }
    }

    Mesh::type_iterator it        = mesh.firstType(spatial_dimension, ghost_type1);
    Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type1);
    UInt nb_quads = 0;
    Real sum_nb_neig = 0;
    UInt max_nb_neig = 0;
    UInt min_nb_neig = std::numeric_limits<UInt>::max();
    for(; it != last_type; ++it) {
      Vector<UInt> & nb_neighor = nb_neighbors(*it, ghost_type1);
      Vector<UInt>::iterator<UInt> nb_neigh = nb_neighor.begin();
      Vector<UInt>::iterator<UInt> end_neigh  = nb_neighor.end();

      for (; nb_neigh != end_neigh; ++nb_neigh, ++nb_quads) {
	UInt nb = *nb_neigh;
	sum_nb_neig += nb;
	max_nb_neig = std::max(max_nb_neig, nb);
	min_nb_neig = std::min(min_nb_neig, nb);
      }
    }


    comm.allReduce(&nb_quads,    1, _so_sum);
    comm.allReduce(&sum_nb_neig, 1, _so_sum);
    comm.allReduce(&max_nb_neig, 1, _so_max);
    comm.allReduce(&min_nb_neig, 1, _so_min);

    if(prank == 0) {
      pout << ghost_type2 << " ";
      pout << "Nb quadrature points: " << nb_quads << std::endl;

      Real mean_nb_neig = sum_nb_neig / Real(nb_quads);
      pout << ghost_type2 << " ";
      pout << "Average nb neighbors: " << mean_nb_neig << "(" << sum_nb_neig << ")" << std::endl;

      pout << ghost_type2 << " ";
      pout << "Max nb neighbors:     " << max_nb_neig << std::endl;

      pout << ghost_type2 << " ";
      pout << "Min nb neighbors:     " << min_nb_neig << std::endl;
    }
  }
  pout.close();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
inline UInt MaterialNonLocal<spatial_dimension, WeightFunction>::getNbDataForElements(const Vector<Element> & elements,
										      SynchronizationTag tag) const {
  UInt nb_quadrature_points = this->getModel().getNbQuadraturePoints(elements);
  UInt size = 0;

  if(tag == _gst_mnl_for_average) {
    typename std::map<ID, NonLocalVariable>::const_iterator it = non_local_variables.begin();
    typename std::map<ID, NonLocalVariable>::const_iterator end = non_local_variables.end();

    for(;it != end; ++it) {
      const NonLocalVariable & non_local_variable = it->second;
      size += non_local_variable.non_local_variable_nb_component * sizeof(Real) * nb_quadrature_points;
    }
  }

  size += weight_func->getNbDataForElements(elements, tag);

  return size;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
inline void MaterialNonLocal<spatial_dimension, WeightFunction>::packElementData(CommunicationBuffer & buffer,
										 const Vector<Element> & elements,
										 SynchronizationTag tag) const {
  if(tag == _gst_mnl_for_average) {
    typename std::map<ID, NonLocalVariable>::const_iterator it = non_local_variables.begin();
    typename std::map<ID, NonLocalVariable>::const_iterator end = non_local_variables.end();

    for(;it != end; ++it) {
      const NonLocalVariable & non_local_variable = it->second;
      this->packElementDataHelper(*non_local_variable.local_variable,
				  buffer, elements);
    }
  }

  weight_func->packElementData(buffer, elements, tag);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
inline void MaterialNonLocal<spatial_dimension, WeightFunction>::unpackElementData(CommunicationBuffer & buffer,
										   const Vector<Element> & elements,
										   SynchronizationTag tag) {
  if(tag == _gst_mnl_for_average) {
    typename std::map<ID, NonLocalVariable>::iterator it = non_local_variables.begin();
    typename std::map<ID, NonLocalVariable>::iterator end = non_local_variables.end();

    for(;it != end; ++it) {
      NonLocalVariable & non_local_variable = it->second;
      this->unpackElementDataHelper(*non_local_variable.local_variable,
				    buffer, elements);
    }
  }

  weight_func->unpackElementData(buffer, elements, tag);
}

/* -------------------------------------------------------------------------- */
// template<UInt spatial_dimension, template <UInt> class WeightFunction>
// inline void MaterialNonLocal<spatial_dimension, WeightFunction>::onElementsAdded(const Vector<Element> & element_list) {
//   AKANTU_DEBUG_IN();

//   Material::onElementsAdded(element_list, event);

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
inline void MaterialNonLocal<spatial_dimension, WeightFunction>::onElementsRemoved(const Vector<Element> & element_list,
										   const ByElementTypeUInt & new_numbering,
										   __attribute__((unused)) const RemovedElementsEvent & event) {
  AKANTU_DEBUG_IN();

  Material::onElementsRemoved(element_list, new_numbering, event);

  pair_type::const_iterator first_pair_types = existing_pairs[1].begin();
  pair_type::const_iterator last_pair_types = existing_pairs[1].end();

  // Renumber element to keep
  for (; first_pair_types != last_pair_types; ++first_pair_types) {
    ElementType type2 = first_pair_types->second;
    GhostType ghost_type2 = _ghost;
    UInt nb_quad2 = this->model->getFEM().getNbQuadraturePoints(type2);

    Vector<UInt> & pairs =
      pair_list(first_pair_types->first, _not_ghost)(first_pair_types->second, ghost_type2);
    Vector<UInt>::iterator< types::Vector<UInt> > first_pair = pairs.begin(2);
    Vector<UInt>::iterator< types::Vector<UInt> > last_pair  = pairs.end(2);
    for(;first_pair != last_pair; ++first_pair) {
      UInt _q2 = (*first_pair)(1);
      const Vector<UInt> & renumbering = new_numbering(type2, ghost_type2);
      UInt el = _q2 / nb_quad2;
      UInt new_el = renumbering(el);
      AKANTU_DEBUG_ASSERT(new_el != UInt(-1), "A local quad as been removed instead f just renumbered");
      (*first_pair)(1) = new_el * nb_quad2 + _q2 % nb_quad2;
    }
  }

  AKANTU_DEBUG_OUT();
}
