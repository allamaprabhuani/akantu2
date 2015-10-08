/**
 * @file   non_local_neighborhood_base.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Mon Sep 21 18:10:49 2015
 *
 * @brief  Implementation of non-local neighborhood base
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
#include "non_local_neighborhood_base.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
NonLocalNeighborhoodBase::NonLocalNeighborhoodBase(const SolidMechanicsModel & model,
						   const ElementTypeMapReal & quad_coordinates,
						   const ID & id,
						   const MemoryID & memory_id)  :
  Memory(id, memory_id),
  Parsable(_st_non_local, id),
  model(model),
  non_local_radius(0.),
  spatial_grid(NULL), 
  is_creating_grid(false), 
  grid_synchronizer(NULL),
  quad_coordinates(quad_coordinates),
  spatial_dimension(this->model.getSpatialDimension()),
  synch_registry(NULL) {

  AKANTU_DEBUG_IN();

  this->registerParam("radius"       , non_local_radius             , 100.,
		      _pat_parsable | _pat_readable  , "Non local radius");


  this->createSynchronizerRegistry(this);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
NonLocalNeighborhoodBase::~NonLocalNeighborhoodBase() {
  AKANTU_DEBUG_IN();

  delete spatial_grid;
  delete grid_synchronizer;
  delete synch_registry;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NonLocalNeighborhoodBase::initNeighborhood() {
  AKANTU_DEBUG_IN();
  //  Material::initMaterial();
  Mesh & mesh = this->model.getMesh();

  ElementTypeMap<UInt> nb_ghost_protected;
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, _ghost);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _ghost);
  for(; it != last_type; ++it)
    nb_ghost_protected(mesh.getNbElement(*it, _ghost), *it, _ghost);

  AKANTU_DEBUG_INFO("Creating the grid");
  this->createGrid();

  /// todo: create pairs, set radius of weight function, intialize weight function, compute weights

  AKANTU_DEBUG_OUT();
}

/* ------------------------------------------------------------------------- */
void NonLocalNeighborhoodBase::createGrid() {
  AKANTU_DEBUG_IN();

  const Real safety_factor = 1.2; // for the cell grid spacing
  Mesh & mesh = this->model.getMesh();
  mesh.computeBoundingBox();

  const Vector<Real> & lower_bounds = mesh.getLocalLowerBounds();
  const Vector<Real> & upper_bounds = mesh.getLocalUpperBounds();
  Vector<Real> center = 0.5 * (upper_bounds + lower_bounds);
  Vector<Real> spacing(spatial_dimension, this->non_local_radius * safety_factor);

  spatial_grid = new SpatialGrid<QuadraturePoint>(spatial_dimension, spacing, center);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NonLocalNeighborhoodBase::updatePairList() {
  AKANTU_DEBUG_IN();

  //// loop over all quads -> all cells
  SpatialGrid<QuadraturePoint>::cells_iterator cell_it = spatial_grid->beginCells();
  SpatialGrid<QuadraturePoint>::cells_iterator cell_end = spatial_grid->endCells();

  Vector<Real> q1_coords(spatial_dimension);
  Vector<Real> q2_coords(spatial_dimension);
  QuadraturePoint q1;
  QuadraturePoint q2;

  UInt counter = 0;
  for (; cell_it != cell_end; ++cell_it) {
    AKANTU_DEBUG_INFO("Looping on next cell");
    SpatialGrid<QuadraturePoint>::Cell::iterator first_quad =
      spatial_grid->beginCell(*cell_it);
    SpatialGrid<QuadraturePoint>::Cell::iterator last_quad =
      spatial_grid->endCell(*cell_it);
  
    for (;first_quad != last_quad; ++first_quad, ++counter){
      q1 = *first_quad;
      Array<Real>::const_vector_iterator coords_type_1_it = this->quad_coordinates(q1.type, q1.ghost_type).begin(spatial_dimension);
      q1_coords = coords_type_1_it[q1.global_num];
      AKANTU_DEBUG_INFO("Current quadrature point in this cell: " << q1);
      SpatialGrid<QuadraturePoint>::CellID cell_id = spatial_grid->getCellID(q1_coords);
      /// loop over all the neighbouring cells of the current quad
      SpatialGrid<QuadraturePoint>::neighbor_cells_iterator first_neigh_cell =
  	spatial_grid->beginNeighborCells(cell_id);
      SpatialGrid<QuadraturePoint>::neighbor_cells_iterator last_neigh_cell =
  	spatial_grid->endNeighborCells(cell_id);

      for (; first_neigh_cell != last_neigh_cell; ++first_neigh_cell) {
      	SpatialGrid<QuadraturePoint>::Cell::iterator first_neigh_quad =
      	  spatial_grid->beginCell(*first_neigh_cell);
      	SpatialGrid<QuadraturePoint>::Cell::iterator last_neigh_quad =
      	  spatial_grid->endCell(*first_neigh_cell);

      	// loop over the quadrature point in the current neighboring cell
      	for (;first_neigh_quad != last_neigh_quad; ++first_neigh_quad){
      	  q2 = *first_neigh_quad;
      Array<Real>::const_vector_iterator coords_type_2_it = this->quad_coordinates(q2.type, q2.ghost_type).begin(spatial_dimension);
	  q2_coords = coords_type_2_it[q2.global_num];

      	  Real distance = q1_coords.distance(q2_coords);

      	  if(distance <= this->non_local_radius + Math::getTolerance()  &&
      	     (q2.ghost_type == _ghost ||
      	      (q2.ghost_type == _not_ghost && q1.global_num <= q2.global_num))) { // storing only half lists
      	    pair_list[q2.ghost_type].push_back(std::make_pair(q1, q2));
      	  }
      	}
      }
    }

  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NonLocalNeighborhoodBase::createSynchronizerRegistry(DataAccessor * data_accessor){
  synch_registry = new SynchronizerRegistry(*data_accessor);
}

/* -------------------------------------------------------------------------- */
void NonLocalNeighborhoodBase::createGridSynchronizer() {
  this->is_creating_grid = true;
  std::set<SynchronizationTag> tags;
  tags.insert(_gst_mnl_for_average);
  tags.insert(_gst_mnl_weight);

  std::stringstream sstr; sstr << getID() << ":grid_synchronizer";
  this->grid_synchronizer = GridSynchronizer::createGridSynchronizer(this->model.getMesh(),
								     *spatial_grid,
								     sstr.str(),
								     synch_registry,
								     tags);
  this->is_creating_grid = false;

}

/* -------------------------------------------------------------------------- */
void NonLocalNeighborhoodBase::savePairs(const std::string & filename) const {
  std::ofstream pout;

  std::stringstream sstr;

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int prank = comm.whoAmI();
  sstr << filename << "." << prank;

  pout.open(sstr.str().c_str());

  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type2 = (GhostType) gt;

    PairList::const_iterator first_pair = pair_list[ghost_type2].begin();
    PairList::const_iterator last_pair  = pair_list[ghost_type2].end();

    for(;first_pair != last_pair; ++first_pair) {

      const QuadraturePoint & q1 = first_pair->first;
      const QuadraturePoint & q2 = first_pair->second;
      pout << q1 << " " << q2 << " " << std::endl;
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalNeighborhoodBase::saveNeighborCoords(const std::string & filename) const {

  /// this function is not optimazed and only used for tests on small meshes
  /// @todo maybe optimize this function for better performance?

  Vector<Real> q1_coords(spatial_dimension);
  Vector<Real> q2_coords(spatial_dimension);
  QuadraturePoint q1;
  QuadraturePoint q2;

  std::ofstream pout;

  std::stringstream sstr;

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int prank = comm.whoAmI();
  sstr << filename << "." << prank;

  pout.open(sstr.str().c_str());

  /// loop over all the quads and write the position of their neighbors
  SpatialGrid<QuadraturePoint>::cells_iterator cell_it = spatial_grid->beginCells();
  SpatialGrid<QuadraturePoint>::cells_iterator cell_end = spatial_grid->endCells();

  for (; cell_it != cell_end; ++cell_it) {
    SpatialGrid<QuadraturePoint>::Cell::iterator first_quad =
      spatial_grid->beginCell(*cell_it);
    SpatialGrid<QuadraturePoint>::Cell::iterator last_quad =
      spatial_grid->endCell(*cell_it);
  
    for (;first_quad != last_quad; ++first_quad){
      q1 = *first_quad;
      Array<Real>::const_vector_iterator coords_type_1_it = this->quad_coordinates(q1.type, q1.ghost_type).begin(spatial_dimension);
      q1_coords = coords_type_1_it[q1.global_num];
      pout << "#neighbors for quad " << q1.global_num << std::endl;
      pout << q1_coords << std::endl;
      for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
	GhostType ghost_type2 = (GhostType) gt;
	PairList::const_iterator first_pair = pair_list[ghost_type2].begin();
	PairList::const_iterator last_pair  = pair_list[ghost_type2].end();
	for(;first_pair != last_pair; ++first_pair) {
	  if (q1 == first_pair->first && first_pair->second != q1) {
	    q2 = first_pair->second;	 
	    Array<Real>::const_vector_iterator coords_type_2_it = 
	    this->quad_coordinates(q2.type, q2.ghost_type).begin(spatial_dimension);
	    q2_coords = coords_type_2_it[q2.global_num]; 
	    pout << q2_coords << std::endl;  	  
	  }
	  if (q1 == first_pair->second && first_pair->first != q1) {
	    q2 = first_pair->first;
	    Array<Real>::const_vector_iterator coords_type_2_it = 
	    this->quad_coordinates(q2.type, q2.ghost_type).begin(spatial_dimension);
	    q2_coords = coords_type_2_it[q2.global_num]; 
	    pout << q2_coords << std::endl;
	  }
	}
      }
    }
  }
}


__END_AKANTU__
