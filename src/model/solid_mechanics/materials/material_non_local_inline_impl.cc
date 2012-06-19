/**
 * @file   material_non_local_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 25 11:59:39 2011
 *
 * @brief  
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
/* -------------------------------------------------------------------------- */
#include <iostream>
#include <fstream>
/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt DIM, template <UInt> class WeightFunction>
MaterialNonLocal<DIM, WeightFunction>::MaterialNonLocal(SolidMechanicsModel & model,
							const ID & id)  :
  Material(model, id), radius(100.), weight_func(NULL), cell_list(NULL),
  update_weigths(0), compute_stress_calls(0) {
  AKANTU_DEBUG_IN();

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

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::initMaterial() {
  AKANTU_DEBUG_IN();
  //  Material::initMaterial();
  ByElementTypeReal quadrature_points_coordinates("quadrature_points_coordinates_tmp_nl", id);
  computeQuadraturePointsCoordinates(quadrature_points_coordinates);

  weight_func->setRadius(radius);
  weight_func->init();

  createCellList(quadrature_points_coordinates);
  updatePairList(quadrature_points_coordinates);

  computeWeights(quadrature_points_coordinates);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::updateResidual(Vector<Real> & displacement, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  // Update the weights for the non local variable averaging
  if(ghost_type == _not_ghost && this->update_weigths && (this->compute_stress_calls % this->update_weigths == 0)) {
    ByElementTypeReal quadrature_points_coordinates("quadrature_points_coordinates", id);
    computeQuadraturePointsCoordinates(quadrature_points_coordinates);
    computeWeights(quadrature_points_coordinates);
  }
  if(ghost_type == _not_ghost) ++this->compute_stress_calls;


  computeAllStresses(displacement, ghost_type);
  computeNonLocalStress(ghost_type);
  assembleResidual(ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
template<typename T>
void MaterialNonLocal<spatial_dimension, WeightFunction>::weightedAvergageOnNeighbours(const ByElementTypeVector<T> & to_accumulate,
								    ByElementTypeVector<T> & accumulated,
								    UInt nb_degree_of_freedom) const {
  AKANTU_DEBUG_IN();

  std::set< std::pair<ElementType, ElementType> >::const_iterator first_pair_types = existing_pairs.begin();
  std::set< std::pair<ElementType, ElementType> >::const_iterator last_pair_types = existing_pairs.end();

  GhostType ghost_type1, ghost_type2;
  ghost_type1 = ghost_type2 = _not_ghost;

  for (; first_pair_types != last_pair_types; ++first_pair_types) {
    const Vector<UInt> & pairs =
      pair_list(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);
    const Vector<Real> & weights =
      pair_weight(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);


    const Vector<T> & to_acc = to_accumulate(first_pair_types->second, ghost_type2);
    Vector<T> & acc = accumulated(first_pair_types->first, ghost_type1);

    acc.resize(to_acc.getSize());
    acc.clear();

    Vector<UInt>::const_iterator< types::Vector<UInt> > first_pair = pairs.begin(2);
    Vector<UInt>::const_iterator< types::Vector<UInt> > last_pair  = pairs.end(2);
    Vector<Real>::const_iterator< types::Vector<Real> > pair_w = weights.begin(2);

    for(;first_pair != last_pair; ++first_pair, ++pair_w) {
      UInt q1 = (*first_pair)(0);
      UInt q2 = (*first_pair)(1);
      for(UInt d = 0; d < nb_degree_of_freedom; ++d){
	acc(q1, d) += (*pair_w)(0) * to_acc(q2, d);
	acc(q2, d) += (*pair_w)(1) * to_acc(q1, d);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::createCellList(const ByElementTypeReal & quadrature_points_coordinates) {
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

  GhostType ghost_type = _not_ghost;

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);

  // first generate the quad points coordinate and count the number of points per cell
  for(; it != last_type; ++it) {
    const Vector<Real> & quads = quadrature_points_coordinates(*it, ghost_type);
    Vector<Real>::const_iterator<types::RVector> first_quad, last_quad;
    first_quad = quads.begin(spatial_dimension);
    last_quad = quads.end(spatial_dimension);

    for(;first_quad != last_quad; ++first_quad) {
      cell_list->count(*first_quad);
    }
  }

  QuadraturePoint q;
  q.ghost_type = ghost_type;

  // second insert the point in the cells
  cell_list->beginInsertions();
  it = mesh.firstType(spatial_dimension, ghost_type);
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

  cell_list->endInsertions();

  SynchronizerRegistry & synch_registry = this->model->getSynchronizerRegistry();
  std::stringstream sstr; sstr << id << ":grid_synchronizer";
  GridSynchronizer * synch = GridSynchronizer::createGridSynchronizer(mesh,
                                                                      *cell_list,
                                                                      sstr.str());
  synch_registry.registerSynchronizer(*synch, _gst_mnl_damage);

  AKANTU_DEBUG_OUT();
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
    //Real * neigh_quad_positions = NULL;
    Vector<UInt> * neighbors = NULL;

    UInt my_num_quad = 0;

    // loop over quad points

    for(;first_quad != last_quad; ++first_quad, ++my_num_quad) {
      RegularGrid<QuadraturePoint>::Cell cell = cell_list->getCell(*first_quad);

      RegularGrid<QuadraturePoint>::neighbor_cells_iterator first_neigh_cell =
      cell_list->beginNeighborCells(cell);
      RegularGrid<QuadraturePoint>::neighbor_cells_iterator last_neigh_cell =
      cell_list->endNeighborCells(cell);

      // loop over neighbor  cells of the one containing  the current quadrature
      // point
      for (; first_neigh_cell != last_neigh_cell; ++first_neigh_cell) {
        RegularGrid<QuadraturePoint>::iterator first_neigh_quad =
        cell_list->beginCell(*first_neigh_cell);
        RegularGrid<QuadraturePoint>::iterator last_neigh_quad =
        cell_list->endCell(*first_neigh_cell);

        // loop over the quadrature point in the current cell of the cell list
        for (;first_neigh_quad != last_neigh_quad; ++first_neigh_quad){
          QuadraturePoint & quad = *first_neigh_quad;
          UInt nb_quad_per_elem =  this->model->getFEM().getNbQuadraturePoints(quad.type, quad.ghost_type);
          UInt neigh_num_quad = quad.element * nb_quad_per_elem + quad.num_point;

          // little optimization to not search in the map at each quad points
          if(quad.type != current_element_type || quad.ghost_type != current_ghost_type) {
            //            neigh_quad_positions = quadrature_points_coordinates(quad.type,
            //                                                                 quad.ghost_type).storage();
            current_element_type = quad.type;
            current_ghost_type   = quad.ghost_type;
            if(!pairs.exists(current_element_type, current_ghost_type)) {
              neighbors = &(pairs.alloc(0, 2,
                                        current_element_type,
                                        current_ghost_type));
            } else {
              neighbors = &(pairs(current_element_type,
                                  current_ghost_type));
            }
            existing_pairs.insert(std::pair<ElementType, ElementType>(*it, current_element_type));
          }

          // types::RVector neigh_quad(neigh_quad_positions + neigh_num_quad * spatial_dimension,
          //                          spatial_dimension);
          const types::RVector & neigh_quad = quad.getPosition();

          Real distance = first_quad->distance(neigh_quad);
          if(distance <= radius && my_num_quad <= neigh_num_quad) { // sotring only half lists
            UInt pair[2];
            pair[0] = my_num_quad;
	    pair[1] = neigh_num_quad;

            neighbors->push_back(pair);
          }
          //      }
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

  std::set< std::pair<ElementType, ElementType> >::iterator first_pair_types;
  std::set< std::pair<ElementType, ElementType> >::iterator last_pair_types = existing_pairs.end();

  GhostType ghost_type1, ghost_type2;
  ghost_type1 = ghost_type2 = _not_ghost;

  ByElementTypeReal quadrature_points_volumes("quadrature_points_volumes", id, memory_id);
  this->model->getFEM().getMesh().initByElementTypeVector(quadrature_points_volumes, 1, 0);

  weight_func->updateInternals(quadrature_points_volumes);

  // Compute the weights
  first_pair_types = existing_pairs.begin();
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

    const Vector<UInt> & elem_filter = element_filter(type1, ghost_type1);
    UInt nb_quad1 = this->model->getFEM().getNbQuadraturePoints(type1);
    UInt nb_quad2 = this->model->getFEM().getNbQuadraturePoints(type2);
    UInt nb_tot_quad =  nb_quad1 * elem_filter.getSize();;

    Vector<Real> & quads_volumes = quadrature_points_volumes(type1, ghost_type1);
    quads_volumes.resize(nb_tot_quad);
    quads_volumes.clear();

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
      (*weight)(0) = this->weight_func->operator()(r, q1, q2);
      if(_q1 != _q2)
	(*weight)(1) = this->weight_func->operator()(r, q2, q1);
      else
	(*weight)(1) = 0;

      quads_volumes(_q1) += (*weight)(0);
      quads_volumes(_q2) += (*weight)(1);
    }
  }

  //normalize the weights
  first_pair_types = existing_pairs.begin();
  for (; first_pair_types != last_pair_types; ++first_pair_types) {
    ElementType type1 = first_pair_types->first;
    ElementType type2 = first_pair_types->second;

    const Vector<UInt> & pairs = pair_list(type1, ghost_type1)(type2, ghost_type2);
    Vector<Real> & weights = pair_weight(type1, ghost_type1)(type2, ghost_type2);

    Vector<Real> & quads_volumes = quadrature_points_volumes(type1, ghost_type1);

    Vector<UInt>::const_iterator< types::Vector<UInt> > first_pair = pairs.begin(2);
    Vector<UInt>::const_iterator< types::Vector<UInt> > last_pair  = pairs.end(2);
    Vector<Real>::iterator<types::RVector> weight  = weights.begin(2);

    for(;first_pair != last_pair; ++first_pair, ++weight) {
      UInt q1 = (*first_pair)(0);
      UInt q2 = (*first_pair)(1);
      (*weight)(0) *= 1. / quads_volumes(q1);
      (*weight)(1) *= 1. / quads_volumes(q2);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
bool MaterialNonLocal<spatial_dimension, WeightFunction>::setParam(const std::string & key, const std::string & value,
						__attribute__((unused)) const ID & id) {
  std::stringstream sstr(value);
  if(key == "radius") { sstr >> radius; }
  else if(key == "UpdateWeights") { sstr >> update_weigths; }
  else if(!weight_func->setParam(key, value)) return false;
  return true;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::savePairs(const std::string & filename) const {
  std::ofstream pout;
  pout.open(filename.c_str());

  std::set< std::pair<ElementType, ElementType> >::const_iterator first_pair_types = existing_pairs.begin();
  std::set< std::pair<ElementType, ElementType> >::const_iterator last_pair_types = existing_pairs.end();

  GhostType ghost_type1, ghost_type2;
  ghost_type1 = ghost_type2 = _not_ghost;

  for (; first_pair_types != last_pair_types; ++first_pair_types) {
    const Vector<UInt> & pairs =
    pair_list(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);
    const Vector<Real> & weights =
    pair_weight(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);

    pout << "Types : " << first_pair_types->first << " " << first_pair_types->second << std::endl;

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

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::neighbourhoodStatistics(const std::string & filename) const {
  std::ofstream pout;
  pout.open(filename.c_str());

  std::set< std::pair<ElementType, ElementType> >::const_iterator first_pair_types = existing_pairs.begin();
  std::set< std::pair<ElementType, ElementType> >::const_iterator last_pair_types = existing_pairs.end();

  const Mesh & mesh = this->model->getFEM().getMesh();

  ByElementTypeUInt nb_neighbors("nb_neighbours", id, memory_id);
  initInternalVector(nb_neighbors, 1);
  resizeInternalVector(nb_neighbors);

  GhostType ghost_type1, ghost_type2;
  ghost_type1 = ghost_type2 = _not_ghost;

  for (; first_pair_types != last_pair_types; ++first_pair_types) {
    const Vector<UInt> & pairs =
    pair_list(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);

    pout << "Types : " << first_pair_types->first << " " << first_pair_types->second << std::endl;
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
  Real mean_nb_neig = 0;
  UInt max_nb_neig = 0;
  UInt min_nb_neig = std::numeric_limits<UInt>::max();
  for(; it != last_type; ++it) {
    Vector<UInt> & nb_neighor = nb_neighbors(*it, ghost_type1);
    Vector<UInt>::iterator<UInt> nb_neigh = nb_neighor.begin();
    Vector<UInt>::iterator<UInt> end_neigh  = nb_neighor.end();

    for (; nb_neigh != end_neigh; ++nb_neigh, ++nb_quads) {
      UInt nb = *nb_neigh;
      mean_nb_neig += nb;
      max_nb_neig = std::max(max_nb_neig, nb);
      min_nb_neig = std::min(min_nb_neig, nb);
    }
  }

  mean_nb_neig /= Real(nb_quads);

  pout << "Nb quadrature points: " << nb_quads << std::endl;
  pout << "Average nb neighbors: " << mean_nb_neig << std::endl;
  pout << "Max nb neighbors:     " << max_nb_neig << std::endl;
  pout << "Min nb neighbors:     " << min_nb_neig << std::endl;

  pout.close();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material<_non_local> [" << std::endl;
  stream << space << " + Radius          : " << radius << std::endl;
  stream << space << " + UpdateWeights   : " << update_weigths << std::endl;
  stream << space << " + Weight Function : " << *weight_func << std::endl;
  stream << space << "]" << std::endl;
}
