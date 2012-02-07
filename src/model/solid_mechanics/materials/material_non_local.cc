/**
 * @file   material_non_local.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Aug  2 11:21:58 2011
 *
 * @brief  implementation of the common parts for Non local materials
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
#include "material_non_local.hh"
#include "solid_mechanics_model.hh"
#include "fem.hh"
#include "grid_synchronizer.hh"


/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MaterialNonLocal::MaterialNonLocal(Model & model, const ID & id)  :
  Material(model, id), cell_list(NULL) {
  //  quadrature_points_coordinates("quadrature_points_coordinates", id, memory_id) {
  AKANTU_DEBUG_IN();

  // this->model->getFEM().getMesh().initByElementTypeVector(quadrature_points_coordinates,
  //                                                         spatial_dimension, 0);

  is_non_local = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MaterialNonLocal::~MaterialNonLocal() {
  AKANTU_DEBUG_IN();

  if(cell_list) { delete cell_list; };

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialNonLocal::initMaterial() {
  AKANTU_DEBUG_IN();
  //  Material::initMaterial();

  Mesh & mesh = model->getFEM().getMesh();
  ByElementTypeReal quadrature_points_coordinates("quadrature_points_coordinates", id);
  mesh.initByElementTypeVector(quadrature_points_coordinates,
                               spatial_dimension, 0);

  computeQuadraturePointsCoordinates(mesh.getNodes(), quadrature_points_coordinates);

  createCellList(quadrature_points_coordinates);

  updatePairList(quadrature_points_coordinates);
  computeWeights<BaseWeightFonction>(quadrature_points_coordinates);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialNonLocal::createCellList(const ByElementTypeReal & quadrature_points_coordinates) {
  const Real safety_factor = 1.2; // for the cell grid spacing
  Mesh & mesh = model->getFEM().getMesh();
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
    UInt nb_quad    = model->getFEM().getNbQuadraturePoints(*it, ghost_type);

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

  SynchronizerRegistry & synch_registry = model->getSynchronizerRegistry();
  std::stringstream sstr; sstr << id << ":grid_synchronizer";
  GridSynchronizer * synch = GridSynchronizer::createGridSynchronizer(mesh,
                                                                      *cell_list,
                                                                      sstr.str());
  synch_registry.registerSynchronizer(*synch, _gst_mnl_damage);
}

/* -------------------------------------------------------------------------- */
void MaterialNonLocal::updatePairList(const ByElementTypeReal & quadrature_points_coordinates) {
  AKANTU_DEBUG_IN();

  Mesh & mesh = model->getFEM().getMesh();

  GhostType ghost_type = _not_ghost;

  // generate the pair of neighbor depending of the cell_list
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
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
          UInt nb_quad_per_elem =  model->getFEM().getNbQuadraturePoints(quad.type, quad.ghost_type);
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
          if(distance <= radius) {
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
void MaterialNonLocal::computeQuadraturePointsNeighborhoudVolumes(ByElementTypeReal & volumes) const {
  AKANTU_DEBUG_IN();
  GhostType ghost_type = _not_ghost;

  ByElementTypeReal per_quadrature_points_volumes("per_quadrature_points_volumes", id, memory_id);
  model->getFEM().getMesh().initByElementTypeVector(per_quadrature_points_volumes, 1, 0);

  const Mesh::ConnectivityTypeList & type_list = model->getFEM().getMesh().getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;
  // first generate the quad points coordinate and count the number of points per cell
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    const Vector<UInt> & elem_filter = element_filter(*it, ghost_type);

    UInt nb_element  = elem_filter.getSize();
    UInt nb_tot_quad = model->getFEM().getNbQuadraturePoints(*it, ghost_type) * nb_element;

    Vector<Real> & quads_volumes = per_quadrature_points_volumes(*it, ghost_type);
    quads_volumes.resize(nb_tot_quad);

    Vector<Real> ones(nb_tot_quad, 1, 1.);
    //    Vector<Real> per_quads_volumes(nb_tot_quad, 1);

    model->getFEM().integrateOnQuadraturePoints(ones,
                                                quads_volumes,
                                                1, *it, ghost_type, &elem_filter);
  }

  accumulateOnNeighbours(per_quadrature_points_volumes, volumes, 1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<class WeightFunction>
void MaterialNonLocal::computeWeights(const ByElementTypeReal & quadrature_points_coordinates) {
  AKANTU_DEBUG_IN();

  std::set< std::pair<ElementType, ElementType> >::iterator first_pair_types;
  std::set< std::pair<ElementType, ElementType> >::iterator last_pair_types = existing_pairs.end();

  GhostType ghost_type1, ghost_type2;
  ghost_type1 = ghost_type2 = _not_ghost;

  ByElementTypeReal quadrature_points_volumes("quadrature_points_volumes", id, memory_id);
  model->getFEM().getMesh().initByElementTypeVector(quadrature_points_volumes, 1, 0);

  // Compute the weights
  first_pair_types = existing_pairs.begin();
  for (; first_pair_types != last_pair_types; ++first_pair_types) {
    ElementType type1 = first_pair_types->first;
    ElementType type2 = first_pair_types->second;

    const Vector<UInt> & pairs = pair_list(type1, ghost_type1)(type2, ghost_type2);

    ByElementTypeReal & weights_type_1 = pair_weigth(type1, ghost_type1);
    Vector<Real> * tmp_weight = NULL;
    if(!weights_type_1.exists(type2, ghost_type2)) {
      tmp_weight = &(weights_type_1.alloc(0, 1, type2, ghost_type2));
    } else {
      tmp_weight = &(weights_type_1(type2, ghost_type2));
    }
    Vector<Real> & weights = *tmp_weight;
    weights.resize(pairs.getSize());
    weights.clear();

    const Vector<UInt> & elem_filter = element_filter(type1, ghost_type1);
    UInt nb_tot_quad = model->getFEM().getNbQuadraturePoints(type1, ghost_type1) * elem_filter.getSize();;

    Vector<Real> & quads_volumes = quadrature_points_volumes(type1, ghost_type1);
    quads_volumes.resize(nb_tot_quad);
    quads_volumes.clear();

    Vector<Real>::const_iterator<types::RVector> iquads1;
    Vector<Real>::const_iterator<types::RVector> iquads2;
    iquads1 = quadrature_points_coordinates(type1, ghost_type1).begin(spatial_dimension);
    iquads2 = quadrature_points_coordinates(type2, ghost_type2).begin(spatial_dimension);

    Vector<UInt>::const_iterator< types::Vector<UInt> > first_pair = pairs.begin(2);
    Vector<UInt>::const_iterator< types::Vector<UInt> > last_pair  = pairs.end(2);
    Vector<Real>::iterator<Real> weight  = weights.begin();

    WeightFunction func(*this, radius, type1, ghost_type1, type2, ghost_type2);

    // Weight function
    for(;first_pair != last_pair; ++first_pair, ++weight) {
      UInt q1 = (*first_pair)(0);
      UInt q2 = (*first_pair)(1);
      Real r = iquads1[q1].distance(iquads2[q2]);
      *weight = func(r, q1, q2);

      quads_volumes(q1) += *weight;
    }
  }

  //normalize the weights
  first_pair_types = existing_pairs.begin();
  for (; first_pair_types != last_pair_types; ++first_pair_types) {
    ElementType type1 = first_pair_types->first;
    ElementType type2 = first_pair_types->second;

    const Vector<UInt> & pairs = pair_list(type1, ghost_type1)(type2, ghost_type2);
    Vector<Real> & weights = pair_weigth(type1, ghost_type1)(type2, ghost_type2);

    Vector<Real> & quads_volumes = quadrature_points_volumes(type1, ghost_type1);

    Vector<UInt>::const_iterator< types::Vector<UInt> > first_pair = pairs.begin(2);
    Vector<UInt>::const_iterator< types::Vector<UInt> > last_pair  = pairs.end(2);
    Vector<Real>::iterator<Real> weight  = weights.begin();

    for(;first_pair != last_pair; ++first_pair, ++weight) {
      UInt q1 = (*first_pair)(0);
      *weight *= 1. / quads_volumes(q1);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialNonLocal::removeDamaged(const ByElementTypeReal & damage, Real thresold) {
  AKANTU_DEBUG_IN();

  std::set< std::pair<ElementType, ElementType> >::iterator first_pair_types = existing_pairs.begin();
  std::set< std::pair<ElementType, ElementType> >::iterator last_pair_types = existing_pairs.end();

  GhostType ghost_type1, ghost_type2;
  ghost_type1 = ghost_type2 = _not_ghost;

  for (; first_pair_types != last_pair_types; ++first_pair_types) {
    Vector<UInt> & pairs =
    pair_list(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);

    Vector<Real> & weights =
    pair_weigth(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);

    const Vector<Real> & dam1 = damage(first_pair_types->first, ghost_type1);
    const Vector<Real> & dam2 = damage(first_pair_types->second, ghost_type2);

    Vector<UInt>::iterator< types::Vector<UInt> > first_pair = pairs.begin(2);
    Vector<Real>::iterator<Real> weight  = weights.begin();

    while(first_pair != pairs.end(2)) {
      UInt q1 = (*first_pair)(0);
      UInt q2 = (*first_pair)(1);

      if(dam1(q1) >= thresold || dam2(q2) >= thresold) {
        pairs.erase(first_pair);
        weights.erase(weight);
      } else {
        ++first_pair;
        ++weight;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialNonLocal::updateResidual(Vector<Real> & displacement, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  computeStress(displacement, ghost_type);
  computeNonLocalStress(ghost_type);

  assembleResidual(ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool MaterialNonLocal::setParam(const std::string & key, const std::string & value,
                                __attribute__((unused)) const ID & id) {
  std::stringstream sstr(value);
  if(key == "radius") { sstr >> radius; }
  else return false;
  return true;
}

/* -------------------------------------------------------------------------- */
void MaterialNonLocal::savePairs(const std::string & filename) const {
  std::ofstream pout;
  pout.open(filename.c_str());

  std::set< std::pair<ElementType, ElementType> >::const_iterator first_pair_types = existing_pairs.begin();
  std::set< std::pair<ElementType, ElementType> >::const_iterator last_pair_types = existing_pairs.end();

  GhostType ghost_type1, ghost_type2;
  ghost_type1 = ghost_type2 = _not_ghost;

  for (; first_pair_types != last_pair_types; ++first_pair_types) {
    const Vector<UInt> & pairs =
    pair_list(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);
    const Vector<Real> & weigths =
    pair_weigth(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);

    pout << "Types : " << first_pair_types->first << " " << first_pair_types->second << std::endl;

    Vector<UInt>::const_iterator< types::Vector<UInt> > first_pair = pairs.begin(2);
    Vector<UInt>::const_iterator< types::Vector<UInt> > last_pair  = pairs.end(2);
    Real * pair_w = weigths.storage();

    for(;first_pair != last_pair; ++first_pair, ++pair_w) {
      UInt q1 = (*first_pair)(0);
      UInt q2 = (*first_pair)(1);
      pout << q1 << " " << q2 << " "<< *pair_w << std::endl;
    }
  }
}

/* -------------------------------------------------------------------------- */
void MaterialNonLocal::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material<_non_local> [" << std::endl;
  stream << space << " + Radius                      : " << radius << std::endl;
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
__END_AKANTU__
