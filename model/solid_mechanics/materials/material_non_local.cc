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
#include "aka_grid.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MaterialNonLocal::MaterialNonLocal(Model & model, const ID & id)  :
  Material(model, id),
  quadrature_points_coordinates("quadrature_points_coordinates", id, memory_id) {
  AKANTU_DEBUG_IN();

  this->model->getFEM().getMesh().initByElementTypeVector(quadrature_points_coordinates,
							  spatial_dimension, 0);

  is_non_local = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MaterialNonLocal::~MaterialNonLocal() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialNonLocal::initMaterial() {
  AKANTU_DEBUG_IN();
  //  Material::initMaterial();

  updatePairList();
  computeWeights();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialNonLocal::updatePairList() {
  AKANTU_DEBUG_IN();

  const Real safety_factor = 1.2; // for the cell grid spacing
  const Mesh & mesh = model->getFEM().getMesh();
  const_cast<Mesh &>(mesh).computeBoundingBox();

  Real lower_bounds[spatial_dimension];
  Real upper_bounds[spatial_dimension];
  mesh.getLowerBounds(lower_bounds);
  mesh.getUpperBounds(upper_bounds);

  Real spacing[spatial_dimension];
  for (UInt i = 0; i < spatial_dimension; ++i) {
    spacing[i] = radius * safety_factor;
  }


  RegularGrid<QuadraturePoint> cell_list(spatial_dimension, lower_bounds, upper_bounds, spacing);

  GhostType ghost_type = _not_ghost;

  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;


  // first generate the quad points coordinate and count the number of points per cell
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    Vector<UInt> & elem_filter = element_filter(*it, ghost_type);

    UInt nb_element  = elem_filter.getSize();
    UInt nb_tot_quad = model->getFEM().getNbQuadraturePoints(*it, ghost_type) * nb_element;

    Vector<Real> & quads = quadrature_points_coordinates(*it, ghost_type);
    quads.resize(nb_tot_quad);

    model->getFEM().interpolateOnQuadraturePoints(mesh.getNodes(),
						  quads, spatial_dimension,
						  *it, ghost_type, &elem_filter);

    Vector<Real>::iterator<types::RVector> first_quad = quads.begin(spatial_dimension);
    Vector<Real>::iterator<types::RVector> last_quad = quads.end(spatial_dimension);
    for(;first_quad != last_quad; ++first_quad) {
      cell_list.count(*first_quad);
    }
  }


  QuadraturePoint q;
  q.ghost_type = ghost_type;

  // second insert the point in the cells
  cell_list.beginInsertions();
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    Vector<UInt> & elem_filter = element_filter(*it, ghost_type);
    UInt nb_element = elem_filter.getSize();
    UInt nb_quad    = model->getFEM().getNbQuadraturePoints(*it, ghost_type);

    Vector<Real> & quads = quadrature_points_coordinates(*it, ghost_type);

    q.type = *it;

    Vector<Real>::iterator<types::RVector> quad = quads.begin(spatial_dimension);
    UInt * elem = elem_filter.storage();

    for (UInt e = 0; e < nb_element; ++e) {
      q.element = *elem;
      for (UInt nq = 0; nq < nb_quad; ++nq) {
	q.num_point = nq;
	q.setPosition(*quad);
	cell_list.insert(q, *quad);
	++quad;
      }
      ++elem;
    }
  }
  cell_list.endInsertions();


  // generate the pair of neighbor depending of the cell_list
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    // Preparing datas
    Vector<Real> & quads = quadrature_points_coordinates(*it, ghost_type);
    Vector<Real>::iterator<types::RVector> first_quad = quads.begin(spatial_dimension);
    Vector<Real>::iterator<types::RVector> last_quad  = quads.end(spatial_dimension);

    ByElementTypeUInt & pairs = pair_list(ByElementTypeUInt("pairs", id, memory_id),
					  *it,
					  ghost_type);

    ElementType current_element_type = _not_defined;
    GhostType current_ghost_type = _casper;
    Real * neigh_quad_positions = NULL;
    Vector<UInt> * neighbors = NULL;

    UInt my_num_quad = 0;

    // loop over quad points

    for(;first_quad != last_quad; ++first_quad, ++my_num_quad) {
      UInt cell = cell_list.getCell(*first_quad);

      RegularGrid<QuadraturePoint>::neighbor_cells_iterator first_neigh_cell =
	cell_list.beginNeighborCells(cell);
      RegularGrid<QuadraturePoint>::neighbor_cells_iterator last_neigh_cell =
	cell_list.endNeighborCells(cell);

      // loop over neighbor  cells of the one containing  the current quadrature
      // point
      for (; first_neigh_cell != last_neigh_cell; ++first_neigh_cell) {
	RegularGrid<QuadraturePoint>::iterator first_neigh_quad =
	  cell_list.beginCell(*first_neigh_cell);
	RegularGrid<QuadraturePoint>::iterator last_neigh_quad =
	  cell_list.endCell(*first_neigh_cell);

	// loop over the quadrature point in the current cell of the cell list
	for (;first_neigh_quad != last_neigh_quad; ++first_neigh_quad){
	  QuadraturePoint & quad = *first_neigh_quad;
	  UInt nb_quad_per_elem =  model->getFEM().getNbQuadraturePoints(quad.type, quad.ghost_type);
	  UInt neigh_num_quad = quad.element * nb_quad_per_elem + quad.num_point;

	  // store only the pair of type (q1, q2) with q1 < q2
	  //	  if(my_num_quad != neigh_num_quad) {

	    // little optimization to not search in the map at each quad points
	    if(quad.type != current_element_type || quad.ghost_type != current_ghost_type) {
	      neigh_quad_positions = quadrature_points_coordinates(quad.type,
								   quad.ghost_type).storage();
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
	    // 			      spatial_dimension);
	    const types::RVector & neigh_quad = quad.getPosition();

	    Real distance = first_quad->distance(neigh_quad);
	    if(distance <= radius) {
	      UInt pair[2];
	      pair[0] = my_num_quad; pair[1] = neigh_num_quad;
	      neighbors->push_back(pair);
	    }
	    //	  }
	}
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialNonLocal::computeQuadraturePointsNeighborhoudVolumes(ByElementTypeReal & volumes) const {
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
}

/* -------------------------------------------------------------------------- */
void MaterialNonLocal::computeWeights() {
  std::set< std::pair<ElementType, ElementType> >::iterator first_pair_types = existing_pairs.begin();
  std::set< std::pair<ElementType, ElementType> >::iterator last_pair_types = existing_pairs.end();

  GhostType ghost_type1, ghost_type2;
  ghost_type1 = ghost_type2 = _not_ghost;

  ByElementTypeReal quadrature_points_volumes("quadrature_points_volumes", id, memory_id);
  model->getFEM().getMesh().initByElementTypeVector(quadrature_points_volumes, 1, 0);

  //  computeQuadraturePointsNeighborhoudVolumes(quadrature_points_volumes);

  Real R_2 = 1. / (radius * radius);
  for (; first_pair_types != last_pair_types; ++first_pair_types) {
    Vector<UInt> & pairs =
      pair_list(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);

    ByElementTypeReal & weights_type_1 = pair_weigth(first_pair_types->first, ghost_type1);
    Vector<Real> * tmp_weight = NULL;
    if(!weights_type_1.exists(first_pair_types->second, ghost_type2)) {
      tmp_weight = &(weights_type_1.alloc(0, 1, first_pair_types->second, ghost_type2));
    } else {
      tmp_weight = &(weights_type_1(first_pair_types->second, ghost_type2));
    }

    Vector<Real> & weights = *tmp_weight;


    weights.resize(pairs.getSize());
    weights.clear();

    Vector<Real>::iterator<types::RVector> iquads1 =
      quadrature_points_coordinates(first_pair_types->first, ghost_type1).begin(spatial_dimension);
    Vector<Real>::iterator<types::RVector> iquads2 =
      quadrature_points_coordinates(first_pair_types->second, ghost_type2).begin(spatial_dimension);
    //    Vector<Real> & ivolumes = quadrature_points_volumes(first_pair_types->first, ghost_type1);

    Vector<UInt>::iterator< types::Vector<UInt> > first_pair = pairs.begin(2);
    Vector<UInt>::iterator< types::Vector<UInt> > last_pair  = pairs.end(2);

    Vector<Real>::iterator<Real> weight  = weights.begin<Real>();

    for(;first_pair != last_pair; ++first_pair, ++weight) {
      UInt q1 = (*first_pair)(0);
      UInt q2 = (*first_pair)(1);
      Real r = iquads1[q1].distance(iquads2[q2]);
      Real alpha = (1. - r*r * R_2);
      // alpha = alpha * alpha;
      *weight = alpha * alpha;
    }
  }

  first_pair_types = existing_pairs.begin();
  for (; first_pair_types != last_pair_types; ++first_pair_types) {
    const Vector<UInt> & pairs =
      pair_list(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);
    const Vector<Real> & weigths =
      pair_weigth(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);

    const Vector<UInt> & elem_filter = element_filter(first_pair_types->first, ghost_type1);
    UInt nb_element  = elem_filter.getSize();
    UInt nb_tot_quad = model->getFEM().getNbQuadraturePoints(first_pair_types->first, ghost_type1) * nb_element;

    Vector<Real> & quads_volumes = quadrature_points_volumes(first_pair_types->first, ghost_type1);
    quads_volumes.resize(nb_tot_quad);
    quads_volumes.clear();

    Vector<UInt>::const_iterator< types::Vector<UInt> > first_pair = pairs.begin(2);
    Vector<UInt>::const_iterator< types::Vector<UInt> > last_pair  = pairs.end(2);

    Real * pair_w = weigths.storage();

    for(;first_pair != last_pair; ++first_pair, ++pair_w) {
      UInt q1 = (*first_pair)(0);
      //      UInt q2 = (*first_pair)(1);
      quads_volumes(q1)  += *pair_w;
    }
  }

  first_pair_types = existing_pairs.begin();
  for (; first_pair_types != last_pair_types; ++first_pair_types) {
    const Vector<UInt> & pairs =
      pair_list(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);
    const Vector<Real> & weigths =
      pair_weigth(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);

    Vector<Real> & quads_volumes = quadrature_points_volumes(first_pair_types->first, ghost_type1);

    Vector<UInt>::const_iterator< types::Vector<UInt> > first_pair = pairs.begin(2);
    Vector<UInt>::const_iterator< types::Vector<UInt> > last_pair  = pairs.end(2);

    Real * pair_w = weigths.storage();

    for(;first_pair != last_pair; ++first_pair, ++pair_w) {
      UInt q1 = (*first_pair)(0);
      //      UInt q2 = (*first_pair)(1);
      *pair_w *= 1. / quads_volumes(q1);
    }
  }
}


/* -------------------------------------------------------------------------- */
bool MaterialNonLocal::setParam(const std::string & key, const std::string & value,
				const ID & id) {
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
