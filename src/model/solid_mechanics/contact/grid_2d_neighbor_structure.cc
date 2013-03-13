/**
 * @file   grid_2d_neighbor_structure.cc
 *
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 *
 * @date   Thu Dec 09 16:55:22 2010
 *
 * @brief  Grid2dNeighborStructure implementation
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
#include "contact_search.hh"
#include "solid_mechanics_model.hh"
#include "grid_2d_neighbor_structure.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
Grid2dNeighborStructure::Grid2dNeighborStructure(const ContactSearch & contact_search,
						 const Surface & master_surface,
						 const ContactNeighborStructureType & type,
						 const ContactNeighborStructureID & id) :
  ContactNeighborStructure(contact_search, master_surface, type, id),
  mesh(contact_search.getContact().getModel().getFEM().getMesh()) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = contact_search.getContact().getModel().getFEM().getMesh().getSpatialDimension();
  if(spatial_dimension != 2)
    AKANTU_DEBUG_ERROR("Wrong ContactType for contact in 2d!");

  /// make sure that the increments are computed
  const_cast<SolidMechanicsModel &>(contact_search.getContact().getModel()).setIncrementFlagOn();

  /// set increment to zero
  max_increment[0] = 0.;
  max_increment[1] = 0.;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Grid2dNeighborStructure::~Grid2dNeighborStructure() {

  AKANTU_DEBUG_IN();

  delete neighbor_list;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Grid2dNeighborStructure::initNeighborStructure() {

  AKANTU_DEBUG_IN();

  std::stringstream sstr; sstr << id << ":neighbor_list";
  neighbor_list = new NeighborList(sstr.str());
  createGrid(true);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool Grid2dNeighborStructure::check() {

  AKANTU_DEBUG_IN();

  UInt nb_surfaces = mesh.getNbSurfaces();
  Real * inc_val = contact_search.getContact().getModel().getIncrement().values;

  Real max[2] = {0. ,0.};
  UInt * surf_nodes_off = contact_search.getContact().getSurfaceToNodesOffset().values;
  UInt * surf_nodes = contact_search.getContact().getSurfaceToNodes().values;

  for (UInt s = 0; s < nb_surfaces; ++s)
    for (UInt n = surf_nodes_off[s]; n < surf_nodes_off[n+1]; ++n) {
      UInt i_node = surf_nodes[n];
      for (UInt i = 0; i < 2; ++i)
	max[i] = std::max(max[i], fabs(inc_val[2*i_node+i]));
    }

  for (UInt i = 0; i < 2; ++i) {
    max_increment[i] += max[i];
    if(max_increment[i] > spacing) {
      AKANTU_DEBUG_OUT();
      return true;
    }
  }

  AKANTU_DEBUG_OUT();
  return false;
}


/* -------------------------------------------------------------------------- */
void Grid2dNeighborStructure::update() {

  AKANTU_DEBUG_IN();

  delete[] this->neighbor_list;
  std::stringstream sstr; sstr << id << ":neighbor_list";
  neighbor_list = new NeighborList(sstr.str());
  createGrid(false);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Grid2dNeighborStructure::createGrid(bool initial_position) {

  AKANTU_DEBUG_IN();

  Real * coord;
  if (initial_position)
    coord = mesh.getNodes().values;
  else
    coord = contact_search.getContact().getModel().getCurrentPosition().values;

  UInt nb_surfaces = mesh.getNbSurfaces();


  /// get the size of the grid
  spacing = getMinSize(coord);

  /// get bounds of master surface
  Real * x_bounds = new Real[2*nb_surfaces];
  Real * y_bounds = new Real[2*nb_surfaces];
  getBounds(coord, x_bounds, y_bounds);


  /// find intersections between mastersurface and the other ones
  Real x_intersection[2];
  Real y_intersection[2];
  bool intersection = getBoundsIntersection(x_bounds, y_bounds, x_intersection, y_intersection);

  /// if every slave surfaces to far away from master surface neighbor list is empty
  if(intersection == false)
    return;

  /// adjust intersection space
  x_intersection[0] -= 1.49999*spacing;
  x_intersection[1] += 1.49999*spacing;
  y_intersection[0] -= 1.49999*spacing;
  y_intersection[1] += 1.49999*spacing;

  x_intersection[0] = std::max(x_intersection[0], x_bounds[2*master_surface]);
  x_intersection[1] = std::min(x_intersection[0+1], x_bounds[2*master_surface+1]);
  y_intersection[0] = std::max(y_intersection[0], y_bounds[2*master_surface]);
  y_intersection[1] = std::min(y_intersection[0+1], y_bounds[2*master_surface+1]);


  /// define grid dimension
  UInt nb_cells[3];
  nb_cells[0] = std::ceil(fabs(x_intersection[1]-x_intersection[0])/spacing);
  nb_cells[1] = std::ceil(fabs(y_intersection[1]-y_intersection[0])/spacing);
  nb_cells[2] = nb_cells[0]*nb_cells[1];

  Real origin[2];
  origin[0] = x_intersection[0] - (nb_cells[0]*spacing-fabs(x_intersection[1]-x_intersection[0]))/2.;
  origin[1] = y_intersection[0] - (nb_cells[1]*spacing-fabs(y_intersection[1]-y_intersection[0]))/2.;


  /// fill grid cells with master segments
  Array<UInt>  cell_to_segments(0, 1);
  UInt * cell_to_segments_offset = new UInt[nb_cells[2]+1];
  memset(cell_to_segments_offset, 0, (nb_cells[2]+1)*sizeof(UInt));
  traceSegments(coord, origin, nb_cells, cell_to_segments_offset, cell_to_segments);


  /// loop over slave nodes to find out impactor ones
  UInt * surf_nodes_off = contact_search.getContact().getSurfaceToNodesOffset().values;
  UInt * surf_nodes     = contact_search.getContact().getSurfaceToNodes().values;
  //  UInt * impactors      = neighbor_list->impactor_nodes.values;
  UInt * cell_seg_val   = cell_to_segments.values;

  ElementType el_type = _segment_2; /* Only linear element at the moment */
  UInt * conn_val = contact_search.getContact().getModel().getFEM().getMesh().getConnectivity(el_type, _not_ghost).values;
  UInt elem_nodes = Mesh::getNbNodesPerElement(el_type);
  //  std::stringstream sstr_fo; sstr_fo << id << ":facets_offset:" << el_type;
  neighbor_list->facets_offset.alloc(0, 1, el_type, _not_ghost);
  //  std::stringstream sstr_f; sstr_f << id << ":facets:" << el_type;
  neighbor_list->facets.alloc(0, 1, el_type, _not_ghost);
  neighbor_list->facets_offset(el_type, _not_ghost).push_back((UInt)0);

  const Int cell_index[9][2] = {{0,-1},{-1,-1},{-1,0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1}};
  Int nb_x = nb_cells[0];
  Int nb_y = nb_cells[1];

  Array<UInt> checked(0,1);

  /// loop over slave surfaces
  for (UInt s = 0; s < nb_surfaces; ++s) {

    if(s == master_surface)
      continue;

    /// check is slave surface inside grid
    if(x_bounds[2*s+1] < origin[0] || y_bounds[2*s+1] < origin[1])
      continue;
    if(x_bounds[2*s] > origin[0]+spacing*nb_cells[0] || y_bounds[2*s] > origin[1]+spacing*nb_cells[1])
      continue;

    /// loop over slave nodes of considered surface
    for (UInt n = surf_nodes_off[s]; n < surf_nodes_off[s+1]; ++n) {
      UInt i_node = surf_nodes[n];
      Real x_node = (coord[2*i_node] - origin[0])/spacing;
      Real y_node = (coord[2*i_node+1] - origin[1])/spacing;

      Int ii = (int)(x_node);
      Int jj = (int)(y_node);

      /// is node in one grid cell?
      if((ii < 0 || jj < 0) || (ii>=nb_x || jj>=nb_y))
	continue;

      Int i_start; /* which cells are to be considered */
      if(std::ceil(x_node-ii-0.5)==1 && std::ceil(y_node-jj-0.5)==0)
	i_start = 3;
      else {
	i_start = std::ceil(x_node-ii-0.5) + std::ceil(y_node-jj-0.5);
      }

      checked.empty();
      bool stored = false;
      /// look at actual cell and 3 adjacent cells
      for (Int k = -1; k < 3; ++k) {

	Int x_index = ii;
	Int y_index = jj;
	if ( k >= 0) {
	  x_index += cell_index[i_start*2+k][0];
	  y_index += cell_index[i_start*2+k][1];

	  if ((x_index < 0 || y_index < 0) || (x_index>=nb_x || y_index>=nb_y))
	    continue; /* cell does not exists */
	}

	UInt i_cell = (UInt)(x_index+y_index*nb_x);

	/// loop over segments related to the considered cell
	for (UInt el = cell_to_segments_offset[i_cell]; el < cell_to_segments_offset[i_cell+1]; ++el) {

	  UInt i_segment = cell_seg_val[el];
	  bool i_checked = false;
	  for (UInt l=0; l<checked.getSize(); l++) {
	    if(i_segment == checked.values[l]) { /* Segment already visited */
	      i_checked = true;
	      break;
	    }
	  }

	  /// check segment-node distance
	  if (i_checked == false) {
	    checked.push_back(i_segment);
	    Real * x1 = &coord[2*conn_val[i_segment*elem_nodes]];
	    Real * x2 = &coord[2*conn_val[i_segment*elem_nodes+1]];
	    Real * x3 = &coord[2*i_node];

	    Real vec_surf[2];
	    Real vec_dist[2];
	    Math::vector_2d(x1, x2, vec_surf);
	    Math::vector_2d(x1, x3, vec_dist);

	    Real length = Math::distance_2d(x1, x2);
	    Real proj = Math::vectorDot2(vec_surf, vec_dist)/(length*length);

	    if(proj < 0. || proj >1.) {
	      Real dist;
	      if(proj < 0.)
		dist =  Math::distance_2d(x1, x3);
	      else
		dist = Math::distance_2d(x2, x3);
	      if(dist < 0.5*spacing) {
		if (stored == false) {
		  stored = true;
		  neighbor_list->impactor_nodes.push_back(i_node);
		  neighbor_list->facets_offset(el_type, _not_ghost).push_back((UInt)0);
		}
		neighbor_list->facets_offset(el_type, _not_ghost)(neighbor_list->impactor_nodes.getSize())++;
		neighbor_list->facets(el_type, _not_ghost).push_back(i_segment);
	      }
	    }

	    else {
	      Real vec_normal[2];
	      Math::normal2(vec_surf, vec_normal);
	      Real gap = Math::vectorDot2(vec_dist, vec_normal);

	      if(fabs(gap) < 0.5*spacing) {
		if (stored == false) {
		  stored = true;
		  neighbor_list->impactor_nodes.push_back(i_node);
		  neighbor_list->facets_offset(el_type, _not_ghost).push_back((UInt)0);
		}
		neighbor_list->facets_offset(el_type, _not_ghost)(neighbor_list->impactor_nodes.getSize())++;
		neighbor_list->facets(el_type, _not_ghost).push_back(i_segment);
	      }
	    }
	  }
	} /* end loop segments */
      } /* end loop over cells */

    } /* end loop over slave nodes */
  } /* end loop over surfaces */

  /// convert occurence array in a csv one
  UInt * facet_off = neighbor_list->facets_offset(el_type, _not_ghost).storage();
  for (UInt i = 1; i < neighbor_list->facets_offset(el_type, _not_ghost).getSize(); ++i) facet_off[i] += facet_off[i-1];

  /// free temporary vectors
  delete[] x_bounds;
  delete[] y_bounds;
  delete[] cell_to_segments_offset;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real Grid2dNeighborStructure::getMinSize(Real * coord) {

  AKANTU_DEBUG_IN();

  ElementType el_type = _segment_2; /* Only linear element at the moment */
  UInt * conn_val = mesh.getConnectivity(el_type, _not_ghost).values;

  UInt nb_elements = mesh.getConnectivity(el_type, _not_ghost).getSize();
  UInt * surface_id_val = mesh.getSurfaceID(el_type, _not_ghost).values;
  UInt elem_nodes = Mesh::getNbNodesPerElement(el_type);
  Real min_size = std::numeric_limits<Real>::max();

  /// loop over master segments
  for (UInt el = 0; el < nb_elements; ++el)
    if (surface_id_val[el] == master_surface) {
      Real * x1 = &coord[2*conn_val[el*elem_nodes]];
      Real * x2 = &coord[2*conn_val[el*elem_nodes+1]];
      Real length = (x1[0]-x2[0])*(x1[0]-x2[0])+(x1[1]-x2[1])*(x1[1]-x2[1]);
      if (length < min_size)
	min_size = length;
    }

  return sqrt(min_size);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Grid2dNeighborStructure::getBounds(Real * coord, Real * x_bounds, Real * y_bounds) {

  AKANTU_DEBUG_IN();

  /// initialize bounds
  UInt nb_surfaces = mesh.getNbSurfaces();
  for (UInt s = 0; s < nb_surfaces; ++s) {
    x_bounds[s*2] = std::numeric_limits<Real>::max();
    x_bounds[s*2+1] = -std::numeric_limits<Real>::max();
    y_bounds[s*2] = std::numeric_limits<Real>::max();
    y_bounds[s*2+1] = -std::numeric_limits<Real>::max();
  }

  UInt * surf_nodes_off = contact_search.getContact().getSurfaceToNodesOffset().values;
  UInt * surf_nodes     = contact_search.getContact().getSurfaceToNodes().values;

  /// find min and max coordinates of each surface
  for (UInt s = 0; s < nb_surfaces; ++s) {
    for (UInt n = surf_nodes_off[s]; n < surf_nodes_off[s+1]; ++n) {

      UInt i_node = surf_nodes[n];

      if(coord[i_node*2] < x_bounds[s*2])
	x_bounds[s*2] = coord[i_node*2];

      if(coord[i_node*2] > x_bounds[2*s+1])
	x_bounds[2*s+1] = coord[i_node*2];

      if(coord[i_node*2+1] < y_bounds[2*s])
	y_bounds[2*s] = coord[i_node*2+1];

      if(coord[i_node*2+1] > y_bounds[2*s+1])
	y_bounds[2*s+1] = coord[i_node*2+1];
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool Grid2dNeighborStructure::getBoundsIntersection(Real * x_bounds, Real * y_bounds, Real * x_int, Real * y_int) {

  AKANTU_DEBUG_IN();

  bool find_intersection = false;
  UInt nb_surfaces = mesh.getNbSurfaces();

  x_int[0] = std::numeric_limits<Real>::max();
  x_int[1] = -std::numeric_limits<Real>::max();
  y_int[0] = std::numeric_limits<Real>::max();
  y_int[1] = -std::numeric_limits<Real>::max();

  /// enlarge bounds of master surface within a tolerance
  x_bounds[2*master_surface]   -= 1.49999*spacing;
  x_bounds[2*master_surface+1] += 1.49999*spacing;
  y_bounds[2*master_surface]   -= 1.49999*spacing;
  y_bounds[2*master_surface+1] += 1.49999*spacing;

  /// loop over surfaces
  for (UInt s = 0; s < nb_surfaces; ++s) {

    if(s == master_surface)
      continue;

    Real x_temp[2] = {REAL_INIT_VALUE, REAL_INIT_VALUE};
    Real y_temp[2] = {REAL_INIT_VALUE, REAL_INIT_VALUE};

    /// find if intersection in x exists
    if(x_bounds[2*master_surface] > x_bounds[2*s]) {

      /* starting point of possible intersection */
      x_temp[0] = x_bounds[2*master_surface];

      /* find ending point of possible intersection */
      if(x_bounds[2*s+1] < x_bounds[2*master_surface])
	continue;
      else if(x_bounds[2*s+1] < x_bounds[2*master_surface+1])
	x_temp[1] = x_bounds[2*s+1];
      else
	x_temp[1] = x_bounds[2*master_surface+1];
    }
    else {

      /* starting point of possible intersection */
      x_temp[0] = x_bounds[2*s];

      /* find ending point of possible intersection */
      if(x_bounds[2*master_surface+1] < x_bounds[2*s])
	continue;
      else if(x_bounds[2*master_surface+1] < x_bounds[2*s+1])
	x_temp[1] = x_bounds[2*master_surface+1];
      else
	x_temp[1] = x_bounds[2*s+1];
    }

    /// find if intersection in y exists
    if(y_bounds[2*master_surface] > y_bounds[2*s]) {

      /* starting point of possible intersection */
      y_temp[0] = y_bounds[2*master_surface];

      /* find ending point of possible intersection */
      if(y_bounds[2*s+1] < y_bounds[2*master_surface])
	continue;
      else if(y_bounds[2*s+1] < y_bounds[2*master_surface+1])
	y_temp[1] = y_bounds[2*s+1];
      else
	y_temp[1] = y_bounds[2*master_surface+1];
    }
    else {

      /* starting point of possible intersection */
      y_temp[0] = y_bounds[2*s];

      /* find ending point of possible intersection */
      if(y_bounds[2*master_surface+1] < y_bounds[2*s])
	continue;
      else if(y_bounds[2*master_surface+1] < y_bounds[2*s+1])
	y_temp[1] = y_bounds[2*master_surface+1];
      else
	y_temp[1] = y_bounds[2*s+1];
    }

    /// intersection exists, find its minimum/maximum
    find_intersection = true;
    x_int[0] = std::min(x_int[0], x_temp[0]);
    x_int[1] = std::max(x_int[1], x_temp[1]);
    y_int[0] = std::min(y_int[0], y_temp[0]);
    y_int[1] = std::max(y_int[1], y_temp[1]);
  }

  return find_intersection;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Grid2dNeighborStructure::traceSegments(Real * coord, Real * origin,
					    UInt * nb_cells, UInt * cell_to_seg_off, Array<UInt> & cell_to_segments) {

  AKANTU_DEBUG_IN();

  Array<UInt> temp(0, 2);
  UInt index[2] = {0, 0};
  ElementType el_type = _segment_2; /* Only linear element at the moment */
  UInt elem_nodes = Mesh::getNbNodesPerElement(el_type);
  UInt * surface_id_val = mesh.getSurfaceID(el_type, _not_ghost).values;
  UInt * conn_val = mesh.getConnectivity(el_type, _not_ghost).values;
  UInt nb_segments = mesh.getConnectivity(el_type, _not_ghost).getSize();

  Int nb_x = nb_cells[0];
  Int nb_y = nb_cells[1];

  /// Detect cells which are crossed by the segments
  for (UInt el = 0; el < nb_segments; ++el)
    if (surface_id_val[el] == master_surface) {

      /* Compute relative coords and index of the staring and ending cell (in respect to the grid) */

      Real x1 = (coord[2*conn_val[elem_nodes*el]] - origin[0])/spacing;
      Real y1 = (coord[2*conn_val[elem_nodes*el]+1] - origin[1])/spacing;
      Real x2 = (coord[2*conn_val[elem_nodes*el+1]] - origin[0])/spacing;
      Real y2 = (coord[2*conn_val[elem_nodes*el+1]+1] - origin[1])/spacing;

      Int ii = (Int)std::floor(x1);
      Int jj = (Int)std::floor(y1);
      Int kk = (Int)std::floor(x2);
      Int ll = (Int)std::floor(y2);

      /// check if segment outside grid
      if((std::max(ii,kk)<0 || std::min(ii,kk)>=nb_x) || (std::max(jj,ll)<0 || std::min(jj,ll)>=nb_y))
      	continue;

      /* */

      Real dx = fabs(x2-x1);
      Real dy = fabs(y2-y1);
      Int n = 0;
      Int ii_inc, jj_inc;
      Real error;


      if(dx == 0.) {
	ii_inc = 0;
	error =  std::numeric_limits<Real>::max();
      }
      else if(x2 > x1) {
	ii_inc = 1;
	n = kk-ii;
	error = (ii+1-x1)*dy;
      }
      else { /* x1 < x2 */
	ii_inc = -1;
	n = ii-kk;
	error = (x1-ii)*dy;
      }

      if(dy == 0.) {
	jj_inc = 0;
	error =  -std::numeric_limits<Real>::max();
      }
      else if(y2 > y1) {
	jj_inc = 1;
	n += ll-jj;
	error -= (jj+1-y1)*dx;
      }
      else { /* y1 < y2 */
	jj_inc = -1;
	n += jj-ll;
	error -= (y1-jj)*dx;
      }

      while(true) {

	/* check if intersected cell belong to the grid */
	if((ii>=0 && ii< nb_x) && (jj>=0 && jj<nb_y)) {
	  /* Increase offset */
	  index[0] = (UInt)(ii+jj*nb_x);
	  index[1] = el;
	  cell_to_seg_off[index[0]]++;
	  temp.push_back(index);
	}

	if(n==0)
	  break;

	/* Move to the next cell */
	if(error > 0) { /* ..move in y-direction */
	  jj += jj_inc;
	  error -= dx;
	}
	else { /* ..move in x-direction */
	  ii += ii_inc;
	  error += dy;
	}

	n--;
      }
    }

  /// convert the occurrence array in a csr one
  for (UInt i = 1; i < nb_cells[2]; ++i) cell_to_seg_off[i] += cell_to_seg_off[i-1];
  for (UInt i = nb_cells[2]; i > 0; --i) cell_to_seg_off[i] = cell_to_seg_off[i-1];
  cell_to_seg_off[0] = 0;

  /// rearrange segments to get the cell-segments list
  cell_to_segments.resize(cell_to_seg_off[nb_cells[2]]);
  UInt * cell_val = cell_to_segments.values;
  UInt * tmp = temp.values;
  UInt nb_traced = temp.getSize();
  for (UInt i = 0; i < nb_traced; i++) {
    cell_val[cell_to_seg_off[tmp[2*i]]] = tmp[2*i+1];
    cell_to_seg_off[tmp[2*i]]++;
  }

  /// reconvert occurence array in a csr one
  for (UInt i = nb_cells[2]; i > 0; --i) cell_to_seg_off[i] = cell_to_seg_off[i-1];
  cell_to_seg_off[0] = 0;

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
