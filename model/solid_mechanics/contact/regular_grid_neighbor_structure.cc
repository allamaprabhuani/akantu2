/**
 * @file   regular_grid_neighbor_structure.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Oct 11 16:03:17 2010
 *
 * @brief  Specialization of the contact neighbor structure for regular grid
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
//#include <math.h>

/* -------------------------------------------------------------------------- */
#include "regular_grid_neighbor_structure.hh"

#include "contact_search.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension> 
RegularGridNeighborStructure<spatial_dimension>::RegularGridNeighborStructure(const ContactSearch & contact_search,
									      const Surface & master_surface,
									      const ContactNeighborStructureID & id) :
  ContactNeighborStructure(contact_search, master_surface, id), mesh(contact_search.getContact().getModel().getFEM().getMesh()) {
  
  AKANTU_DEBUG_IN();
  
  std::cout << "call of RegularGridNeighborStructure with master surface: " << master_surface << std::endl;

  //  spatial_dimension = mesh.getSpatialDimension();

  /// make sure that the increments are computed
  const_cast<SolidMechanicsModel &>(contact_search.getContact().getModel()).setIncrementFlagOn();

  /// arbitrary initial value
  grid_spacing[0] = 0.1;
  grid_spacing[1] = 0.1;
  grid_spacing[2] = 0.1;

  /// securty factor of max 0.5 is needed for a neighborlist that is always complete 
  security_factor[0] = 0.5;
  security_factor[1] = 0.5;
  security_factor[2] = 0.5;
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
RegularGridNeighborStructure<spatial_dimension>::~RegularGridNeighborStructure() {
  AKANTU_DEBUG_IN();

  

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension> 
void RegularGridNeighborStructure<spatial_dimension>::init() {
  AKANTU_DEBUG_IN();
  
  std::cout << "RegularGridNeighborStructure initialization with master surface: " << this->master_surface << std::endl;
  
  Real * node_coordinates = mesh.getNodes().values;
  this->update(node_coordinates);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension> 
void RegularGridNeighborStructure<spatial_dimension>::update() {
  AKANTU_DEBUG_IN();

  std::cout << "RegularGridNeighborStructure update with master surface: " << this->master_surface << std::endl;

  Real * node_current_position = contact_search.getContact().getModel().getCurrentPosition().values;
  this->update(node_current_position);
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension> 
void RegularGridNeighborStructure<spatial_dimension>::update(Real * node_position) {
  AKANTU_DEBUG_IN();

  UInt nb_surfaces = mesh.getNbSurfaces();

  AKANTU_DEBUG_ASSERT(master_surface < nb_surfaces, "Master surface (" << 
		      master_surface << ") out of surface range (number of surfaces: " << 
		      nb_surfaces << ") !!");
  
  //Real * node_position = contact_search.getContact().getModel().getCurrentPosition().values;

  Real max[nb_surfaces][spatial_dimension];
  Real min[nb_surfaces][spatial_dimension];

  /// initialize max and min table with extrem values
  for(UInt surf = 0; surf < nb_surfaces; ++surf) {
    for(UInt dim = 0; dim < spatial_dimension; ++dim) {
      max[surf][dim] = std::numeric_limits<Real>::min();
      min[surf][dim] = std::numeric_limits<Real>::max();
    }
  }

  // get nodes that are on a given surface
  UInt * surface_to_nodes_offset = contact_search.getContact().getSurfaceToNodesOffset().values;
  UInt * surface_to_nodes        = contact_search.getContact().getSurfaceToNodes().values;

  /// find max and min values of current position for each surface
  for(UInt surf = 0; surf < nb_surfaces; ++surf) {
    UInt min_surf_offset = surface_to_nodes_offset[surf];
    UInt max_surf_offset = surface_to_nodes_offset[surf + 1];
    std::cout << "surface: " << surf << std::endl;
    for(UInt n = min_surf_offset; n < max_surf_offset; ++n) {
      UInt cur_node = surface_to_nodes[n];
      std::cout << "current nodes: " << cur_node << " : ";
      for(UInt dim = 0; dim < spatial_dimension; ++dim) {
	Real cur_position = node_position[cur_node + dim];
	std::cout << cur_position << " ";
	max[surf][dim] = std::max(max[surf][dim], cur_position);
	min[surf][dim] = std::min(min[surf][dim], cur_position);
      }
      std::cout << std::endl;
    } 
  }

  // print for error search
  for(UInt surf = 0; surf < nb_surfaces; ++surf) {
    for(UInt dim = 0; dim < spatial_dimension; ++dim) {
      std::cout << "surface " << surf << " has min: " << min[surf][dim] << " and max: " << max[surf][dim] << " for direction " << dim << std::endl;
    }
  }

  /// define grid geometry around the master surface
  Real grid_min[spatial_dimension];
  Real grid_max[spatial_dimension];
  UInt directional_nb_cells[spatial_dimension];
  UInt nb_cells = 0;
  
  for(UInt dim = 0; dim < spatial_dimension; ++dim) {
    Real grid_length = max[master_surface][dim] - min[master_surface][dim];

    /// get nb of cells needed to cover total length and add a cell to each side (start, end)
    directional_nb_cells[dim] = static_cast<UInt>(ceil(grid_length / grid_spacing[dim])) + 2;
    nb_cells += directional_nb_cells[dim];

    Real additional_grid_length = (directional_nb_cells[dim]*grid_spacing[dim] - grid_length) / 2.;

    /// get minimal and maximal coordinates of the grid
    grid_min[dim] = min[master_surface][dim] - additional_grid_length;
    grid_max[dim] = max[master_surface][dim] + additional_grid_length;
  }

  /// find surfaces being in the grid space
  UInt nb_grid_surfaces = 0;
  UInt grid_surfaces[nb_surfaces];
  UInt not_grid_space = 0;

  for(UInt surf = 0; surf < nb_surfaces; ++surf) {
    for(UInt dim = 0; dim < spatial_dimension; ++dim) {
      if(max[surf][dim] < grid_min[dim] || min[surf][dim] > grid_max[dim])
  	not_grid_space = 1;
    }
    if(not_grid_space == 0 || surf == master_surface) {
      grid_surfaces[nb_grid_surfaces++] = surf;
    }
    not_grid_space = 0;
  }

  /// if number of grid surfaces is equal to 1 we do not need to consider any slave surface
  /// @todo exit with empty neighbor list


  /// assign cell number to all surface nodes (put -1 if out of grid space) (cell numbers start with zero)
  Int not_grid_space_node    = -1;  // should not be same as not_grid_space_surface and be < 0
  Int not_grid_space_surface = -2;  // should not be same as not_grid_space_node and be < 0

  Int directional_cell[spatial_dimension];
  UInt nb_surface_nodes = surface_to_nodes_offset[nb_surfaces];
  Vector<Int> * cell = new Vector<Int>(nb_surface_nodes, 1, not_grid_space_surface);
  Int * cell_val = cell->values;

  /// define the cell number for all surface nodes
  for(UInt surf = 0; surf < nb_grid_surfaces; ++surf) {
    UInt current_surface = grid_surfaces[surf];
    UInt min_surf_offset = surface_to_nodes_offset[current_surface];
    UInt max_surf_offset = surface_to_nodes_offset[current_surface + 1];

    for(UInt n = min_surf_offset; n < max_surf_offset; ++n) {
      UInt cur_node = surface_to_nodes[n];

      /// compute cell index for all directions
      for(UInt dim = 0; dim < spatial_dimension; ++dim) {
	Real cur_position = node_position[cur_node + dim];
	directional_cell[dim] = static_cast<UInt>(floor((cur_position - grid_min[dim])/grid_spacing[dim]));
      }
      
      /// test if out of grid space
      for(UInt dim = 0; dim < spatial_dimension; ++dim) {
	if(directional_cell[dim] < 0 || directional_cell[dim] > directional_nb_cells[dim])
	  cell_val[n] = not_grid_space_node;
      }
    
      /// compute global cell index
      if(cell_val[n] != not_grid_space_node) {
	cell_val[n] = computeCellNb(directional_nb_cells, directional_cell);
	// cell_val[sn] = directional_cell[spatial_dimension - 1];
	// for(Int dim = spatial_dimension - 2; dim > 0; --dim) {
	//   cell_val[sn] *= directional_nb_cells[dim];
	//   cell_val[sn] += directional_cell[dim];
	// }
      }
    }
  } 


  /// define offset arrays for nodes per cell which will be computed below
  UInt * impactor_nodes_cell_offset = new UInt[nb_cells + 1];
  memset(impactor_nodes_cell_offset, 0, nb_cells * sizeof(UInt));
  UInt * master_nodes_cell_offset = new UInt[nb_cells + 1];
  memset(master_nodes_cell_offset, 0, nb_cells * sizeof(UInt));
  
  /// count number of nodes per cell for impactors and master
  for(UInt surf = 0; surf < nb_grid_surfaces; ++surf) {
    UInt current_surface = grid_surfaces[surf];
    UInt min_surf_offset = surface_to_nodes_offset[current_surface];
    UInt max_surf_offset = surface_to_nodes_offset[current_surface + 1];

    /// define temporary pointers
    UInt * nodes_cell_offset;
    if(current_surface == master_surface) {
      nodes_cell_offset = master_nodes_cell_offset;
    } else {
      nodes_cell_offset = impactor_nodes_cell_offset;
    }

    for(UInt n = min_surf_offset; n < max_surf_offset; ++n) {
      UInt node = surface_to_nodes[n];
      UInt cell = cell_val[node];
      if(cell != not_grid_space_node)
	nodes_cell_offset[cell]++;
    }
  }

  /// create two separate offset arrays for impactor nodes and master nodes
  for (UInt i = 1; i < nb_cells; ++i) {
    impactor_nodes_cell_offset[i] += impactor_nodes_cell_offset[i - 1];
    master_nodes_cell_offset  [i] += master_nodes_cell_offset  [i - 1];
  }
  for (UInt i = nb_cells; i > 0; --i) {
    impactor_nodes_cell_offset[i] = impactor_nodes_cell_offset[i - 1];
    master_nodes_cell_offset  [i] = master_nodes_cell_offset  [i - 1];
  }
  impactor_nodes_cell_offset[0] = 0;
  master_nodes_cell_offset  [0] = 0;


  /// find all impactor and master nodes in a cell
  UInt * impactor_nodes_cell = new UInt[impactor_nodes_cell_offset[nb_cells]];
  UInt * master_nodes_cell   = new UInt[master_nodes_cell_offset  [nb_cells]];
  cell_val = cell->values;

  for(UInt surf = 0; surf < nb_grid_surfaces; ++surf) {
    UInt current_surface = grid_surfaces[surf];
    UInt min_surf_offset = surface_to_nodes_offset[current_surface];
    UInt max_surf_offset = surface_to_nodes_offset[current_surface + 1];

    /// define temporary variables
    UInt * nodes_cell;
    UInt * nodes_cell_offset;
    if(current_surface == master_surface) {
      nodes_cell        = master_nodes_cell;
      nodes_cell_offset = master_nodes_cell_offset;
    } else {
      nodes_cell        = impactor_nodes_cell;
      nodes_cell_offset = impactor_nodes_cell_offset;
    }

    /// loop over the nodes of surf and create nodes_cell
    for(UInt n = min_surf_offset; n < max_surf_offset; ++n) {
      UInt node = surface_to_nodes[n];
      UInt cell = cell_val[node];
      if(cell != not_grid_space_node)
	nodes_cell[nodes_cell_offset[cell]++] = node;
    }
  }

  for (UInt i = nb_cells; i > 0; --i) {
    impactor_nodes_cell_offset[i] = impactor_nodes_cell_offset[i - 1];
    master_nodes_cell_offset  [i] =   master_nodes_cell_offset[i - 1];
  }
  impactor_nodes_cell_offset[0] = 0;
  master_nodes_cell_offset  [0] = 0;

  

  cell_val = cell->values;
  UInt nb_impactor_nodes = impactor_nodes_cell_offset[nb_cells];
  
  /// declare neighbor list
  //  NeighborList neighbor_list;
  //neighbor_list.nb_nodes = 0;
  //neighbor_list.impactor_nodes = new Vector<UInt>(nb_impactor_nodes, 1);
  neighbor_list.impactor_nodes.resize(nb_impactor_nodes);
  UInt * impactor_nodes_val = neighbor_list.impactor_nodes.values;
  // for (UInt i = 0; i < _max_element_type; ++i) {
  //   neighbor_list.facets_offset[i] = NULL;
  //   neighbor_list.facets       [i] = NULL;
  // }
  
  /// define maximal number of neighbor cells and include it-self
  UInt max_nb_neighbor_cells;
  if(spatial_dimension == 2) {
    max_nb_neighbor_cells = 9;
  }
  else if(spatial_dimension == 3) {
    max_nb_neighbor_cells = 27;
  }

  //const ByElementTypeUInt & node_to_elements_offset = contact_search.getContact().getNodeToElementsOffset();
  //ByElementTypeUInt node_to_elements = contact_search.getContact().getNodeToElements();

  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  /// find existing surface element types
  UInt nb_types = type_list.size();
  UInt nb_facet_types = 0;
  ElementType facet_type[_max_element_type];
 
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    const Vector<UInt> & nte = contact_search.getContact().getNodeToElements(type);
    if(type == (spatial_dimension-1) && nte.getSize() != 0) {
      facet_type[nb_facet_types++] = mesh.getFacetElementType(type);
    }
  }

  /// loop over all existing surface element types
  for (UInt el_type = 0; el_type < nb_facet_types; ++el_type) {
    ElementType type = facet_type[el_type];

    const Vector<UInt> & node_to_elements_offset = contact_search.getContact().getNodeToElementsOffset(type);
    const Vector<UInt> & node_to_elements = contact_search.getContact().getNodeToElements(type);

    UInt * node_to_elements_offset_val = node_to_elements_offset.values;
    UInt * node_to_elements_val        = node_to_elements.values;;

    Vector<bool> * visited_node = new Vector<bool>(nb_impactor_nodes, 1, false); // does it need a delete at the end ?
    bool * visited_node_val = visited_node->values;
    UInt neighbor_cells[max_nb_neighbor_cells];
    
    std::stringstream sstr_name_offset; 
    sstr_name_offset << id << ":facets_offset:" << type;
    neighbor_list.facets_offset[type] = new Vector<UInt>(0, 1, sstr_name_offset.str());
    Vector<UInt> & tmp_facets_offset = *(neighbor_list.facets_offset[type]);
    UInt * tmp_facets_offset_val = tmp_facets_offset.values;
    std::stringstream sstr_name;
    sstr_name << id << ":facets:" << type;
    neighbor_list.facets[type] = new Vector<UInt>(0, 1, sstr_name.str());// declare vector that is ??
    Vector<UInt> & tmp_facets = *(neighbor_list.facets[type]);

    for(UInt in = 0; in < nb_impactor_nodes; ++in) {
      UInt current_impactor_node = impactor_nodes_cell[in];
      /// test if nodes has not already been visited
      if(!visited_node_val[in]) {
	
	/// find and store cell numbers of neighbor cells and it-self
	AKANTU_DEBUG_ASSERT(cell_val[current_impactor_node] >= 0, 
			    "Bad cell index. This case normally should not happen !!");
	UInt current_cell = cell_val[current_impactor_node];
	UInt nb_neighbor_cells = computeNeighborCells(current_cell, neighbor_cells, directional_nb_cells);
	neighbor_cells[nb_neighbor_cells++] = current_cell;
	
	/// define a set in which the found master surface elements are stored
	std::set<UInt> master_surface_elements;
	std::set<UInt>::iterator it_set;
	
	/// find all master elements that are in the considered region
	for(UInt cl = 0; cl < nb_neighbor_cells; ++cl) {
	  
	  /// get cell number and offset range
	  UInt considered_cell = neighbor_cells[cl];
	  UInt min_master_offset = master_nodes_cell_offset[considered_cell];
	  UInt max_master_offset = master_nodes_cell_offset[considered_cell + 1];
	  
	  for(UInt mn = min_master_offset; mn < max_master_offset; ++mn) {
	    UInt master_node = master_nodes_cell[mn];
	    UInt min_element_offset = node_to_elements_offset_val[master_node];
	    UInt max_element_offset = node_to_elements_offset_val[master_node + 1];
	    
	    for(UInt el = min_element_offset; el < max_element_offset; ++el)
	      master_surface_elements.insert(node_to_elements_val[el]);
	  }
	}
	UInt min_impactor_offset = impactor_nodes_cell_offset[current_cell];
	UInt max_impactor_offset = impactor_nodes_cell_offset[current_cell + 1];
	
	for(UInt imp = min_impactor_offset; imp < max_impactor_offset; ++imp) {
	  UInt impactor_node = impactor_nodes_cell[imp];
	  impactor_nodes_val[neighbor_list.nb_nodes++] = impactor_node;
	  for(it_set = master_surface_elements.begin(); it_set != master_surface_elements.end(); it_set++) {
	    tmp_facets.push_back(*it_set);
	  }
	  tmp_facets_offset.push_back(master_surface_elements.size());
	  for(UInt inn = 0; inn < nb_impactor_nodes; ++inn) {
	    if(impactor_nodes_cell[inn] == impactor_node) {
	      visited_node_val[inn] = true;
	      break;
	    }
	  }
	}
      }
    }
    
    for (UInt i = 1; i < neighbor_list.nb_nodes; ++i) tmp_facets_offset_val[i] += tmp_facets_offset_val[i - 1];
    for (UInt i = neighbor_list.nb_nodes; i > 0; --i) tmp_facets_offset_val[i]  = tmp_facets_offset_val[i - 1];
    tmp_facets_offset_val[0] = 0;
    
    delete visited_node;
    
  }
  
  delete cell;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension> 
bool RegularGridNeighborStructure<spatial_dimension>::check() {
  
  /// test if increment flag is on
  /// @todo compute test if increment flag is on and give warning if not

  AKANTU_DEBUG_IN();

  bool need_update = false;
  UInt nb_surfaces = mesh.getNbSurfaces();
  Real * current_increment = contact_search.getContact().getModel().getIncrement().values;

  Real max[spatial_dimension];
  
  /// initialize max table with extrem values
  for(UInt dim = 0; dim < spatial_dimension; ++dim)
    max[dim] = std::numeric_limits<Real>::min();

  // get the nodes that are on the surfaces
  UInt * surface_to_nodes_offset = contact_search.getContact().getSurfaceToNodesOffset().values;
  UInt * surface_to_nodes        = contact_search.getContact().getSurfaceToNodes().values;

  /// find maximal increment of surface nodes in all directions
  for(UInt surf = 0; surf < nb_surfaces; ++surf) {
    UInt min_surf_offset = surface_to_nodes_offset[surf];
    UInt max_surf_offset = surface_to_nodes_offset[surf + 1];
    for(UInt n = min_surf_offset; n < max_surf_offset; ++n) {
      UInt cur_node = surface_to_nodes[n];
      for(UInt dim = 0; dim < spatial_dimension; ++dim) {
  	Real cur_increment = current_increment[cur_node + dim];
  	max[dim] = std::max(max[dim], fabs(cur_increment));
      }
    } 
  }
  
  /// test if the total increment is larger than a critical distance
  for(UInt dim = 0; dim < spatial_dimension; ++dim) {
    max_increment[dim] += max[dim];
    if(max_increment[dim] > security_factor[dim] * grid_spacing[dim]) {
      need_update = true;
      break;
    }
  }

  AKANTU_DEBUG_OUT();
  return need_update;
}

template class RegularGridNeighborStructure<2>;
template class RegularGridNeighborStructure<3>;

__END_AKANTU__
