/**
 * @file   boundary.cc
 *
 * @author Dana Christen <dana.christen@gmail.com>
 *
 * @date   Wed Mar 06 09:30:00 2013
 *
 * @brief  Stores information relevent to the notion of domain boundary and surfaces.
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
#include <sstream>
#include <algorithm>
#include <iterator>
#include "boundary.hh"
#include "mesh.hh"
#include "aka_csr.hh"
#include "mesh_utils.hh"
#include "sub_boundary.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
Boundary::Boundary(const Mesh & mesh, const ID & id, const ID & parent_id, const MemoryID & mem_id)
:memory_id(mem_id), mesh(mesh)
{
  AKANTU_DEBUG_IN();

  std::stringstream sstr;
  if(parent_id != "") sstr << parent_id << ":";
  sstr << id;

  this->id = sstr.str();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Boundary::~Boundary() {
  for(BoundaryList::iterator boundaries_iter = boundaries.begin(); boundaries_iter != boundaries.end(); boundaries_iter++)
  {
    delete boundaries_iter->second;
  }
}

/* -------------------------------------------------------------------------- */
Boundary::BoundaryTypeSet Boundary::getBoundaryElementTypes() {

  // find surface element type that exists in this mesh (from David's code)
  UInt dim = mesh.getSpatialDimension();
  BoundaryTypeSet surface_types;
  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  for (Mesh::ConnectivityTypeList::const_iterator it = type_list.begin(); it != type_list.end(); ++it)
  {
    ElementType surface_type = mesh.getFacetType(*it);
    if (mesh.getSpatialDimension(*it) == dim &&	mesh.getNbElement(surface_type) != 0)
    {
      surface_types.insert(surface_type);
    }
  }
  return surface_types;
}

/* -------------------------------------------------------------------------- */
void Boundary::addElementAndNodesToBoundaryAlloc(const std::string & boundary_name,
						 const ElementType & elem_type,
						 UInt elem_id,
						 const GhostType & ghost_type) {
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(elem_type);
  BoundaryList::iterator boundaries_iter = boundaries.find(boundary_name);

  const Array<UInt> & connectivity = mesh.getConnectivity(elem_type, ghost_type);
  if(boundaries_iter == boundaries.end()) {
    SubBoundary * sub_b = new SubBoundary(boundary_name, std::string(id + ":sub_boundary:" + boundary_name), memory_id);
    boundaries_iter = boundaries.insert(boundaries_iter, std::pair<std::string, SubBoundary*>(boundary_name, sub_b));
  }

  boundaries_iter->second->addElement(elem_type, elem_id, ghost_type);
  for (UInt n(0); n < nb_nodes_per_element; n++) {
    boundaries_iter->second->addNode(connectivity(elem_id, n));
  }
}

/* -------------------------------------------------------------------------- */
void Boundary::createBoundariesFromGeometry() {
  AKANTU_DEBUG_IN();

  CSR<UInt> node_to_elem;

  /// Get list of surface elements
  UInt spatial_dimension = mesh.getSpatialDimension();

  //  buildNode2Elements(mesh, node_offset, node_to_elem, spatial_dimension-1);
  MeshUtils::buildNode2Elements(mesh, node_to_elem, spatial_dimension-1);

  /// Find which types of elements have been linearized
  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  UInt nb_types = type_list.size();
  ElementType lin_element_type[nb_types];
  UInt nb_lin_types = 0;

  UInt nb_nodes_per_element[nb_types];
  UInt nb_nodes_per_element_p1[nb_types];

  UInt * conn_val[nb_types];
  UInt nb_element[nb_types+1];

  ElementType type_p1;

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    //std::cout <<"Main type: " << *it << std::endl;

    ElementType type = *it;
    if(Mesh::getSpatialDimension(type) != spatial_dimension) {
      continue;
    }

    ElementType facet_type = mesh.getFacetType(type);
    //std::cout <<"Facet type: " << facet_type << std::endl;
    lin_element_type[nb_lin_types] = facet_type;
    nb_nodes_per_element[nb_lin_types]    = Mesh::getNbNodesPerElement(facet_type);
    type_p1 = Mesh::getP1ElementType(facet_type);
    nb_nodes_per_element_p1[nb_lin_types] = Mesh::getNbNodesPerElement(type_p1);

    conn_val[nb_lin_types] = mesh.getConnectivity(facet_type, _not_ghost).values;
    nb_element[nb_lin_types] = mesh.getNbElement(facet_type, _not_ghost);
    nb_lin_types++;
  }

  for (UInt i = 1; i < nb_lin_types; ++i) {
    nb_element[i] += nb_element[i+1];
  }

  for (UInt i = nb_lin_types; i > 0; --i) {
    nb_element[i] = nb_element[i-1];
  }

  nb_element[0] = 0;

  /// Find close surfaces
  Array<Int> surface_value_id(1, nb_element[nb_lin_types], -1);
  Int * surf_val = surface_value_id.storage();
  UInt nb_boundaries = 0;

  UInt nb_cecked_elements;
  UInt nb_elements_to_ceck;
  UInt * elements_to_ceck = new UInt [nb_element[nb_lin_types]];
  memset(elements_to_ceck, 0, nb_element[nb_lin_types]*sizeof(UInt));

  for (UInt lin_el = 0; lin_el < nb_element[nb_lin_types]; ++lin_el) {

    if(surf_val[lin_el] != -1) continue; /* Surface id already assigned */

    /* First element of new surface */
    surf_val[lin_el] = nb_boundaries;
    nb_cecked_elements = 0;
    nb_elements_to_ceck = 1;
    memset(elements_to_ceck, 0, nb_element[nb_lin_types]*sizeof(UInt));
    elements_to_ceck[0] = lin_el;

    // Find others elements belonging to this surface
    while(nb_cecked_elements < nb_elements_to_ceck) {

      UInt ceck_lin_el = elements_to_ceck[nb_cecked_elements];

      // Transform linearized index of element into ElementType one
      UInt lin_type_nb = 0;
      while (ceck_lin_el >= nb_element[lin_type_nb+1]) {
	      lin_type_nb++;
	    }
      UInt ceck_el = ceck_lin_el - nb_element[lin_type_nb];

      // Get connected elements
      UInt el_offset = ceck_el*nb_nodes_per_element[lin_type_nb];
      for (UInt n = 0; n < nb_nodes_per_element_p1[lin_type_nb]; ++n) {
	      UInt node_id = conn_val[lin_type_nb][el_offset + n];
	      CSR<UInt>::iterator it_n;
	      for (it_n = node_to_elem.begin(node_id); it_n != node_to_elem.end(node_id); ++it_n) {
	        if(surf_val[*it_n] == -1) { /* Found new surface element */
	          surf_val[*it_n] = nb_boundaries;
	          elements_to_ceck[nb_elements_to_ceck] = *it_n;
	          nb_elements_to_ceck++;
          }
        }
      }
      nb_cecked_elements++;
    }
    nb_boundaries++;
  }

  delete [] elements_to_ceck;

  /// Transform local linearized element index in the global one
  for (UInt i = 0; i < nb_lin_types; ++i) {
    nb_element[i] = nb_element[i+1] - nb_element[i];
  }

  UInt el_offset = 0;

  for(UInt type_it = 0; type_it < nb_lin_types; ++type_it) {
    ElementType type = lin_element_type[type_it];
    //std::cout << "Second type: " << type << std::endl;
    for (UInt el = 0; el < nb_element[type_it]; ++el) {
      std::stringstream sstr;
      sstr << surf_val[el+el_offset];
      std::string surf_name(sstr.str());
      //std::cout << "Inserting elements and nodes to surface " << surf_name << std::endl;
      addElementAndNodesToBoundaryAlloc(surf_name, type, el);
      }
    //std::cout << "Offset +" << nb_element[type_it] << std::endl;
    el_offset += nb_element[type_it];
  }

  BoundaryList::iterator subB_iter = boundaries.begin();
  BoundaryList::iterator subB_iter_end = boundaries.end();

  for(;subB_iter != subB_iter_end; ++subB_iter) {
    SubBoundary & sub = *subB_iter->second;
    sub.cleanUpNodeList();
    sub.addDumpFilteredMesh(mesh,
			    sub.elements,
			    sub.nodes,
			    mesh.getSpatialDimension() - 1,
			    _not_ghost,
			    _ek_regular);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<typename T>
void Boundary::createBoundariesFromMeshData(const std::string & dataset_name)
{
  UInt spatial_dimension = mesh.getSpatialDimension();
  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();  gt != ghost_type_t::end(); ++gt) {
    Mesh::type_iterator type_it = mesh.firstType(spatial_dimension - 1, *gt);
    Mesh::type_iterator type_end  = mesh.lastType(spatial_dimension - 1, *gt);
    for (; type_it != type_end; ++type_it) {
      // Dirty check to assert that the processor has indeed this type of element
      // associated to the dataset name
      try {
        mesh.getData<T>(*type_it, dataset_name, *gt);
      } catch(akantu::debug::Exception e) {
        if(psize == 1) {
          AKANTU_DEBUG_ERROR("Error building the boundaries. Type " << *type_it << "/" << *gt << "not registered in dataset  \"" << dataset_name << "\".");
        }
        else {
          AKANTU_DEBUG_INFO("Rank " << comm.whoAmI(); << " does not have any elements of type " << *type_it << "/" << *gt << "in dataset \"" << dataset_name << "\".");
          continue;
        }
      }

      const Array<T> & dataset = mesh.getData<T>(*type_it, dataset_name, *gt);
      UInt nb_element = mesh.getNbElement(*type_it, *gt);

      AKANTU_DEBUG_ASSERT(dataset.getSize() == nb_element,
			  "Not the same number of elements in the map from MeshData and in the mesh!");

      // FIXME This could be improved (performance...)
      for(UInt i(0); i < nb_element; ++i) {
        std::stringstream sstr;
        sstr << dataset(i);
        std::string boundary_name = sstr.str();
        addElementAndNodesToBoundaryAlloc(boundary_name, *type_it, i, *gt);
      }
    }
  }

  BoundaryList::iterator subB_iter = boundaries.begin();
  BoundaryList::iterator subB_iter_end = boundaries.end();

  for(;subB_iter != subB_iter_end; ++subB_iter) {
    SubBoundary & sub = *subB_iter->second;
    sub.cleanUpNodeList();
    sub.addDumpFilteredMesh(mesh,
			    sub.elements,
			    sub.nodes,
			    mesh.getSpatialDimension() - 1,
			    _not_ghost,
			    _ek_regular);
  }
}

template void Boundary::createBoundariesFromMeshData<std::string>(const std::string & dataset_name);
template void Boundary::createBoundariesFromMeshData<UInt>(const std::string & dataset_name);

/* -------------------------------------------------------------------------- */
void Boundary::createSubBoundaryFromNodeGroup(const std::string & name,
					      const Array<UInt> & node_group) {
  CSR<Element> node_to_elem;
  MeshUtils::buildNode2Elements(mesh, node_to_elem, mesh.getSpatialDimension() - 1);

  std::set<Element> seen;

  BoundaryList::iterator boundaries_iter = boundaries.find(name);
  if(boundaries_iter == boundaries.end()) {
    // Manual construction of the sub-boundary ID
    SubBoundary * sub_b = new SubBoundary(name, std::string(id + "_sub_boundary_" + name), memory_id);
    boundaries_iter = boundaries.insert(boundaries_iter, std::pair<std::string, SubBoundary*>(name, sub_b));
  }

  SubBoundary & sub_bound = *boundaries_iter->second;

  Array<UInt>::const_iterator<> itn  = node_group.begin();
  Array<UInt>::const_iterator<> endn = node_group.end();
  for (;itn != endn; ++itn) {
    CSR<Element>::iterator ite = node_to_elem.begin(*itn);
    CSR<Element>::iterator ende = node_to_elem.end(*itn);
    for (;ite != ende; ++ite) {
      const Element & elem = *ite;
      if(seen.find(elem) != seen.end()) continue;

      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(elem.type);
      Array<UInt>::const_iterator< Vector<UInt> > conn_it =
      mesh.getConnectivity(elem.type, elem.ghost_type).begin(nb_nodes_per_element);

      const Vector<UInt> & conn = conn_it[elem.element];
      UInt count = 0;
      for (UInt n = 0; n < conn.size(); ++n) {
        count += (node_group.find(conn(n)) != -1 ? 1 : 0);
      }

      if(count == nb_nodes_per_element) {
        sub_bound.addElement(elem.type, elem.element, elem.ghost_type);
        for (UInt n(0); n < nb_nodes_per_element; n++) {
          sub_bound.addNode(conn(n));
        }
      }

      seen.insert(elem);
    }
  }

  sub_bound.cleanUpNodeList();
  sub_bound.addDumpFilteredMesh(mesh,
				sub_bound.elements,
				sub_bound.nodes,
				mesh.getSpatialDimension() - 1,
				_not_ghost,
				_ek_regular);
}


/* -------------------------------------------------------------------------- */
void Boundary::createBoundariesFromMeshData(const std::string & dataset_name) {
  createBoundariesFromMeshData<std::string>(dataset_name);
}

/* -------------------------------------------------------------------------- */
void Boundary::printself(std::ostream & stream) const {
  for(const_iterator it(begin()); it != end(); ++it) {
    stream << "-- Boundary \"" << it->getName() << "\":" << std::endl;

    // Loop over the nodes
    stream << "---- " << it->getNbNodes() << " nodes";
    if(AKANTU_DEBUG_TEST(dblDump)) {
        stream << ":" << std::endl;
      stream << "------ ";
      for(SubBoundary::nodes_const_iterator nodes_it(it->nodes_begin()); nodes_it!= it->nodes_end(); ++nodes_it) {
        stream << *nodes_it << " ";
      }
    }
    stream << std::endl;

    // Loop over the element types
    stream << "---- Elements:" << std::endl;
    Mesh::type_iterator type_it = mesh.firstType(mesh.getSpatialDimension()-1);
    Mesh::type_iterator type_end = mesh.lastType(mesh.getSpatialDimension()-1);
    for(; type_it != type_end; ++type_it) {
      stream << "------ Type " << *type_it;
      const Array<UInt> & element_ids = it->getElements(*type_it);
      if(AKANTU_DEBUG_TEST(dblDump)) {
        stream << ":" << std::endl;
        stream << "-------- ";
        // Loop over the corresponding elements in the SubBoundary
        for(UInt i(0); i<element_ids.getSize(); ++i) {
          stream << element_ids(i) << " ";
        }
      }
      stream << std::endl;
      //std::cout << it->getElements(*type_it);
    }
  }
  stream << std::endl;
}

/* -------------------------------------------------------------------------- */
void Boundary::dump() {
  for(iterator it(begin()); it != end(); ++it) {
    it->dump();
  }
}

__END_AKANTU__

