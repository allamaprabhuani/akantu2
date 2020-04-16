/**
 * @file   surface_selector.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Fri Jun 21 2019
 * @date last modification: Fri Jun 21 2019
 *
 * @brief  Surface selector for contact detector
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "surface_selector.hh"
#include "model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SurfaceSelector::SurfaceSelector(Mesh & mesh)
    : Parsable(ParserType::_contact_detector), mesh(mesh) {}


/* -------------------------------------------------------------------------- */
/**
 * class that selects contact surface from physical names 
 */
PhysicalSurfaceSelector::PhysicalSurfaceSelector(Mesh & mesh)
    : SurfaceSelector(mesh) {

  const Parser & parser = getStaticParser();

  const ParserSection & section =
      *(parser.getSubSections(ParserType::_contact_detector).first);
  
  master = section.getParameterValue<std::string>("master");
  slave = section.getParameterValue<std::string>("slave");

  UInt surface_dimension = mesh.getSpatialDimension() - 1;
  auto & group = mesh.createElementGroup("contact_surface",
					 surface_dimension);
  group.append(mesh.getElementGroup(master));
  group.append(mesh.getElementGroup(slave));

  group.optimize();
}

/* -------------------------------------------------------------------------- */
Array<UInt> & PhysicalSurfaceSelector::getMasterList() {
  return mesh.getElementGroup(master).getNodeGroup().getNodes();
}

/* -------------------------------------------------------------------------- */
Array<UInt> & PhysicalSurfaceSelector::getSlaveList() {
  return mesh.getElementGroup(slave).getNodeGroup().getNodes();
}


/* -------------------------------------------------------------------------- */
/**
 * class that selects contact surface from cohesive elements
 */
#if defined(AKANTU_COHESIVE_ELEMENT)
/* -------------------------------------------------------------------------- */
CohesiveSurfaceSelector::CohesiveSurfaceSelector(Mesh & mesh)
    : SurfaceSelector(mesh), mesh_facets(mesh.getMeshFacets()) {
  this->mesh.registerEventHandler(*this, _ehp_lowest);
}

/* -------------------------------------------------------------------------- */
void CohesiveSurfaceSelector::onNodesAdded(const Array<UInt> & new_nodes,
                                           const NewNodesEvent & event) {

  if (not aka::is_of_type<CohesiveNewNodesEvent>(event))
    return;

  const auto & cohesive_event = aka::as_type<CohesiveNewNodesEvent>(event);
  const auto & old_nodes = cohesive_event.getOldNodesList();

  UInt nb_new_nodes = new_nodes.size();
  UInt nb_old_nodes = old_nodes.size();

  for (auto n : arange(nb_new_nodes)) {
    new_nodes_list.push_back(new_nodes(n));
  }

  for (auto n : arange(nb_old_nodes)) {
    new_nodes_list.push_back(old_nodes(n));
  }

  mesh_facets.fillNodesToElements(mesh.getSpatialDimension() - 1);

  UInt surface_dimension = mesh.getSpatialDimension() - 1;
  auto & group = mesh_facets.createElementGroup("contact_surface",
						surface_dimension, true);

  for (auto node : new_nodes_list) {
    Array<Element> all_elements;
    mesh_facets.getAssociatedElements(node, all_elements);

    Array<Element> mesh_facet_elements;
    this->filterBoundaryElements(all_elements, mesh_facet_elements);

    for (auto nb_elem : arange(mesh_facet_elements.size()))
      group.add(mesh_facet_elements[nb_elem], true);
  }
  group.optimize();
}

/* -------------------------------------------------------------------------- */
void CohesiveSurfaceSelector::filterBoundaryElements(
    Array<Element> & elements, Array<Element> & boundary_elements) {

  for (auto elem : elements) {
    
    const auto & element_to_subelement =
        mesh_facets.getElementToSubelement(elem.type)(elem.element);

    UInt nb_subelements_regular = 0;
    for (auto subelem : element_to_subelement) {
      if (subelem == ElementNull)
        continue;

      if (subelem.kind() == _ek_regular) {
        ++nb_subelements_regular;
      }
    }

    auto nb_subelements = element_to_subelement.size();

    if (nb_subelements_regular < nb_subelements and nb_subelements != 1) {
      boundary_elements.push_back(elem);
    }
  }
}

/* -------------------------------------------------------------------------- */
Array<UInt> & CohesiveSurfaceSelector::getMasterList() {
  return this->getNewNodesList();
}

/* -------------------------------------------------------------------------- */
Array<UInt> & CohesiveSurfaceSelector::getSlaveList() {
  return this->getNewNodesList();
}

/* -------------------------------------------------------------------------- */
/**
 * class that selects contact surface from both cohesive elements and
 * physical names
 */
AllSurfaceSelector::AllSurfaceSelector(Mesh & mesh)
    : SurfaceSelector(mesh), mesh_facets(mesh.getMeshFacets()) {
  this->mesh.registerEventHandler(*this, _ehp_lowest);

  const Parser & parser = getStaticParser();

  const ParserSection & section =
      *(parser.getSubSections(ParserType::_contact_detector).first);

  master = section.getParameterValue<std::string>("master");
  slave = section.getParameterValue<std::string>("slave");

  UInt surface_dimension = this->mesh.getSpatialDimension() - 1;
  auto & group = mesh_facets.createElementGroup("contact_surface",
						surface_dimension);
  group.append(mesh_facets.getElementGroup(master));
  group.append(mesh_facets.getElementGroup(slave));

  group.optimize();
}

/* -------------------------------------------------------------------------- */
void AllSurfaceSelector::onNodesAdded(const Array<UInt> & new_nodes,
                                      const NewNodesEvent & event) {

  if (not aka::is_of_type<CohesiveNewNodesEvent>(event))
    return;

  const auto & cohesive_event = aka::as_type<CohesiveNewNodesEvent>(event);
  const auto & old_nodes = cohesive_event.getOldNodesList();

  UInt nb_new_nodes = new_nodes.size();
  UInt nb_old_nodes = old_nodes.size();

  auto & slave_group = mesh_facets.getElementGroup(slave).getNodeGroup();
  auto & master_group = mesh_facets.getElementGroup(master).getNodeGroup();

  for (auto n : arange(nb_new_nodes)) {
    new_nodes_list.push_back(new_nodes(n));
    slave_group.add(new_nodes(n));
    master_group.add(new_nodes(n));
  }

  for (auto n : arange(nb_old_nodes)) {
    new_nodes_list.push_back(old_nodes(n));
    slave_group.add(old_nodes(n));
    master_group.add(old_nodes(n));
  }

  mesh_facets.fillNodesToElements(mesh.getSpatialDimension() - 1);

  auto & group = mesh_facets.getElementGroup("contact_surface");

  for (auto node : new_nodes_list) {
    Array<Element> all_elements;
    mesh_facets.getAssociatedElements(node, all_elements);

    Array<Element> mesh_facet_elements;
    this->filterBoundaryElements(all_elements, mesh_facet_elements);

    for (auto nb_elem : arange(mesh_facet_elements.size()))
      group.add(mesh_facet_elements[nb_elem], true);
  }

  group.optimize();
}

/* -------------------------------------------------------------------------- */
void AllSurfaceSelector::filterBoundaryElements(
    Array<Element> & elements, Array<Element> & boundary_elements) {

  for (auto elem : elements) {

    // to ensure that normal is always outwards from master surface
    const auto & element_to_subelement =
        mesh_facets.getElementToSubelement(elem.type)(elem.element);

    UInt nb_subelements_regular = 0;
    for (auto subelem : element_to_subelement) {
      if (subelem == ElementNull)
        continue;

      if (subelem.kind() == _ek_regular) {
        ++nb_subelements_regular;
      }
    }

    auto nb_subelements = element_to_subelement.size();

    if (nb_subelements_regular < nb_subelements) {
      boundary_elements.push_back(elem);
    }
  }
}

/* -------------------------------------------------------------------------- */
Array<UInt> & AllSurfaceSelector::getMasterList() {
  return this->getNewNodesList();
}

/* -------------------------------------------------------------------------- */
Array<UInt> & AllSurfaceSelector::getSlaveList() {
  return this->getNewNodesList();
}

#endif

} // namespace akantu
