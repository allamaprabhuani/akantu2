/**
 * @file   surface_selector.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Sun Jun 30 2019
 * @date last modification: Fri Sep 18 2020
 *
 * @brief  Surface selector for contact detector
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "surface_selector.hh"
#include "geometry_utils.hh"
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
  auto & group = mesh.createElementGroup("contact_surface", surface_dimension);
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

  UInt surface_dimension = mesh.getSpatialDimension() - 1;
  mesh_facets.createElementGroup("contact_surface", surface_dimension, true);
}

/* -------------------------------------------------------------------------- */
void CohesiveSurfaceSelector::onElementsAdded(
    const Array<Element> & element_list,
    __attribute__((unused)) const NewElementsEvent & event) {

  auto & group = mesh_facets.getElementGroup("contact_surface");

  for (auto elem : element_list) {
    if (elem.kind() != _ek_cohesive) {
      continue;
    }

    const auto & subelement_to_element =
        mesh_facets.getSubelementToElement(elem.type);

    auto && facets = Vector<Element>(
        make_view(subelement_to_element, subelement_to_element.getNbComponent())
            .begin()[elem.element]);

    for (auto facet : facets) {
      group.add(facet, true);
    }
  }

  group.optimize();
}

/* -------------------------------------------------------------------------- */
void CohesiveSurfaceSelector::onNodesAdded(const Array<UInt> & /*nodes_list*/,
                                           const NewNodesEvent & event) {

  if (not aka::is_of_type<CohesiveNewNodesEvent>(event)) {
    return;
  }

  mesh_facets.fillNodesToElements(mesh.getSpatialDimension() - 1);
}

/* -------------------------------------------------------------------------- */
Array<UInt> & CohesiveSurfaceSelector::getMasterList() {
  return mesh_facets.getElementGroup("contact_surface")
      .getNodeGroup()
      .getNodes();
}

/* -------------------------------------------------------------------------- */
Array<UInt> & CohesiveSurfaceSelector::getSlaveList() {
  return mesh_facets.getElementGroup("contact_surface")
      .getNodeGroup()
      .getNodes();
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
  auto & group =
      mesh_facets.createElementGroup("contact_surface", surface_dimension);
  group.append(mesh_facets.getElementGroup(master));
  group.append(mesh_facets.getElementGroup(slave));

  group.optimize();
}

/* -------------------------------------------------------------------------- */
void AllSurfaceSelector::onElementsAdded(const Array<Element> & element_list,
                                         __attribute__((unused))
                                         const NewElementsEvent & event) {

  auto & group = mesh_facets.getElementGroup("contact_surface");

  for (auto elem : element_list) {
    if (elem.kind() != _ek_cohesive) {
      continue;
    }

    const auto & subelement_to_element =
        mesh_facets.getSubelementToElement(elem.type);

    auto && facets = Vector<Element>(
        make_view(subelement_to_element, subelement_to_element.getNbComponent())
            .begin()[elem.element]);

    for (auto facet : facets) {
      group.add(facet, true);
    }
  }

  group.optimize();
}

/* -------------------------------------------------------------------------- */
void AllSurfaceSelector::onNodesAdded(const Array<UInt> & /* nodes_list */,
                                      const NewNodesEvent & event) {

  if (not aka::is_of_type<CohesiveNewNodesEvent>(event)) {
    return;
  }

  mesh_facets.fillNodesToElements(mesh.getSpatialDimension() - 1);
}

/* -------------------------------------------------------------------------- */
Array<UInt> & AllSurfaceSelector::getMasterList() {
  return mesh_facets.getElementGroup("contact_surface")
      .getNodeGroup()
      .getNodes();
}

/* -------------------------------------------------------------------------- */
Array<UInt> & AllSurfaceSelector::getSlaveList() {
  return mesh_facets.getElementGroup("contact_surface")
      .getNodeGroup()
      .getNodes();
}

#endif

} // namespace akantu
