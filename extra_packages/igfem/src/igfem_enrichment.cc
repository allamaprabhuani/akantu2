/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
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
 */

/* -------------------------------------------------------------------------- */
#include "igfem_enrichment.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
IGFEMEnrichment::IGFEMEnrichment(Mesh & mesh) : intersector_sphere(mesh) {}

/* -------------------------------------------------------------------------- */
void IGFEMEnrichment::initialize() { intersector_sphere.init(); }

/* -------------------------------------------------------------------------- */
void IGFEMEnrichment::update(ID domain) {
  if (domain == "")
    domain = default_geometry;
  Geometry & geometry = getGeometry(domain);
  intersector_sphere.buildIGFEMMeshFromSpheres(geometry);
}

/* -------------------------------------------------------------------------- */
void IGFEMEnrichment::unRegisterGeometryObject(const ID & domain) {

  GeometryMap::iterator it = geometries.find(domain);
  AKANTU_DEBUG_ASSERT(it != geometries.end(), "Geometry object with domain "
                                                  << domain
                                                  << " was not found");
  geometries.erase(it);
  if (!geometries.empty())
    default_geometry = (*geometries.begin()).first;
}

/* -------------------------------------------------------------------------- */
void IGFEMEnrichment::registerGeometryObject(Geometry & geometry,
                                             const ID & domain) {
  if (geometries.size() == 0)
    default_geometry = domain;

#ifndef AKANTU_NDEBUG
  GeometryMap::iterator it = geometries.find(domain);
  AKANTU_DEBUG_ASSERT(it == geometries.end(), "Geometry object with domain "
                                                  << domain
                                                  << " was already created");
#endif

  std::stringstream sstr;
  sstr << "geometry:" << domain;
  geometries[domain] = &geometry;
}

/* -------------------------------------------------------------------------- */
IGFEMEnrichment::Geometry & IGFEMEnrichment::getGeometry(ID & domain) const {
  AKANTU_DEBUG_IN();

  if (domain == "")
    domain = default_geometry;

  GeometryMap::const_iterator it = geometries.find(domain);
  AKANTU_DEBUG_ASSERT(it != geometries.end(),
                      "The geometry " << domain << " is not registered");

  AKANTU_DEBUG_OUT();
  return *(it->second);
}

/* -------------------------------------------------------------------------- */
void IGFEMEnrichment::moveInterface(Real new_position, ID domain) {
  if (domain == "")
    domain = default_geometry;
  Geometry & geometry = getGeometry(domain);

  /// for this type of IGFEM enrichment the geometry consists of a list of
  /// spheres
  /// -> need to loop over spheres and change their radius,
  /// which specifies the position of interfaces
  Geometry::const_iterator query_it = geometry.begin();
  Geometry sphere_list;
  for (; query_it != geometry.end(); ++query_it) {
    SK::Sphere_3 sphere(query_it->center(), new_position * new_position);
    sphere_list.push_back(sphere);
  }
  geometry.zero();
  geometry = sphere_list;

  this->update(domain);
}

} // namespace akantu
