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

} // namespace akantu

/// place here includes

namespace akantu {

/* -------------------------------------------------------------------------- */
inline bool IGFEMEnrichment::isInside(const Vector<Real> & point,
                                      ID domain) const {
  if (domain == "")
    domain = default_geometry;
  Geometry & spheres = this->getGeometry(domain);
  SK::Point_3 p(point(0), point(1), 0.);
  std::list<Spherical::Sphere_3>::const_iterator begin = spheres.begin();
  std::list<Spherical::Sphere_3>::const_iterator end = spheres.end();
  for (; begin != end; ++begin) {
    if (!(begin->has_on_unbounded_side(p)))
      return true;
  }
  return false;
}
