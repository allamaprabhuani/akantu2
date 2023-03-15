/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef AKANTU_TRIANGLE_HH_
#define AKANTU_TRIANGLE_HH_

#include "aka_common.hh"

#include "mesh_geom_common.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */

/// Class used for substitution of CGAL::Triangle_3 primitive
template <typename K> class Triangle : public CGAL::Triangle_3<K> {
  using parent = CGAL::Triangle_3<K>;

public:
  /// Default constructor
  Triangle() = default;

  /// Copy constructor
  Triangle(const Triangle & other) = default;
  Triangle(Triangle && other) noexcept = default;

  Triangle & operator=(const Triangle & other) = default;
  Triangle & operator=(Triangle && other) noexcept = default;

  /// Construct from 3 points
  Triangle(const CGAL::Point_3<K> & a, const CGAL::Point_3<K> & b,
           const CGAL::Point_3<K> & c)
      : parent(a, b, c) {}

public:
  UInt id() const { return meshId; }
  void setId(UInt newId) { meshId = newId; }

protected:
  /// Id of the element represented by the primitive
  UInt meshId{0};
};

} // namespace akantu

#endif // AKANTU_TRIANGLE_HH_
