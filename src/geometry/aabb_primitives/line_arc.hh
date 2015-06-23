/**
 * @file   segment.hh
 *
 * @author Clément Roux-Langlois <clement.roux@epfl.ch>
 *
 * @date creation: Wed Mar 13 2015
 * @date last modification: Wed Mar 13 2015
 *
 * @brief  Segment classe (geometry) for AABB CGAL algos
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_LINE_ARC_HH__
#define __AKANTU_LINE_ARC_HH__

#include "aka_common.hh"

#include "mesh_geom_common.hh"

__BEGIN_AKANTU__
  
/* -------------------------------------------------------------------------- */

/// Class used for substitution of CGAL::Triangle_3 primitive
template<typename K>
class Line_arc : public CGAL::Line_arc_3<K> {
public:
  /// Default constructor
  Line_arc() :
    CGAL::Line_arc_3<K>(), meshId(0) {}

  /// Copy constructor
  Line_arc(const Line_arc & other) :
    CGAL::Line_arc_3<K>(other), meshId(other.meshId) {}

  /// Construct from 3 points
  // "CGAL-4.5/doc_html/Circular_kernel_3/classCGAL_1_1Line__arc__3.html"
  Line_arc(const CGAL::Line_3<K> & l, const CGAL::Circular_arc_point_3<K> & a,
	  const CGAL::Circular_arc_point_3<K> & b):
    CGAL::Line_arc_3<K>(l, a, b), meshId(0) {}

public:
  UInt id() const { return meshId; }
  void setId(UInt newId) { meshId = newId; }

protected:
  /// Id of the element represented by the primitive
  UInt meshId;
};

__END_AKANTU__

#endif // __AKANTU_LINE_ARC_HH__
