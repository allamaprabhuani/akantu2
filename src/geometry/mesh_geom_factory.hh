/**
 * @file   mesh_geom_factory.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Feb 26 2015
 * @date last modification: Fri Mar 6 2015
 *
 * @brief  Class for constructing the CGAL primitives of a mesh
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

#ifndef __AKANTU_MESH_GEOM_FACTORY_HH__
#define __AKANTU_MESH_GEOM_FACTORY_HH__

#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_geom_abstract.hh"
#include "tree_type_helper.hh"

#include <CGAL/Cartesian.h>

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

typedef CGAL::Cartesian<Real> K;

/**
 * @brief Class used to construct AABB tree for intersection computations
 *
 * This class constructs a CGAL AABB tree of one type of element in a mesh
 * for fast intersection computations.
 *
 * It also responds to queries on this tree (intersection of the tree with primitives)
 * and construct the resulting mesh.
 */
template<UInt dim, ElementType el_type>
class MeshGeomFactory : public MeshGeomAbstract {

public:
  /// Construct from mesh
  explicit MeshGeomFactory(const Mesh & mesh);

  /// Desctructor
  virtual ~MeshGeomFactory();

public:
  /// Construct AABB tree for fast intersection computing
  virtual void constructData();

  /// Compute the number of intersections with primitive
  virtual UInt numberOfIntersectionsWithInterface(const K::Segment_3 & interface) const;

  /// Create a mesh of the intersection with a linear interface
  virtual void meshOfLinearInterface(const Interface & interface, Mesh & interface_mesh);

  /// Construct a primitive and add it to the list
  /**
   * This function needs to be specialized for every type that is wished to be supported.
   * @param node_coordinates coordinates of the nodes making up the element
   * @param id element number
   * @see MeshGeomContainer
   */
  void addPrimitive(const Matrix<Real> & node_coordinates, UInt id);

  /// Construct segment list from intersections and remove duplicates
  void constructSegments(
      const std::list<typename TreeTypeHelper<dim, el_type>::linear_intersection> & intersections,
      std::list<std::pair<K::Segment_3, UInt> > & segments,
      const K::Segment_3 & interface);

  /// Getter for the AABB tree
  const typename TreeTypeHelper<dim, el_type>::tree & getTree() const { return *data_tree; }

protected:
  /// AABB data tree
  typename TreeTypeHelper<dim, el_type>::tree * data_tree;

  /// Primitive list
  typename TreeTypeHelper<dim, el_type>::container_type primitive_list;
};


__END_AKANTU__

#include "mesh_geom_factory_tmpl.hh"

#endif // __AKANTU_MESH_GEOM_FACTORY_HH__
