/**
 * @file   mesh_geom_abstract.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Thu Jan 14 2016
 *
 * @brief  Class for constructing the CGAL primitives of a mesh
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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

#ifndef __AKANTU_MESH_GEOM_ABSTRACT_HH__
#define __AKANTU_MESH_GEOM_ABSTRACT_HH__

#include "aka_common.hh"
#include "mesh.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/// Abstract class for mesh geometry operations
class MeshGeomAbstract {

public:
  /// Construct from mesh
  explicit MeshGeomAbstract(Mesh & mesh) : mesh(mesh) {};

  /// Destructor
  virtual ~MeshGeomAbstract() {};

public:
  /// Construct geometric data for computational geometry algorithms
  virtual void constructData(GhostType ghost_type = _not_ghost) = 0;

protected:
  /// Mesh used to construct the primitives
  Mesh & mesh;
};

__END_AKANTU__

#endif // __AKANTU_MESH_GEOM_ABSTRACT_HH__
