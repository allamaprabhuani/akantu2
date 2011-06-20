/**
 * @file   mesh_utils.hh
 * @author Guillaume ANCIAUX <guillaume.anciaux@epfl.ch>
 * @date   Wed Aug 18 14:03:39 2010
 *
 * @brief All mesh utils necessary for various tasks
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

#ifndef __AKANTU_MESH_UTILS_HH__
#define __AKANTU_MESH_UTILS_HH__

#include "aka_common.hh"
#include "mesh.hh"
#include "aka_csr.hh"


__BEGIN_AKANTU__

class MeshUtils {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MeshUtils();
  virtual ~MeshUtils();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// build map from nodes to elements
  static void buildNode2Elements(const Mesh & mesh, CSR<UInt> & node_to_elem, UInt spatial_dimension = 0);
  //  static void buildNode2Elements(const Mesh & mesh, Vector<UInt> & node_offset, Vector<UInt> & node_to_elem, UInt spatial_dimension = 0);
  /// build map from nodes to elements for a specific element type
  static void buildNode2ElementsByElementType(const Mesh & mesh, ElementType type, CSR<UInt> & node_to_elem);
  //  static void buildNode2ElementsByElementType(const Mesh & mesh, ElementType type, Vector<UInt> & node_offset, Vector<UInt> & node_to_elem);
  /// build facets elements : boundary and/or internals
  static void buildFacets(Mesh & mesh, bool boundary_flag = true, bool internal_flag = false);
  /// build normal to some elements
  //  static void buildNormals(Mesh & mesh, UInt spatial_dimension=0);

  /// take  the local_connectivity  array  as  the array  of  local and  ghost
  /// connectivity, renumber the nodes and set the connectivity of the mesh
  static void renumberMeshNodes(Mesh & mesh,
				UInt * local_connectivities,
				UInt nb_local_element,
				UInt nb_ghost_element,
				ElementType type,
				Vector<UInt> & old_nodes);

  static void setUIntData(Mesh & mesh, UInt * data, UInt nb_tags, const ElementType & type);

  /// Detect closed surfaces of the mesh and save the surface id
  /// of the surface elements in the array surface_id
  static void buildSurfaceID(Mesh & mesh);

  /// compute pbc pair for on given direction
  static void computePBCMap(const Mesh & mymesh,const UInt dir,
		     std::map<UInt,UInt> & pbc_pair);


  // /// tweak mesh connectivity to activate pbc
  // static void tweakConnectivityForPBC(Mesh & mesh, 
  // 				      bool flag_x,
  // 				      bool flag_y = false,
  // 				      bool flag_z = false);

  /// create a multimap of nodes per surfaces
  static void buildNodesPerSurface(const Mesh & mesh, CSR<UInt> & nodes_per_surface);

  /// function to print the contain of the class
  //  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "mesh_utils_inline_impl.cc"

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const MeshUtils & _this)
// {
//   _this.printself(stream);
//   return stream;
// }

__END_AKANTU__
#endif /* __AKANTU_MESH_UTILS_HH__ */
