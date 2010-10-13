/**
 * @file   mesh_utils.hh
 * @author Guillaume ANCIAUX <anciaux@lsmscluster1.epfl.ch>
 * @date   Wed Aug 18 14:03:39 2010
 *
 * @brief All mesh utils necessary for various tasks
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MESH_UTILS_HH__
#define __AKANTU_MESH_UTILS_HH__

#include "aka_common.hh"
#include "mesh.hh"


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
  static void buildNode2Elements(const Mesh & mesh, Vector<UInt> & node_offset, Vector<UInt> & node_to_elem);
  /// build map from nodes to elements for a specific element type
  static void buildNode2ElementsByElementType(const Mesh & mesh, ElementType type, Vector<UInt> & node_offset, Vector<UInt> & node_to_elem);
  /// build facets elements : boundary and/or internals
  static void buildFacets(Mesh & mesh, bool boundary_flag=1, bool internal_flag=0);
  /// build normal to some elements
  static void buildNormals(Mesh & mesh, UInt spatial_dimension=0);

  /// take  the local_connectivity  array  as  the array  of  local and  ghost
  /// connectivity, renumber the nodes and set the connectivity of the mesh
  static void renumberMeshNodes(Mesh & mesh,
				UInt * local_connectivities,
				UInt nb_local_element,
				UInt nb_ghost_element,
				ElementType type,
				Vector<UInt> & old_nodes);

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
