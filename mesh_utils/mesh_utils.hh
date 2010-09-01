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
  void static buildNode2Elements(const Mesh & mesh, Vector<UInt> & node_offset, Vector<UInt> & node_to_elem);
  /// build facets elements : boundary and/or internals
  void static buildFacets(Mesh & mesh, bool boundary_flag=1, bool internal_flag=0);
  /// build normal to some elements
  void static buildNormals(Mesh & mesh, UInt spatial_dimension=0);

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
