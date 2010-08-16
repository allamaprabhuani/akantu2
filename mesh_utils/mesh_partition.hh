/**
 * @file   mesh_partition.h
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 12 16:24:40 2010
 *
 * @brief  tools to partitionate a mesh
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MESH_PARTITION_H__
#define __AKANTU_MESH_PARTITION_H__

/* -------------------------------------------------------------------------- */
#include "aka_memory.hh"
#include "mesh.hh"

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

class MeshPartition : public Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  MeshPartition(const Mesh & mesh, const MemoryID & memory_id = 0) :
    Memory(memory_id), mesh(mesh) {};
  virtual ~MeshPartition();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  
  virtual partitionate() = 0;

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const = 0;
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  const Mesh & mesh;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "mesh_partition_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MeshPartition & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_MESH_PARTITION_H__ */
