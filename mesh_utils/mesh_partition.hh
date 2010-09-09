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

  MeshPartition(const Mesh & mesh, UInt spatial_dimension,
		const MemoryID & memory_id = 0);

  virtual ~MeshPartition();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void partitionate(UInt nb_part) = 0;

  /// function to print the contain of the class
  //virtual void printself(std::ostream & stream, int indent = 0) const = 0;

protected:

  /// build the dual graph of the mesh, for all element of spatial_dimension
  void buildDualGraph(Vector<Int> & dxadj, Vector<Int> & dadjncy);

  /// fill the partitions array with a given linearized partition information
  void fillPartitionInformations(const Mesh & mesh, const Int * linearized_partitions);
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Partition, partitions, const Vector<UInt> &);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(GhostPartition,
				   ghost_partitions, const Vector<UInt> &);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(GhostPartitionOffset,
				   ghost_partitions_offset, const Vector<UInt> &);

  AKANTU_GET_MACRO(NbPartition, nb_partitions, UInt);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id
  std::string id;

  /// the mesh to partition
  const Mesh & mesh;

  /// dimension of the elements to consider in the mesh
  UInt spatial_dimension;

  /// number of partitions
  UInt nb_partitions;

  /// partition numbers
  ByElementTypeUInt partitions;

  ByElementTypeUInt ghost_partitions;

  ByElementTypeUInt ghost_partitions_offset;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "mesh_partition_inline_impl.cc"

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const MeshPartition & _this)
// {
//   _this.printself(stream);
//   return stream;
// }


__END_AKANTU__

#endif /* __AKANTU_MESH_PARTITION_H__ */
