/**
 * @file   mesh_partition.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 12 16:24:40 2010
 *
 * @brief  tools to partitionate a mesh
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

#ifndef __AKANTU_MESH_PARTITION_HH__
#define __AKANTU_MESH_PARTITION_HH__

/* -------------------------------------------------------------------------- */
#include "aka_memory.hh"
#include "mesh.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class MeshPartition : protected Memory {
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

  virtual void reorder() = 0;

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

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Partition, partitions, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(GhostPartition, ghost_partitions, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(GhostPartitionOffset, ghost_partitions_offset, UInt);

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

  Vector<UInt> * permutation;
};


/* -------------------------------------------------------------------------- */
/* __aka_inline__ functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "mesh_partition_inline_impl.cc"

/// standard output stream operator
// __aka_inline__ std::ostream & operator <<(std::ostream & stream, const MeshPartition & _this)
// {
//   _this.printself(stream);
//   return stream;
// }

__END_AKANTU__

#ifdef AKANTU_USE_SCOTCH
#include "mesh_partition_scotch.hh"
#endif

#endif /* __AKANTU_MESH_PARTITION_HH__ */
