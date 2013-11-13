/**
 * @file   mesh_partition.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Aug 16 13:17:22 2010
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
#include "aka_csr.hh"
#include "mesh.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class MeshPartition : protected Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MeshPartition(const Mesh & mesh, UInt spatial_dimension,
                const ID & id = "MeshPartitioner",
		const MemoryID & memory_id = 0);

  virtual ~MeshPartition();

  class EdgeLoadFunctor {
  public:
    virtual Int operator()(__attribute__((unused)) const Element & el1,
			   __attribute__((unused)) const Element & el2) const = 0;
  };

  class ConstEdgeLoadFunctor : public EdgeLoadFunctor {
  public:
    virtual inline Int operator()(__attribute__((unused)) const Element & el1,
				  __attribute__((unused)) const Element & el2) const {
      return 1;
    }
  };


  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void partitionate(UInt nb_part,
			    const EdgeLoadFunctor & edge_load_func = ConstEdgeLoadFunctor(),
			    const Array<UInt> & pairs = Array<UInt>()) = 0;

  virtual void reorder() = 0;

  /// function to print the contain of the class
  //virtual void printself(std::ostream & stream, int indent = 0) const = 0;

  /// fill the partitions array with a given linearized partition information
  void fillPartitionInformation(const Mesh & mesh, const Int * linearized_partitions);

protected:

  /// build the dual graph of the mesh, for all element of spatial_dimension
  void buildDualGraph(Array<Int> & dxadj, Array<Int> & dadjncy,
		      Array<Int> & edge_loads,
		      const EdgeLoadFunctor & edge_load_func);

  /// tweak the mesh to handle the PBC pairs
  void tweakConnectivity(const Array<UInt> & pairs);
  /// restore the mesh that has been tweaked
  void restoreConnectivity();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(Partitions, partitions, const ByElementTypeUInt &);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Partition, partitions, UInt);
  // AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(GhostPartition, ghost_partitions, UInt);
  // AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(GhostPartitionOffset, ghost_partitions_offset, UInt);

  AKANTU_GET_MACRO(GhostPartitionCSR, ghost_partitions_csr, const ByElementType< CSR<UInt> > &);

  AKANTU_GET_MACRO(NbPartition, nb_partitions, UInt);
  AKANTU_SET_MACRO(NbPartition, nb_partitions, UInt);

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

  ByElementType< CSR<UInt> > ghost_partitions_csr;
  ByElementTypeUInt ghost_partitions;
  ByElementTypeUInt ghost_partitions_offset;

  Array<UInt> * permutation;

  ByElementTypeUInt saved_connectivity;

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

#ifdef AKANTU_USE_SCOTCH
#include "mesh_partition_scotch.hh"
#endif

#endif /* __AKANTU_MESH_PARTITION_HH__ */
