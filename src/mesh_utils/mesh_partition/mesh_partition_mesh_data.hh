/**
 * @file   mesh_partition_mesh_data.hh
 *
 * @author Dana Christen <dana.christen@epfl.ch>
 *
 * @date   Mon May 01 10:22:00 2013 (On Labor Day, such a shame!)
 *
 * @brief  mesh partitioning based on data provided in the mesh
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

#ifndef __AKANTU_MESH_PARTITION_MESH_DATA_HH__
#define __AKANTU_MESH_PARTITION_MESH_DATA_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh_partition.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class MeshPartitionMeshData : public MeshPartition {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MeshPartitionMeshData(const Mesh & mesh, UInt spatial_dimension,
                        const MemoryID & memory_id = 0);

  MeshPartitionMeshData(const Mesh & mesh,
                        const ByElementTypeArray<UInt> & mapping,
                        UInt spatial_dimension,
                        const MemoryID & memory_id = 0);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void partitionate(UInt nb_part,
                            const EdgeLoadFunctor & edge_load_func = ConstEdgeLoadFunctor(),
                            const Array<UInt> & pairs = Array<UInt>());

  virtual void reorder();

  void setPartitionMapping(const ByElementTypeArray<UInt> & mapping);

  void setPartitionMappingFromMeshData(const std::string & data_name);

  private:

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  const ByElementTypeArray<UInt> * partition_mapping;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */


__END_AKANTU__

#endif /* __AKANTU_MESH_PARTITION_MESH_DATA_HH__ */
