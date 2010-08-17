/**
 * @file   mesh_partition_scotch.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Aug 13 11:54:11 2010
 *
 * @brief  implementation of the MeshPartitionScotch class
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "mesh_partition_scotch.hh"

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MeshPartitionScotch::MeshPartitionScotch(const Mesh & mesh) : MeshPartition(mesh) {

}

/* -------------------------------------------------------------------------- */
void MeshPartitionScotch::partitionate();

__END_AKANTU__

