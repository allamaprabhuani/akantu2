/**
 * @file   mesh.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jun 16 12:02:26 2010
 *
 * @brief  class handling meshes
 *
 * @section LICENSE
 *
 * <insert lisence here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "mesh.hh"

/* -------------------------------------------------------------------------- */

Mesh::Mesh(MemoryID memory_id, unsigned int spatial_dimension) :
  Memory(memory_id) {

  nodes = NULL;

}
