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

__BEGIN_AKANTU__

Mesh::Mesh(unsigned int spatial_dimension, MemoryID memory_id) :
  Memory(memory_id) {

  this->spatial_dimension = spatial_dimension;
  nodes = NULL;
}

__END_AKANTU__
