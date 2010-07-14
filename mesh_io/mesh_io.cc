/**
 * @file   mesh_io.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jul 14 16:51:22 2010
 *
 * @brief  common part for all mesh io classes
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "common.hh"
#include "mesh_io.hh"

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MeshIO::MeshIO() {
  canReadSurface      = false;
  canReadExtendedData = false;
}

/* -------------------------------------------------------------------------- */
MeshIO::~MeshIO() {

}
/* -------------------------------------------------------------------------- */

__END_AKANTU__
