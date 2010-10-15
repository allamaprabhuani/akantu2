/**
 * @file   mesh_io_msh.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jul 14 14:27:11 2010
 *
 * @brief  unit test for the MeshIOMSH class
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
#include <cstdlib>

/* -------------------------------------------------------------------------- */
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"

/* -------------------------------------------------------------------------- */


int main(int argc, char *argv[]) {
  akantu::MeshIOMSH mesh_io;
  akantu::Mesh mesh(3);

  mesh_io.read("./cube.msh", mesh);

  std::cout << mesh << std::endl;

  return EXIT_SUCCESS;
}
