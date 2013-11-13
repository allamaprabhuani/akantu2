#include <iostream>
#include <sstream>
#include "aka_common.hh"
#include "mesh.hh"

using namespace akantu;

/* -------------------------------------------------------------------------- */

int main(int argc, char* argv[]) {
  UInt spatialDimension(3);

  akantu::initialize(argc, argv);

  Mesh mesh(spatialDimension, "mesh_names");

  std::cout << "Loading the mesh." << std::endl;

  //    mesh.read("./cube_physical_names.msh");
  mesh.read("./cube_physical_names.msh");
  std::stringstream sstr;

  std::cout << "Examining mesh:" << std::endl;

  // Inspection of the number of boundaries
  UInt nb_boundaries= mesh.getNbElementGroups();
  AKANTU_DEBUG_INFO(nb_boundaries << " boundaries advertised initially by Mesh.");
  AKANTU_DEBUG_INFO("Building boundaries");

  // Two methods: either building using data loaded from the mesh file in MeshData
  // or build with automatic numbering
  mesh.createGroupsFromMeshData<std::string>("physical_names");

  // Second inspection of the number of boundaries (should not be 0)
  nb_boundaries = mesh.getNbElementGroups();

  AKANTU_DEBUG_INFO(nb_boundaries << " boundaries advertised by Mesh.");
  AKANTU_DEBUG_ASSERT(nb_boundaries != 0, "No boundary detected!");

  std::cout << (*dynamic_cast<GroupManager*>(&mesh)) << std::endl;

  akantu::finalize();

  return EXIT_SUCCESS;
}


