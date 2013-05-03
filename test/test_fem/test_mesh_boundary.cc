#include <iostream>
#include <sstream>
#include "aka_common.hh"
#include "mesh.hh"
#include "sub_boundary.hh"

using namespace akantu;

/* -------------------------------------------------------------------------- */

int main(int argc, char* argv[])
{
    UInt spatialDimension(3);

    akantu::initialize(argc, argv);

    Mesh mesh(spatialDimension, "mesh_names");

    std::cout << "Loading the mesh." << std::endl;

//    mesh.read("./cube_physical_names.msh");
    mesh.read("./cube_physical_names.msh");
    std::stringstream sstr;

    // References to Boud
    const Boundary & boundary = mesh.getBoundary();
    Boundary & boundary_not_const = mesh.getBoundary();

    std::cout << "Examining mesh:" << std::endl;

    // Inspection of the number of boundaries
    UInt mesh_nb_boundaries_mesh;
    UInt mesh_nb_boundaries_boundary;
    mesh_nb_boundaries_mesh     = mesh.getNbBoundaries();
    mesh_nb_boundaries_boundary = boundary.getNbBoundaries();
    std::cout << "-- " << mesh_nb_boundaries_mesh << " boundaries advertised initially by Mesh." << std::endl;
    std::cout << "-- " << mesh_nb_boundaries_boundary << " boundaries advertised initially by Boundary." << std::endl;
    AKANTU_DEBUG_ASSERT(mesh_nb_boundaries_mesh == mesh_nb_boundaries_boundary,
      "Not the same number of boundaries advertised by Mesh and Boundary (should be 0 both)!");

    std::cout << "-- Building boundaries" <<  std::endl;

    // Two methods: either building using data loaded from the mesh file in MeshData
    // or build with automatic numbering
//    boundary_not_const.createBoundariesFromGeometry();
    boundary_not_const.createBoundariesFromMeshData("physical_names");

    // Second inspection of the number of boundaries (should not be 0)
    mesh_nb_boundaries_mesh     = mesh.getNbBoundaries();
    mesh_nb_boundaries_boundary = boundary.getNbBoundaries();

    std::cout << "-- " << mesh_nb_boundaries_mesh << " boundaries advertised by Mesh." << std::endl;
    std::cout << "-- " << mesh_nb_boundaries_boundary << " boundaries advertised by Boundary." << std::endl;
    AKANTU_DEBUG_ASSERT(mesh_nb_boundaries_mesh == mesh_nb_boundaries_boundary,
      "Not the same number of boundaries advertised by Mesh and Boundary!");
    AKANTU_DEBUG_ASSERT(mesh_nb_boundaries_mesh != 0,
      "No boundary detected!");

    UInt mesh_count = 0;

    // Loop over (Sub)Boundar(ies)
    for(Boundary::const_iterator it(boundary.begin()); it != boundary.end(); ++it) {
      it->printself(std::cout);
      mesh_count++;
    }

    // Check if the right number of boundaries were retrieved
    AKANTU_DEBUG_ASSERT(mesh_nb_boundaries_boundary == mesh_count,
      "Iterators do not yield the same number of boundaries as the corresponding getter!");

    akantu::finalize();

    return EXIT_SUCCESS;
}


