#include <iostream>
#include <sstream>
#include "aka_common.hh"
#include "aka_error.hh"
#include "mesh.hh"
#include "solid_mechanics_model.hh"
#include "boundary_condition_functor.hh"

using namespace akantu;

/* -------------------------------------------------------------------------- */

int main(int argc, char* argv[])
{
    UInt spatial_dimension(3);

    akantu::initialize(argc, argv);

    Mesh mesh(spatial_dimension, "mesh_names");

    std::cout << "Loading the mesh." << std::endl;

    mesh.read("./cube_physical_names.msh");

    Boundary & boundary = mesh.getBoundary();
    boundary.createBoundariesFromMeshData("physical_names");
    std::stringstream sstr;

    SolidMechanicsModel model(mesh);
    model.initFull("material.dat");
    std::cout << model.getMaterial(0) << std::endl;

    Vector<Real> surface_traction(3);
    surface_traction(0)=0.0;
    surface_traction(1)=0.0;
    surface_traction(2)=-1.0;

    Matrix<Real> surface_stress(3, 3, 0.0);
    surface_stress(0,0)=0.0;
    surface_stress(1,1)=0.0;
    surface_stress(2,2)=-1.0;

    model.applyBC(BC::Dirichlet::FixedValue(13.0, BC::_x), "Bottom");
    model.applyBC(BC::Dirichlet::FixedValue(13.0, BC::_y), "Bottom");
    model.applyBC(BC::Dirichlet::FixedValue(13.0, BC::_z), "Bottom");

    //model.applyBC(BC::Neumann::FromSameDim(surface_traction), "Top");
    model.applyBC(BC::Neumann::FromHigherDim(surface_stress), "Top");

    debug::setDebugLevel(dblTest);
    boundary.printself(std::cout);
    std::cout << model.getDisplacement();
    std::cout << model.getForce();
    debug::setDebugLevel(dblInfo);

    akantu::finalize();

    return EXIT_SUCCESS;
}


