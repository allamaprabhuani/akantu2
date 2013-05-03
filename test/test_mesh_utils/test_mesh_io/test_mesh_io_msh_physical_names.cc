#include <iostream>
#include <sstream>
#include "aka_common.hh"
#include "mesh.hh"

using namespace akantu;

/* -------------------------------------------------------------------------- */

int main(int argc, char* argv[])
{
    UInt spatialDimension(3);
    
    akantu::initialize(argc, argv);
    
    Mesh mesh(spatialDimension);

    mesh.read("./cube_physical_names.msh");
    std::stringstream sstr;
    
    for(Mesh::type_iterator type_it = mesh.firstType(); type_it != mesh.lastType(); ++type_it) {
      const Array<std::string> & name_vec = mesh.getData<std::string>(*type_it, std::string("physical_names") );
      for(UInt i(0); i < name_vec.getSize(); i++) {
        std::cout << "Element " << i << " (of type " << *type_it << ") has physical name " << name_vec(i) << "." << std::endl;
      }
    }
    
    akantu::finalize();
    
    return EXIT_SUCCESS;
}



