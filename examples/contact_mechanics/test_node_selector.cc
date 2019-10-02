#include "aka_common.hh"
#include "surface_selector.hh"
#include "contact_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {

  const UInt spatial_dimension = 3;

  initialize("material.dat", argc, argv);

  Mesh mesh(spatial_dimension);
  mesh.read("contact_hertz_2d.msh");

  ContactMechanicsModel model(mesh);

  PhysicalSurfaceSelector selector(mesh);
  auto & slave = selector.getSlaveList();
  auto & master = selector.getMasterList();
  
  for (auto & s : slave) {
    std::cerr << s << std::endl;
  }

  for (auto m : master) {
    std::cerr << m << std::endl;
  }
  
  return 0;
}
