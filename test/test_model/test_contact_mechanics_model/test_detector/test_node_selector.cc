#include "aka_common.hh"
#include "node_selector.hh"
#include "contact_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {

  const UInt spatial_dimension = 3;

  initialize("options.dat", argc, argv);

  Mesh mesh(spatial_dimension);
  mesh.read("node.msh");

  ContactMechanicsModel model(mesh);

  PhysicalSurfaceNodeSelector selector(model);
  auto & slave = selector.getSlaveList();

  for (auto & s : slave) {
    std::cerr << s << std::endl;
  }
  
  return 0;
}
