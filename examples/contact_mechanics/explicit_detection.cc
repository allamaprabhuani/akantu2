/* -------------------------------------------------------------------------- */
#include "contact_mechanics_model.hh"

#include "dumper_elemental_field.hh"
#include "dumper_nodal_field.hh"

#include "dumper_iohelper_paraview.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {

  std::string mesh_file = "explicit.msh";
  std::string material_file = "options_explicit.dat";
    
  const UInt spatial_dimension = 2;
  initialize(material_file, argc, argv);
  
  Mesh mesh(spatial_dimension);
  mesh.read(mesh_file);

  ContactDetector detector(mesh);

  auto &&surface_selector = std::make_shared<PhysicalSurfaceSelector>(mesh);
  detector.setSurfaceSelector(surface_selector);
   
  std::map<UInt, ContactElement> contact_map;

  detector.search(contact_map);

  UInt nb_nodes = mesh.getNbNodes();

  Array<Real> gaps(nb_nodes, 1);
  Array<Real> normals(nb_nodes, spatial_dimension);
  
  for (auto & entry: contact_map) {
    const auto & node = entry.first;
    const auto & element = entry.second;

    gaps[node] = element.gap;
    for (auto i : arange(spatial_dimension)) {
      normals(node, i) = element.normal[i];
    }
  }

  DumperParaview dumper("explicit_detection", "./paraview", false);
  dumper.registerMesh(mesh);

  auto gap_field = std::make_shared<dumper::NodalField<Real>>(gaps);
  auto normal_field = std::make_shared<dumper::NodalField<Real>>(normals);
  
  dumper.registerField("gaps", gap_field);
  dumper.registerField("normals", normal_field);
  dumper.dump();
  
  finalize();
  return EXIT_SUCCESS;
}

