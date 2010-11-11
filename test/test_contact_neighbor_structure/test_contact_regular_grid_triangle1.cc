/**
 * @file   test_contact_regular_grid.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Oct 26 16:58:42 2010
 *
 * @brief  test regular grid neighbor structure for 3d case
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "contact.hh"
#include "contact_neighbor_structure.hh"



#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  int dim = 2;

  /// load mesh
  Mesh my_mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("squares.msh", my_mesh);

  /// build facet connectivity and surface id
  MeshUtils::buildFacets(my_mesh,1,0);
  MeshUtils::buildSurfaceID(my_mesh);
 
  unsigned int nb_nodes = my_mesh.getNbNodes();

  
  /// dump facet and surface information to paraview
#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  dumper.SetMode(TEXT);
  
  dumper.SetPoints(my_mesh.getNodes().values, dim, nb_nodes, "test-surface-extraction_2d");
  dumper.SetConnectivity((int*)my_mesh.getConnectivity(_triangle_1).values,
   			 TRIANGLE1, my_mesh.getNbElement(_triangle_1), C_MODE);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
  
  DumperParaview dumper_surface;
  dumper_surface.SetMode(TEXT);

  dumper_surface.SetPoints(my_mesh.getNodes().values, dim, nb_nodes, "test-surface-extraction_boundary_2d");
  
  dumper_surface.SetConnectivity((int *)my_mesh.getConnectivity(_line_1).values,
				 LINE1, my_mesh.getNbElement(_line_1), C_MODE);
  double * surf_id = new double [my_mesh.getSurfaceId(_line_1).getSize()];
  for (UInt i = 0; i < my_mesh.getSurfaceId(_line_1).getSize(); ++i)
    surf_id[i] = (double)my_mesh.getSurfaceId(_line_1).values[i];
  dumper_surface.AddElemDataField(surf_id, 1, "surface_id");
  delete [] surf_id;
  dumper_surface.SetPrefix("paraview/");
  dumper_surface.Init();
  dumper_surface.Dump();
#endif //AKANTU_USE_IOHELPER

  /// declaration of model
  SolidMechanicsModel  my_model(my_mesh);
  /// model initialization
  my_model.initVectors();
  /// initialize the vectors
  //UInt nb_nodes = my_model.getFEM().getMesh().getNbNodes();
  memset(my_model.getForce().values,        0, dim*nb_nodes*sizeof(Real));
  memset(my_model.getVelocity().values,     0, dim*nb_nodes*sizeof(Real));
  memset(my_model.getAcceleration().values, 0, dim*nb_nodes*sizeof(Real));
  memset(my_model.getDisplacement().values, 0, dim*nb_nodes*sizeof(Real));
  
  my_model.readMaterials("material.dat");
  my_model.initMaterials();
  my_model.initModel();

  Real time_step = my_model.getStableTimeStep();
  my_model.setTimeStep(time_step/10.);

  my_model.assembleMassLumped();

  //std::cout << *my_model << std::endl;

  /// contact declaration
  Contact * my_contact = Contact::newContact(my_model, 
					     _ct_3d_expli, 
					     _cst_3d_expli, 
					     _cnst_regular_grid);

  my_contact->initContact(false);

  Surface master = 0;

  my_contact->addMasterSurface(master);
  
  NeighborList * my_neighbor_list = const_cast<ContactNeighborStructure&>(my_contact->getContactSearch().getContactNeighborStructure(master)).getNeighborList();

#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper_neighbor;
  dumper_neighbor.SetMode(TEXT);

  dumper_neighbor.SetPoints(my_mesh.getNodes().values, dim, nb_nodes, "test-neighbor-elements_2d");
  
  dumper_neighbor.SetConnectivity((int *)my_mesh.getConnectivity(_line_1).values,
				 LINE1, my_mesh.getNbElement(_line_1), C_MODE);

  UInt nb_nodes_neigh = my_neighbor_list->nb_nodes;
  Vector<UInt> impact_nodes = my_neighbor_list->impactor_nodes;
  UInt * impact_nodes_val = impact_nodes.values;

  /// print impactor nodes
  std::cout << "we have " << nb_nodes_neigh << " impactor nodes:";
  for (UInt i = 0; i < nb_nodes_neigh; ++i) {
    std::cout << " " << impact_nodes_val[i];
  }
  std::cout << std::endl;

  UInt * node_to_elem_offset_val = my_neighbor_list->facets_offset[_line_1]->values;
  UInt * node_to_elem_val = my_neighbor_list->facets[_line_1]->values;

  double * neigh_elem = new double [my_mesh.getNbElement(_line_1)];
  for (UInt i = 0; i < my_mesh.getNbElement(_line_1); ++i)
    neigh_elem[i] = 0.0; 
  
  UInt visualize_node = 1;
  UInt n = impact_nodes_val[visualize_node];
  std::cout << "plot for node: " << n << std::endl;
  for (UInt i = node_to_elem_offset_val[visualize_node]; i < node_to_elem_offset_val[visualize_node+1]; ++i)
    neigh_elem[node_to_elem_val[i]] = 1.;

  dumper_neighbor.AddElemDataField(neigh_elem, 1, "neighbor id");
 
  dumper_neighbor.SetPrefix("paraview/");
  dumper_neighbor.Init();
  dumper_neighbor.Dump();

  delete [] neigh_elem;

#endif //AKANTU_USE_IOHELPER
  
  finalize();

  return EXIT_SUCCESS;
}
