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



#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  int dim = 3;

  /// load mesh
  Mesh my_mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("cubes.msh", my_mesh);
  
  /// build facet connectivity and surface id
  MeshUtils::buildFacets(my_mesh,1,0);
  MeshUtils::buildSurfaceID(my_mesh);

  unsigned int nb_nodes = my_mesh.getNbNodes();

  /// declaration of model
  SolidMechanicsModel  my_model(my_mesh);
  /// model initialization
  my_model.initVectors();
  /// initialize the vectors
  //UInt nb_nodes = my_model.getFEM().getMesh().getNbNodes();
  memset(my_model.getForce().values,        0, 3*nb_nodes*sizeof(Real));
  memset(my_model.getVelocity().values,     0, 3*nb_nodes*sizeof(Real));
  memset(my_model.getAcceleration().values, 0, 3*nb_nodes*sizeof(Real));
  memset(my_model.getDisplacement().values, 0, 3*nb_nodes*sizeof(Real));


  /// dump facet and surface information to paraview
#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  dumper.SetMode(TEXT);
  
  dumper.SetPoints(my_mesh.getNodes().values, dim, nb_nodes, "test-surface-extraction");
  dumper.SetConnectivity((int*)my_mesh.getConnectivity(_tetrahedra_1).values,
   			 TETRA1, my_mesh.getNbElement(_tetrahedra_1), C_MODE);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
  
  DumperParaview dumper_surface;
  dumper_surface.SetMode(TEXT);

  dumper_surface.SetPoints(my_mesh.getNodes().values, dim, nb_nodes, "test-surface-extraction_boundary");
  
  dumper_surface.SetConnectivity((int *)my_mesh.getConnectivity(_triangle_1).values,
				 TRIANGLE1, my_mesh.getNbElement(_triangle_1), C_MODE);
  double * surf_id = new double [my_mesh.getSurfaceId(_triangle_1).getSize()];
  for (UInt i = 0; i < my_mesh.getSurfaceId(_triangle_1).getSize(); ++i)
    surf_id[i] = (double)my_mesh.getSurfaceId(_triangle_1).values[i];
  dumper_surface.AddElemDataField(surf_id, 1, "surface_id");
  delete [] surf_id;
  dumper_surface.SetPrefix("paraview/");
  dumper_surface.Init();
  dumper_surface.Dump();
#endif //AKANTU_USE_IOHELPER



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

  finalize();

  return EXIT_SUCCESS;
}
