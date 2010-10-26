/**
 * @file   test_surface_extraction_2d.cc
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @date   Mon Oct 25 09:47:15 2010
 *
 * @brief  
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

#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  int dim = 2;

  Mesh mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("squares.msh", mesh);
  
  MeshUtils::buildFacets(mesh,1,0);
  MeshUtils::buildSurfaceID(mesh);

  unsigned int nb_nodes = mesh.getNbNodes();
#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  dumper.SetMode(TEXT);

  dumper.SetPoints(mesh.getNodes().values, dim, nb_nodes, "test-surface-extraction");
  dumper.SetConnectivity((int*)mesh.getConnectivity(_triangle_1).values,
   			 TRIANGLE1, mesh.getNbElement(_triangle_1), C_MODE);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();

  DumperParaview dumper_surface;
  dumper_surface.SetMode(TEXT);

  dumper_surface.SetPoints(mesh.getNodes().values, dim, nb_nodes, "test-surface-extraction_boundary");
  
  dumper_surface.SetConnectivity((int *)mesh.getConnectivity(_line_1).values,
			       LINE1, mesh.getNbElement(_line_1), C_MODE);
  double * surf_id = new double [mesh.getSurfaceId(_line_1).getSize()];
  for (UInt i = 0; i < mesh.getSurfaceId(_line_1).getSize(); ++i)
    surf_id[i] = (double)mesh.getSurfaceId(_line_1).values[i];
  dumper_surface.AddElemDataField(surf_id, 1, "surface_id");
  delete [] surf_id;
  dumper_surface.SetPrefix("paraview/");
  dumper_surface.Init();
  dumper_surface.Dump();

#endif //AKANTU_USE_IOHELPER

  return EXIT_SUCCESS;
}
