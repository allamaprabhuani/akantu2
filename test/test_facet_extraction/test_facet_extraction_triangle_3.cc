/**
 * @file   test_facet_extraction.cc
 * @author Guillaume ANCIAUX <anciaux@lsmscluster1.epfl.ch>
 * @date   Thu Aug 19 13:05:27 2010
 *
 * @brief  test of internal facet extraction
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  int dim = 2;

  Mesh mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("square.msh", mesh);
  
  MeshUtils::buildFacets(mesh,1,1);

  unsigned int nb_nodes = mesh.getNbNodes();
#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  dumper.SetMode(TEXT);

  dumper.SetPoints(mesh.getNodes().values, dim, nb_nodes, "test-facet-extraction");
  dumper.SetConnectivity((int*)mesh.getConnectivity(_triangle_3).values,
   			 TRIANGLE1, mesh.getNbElement(_triangle_3), C_MODE);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();

  DumperParaview dumper_facet;
  dumper_facet.SetMode(TEXT);

  dumper_facet.SetPoints(mesh.getNodes().values, dim, nb_nodes, "test-facet-extraction_boundary");
  dumper_facet.SetConnectivity((int*)mesh.getConnectivity(_segment_2).values,
			       LINE1, mesh.getNbElement(_segment_2), C_MODE);
  dumper_facet.SetPrefix("paraview/");
  dumper_facet.Init();
  dumper_facet.Dump();

  dumper_facet.SetPoints(mesh.getNodes().values, dim, nb_nodes, "test-facet-extraction_internal");
  dumper_facet.SetConnectivity((int*)mesh.getInternalFacetsMesh().getConnectivity(_segment_2).values,
  			       LINE1, mesh.getInternalFacetsMesh().getNbElement(_segment_2), C_MODE);
  dumper_facet.Init();
  dumper_facet.Dump();

#endif //AKANTU_USE_IOHELPER

  return EXIT_SUCCESS;
}
