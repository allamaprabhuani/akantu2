/**
 * @file   test_mesh_partition_scotch.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
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
#include "fem.hh"
#include "mesh_io_msh.hh"
#include "mesh_partition_scotch.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER


/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  akantu::initialize(&argc, &argv);

  akantu::debug::setDebugLevel(akantu::dblDump);

  int dim = 2;
#ifdef AKANTU_USE_IOHELPER
  akantu::ElementType type = akantu::_triangle_6;
#endif //AKANTU_USE_IOHELPER
  akantu::Mesh mesh(dim);

  akantu::MeshIOMSH mesh_io;
  mesh_io.read("triangle.msh", mesh);
  akantu::MeshPartition * partition = new akantu::MeshPartitionScotch(mesh, dim);
  partition->partitionate(8);


#ifdef AKANTU_USE_IOHELPER
  unsigned int nb_nodes = mesh.getNbNodes();
  unsigned int nb_element = mesh.getNbElement(type);

  DumperParaview dumper;
  dumper.SetMode(TEXT);
  dumper.SetPoints(mesh.getNodes().values, dim, nb_nodes, "test-scotch-partition");
  dumper.SetConnectivity((int*) mesh.getConnectivity(type).values,
   			 TRIANGLE2, nb_element, C_MODE);

  akantu::UInt  nb_quadrature_points = akantu::FEM::getNbQuadraturePoints(type);
  double * part = new double[nb_element*nb_quadrature_points];
  akantu::UInt * part_val = partition->getPartition(type).values;
  for (unsigned int i = 0; i < nb_element; ++i)
    for (unsigned int q = 0; q < nb_quadrature_points; ++q)
      part[i*nb_quadrature_points + q] = part_val[i];
  dumper.AddElemDataField(part, 1, "partitions");
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();

  delete [] part;

#endif //AKANTU_USE_IOHELPER

  delete partition;

  akantu::finalize();

  return EXIT_SUCCESS;
}
