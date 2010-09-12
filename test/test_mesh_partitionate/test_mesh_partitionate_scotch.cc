/**
 * @file   test_mesh_partition_scotch.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 19 13:05:27 2010
 *
 * @brief  test of internal facet extraction
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
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

  int dim = 2;
  akantu::ElementType type = akantu::_triangle_1;
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
   			 TRIANGLE1, nb_element, C_MODE);
  double * part = new double[nb_element];
  for (unsigned int i = 0; i < nb_element; ++i)
    part[i] = partition->getPartition(type).values[i];
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
