/**
 * @file   test_facet_extraction.cc
 * @author Guillaume ANCIAUX <anciaux@lsmscluster1.epfl.ch>
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


int main(int argc, char *argv[])
{
  int dim = 3;
  akantu::ElementType type = akantu::_tetrahedra_1;

  akantu::Mesh mesh(dim);

  akantu::MeshIOMSH mesh_io;
  mesh_io.read("cube.msh", mesh);

  akantu::MeshPartitionScotch partitioner(mesh, dim);

  partitioner.partitionate(5);

#ifdef AKANTU_USE_IOHELPER
  unsigned int nb_nodes = mesh.getNbNodes();
  unsigned int nb_element = mesh.getNbElement(type);

  DumperParaview dumper;
  dumper.SetMode(TEXT);
  dumper.SetPoints(mesh.getNodes().values, dim, nb_nodes, "test-scotch-partition");
  dumper.SetConnectivity((int*)mesh.getConnectivity(type).values,
   			 TETRA1, nb_element, C_MODE);
  double * part = new double[nb_element];
  for (unsigned int i = 0; i < nb_element; ++i)
    part[i] = partitioner.getPartition(type).values[i];
  dumper.AddElemDataField(part, 1, "partitions");
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER

  return EXIT_SUCCESS;
}
