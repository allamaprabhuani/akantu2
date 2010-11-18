/**
 * @file   test_synchronizer_communication.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 19 13:05:27 2010
 *
 * @brief  test to synchronize barycenters
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "static_communicator.hh"
#include "communicator.hh"
#include "ghost_synchronizer.hh"
#include "mesh.hh"
#include "mesh_io_msh.hh"
#include "mesh_partition_scotch.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER


class TestSynchronizer : public akantu::GhostSynchronizer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  TestSynchronizer(const akantu::Mesh & mesh);
  ~TestSynchronizer();

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(GhostBarycenter, ghost_barycenter, const akantu::Vector<akantu::Real>);

  /* ------------------------------------------------------------------------ */
  /* Ghost Synchronizer inherited members                                     */
  /* ------------------------------------------------------------------------ */
protected:
  virtual akantu::UInt getNbDataToPack(const akantu::Element & element,
				       akantu::GhostSynchronizationTag tag) const;
  virtual akantu::UInt getNbDataToUnpack(const akantu::Element & element,
					 akantu::GhostSynchronizationTag tag) const;
  virtual void packData(akantu::Real ** buffer,
			const akantu::Element & element,
			akantu::GhostSynchronizationTag tag) const;
  virtual void unpackData(akantu::Real ** buffer,
			  const akantu::Element & element,
			  akantu::GhostSynchronizationTag tag) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  std::string id;

  akantu::ByElementTypeReal ghost_barycenter;

  const akantu::Mesh & mesh;
};

/* -------------------------------------------------------------------------- */
/* TestSynchronizer implementation                                            */
/* -------------------------------------------------------------------------- */
TestSynchronizer::TestSynchronizer(const akantu::Mesh & mesh) : mesh(mesh) {
  akantu::UInt spatial_dimension = mesh.getSpatialDimension();

  id = "test_synchronizer";

  akantu::Mesh::ConnectivityTypeList::const_iterator it;
  const akantu::Mesh::ConnectivityTypeList & ghost_type_list = mesh.getConnectivityTypeList(akantu::_ghost);
  for(it = ghost_type_list.begin(); it != ghost_type_list.end(); ++it) {
    if(akantu::Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    akantu::UInt nb_ghost_element = mesh.getNbGhostElement(*it);
    ghost_barycenter[*it] = new akantu::Vector<akantu::Real>(nb_ghost_element,
							     spatial_dimension,
							     std::numeric_limits<akantu::Real>::quiet_NaN(),
							     "ghost_barycenter");

  }
}

TestSynchronizer::~TestSynchronizer() {
  akantu::UInt spatial_dimension = mesh.getSpatialDimension();

  akantu::Mesh::ConnectivityTypeList::const_iterator it;
  const akantu::Mesh::ConnectivityTypeList & ghost_type_list = mesh.getConnectivityTypeList(akantu::_ghost);
  for(it = ghost_type_list.begin(); it != ghost_type_list.end(); ++it) {
    if(akantu::Mesh::getSpatialDimension(*it) != spatial_dimension) continue;
    delete ghost_barycenter[*it];
  }
}

akantu::UInt TestSynchronizer::getNbDataToPack(const akantu::Element & element,
					       akantu::GhostSynchronizationTag tag) const {
  return akantu::Mesh::getSpatialDimension(element.type);
}

akantu::UInt TestSynchronizer::getNbDataToUnpack(const akantu::Element & element,
						 akantu::GhostSynchronizationTag tag) const {
  return akantu::Mesh::getSpatialDimension(element.type);
}

void TestSynchronizer::packData(akantu::Real ** buffer,
				const akantu::Element & element,
				akantu::GhostSynchronizationTag tag) const {
  akantu::UInt spatial_dimension = akantu::Mesh::getSpatialDimension(element.type);
  mesh.getBarycenter(element.element, element.type, *buffer);

  *buffer += spatial_dimension;
}

void TestSynchronizer::unpackData(akantu::Real ** buffer,
				  const akantu::Element & element,
				  akantu::GhostSynchronizationTag tag) const {
  akantu::UInt spatial_dimension = akantu::Mesh::getSpatialDimension(element.type);
  memcpy(ghost_barycenter[element.type]->values + element.element * spatial_dimension,
	 *buffer, spatial_dimension * sizeof(akantu::Real));

  *buffer += spatial_dimension;
}




/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  akantu::initialize(&argc, &argv);

  int dim = 2;
  akantu::ElementType type = akantu::_triangle_3;
  akantu::Mesh mesh(dim);

  akantu::StaticCommunicator * comm = akantu::StaticCommunicator::getStaticCommunicator();
  akantu::Int psize = comm->getNbProc();
  akantu::Int prank = comm->whoAmI();

  akantu::Communicator * communicator;
  if(prank == 0) {
    akantu::MeshIOMSH mesh_io;
    mesh_io.read("triangle.msh", mesh);
    akantu::MeshPartition * partition = new akantu::MeshPartitionScotch(mesh, dim);
    partition->partitionate(psize);
    communicator = akantu::Communicator::createCommunicatorDistributeMesh(mesh, partition);
    delete partition;
  } else {
    communicator = akantu::Communicator::createCommunicatorDistributeMesh(mesh, NULL);
  }

  comm->barrier();

  AKANTU_DEBUG_INFO("Creating TestSynchronizer");
  TestSynchronizer test_synchronizer(mesh);
  test_synchronizer.registerSynchronizer(*communicator);

  AKANTU_DEBUG_INFO("Registering tag");
  test_synchronizer.registerTag(akantu::_gst_test, "barycenter");

  AKANTU_DEBUG_INFO("Synchronizing tag");
  test_synchronizer.synchronize(akantu::_gst_test);


  akantu::Mesh::ConnectivityTypeList::const_iterator it;
  const akantu::Mesh::ConnectivityTypeList & ghost_type_list = mesh.getConnectivityTypeList(akantu::_ghost);
  for(it = ghost_type_list.begin(); it != ghost_type_list.end(); ++it) {
    if(akantu::Mesh::getSpatialDimension(*it) != mesh.getSpatialDimension()) continue;

    akantu::UInt spatial_dimension = akantu::Mesh::getSpatialDimension(*it);
    akantu::UInt nb_ghost_element = mesh.getNbGhostElement(*it);

    for (akantu::UInt el = 0; el < nb_ghost_element; ++el) {
      akantu::Real barycenter[spatial_dimension];
      mesh.getBarycenter(el, *it, barycenter, akantu::_ghost);

      for (akantu::UInt i = 0; i < spatial_dimension; ++i) {
	if(fabs(barycenter[i] - test_synchronizer.getGhostBarycenter(*it).values[el * spatial_dimension + i])
	   > std::numeric_limits<akantu::Real>::epsilon())
	  AKANTU_DEBUG_ERROR("The barycenter of ghost element " << el
			     << " on proc " << prank
			     << " does not match the one get during synchronisation" );
      }

    }
  }


#ifdef AKANTU_USE_IOHELPER
  unsigned int nb_nodes = mesh.getNbNodes();
  unsigned int nb_element = mesh.getNbElement(type);

  DumperParaview dumper;
  dumper.SetMode(TEXT);
  dumper.SetParallelContext(prank, psize);
  dumper.SetPoints(mesh.getNodes().values, dim, nb_nodes, "test-scotch-partition");
  dumper.SetConnectivity((int*) mesh.getConnectivity(type).values,
   			 TRIANGLE1, nb_element, C_MODE);
  double * part = new double[nb_element];
  for (unsigned int i = 0; i < nb_element; ++i)
    part[i] = prank;
  dumper.AddElemDataField(part, 1, "partitions");
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();

  delete [] part;

  unsigned int nb_ghost_element = mesh.getNbGhostElement(type);
  DumperParaview dumper_ghost;
  dumper_ghost.SetMode(TEXT);
  dumper_ghost.SetParallelContext(prank, psize);
  dumper_ghost.SetPoints(mesh.getNodes().values, dim, nb_nodes, "test-scotch-partition_ghost");
  dumper_ghost.SetConnectivity((int*) mesh.getGhostConnectivity(type).values,
   			 TRIANGLE1, nb_ghost_element, C_MODE);
  part = new double[nb_ghost_element];
  for (unsigned int i = 0; i < nb_ghost_element; ++i)
    part[i] = prank;
  dumper_ghost.AddElemDataField(part, 1, "partitions");
  dumper_ghost.SetPrefix("paraview/");
  dumper_ghost.Init();
  dumper_ghost.Dump();

  delete [] part;


#endif //AKANTU_USE_IOHELPER

  akantu::finalize();

  return EXIT_SUCCESS;
}
