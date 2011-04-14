/**
 * @file   test_synchronizer_communication.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 19 13:05:27 2010
 *
 * @brief  test to synchronize barycenters
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
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
#include "communication_buffer.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

class TestSynchronizer : public GhostSynchronizer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  TestSynchronizer(const Mesh & mesh);
  ~TestSynchronizer();

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(GhostBarycenter, ghost_barycenter, const Vector<Real>);

  /* ------------------------------------------------------------------------ */
  /* Ghost Synchronizer inherited members                                     */
  /* ------------------------------------------------------------------------ */
protected:
  virtual UInt getNbDataToPack(const Element & element,
				       GhostSynchronizationTag tag) const;
  virtual UInt getNbDataToUnpack(const Element & element,
					 GhostSynchronizationTag tag) const;
  virtual void packData(CommunicationBuffer & buffer,
			const Element & element,
			GhostSynchronizationTag tag) const;
  virtual void unpackData(CommunicationBuffer & buffer,
			  const Element & element,
			  GhostSynchronizationTag tag) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  std::string id;

  ByElementTypeReal ghost_barycenter;

  const Mesh & mesh;
};

/* -------------------------------------------------------------------------- */
/* TestSynchronizer implementation                                            */
/* -------------------------------------------------------------------------- */
TestSynchronizer::TestSynchronizer(const Mesh & mesh) : mesh(mesh) {
  UInt spatial_dimension = mesh.getSpatialDimension();

  id = "test_synchronizer";

  Mesh::ConnectivityTypeList::const_iterator it;
  const Mesh::ConnectivityTypeList & ghost_type_list = mesh.getConnectivityTypeList(_ghost);
  for(it = ghost_type_list.begin(); it != ghost_type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    UInt nb_ghost_element = mesh.getNbGhostElement(*it);
    ghost_barycenter[*it] = new Vector<Real>(nb_ghost_element,
							     spatial_dimension,
							     std::numeric_limits<Real>::quiet_NaN(),
							     "ghost_barycenter");

  }
}

TestSynchronizer::~TestSynchronizer() {
  UInt spatial_dimension = mesh.getSpatialDimension();

  Mesh::ConnectivityTypeList::const_iterator it;
  const Mesh::ConnectivityTypeList & ghost_type_list = mesh.getConnectivityTypeList(_ghost);
  for(it = ghost_type_list.begin(); it != ghost_type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;
    delete ghost_barycenter[*it];
  }
}

UInt TestSynchronizer::getNbDataToPack(const Element & element,
					       __attribute__ ((unused)) GhostSynchronizationTag tag) const {
  return Mesh::getSpatialDimension(element.type) * sizeof(Real);
}

UInt TestSynchronizer::getNbDataToUnpack(const Element & element,
						 __attribute__ ((unused)) GhostSynchronizationTag tag) const {
  return Mesh::getSpatialDimension(element.type) * sizeof(Real);
}

void TestSynchronizer::packData(CommunicationBuffer & buffer,
				const Element & element,
				__attribute__ ((unused)) GhostSynchronizationTag tag) const {
  UInt spatial_dimension = Mesh::getSpatialDimension(element.type);
  types::RVector bary(spatial_dimension);
  mesh.getBarycenter(element.element, element.type, bary.storage());

  buffer << bary;
}

void TestSynchronizer::unpackData(CommunicationBuffer & buffer,
				  const Element & element,
				  __attribute__ ((unused)) GhostSynchronizationTag tag) const {
  UInt spatial_dimension = Mesh::getSpatialDimension(element.type);
  Vector<Real>::iterator<types::RVector> bary = ghost_barycenter[element.type]->begin(spatial_dimension);
  buffer >> bary[element.element];
}

/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  initialize(&argc, &argv);

  int dim = 2;
  ElementType type = _triangle_3;

  Mesh mesh(dim);

  StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm->getNbProc();
  Int prank = comm->whoAmI();

  Communicator * communicator;
  if(prank == 0) {
    MeshIOMSH mesh_io;
    mesh_io.read("triangle.msh", mesh);
    MeshPartition * partition = new MeshPartitionScotch(mesh, dim);
    partition->partitionate(psize);
    communicator = Communicator::createCommunicatorDistributeMesh(mesh, partition);
    delete partition;
  } else {
    communicator = Communicator::createCommunicatorDistributeMesh(mesh, NULL);
  }

  comm->barrier();

  AKANTU_DEBUG_INFO("Creating TestSynchronizer");
  TestSynchronizer test_synchronizer(mesh);
  test_synchronizer.registerSynchronizer(*communicator);

  AKANTU_DEBUG_INFO("Registering tag");
  test_synchronizer.registerTag(_gst_test, "barycenter");

  AKANTU_DEBUG_INFO("Synchronizing tag");
  test_synchronizer.synchronize(_gst_test);


  Mesh::ConnectivityTypeList::const_iterator it;
  const Mesh::ConnectivityTypeList & ghost_type_list = mesh.getConnectivityTypeList(_ghost);
  for(it = ghost_type_list.begin(); it != ghost_type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != mesh.getSpatialDimension()) continue;

    UInt spatial_dimension = Mesh::getSpatialDimension(*it);
    UInt nb_ghost_element = mesh.getNbGhostElement(*it);

    for (UInt el = 0; el < nb_ghost_element; ++el) {
      Real barycenter[spatial_dimension];
      mesh.getBarycenter(el, *it, barycenter, _ghost);

      for (UInt i = 0; i < spatial_dimension; ++i) {
	if(fabs(barycenter[i] - test_synchronizer.getGhostBarycenter(*it).values[el * spatial_dimension + i])
	   > std::numeric_limits<Real>::epsilon())
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

  finalize();

  return EXIT_SUCCESS;
}
