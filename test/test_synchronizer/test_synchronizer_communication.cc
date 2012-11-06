/**
 * @file   test_synchronizer_communication.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Sun Sep 12 23:37:43 2010
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
#include "distributed_synchronizer.hh"
#include "synchronizer_registry.hh"
#include "mesh.hh"
#include "mesh_io_msh.hh"
#include "mesh_partition_scotch.hh"
#include "communication_buffer.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"

#endif //AKANTU_USE_IOHELPER

using namespace akantu;

class TestAccessor : public DataAccessor {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  TestAccessor(const Mesh & mesh);
  ~TestAccessor();

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(GhostBarycenter, ghost_barycenter, Real);

  /* ------------------------------------------------------------------------ */
  /* Ghost Synchronizer inherited members                                     */
  /* ------------------------------------------------------------------------ */
protected:
  virtual UInt getNbDataForElements(const Vector<Element> & elements,
				    SynchronizationTag tag) const;
  virtual void packElementData(CommunicationBuffer & buffer,
			       const Vector<Element> & elements,
			       SynchronizationTag tag) const;
  virtual void unpackElementData(CommunicationBuffer & buffer,
				 const Vector<Element> & elements,
				 SynchronizationTag tag);

  virtual UInt getNbDataToPack(SynchronizationTag tag) const;
  virtual UInt getNbDataToUnpack(SynchronizationTag tag) const;
  virtual void packData(CommunicationBuffer & buffer,
			const UInt index,
			SynchronizationTag tag) const;
  virtual void unpackData(CommunicationBuffer & buffer,
			  const UInt index,
			  SynchronizationTag tag);


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
TestAccessor::TestAccessor(const Mesh & mesh) : ghost_barycenter("ghost_barycenter", id), mesh(mesh) {
  UInt spatial_dimension = mesh.getSpatialDimension();

  id = "test_synchronizer";

  Mesh::ConnectivityTypeList::const_iterator it;
  const Mesh::ConnectivityTypeList & ghost_type_list = mesh.getConnectivityTypeList(_ghost);
  for(it = ghost_type_list.begin(); it != ghost_type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    UInt nb_ghost_element = mesh.getNbElement(*it,_ghost);
    ghost_barycenter.alloc(nb_ghost_element, spatial_dimension, *it);

  }
}

TestAccessor::~TestAccessor() {

}

UInt TestAccessor::getNbDataForElements(const Vector<Element> & elements,
					__attribute__ ((unused)) SynchronizationTag tag) const {
  return Mesh::getSpatialDimension(elements(0).type) * sizeof(Real) * elements.getSize();
}

void TestAccessor::packElementData(CommunicationBuffer & buffer,
				   const Vector<Element> & elements,
				   __attribute__ ((unused)) SynchronizationTag tag) const {
  Vector<Element>::const_iterator<Element> bit  = elements.begin();
  Vector<Element>::const_iterator<Element> bend = elements.end();
  for (; bit != bend; ++bit) {
    const Element & element = *bit;
    UInt spatial_dimension = Mesh::getSpatialDimension(element.type);
    types::RVector bary(spatial_dimension);
    mesh.getBarycenter(element.element, element.type, bary.storage());

    buffer << bary;
  }
}

void TestAccessor::unpackElementData(CommunicationBuffer & buffer,
				     const Vector<Element> & elements,
				     __attribute__ ((unused)) SynchronizationTag tag) {
  Vector<Element>::const_iterator<Element> bit  = elements.begin();
  Vector<Element>::const_iterator<Element> bend = elements.end();
  for (; bit != bend; ++bit) {
    const Element & element = *bit;
    UInt spatial_dimension = Mesh::getSpatialDimension(element.type);
    Vector<Real>::iterator<types::RVector> bary =
      ghost_barycenter(element.type).begin(spatial_dimension);
    buffer >> bary[element.element];
  }
}

UInt TestAccessor::getNbDataToPack(__attribute__ ((unused)) SynchronizationTag tag) const {
  return 0;
}

UInt TestAccessor::getNbDataToUnpack(__attribute__ ((unused)) SynchronizationTag tag) const {
  return 0;
}

void TestAccessor::packData(__attribute__ ((unused)) CommunicationBuffer & buffer,
			    __attribute__ ((unused)) const UInt index,
			    __attribute__ ((unused)) SynchronizationTag tag) const {
}

void TestAccessor::unpackData(__attribute__ ((unused)) CommunicationBuffer & buffer,
			      __attribute__ ((unused)) const UInt index,
			      __attribute__ ((unused)) SynchronizationTag tag) {
}


/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  initialize(argc, argv);

  int dim = 2;
  ElementType type = _triangle_3;

  Mesh mesh(dim);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  DistributedSynchronizer * communicator;
  if(prank == 0) {
    MeshIOMSH mesh_io;
    mesh_io.read("triangle.msh", mesh);
    MeshPartition * partition = new MeshPartitionScotch(mesh, dim);
    partition->partitionate(psize);
    communicator = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, partition);
    delete partition;
  } else {
    communicator = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, NULL);
  }

  comm.barrier();

  AKANTU_DEBUG_INFO("Creating TestAccessor");
  TestAccessor test_accessor(mesh);
  SynchronizerRegistry synch_registry(test_accessor);
  synch_registry.registerSynchronizer(*communicator,_gst_test);

  AKANTU_DEBUG_INFO("Synchronizing tag");
  synch_registry.synchronize(_gst_test);


  Mesh::ConnectivityTypeList::const_iterator it;
  const Mesh::ConnectivityTypeList & ghost_type_list = mesh.getConnectivityTypeList(_ghost);
  for(it = ghost_type_list.begin(); it != ghost_type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != mesh.getSpatialDimension()) continue;

    UInt spatial_dimension = Mesh::getSpatialDimension(*it);
    UInt nb_ghost_element = mesh.getNbElement(*it,_ghost);

    for (UInt el = 0; el < nb_ghost_element; ++el) {
      Real barycenter[spatial_dimension];
      mesh.getBarycenter(el, *it, barycenter, _ghost);

      for (UInt i = 0; i < spatial_dimension; ++i) {
	if(fabs(barycenter[i] - test_accessor.getGhostBarycenter(*it).values[el * spatial_dimension + i])
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

  iohelper::DumperParaview dumper;
  dumper.setMode(iohelper::TEXT);
  dumper.setParallelContext(prank, psize);
  dumper.setPoints(mesh.getNodes().values, dim, nb_nodes, "test-scotch-partition");
  dumper.setConnectivity((int*) mesh.getConnectivity(type).values,
   			 iohelper::TRIANGLE1, nb_element, iohelper::C_MODE);
  double * part = new double[nb_element];
  for (unsigned int i = 0; i < nb_element; ++i)
    part[i] = prank;
  dumper.addElemDataField("partitions", part, 1, nb_element);
  dumper.setPrefix("paraview/");
  dumper.init();
  dumper.dump();

  delete [] part;

  unsigned int nb_ghost_element = mesh.getNbElement(type,_ghost);
  iohelper::DumperParaview dumper_ghost;
  dumper_ghost.setMode(iohelper::TEXT);
  dumper_ghost.setParallelContext(prank, psize);
  dumper_ghost.setPoints(mesh.getNodes().values, dim, nb_nodes, "test-scotch-partition_ghost");
  dumper_ghost.setConnectivity((int*) mesh.getConnectivity(type,_ghost).values,
   			 iohelper::TRIANGLE1, nb_ghost_element, iohelper::C_MODE);
  part = new double[nb_ghost_element];
  for (unsigned int i = 0; i < nb_ghost_element; ++i)
    part[i] = prank;

  dumper_ghost.addElemDataField("partitions", part, 1, nb_ghost_element);
  dumper_ghost.setPrefix("paraview/");
  dumper_ghost.init();
  dumper_ghost.dump();

  delete [] part;


#endif //AKANTU_USE_IOHELPER

  finalize();

  return EXIT_SUCCESS;
}
