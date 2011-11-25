/**
 * @file   test_grid_synchronizer.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Nov  7 11:58:02 2011
 *
 * @brief  test the GridSynchronizer object
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
#include "aka_common.hh"
#include "aka_grid.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "grid_synchronizer.hh"
#include "mesh_partition.hh"
#include "synchronizer_registry.hh"

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
  virtual UInt getNbDataToPack(const Element & element,
				       SynchronizationTag tag) const;
  virtual UInt getNbDataToUnpack(const Element & element,
					 SynchronizationTag tag) const;
  virtual void packData(CommunicationBuffer & buffer,
			const Element & element,
			SynchronizationTag tag) const;
  virtual void unpackData(CommunicationBuffer & buffer,
			  const Element & element,
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

UInt TestAccessor::getNbDataToPack(const Element & element,
					       __attribute__ ((unused)) SynchronizationTag tag) const {
  return Mesh::getSpatialDimension(element.type) * sizeof(Real);
}

UInt TestAccessor::getNbDataToUnpack(const Element & element,
						 __attribute__ ((unused)) SynchronizationTag tag) const {
  return Mesh::getSpatialDimension(element.type) * sizeof(Real);
}

void TestAccessor::packData(CommunicationBuffer & buffer,
				const Element & element,
				__attribute__ ((unused)) SynchronizationTag tag) const {
  UInt spatial_dimension = Mesh::getSpatialDimension(element.type);
  types::RVector bary(spatial_dimension);
  mesh.getBarycenter(element.element, element.type, bary.storage());

  buffer << bary;
}

void TestAccessor::unpackData(CommunicationBuffer & buffer,
				  const Element & element,
				  __attribute__ ((unused)) SynchronizationTag tag) {
  UInt spatial_dimension = Mesh::getSpatialDimension(element.type);
  Vector<Real>::iterator<types::RVector> bary =
    ghost_barycenter(element.type).begin(spatial_dimension);
  buffer >> bary[element.element];
}

int main(int argc, char *argv[]) {
  akantu::initialize(&argc, &argv);

  UInt spatial_dimension = 2;
  //  ElementType type = _triangle_3;

  Mesh mesh(spatial_dimension);

  StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm->getNbProc();
  Int prank = comm->whoAmI();

  DistributedSynchronizer * communicator;
  if(prank == 0) {
    MeshIOMSH mesh_io;
    mesh_io.read("triangle.msh", mesh);
    MeshPartition * partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
    communicator = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, partition);
    delete partition;
  } else {
    communicator = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, NULL);
  }

  comm->barrier();

  AKANTU_DEBUG_INFO("Creating TestAccessor");
  TestAccessor test_accessor(mesh);
  SynchronizerRegistry synch_registry(test_accessor);

  mesh.computeBoundingBox();

  Real lower_bounds[spatial_dimension];
  Real upper_bounds[spatial_dimension];

  mesh.getLocalLowerBounds(lower_bounds);
  mesh.getLocalUpperBounds(upper_bounds);

  Real spacing[spatial_dimension];

  for (UInt i = 0; i < spatial_dimension; ++i) {
    spacing[i] = (upper_bounds[i] - lower_bounds[i]) / 100.;
  }


  RegularGrid<Element> grid(spatial_dimension, lower_bounds, upper_bounds, spacing);


  GhostType ghost_type = _not_ghost;

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);

  ByElementTypeReal barycenters("", "", 0);
  mesh.initByElementTypeVector(barycenters, spatial_dimension, 0);
  // first generate the quad points coordinate and count the number of points per cell
  for(; it != last_type; ++it) {
    UInt nb_element = mesh.getNbElement(*it);
    Vector<Real> & barycenter = barycenters(*it);
    barycenter.resize(nb_element);
    Vector<Real>::iterator<types::RVector> bary_it = barycenter.begin(spatial_dimension);
    for (UInt elem = 0; elem < nb_element; ++elem) {
      mesh.getBarycenter(elem, *it, bary_it->storage());
      grid.count(*bary_it);
      ++bary_it;
    }
  }


  Element e;
  e.ghost_type = ghost_type;

  // second insert the point in the cells
  grid.beginInsertions();
  it = mesh.firstType(spatial_dimension, ghost_type);
  for(; it != last_type; ++it) {
    UInt nb_element = mesh.getNbElement(*it);
    e.type = *it;
    Vector<Real> & barycenter = barycenters(*it);
    Vector<Real>::iterator<types::RVector> bary_it = barycenter.begin(spatial_dimension);
    for (UInt elem = 0; elem < nb_element; ++elem) {
      e.element = elem;
      grid.insert(e, *bary_it);
      ++bary_it;
    }
  }
  grid.endInsertions();

  GridSynchronizer * grid_communicator = GridSynchronizer::createGridSynchronizer(mesh, grid);

  synch_registry.registerSynchronizer(*grid_communicator,_gst_test);

  AKANTU_DEBUG_INFO("Synchronizing tag");
  synch_registry.synchronize(_gst_test);




  akantu::finalize();

  return EXIT_SUCCESS;
}
