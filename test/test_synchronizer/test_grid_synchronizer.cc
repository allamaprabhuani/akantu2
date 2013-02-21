/**
 * @file   test_grid_synchronizer.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Nov 25 17:00:17 2011
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
#include "aka_grid_dynamic.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "grid_synchronizer.hh"
#include "mesh_partition.hh"
#include "synchronizer_registry.hh"
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

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Barycenter, barycenter, Real);

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


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  std::string id;

  ByElementTypeReal barycenter;

  const Mesh & mesh;
};


/* -------------------------------------------------------------------------- */
/* TestSynchronizer implementation                                            */
/* -------------------------------------------------------------------------- */
TestAccessor::TestAccessor(const Mesh & mesh) : barycenter("barycenter", id), mesh(mesh) {
  UInt spatial_dimension = mesh.getSpatialDimension();

  id = "test_synchronizer";

  Mesh::ConnectivityTypeList::const_iterator it;
  try {
    const Mesh::ConnectivityTypeList & ghost_type_list = mesh.getConnectivityTypeList(_ghost);
    for(it = ghost_type_list.begin(); it != ghost_type_list.end(); ++it) {
      if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

      UInt nb_ghost_element = mesh.getNbElement(*it,_ghost);
      barycenter.alloc(nb_ghost_element, spatial_dimension, *it);
    }
  } catch (debug::Exception & e) { std::cout << e << std::endl; }
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
  UInt spatial_dimension = mesh.getSpatialDimension();
  Vector<Element>::const_iterator<Element> bit  = elements.begin();
  Vector<Element>::const_iterator<Element> bend = elements.end();
  for (; bit != bend; ++bit) {
    const Element & element = *bit;

    types::RVector bary(barycenter(element.type,
				   element.ghost_type).storage() + element.element * spatial_dimension,
			spatial_dimension);;
    buffer << barycenter;
  }
}

void TestAccessor::unpackElementData(CommunicationBuffer & buffer,
				     const Vector<Element> & elements,
				     __attribute__ ((unused)) SynchronizationTag tag) {
  UInt spatial_dimension = mesh.getSpatialDimension();
  Vector<Element>::const_iterator<Element> bit  = elements.begin();
  Vector<Element>::const_iterator<Element> bend = elements.end();
  for (; bit != bend; ++bit) {
    const Element & element = *bit;

    types::RVector barycenter_loc(spatial_dimension);
    mesh.getBarycenter(element.element, element.type, barycenter_loc.storage(), element.ghost_type);

    types::RVector barycenter(spatial_dimension);
    buffer >> barycenter;
    Real tolerance = 1e-15;
    for (UInt i = 0; i < spatial_dimension; ++i) {
      if(!(std::abs(barycenter(i) - barycenter_loc(i)) <= tolerance))
	AKANTU_DEBUG_ERROR("Unpacking an unknown value for the element: "
			   << element
			   << "(barycenter[" << i << "] = " << barycenter_loc(i)
			   << " and buffer[" << i << "] = " << barycenter(i) << ") - tag: " << tag);
    }
  }
}


/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
  akantu::initialize(argc, argv);

  UInt spatial_dimension = 2;
  ElementType type = _triangle_3;

  Mesh mesh(spatial_dimension);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  if(prank == 0) {
    MeshIOMSH mesh_io;
    mesh_io.read("triangle.msh", mesh);
    MeshPartition * partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
    DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, partition);
    delete partition;
  } else {
    DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, NULL);
  }

#ifdef AKANTU_USE_IOHELPER
  double * part;
  unsigned int nb_nodes = mesh.getNbNodes();
  unsigned int nb_element = mesh.getNbElement(type);

  iohelper::DumperParaview dumper;
  dumper.setMode(iohelper::TEXT);
  dumper.setParallelContext(prank, psize);
  dumper.setPoints(mesh.getNodes().values, spatial_dimension, nb_nodes, "test-grid-synchronizer");
  dumper.setConnectivity((int*) mesh.getConnectivity(type).values,
   			 iohelper::TRIANGLE1, nb_element, iohelper::C_MODE);
  part = new double[nb_element];
  for (unsigned int i = 0; i < nb_element; ++i)
    part[i] = prank;
  dumper.addElemDataField("partitions", part, 1, nb_element);
  dumper.setPrefix("paraview/");
  dumper.init();
  dumper.dump();
  delete [] part;

  iohelper::DumperParaview dumper_ghost;
  dumper_ghost.setMode(iohelper::TEXT);
  dumper_ghost.setParallelContext(prank, psize);

  try {
    unsigned int nb_ghost_element = mesh.getNbElement(type,_ghost);
    dumper_ghost.setPoints(mesh.getNodes().values, spatial_dimension, nb_nodes, "test-grid-synchronizer-ghost");
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
  } catch (debug::Exception & e) { std::cout << e << std::endl; }
#endif


  comm.barrier();

  mesh.computeBoundingBox();

  Real lower_bounds[spatial_dimension];
  Real upper_bounds[spatial_dimension];

  mesh.getLowerBounds(lower_bounds);
  mesh.getUpperBounds(upper_bounds);

  types::Vector<Real> spacing(spatial_dimension);
  types::Vector<Real> center(spatial_dimension);

  for (UInt i = 0; i < spatial_dimension; ++i) {
    spacing[i] = (upper_bounds[i] - lower_bounds[i]) / 10.;
    center[i] = (upper_bounds[i] + lower_bounds[i]) / 2.;
  }

  SpacingGrid<Element> grid(spatial_dimension, spacing, center);

  GhostType ghost_type = _not_ghost;

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);

  ByElementTypeReal barycenters("", "", 0);
  mesh.initByElementTypeVector(barycenters, spatial_dimension, 0);

  Element e;
  e.ghost_type = ghost_type;

  it = mesh.firstType(spatial_dimension, ghost_type);
  for(; it != last_type; ++it) {
    UInt nb_element = mesh.getNbElement(*it);
    e.type = *it;
    Vector<Real> & barycenter = barycenters(*it);
    barycenter.resize(nb_element);

    Vector<Real>::iterator<types::RVector> bary_it = barycenter.begin(spatial_dimension);
    for (UInt elem = 0; elem < nb_element; ++elem) {
      mesh.getBarycenter(elem, *it, bary_it->storage());
      e.element = elem;
      grid.insert(e, *bary_it);
      ++bary_it;
    }
  }

  MeshIOMSH mesh_io;

  std::stringstream sstr; sstr << "mesh_" << prank << ".msh";
  mesh_io.write(sstr.str(), mesh);

  Mesh grid_mesh(spatial_dimension, "grid_mesh", 0);
  std::stringstream sstr_grid; sstr_grid << "grid_mesh_" << prank << ".msh";
  grid.saveAsMesh(grid_mesh);
  mesh_io.write(sstr_grid.str(), grid_mesh);

  GridSynchronizer * grid_communicator = GridSynchronizer::createGridSynchronizer(mesh, grid);

  AKANTU_DEBUG_INFO("Creating TestAccessor");
  TestAccessor test_accessor(mesh);
  SynchronizerRegistry synch_registry(test_accessor);

  synch_registry.registerSynchronizer(*grid_communicator, _gst_test);

  AKANTU_DEBUG_INFO("Synchronizing tag");
  synch_registry.synchronize(_gst_test);

#ifdef AKANTU_USE_IOHELPER
  try {
    UInt nb_ghost_element = mesh.getNbElement(type, _ghost);
    UInt nb_nodes = mesh.getNbNodes();

    iohelper::DumperParaview dumper_ghost;
    dumper_ghost.setMode(iohelper::TEXT);
    dumper_ghost.setParallelContext(prank, psize);

    dumper_ghost.setMode(iohelper::TEXT);
    dumper_ghost.setParallelContext(prank, psize);
    dumper_ghost.setPoints(mesh.getNodes().values, spatial_dimension, nb_nodes, "test-grid-synchronizer-ghost");
    dumper_ghost.setConnectivity((int*) mesh.getConnectivity(type,_ghost).values,
				 iohelper::TRIANGLE1, nb_ghost_element, iohelper::C_MODE);
    part = new double[nb_ghost_element];
    for (unsigned int i = 0; i < nb_ghost_element; ++i)
      part[i] = prank;
    dumper_ghost.addElemDataField("partitions", part, 1, nb_ghost_element);
    dumper_ghost.dump();
    delete [] part;

    unsigned int nb_grid_element = grid_mesh.getNbElement(_quadrangle_4);
    iohelper::DumperParaview dumper_grid;
    dumper_grid.setMode(iohelper::TEXT);
    dumper_grid.setParallelContext(prank, psize);
    dumper_grid.setPoints(grid_mesh.getNodes().values, spatial_dimension,
			  grid_mesh.getNbNodes(), "test-grid-synchronizer-grid");
    dumper_grid.setConnectivity((int*) grid_mesh.getConnectivity(_quadrangle_4).values,
				iohelper::QUAD1, nb_grid_element, iohelper::C_MODE);

    part = new double[nb_grid_element];
    std::fill_n(part, nb_grid_element, prank);

    dumper_grid.addElemDataField("partitions", part, 1, nb_grid_element);
    dumper_grid.setPrefix("paraview/");
    dumper_grid.init();
    dumper_grid.dump();
    delete [] part;
  } catch(debug::Exception & e) { std::cout << e << std::endl; }
#endif //AKANTU_USE_IOHELPER


  akantu::finalize();

  return EXIT_SUCCESS;
}
