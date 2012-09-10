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
    barycenter(element.type).begin(spatial_dimension);
  buffer >> bary[element.element];
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
  dumper.SetMode(iohelper::TEXT);
  dumper.SetParallelContext(prank, psize);
  dumper.SetPoints(mesh.getNodes().values, spatial_dimension, nb_nodes, "test-grid-synchronizer");
  dumper.SetConnectivity((int*) mesh.getConnectivity(type).values,
   			 iohelper::TRIANGLE1, nb_element, iohelper::C_MODE);
  part = new double[nb_element];
  for (unsigned int i = 0; i < nb_element; ++i)
    part[i] = prank;
  dumper.AddElemDataField(part, 1, "partitions");
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
  delete [] part;

  iohelper::DumperParaview dumper_ghost;
  dumper_ghost.SetMode(iohelper::TEXT);
  dumper_ghost.SetParallelContext(prank, psize);

  try {
    unsigned int nb_ghost_element = mesh.getNbElement(type,_ghost);
    dumper_ghost.SetPoints(mesh.getNodes().values, spatial_dimension, nb_nodes, "test-grid-synchronizer-ghost");
    dumper_ghost.SetConnectivity((int*) mesh.getConnectivity(type,_ghost).values,
				 iohelper::TRIANGLE1, nb_ghost_element, iohelper::C_MODE);
    part = new double[nb_ghost_element];
    for (unsigned int i = 0; i < nb_ghost_element; ++i)
      part[i] = prank;
    dumper_ghost.AddElemDataField(part, 1, "partitions");
    dumper_ghost.SetPrefix("paraview/");
    dumper_ghost.Init();
    dumper_ghost.Dump();
    delete [] part;
  } catch (debug::Exception & e) { std::cout << e << std::endl; }
#endif


  comm.barrier();

  mesh.computeBoundingBox();

  Real lower_bounds[spatial_dimension];
  Real upper_bounds[spatial_dimension];

  mesh.getLowerBounds(lower_bounds);
  mesh.getUpperBounds(upper_bounds);

  Real spacing[spatial_dimension];

  for (UInt i = 0; i < spatial_dimension; ++i) {
    spacing[i] = (upper_bounds[i] - lower_bounds[i]) / 10.;
  }

  Real local_lower_bounds[spatial_dimension];
  Real local_upper_bounds[spatial_dimension];

  mesh.getLocalLowerBounds(local_lower_bounds);
  mesh.getLocalUpperBounds(local_upper_bounds);

  RegularGrid<Element> grid(spatial_dimension, local_lower_bounds, local_upper_bounds, spacing);

  GhostType ghost_type = _not_ghost;

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);

  std::cout << prank <<  " TTTTTTTOOOOOOOOOOOOOTTTTTTTTTTTTTTTTOOOOOOOOOOOOO" << grid <<std::endl;

  ByElementTypeReal barycenters("", "", 0);
  mesh.initByElementTypeVector(barycenters, spatial_dimension, 0);
  // first generate the quad points coordinate and count the number of points per cell
  for(; it != last_type; ++it) {
    UInt nb_element = mesh.getNbElement(*it);
    Vector<Real> & barycenter = barycenters(*it);
    barycenter.resize(nb_element);

    //     Vector<Real>::iterator<types::RVector> bary_it = barycenter.begin(spatial_dimension);
    // for (UInt elem = 0; elem < nb_element; ++elem) {
    //   
    //   grid.count(*bary_it);
    //   ++bary_it;
    // }
  }

  Element e;
  e.ghost_type = ghost_type;

  std::cout << prank <<  " TTTTTTTOOOOOOOOOOOOOTTTTTTTTTTTTTTTTOOOOOOOOOOOOO" << std::endl;

  // second insert the point in the cells
  //  grid.beginInsertions();
  it = mesh.firstType(spatial_dimension, ghost_type);
  for(; it != last_type; ++it) {
    UInt nb_element = mesh.getNbElement(*it);
    e.type = *it;
    Vector<Real> & barycenter = barycenters(*it);
    Vector<Real>::iterator<types::RVector> bary_it = barycenter.begin(spatial_dimension);
    for (UInt elem = 0; elem < nb_element; ++elem) {
      mesh.getBarycenter(elem, *it, bary_it->storage());
      e.element = elem;
      grid.insert(e, *bary_it);
      ++bary_it;
    }
  }
  //  grid.endInsertions();

  MeshIOMSH mesh_io;

  std::stringstream sstr; sstr << "mesh_" << prank << ".msh";
  mesh_io.write(sstr.str(), mesh);

  Mesh grid_mesh(spatial_dimension, "grid_mesh", 0);
  std::stringstream sstr_grid; sstr_grid << "grid_mesh_" << prank << ".msh";
  grid.saveAsMesh(grid_mesh);
  mesh_io.write(sstr_grid.str(), grid_mesh);

  std::cout << "TTTTTTTOOOOOOOOOOOOOTTTTTTTTTTTTTTTTOOOOOOOOOOOOO" << std::endl;
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
    dumper_ghost.SetMode(iohelper::TEXT);
    dumper_ghost.SetParallelContext(prank, psize);

    dumper_ghost.SetMode(iohelper::TEXT);
    dumper_ghost.SetParallelContext(prank, psize);
    dumper_ghost.SetPoints(mesh.getNodes().values, spatial_dimension, nb_nodes, "test-grid-synchronizer-ghost");
    dumper_ghost.SetConnectivity((int*) mesh.getConnectivity(type,_ghost).values,
				 iohelper::TRIANGLE1, nb_ghost_element, iohelper::C_MODE);
    part = new double[nb_ghost_element];
    for (unsigned int i = 0; i < nb_ghost_element; ++i)
      part[i] = prank;
    dumper_ghost.AddElemDataField(part, 1, "partitions");
    dumper_ghost.Dump();
    delete [] part;

    unsigned int nb_grid_element = grid_mesh.getNbElement(_quadrangle_4);
    unsigned int nb_quadrature_points = 4;
    iohelper::DumperParaview dumper_grid;
    dumper_grid.SetMode(iohelper::TEXT);
    dumper_grid.SetParallelContext(prank, psize);
    dumper_grid.SetPoints(grid_mesh.getNodes().values, spatial_dimension,
			  grid_mesh.getNbNodes(), "test-grid-synchronizer-grid");
    dumper_grid.SetConnectivity((int*) grid_mesh.getConnectivity(_quadrangle_4).values,
				iohelper::QUAD1, nb_grid_element, iohelper::C_MODE);

    part = new double[nb_grid_element * nb_quadrature_points];
    std::fill_n(part, nb_grid_element * nb_quadrature_points, prank);

    dumper_grid.AddElemDataField(part, 1, "partitions");
    dumper_grid.SetPrefix("paraview/");
    dumper_grid.Init();
    dumper_grid.Dump();
    delete [] part;
  } catch(debug::Exception & e) { std::cout << e << std::endl; }
#endif //AKANTU_USE_IOHELPER


  akantu::finalize();

  return EXIT_SUCCESS;
}
