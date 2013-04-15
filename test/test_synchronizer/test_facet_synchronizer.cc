/**
 * @file   test_facet_synchronizer.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Tue Apr  2 15:05:49 2013
 *
 * @brief  Facet synchronizer test
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

#include "test_data_accessor.hh"
#include "facet_synchronizer.hh"
#include "mesh_utils.hh"
#include "synchronizer_registry.hh"


/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  akantu::initialize(argc, argv);

  UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();
  DistributedSynchronizer * dist = NULL;

  /// partition the mesh
  if(prank == 0) {
    mesh.read("facet.msh");
    MeshPartition * partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
    dist = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, partition);
    delete partition;
  } else {
    dist = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, NULL);
  }

  /// create facets
  Mesh mesh_facets(spatial_dimension, mesh.getNodes(), "mesh_facets");

  ByElementTypeUInt prank_to_element("prank_to_element", "prank_to_el");
  dist->buildPrankToElement(prank_to_element);

  MeshUtils::buildAllFacetsParallel(mesh, mesh_facets, prank_to_element);

  debug::setDebugLevel(dblDump);
  std::cout << mesh << std::endl;
  debug::setDebugLevel(dblInfo);

  debug::setDebugLevel(dblDump);
  std::cout << mesh_facets << std::endl;
  debug::setDebugLevel(dblInfo);

  /// compute barycenter for each facet
  ByElementTypeReal barycenters("", "", 0);
  mesh_facets.initByElementTypeArray(barycenters, spatial_dimension, spatial_dimension - 1);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;
    Mesh::type_iterator it = mesh_facets.firstType(spatial_dimension - 1, ghost_type);
    Mesh::type_iterator last_type = mesh_facets.lastType(spatial_dimension - 1, ghost_type);

    for(; it != last_type; ++it) {
      UInt nb_element = mesh_facets.getNbElement(*it, ghost_type);
      Array<Real> & barycenter = barycenters(*it, ghost_type);
      barycenter.resize(nb_element);

      Array<Real>::iterator< Vector<Real> > bary_it = barycenter.begin(spatial_dimension);
      for (UInt elem = 0; elem < nb_element; ++elem, ++bary_it)
	mesh_facets.getBarycenter(elem, *it, bary_it->storage(), ghost_type);
    }
  }

  /// setup facet communications
  FacetSynchronizer * facet_synchronizer =
    FacetSynchronizer::createFacetSynchronizer(*dist, mesh_facets);

  AKANTU_DEBUG_INFO("Creating TestAccessor");
  TestAccessor test_accessor(mesh, barycenters);
  SynchronizerRegistry synch_registry(test_accessor);

  synch_registry.registerSynchronizer(*facet_synchronizer, _gst_test);

  /// synchronize facets and check results
  AKANTU_DEBUG_INFO("Synchronizing tag");
  synch_registry.synchronize(_gst_test);

  delete facet_synchronizer;
  delete dist;
  akantu::finalize();

  return EXIT_SUCCESS;
}
