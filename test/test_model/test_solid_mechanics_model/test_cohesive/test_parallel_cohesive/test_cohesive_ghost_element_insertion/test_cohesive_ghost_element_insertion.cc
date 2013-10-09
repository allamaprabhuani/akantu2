/**
 * @file   test_cohesive_ghost_element_insertion.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Thu Apr  4 16:59:43 2013
 *
 * @brief  test for ghost cohesive elements insertion
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

#include "mesh_utils.hh"
#include "mesh_partition.hh"
#include "solid_mechanics_model_cohesive.hh"

/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  initialize(argc, argv);
  debug::setDebugLevel(dblInfo);

  UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  akantu::MeshPartition * partition = NULL;
  if(prank == 0) {
    // Read the mesh
    mesh.read("mesh.msh");

    /// partition the mesh
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    debug::setDebugLevel(dblDump);
    partition->partitionate(psize);
    debug::setDebugLevel(dblInfo);
  }

  SolidMechanicsModelCohesive model(mesh);
  model.initParallel(partition, NULL, true);
  model.initFull("material.dat", _explicit_lumped_mass, true);

  Mesh & mesh_facets = const_cast<Mesh&>(model.getMeshFacets());

  ByElementTypeArray<bool> facet_insertion("facet_insertion", "");
  mesh_facets.initByElementTypeArray(facet_insertion, 1, spatial_dimension - 1);

  if(prank == 0) {
    for (ghost_type_t::iterator gt = ghost_type_t::begin();
	 gt != ghost_type_t::end(); ++gt) {

      GhostType type_ghost = *gt;

      Mesh::type_iterator it = mesh_facets.firstType(spatial_dimension - 1, type_ghost);
      Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1, type_ghost);

      for (; it != end; ++it) {
	UInt nb_facet = mesh_facets.getNbElement(*it, type_ghost);
	if (nb_facet != 0) {
	  Array<bool> & facet_ins = facet_insertion(*it, type_ghost);
	  facet_ins.resize(nb_facet);
	  facet_ins.clear();

	  Real * bary_facet = new Real[spatial_dimension];
	  for (UInt f = 0; f < nb_facet; ++f) {
	    mesh_facets.getBarycenter(f, *it, bary_facet, type_ghost);
	    if ((bary_facet[0] > -0.26 && bary_facet[0] < 0.26) &&
		(bary_facet[1] < 0.51)) facet_ins(f) = true;
	  }
	  delete[] bary_facet;

	}
      }
    }
  }

  if(prank == 0) {
    debug::setDebugLevel(dblDump);
    std::cout << mesh << std::endl;
    std::cout << mesh_facets << std::endl;
    debug::setDebugLevel(dblInfo);
  }

  MeshUtils::insertCohesiveElements(mesh, mesh_facets, facet_insertion, true);

  if(prank == 0) {
    debug::setDebugLevel(dblDump);
    std::cout << mesh << std::endl;
    std::cout << mesh_facets << std::endl;
    debug::setDebugLevel(dblInfo);
  }

  finalize();
  if(prank == 0) std::cout << "OK: Test passed!" << std::endl;
  return EXIT_SUCCESS;
}
