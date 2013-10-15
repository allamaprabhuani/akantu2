/**
 * @file   solid_mechanics_model_cohesive_parallel.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Fri Jun 14 16:59:55 2013
 *
 * @brief  Functions for parallel cohesive elements
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
#include "solid_mechanics_model_cohesive.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initParallel(MeshPartition * partition,
                                               DataAccessor * data_accessor,
                                               bool extrinsic) {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel::initParallel(partition, data_accessor);

  /// create the distributed synchronizer for cohesive elements
  cohesive_distributed_synchronizer =
    new DistributedSynchronizer(mesh, "cohesive_distributed_synchronizer");

  synch_registry->registerSynchronizer(*cohesive_distributed_synchronizer,
  				       _gst_material_id);

  synch_registry->registerSynchronizer(*cohesive_distributed_synchronizer,
  				       _gst_smm_stress);

  synch_registry->registerSynchronizer(*cohesive_distributed_synchronizer,
  				       _gst_smm_boundary);

  DistributedSynchronizer & distributed_synchronizer =
    dynamic_cast<DistributedSynchronizer &>(*synch_parallel);

  distributed_synchronizer.filterElementsByKind(cohesive_distributed_synchronizer,
						_ek_cohesive);

  /// create the facet synchronizer for extrinsic simulations
  if (extrinsic) {
    rank_to_element = new ByElementTypeUInt("prank_to_element", id);
    ByElementTypeUInt & prank_to_element = *rank_to_element;

    distributed_synchronizer.buildPrankToElement(prank_to_element);

    MeshUtils::buildAllFacetsParallel(mesh, mesh_facets, prank_to_element);
    facet_generated = true;

    facet_synchronizer =
      FacetSynchronizer::createFacetSynchronizer(distributed_synchronizer,
                                                 mesh_facets);

    synch_registry->registerSynchronizer(*facet_synchronizer, _gst_smmc_facets);
    synch_registry->registerSynchronizer(*facet_synchronizer, _gst_smmc_facets_conn);

    synchronizeGhostFacets();

    facet_stress_synchronizer =
      FacetStressSynchronizer::createFacetStressSynchronizer(*facet_synchronizer,
							     mesh_facets);

    synch_registry->registerSynchronizer(*facet_stress_synchronizer,
					 _gst_smmc_facets_stress);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::fillGlobalConnectivity(ByElementTypeUInt & global_connectivity,
							 GhostType ghost_type,
							 bool just_init) {
  AKANTU_DEBUG_IN();

  Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1, ghost_type);
  Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1, ghost_type);

  for(; it != end; ++it) {
    ElementType type_facet = *it;

    Array<UInt> & connectivity = mesh_facets.getConnectivity(type_facet, ghost_type);
    Array<UInt> & g_connectivity = global_connectivity(type_facet, ghost_type);

    UInt nb_nodes_per_facet = connectivity.getNbComponent();
    UInt nb_facet = connectivity.getSize();

    g_connectivity.extendComponentsInterlaced(nb_nodes_per_facet, 1);
    g_connectivity.resize(nb_facet);

    if (just_init) continue;

    mesh_facets.getGlobalConnectivity(g_connectivity, type_facet, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::synchronizeGhostFacets() {
  AKANTU_DEBUG_IN();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();

  if (psize > 1) {

    /// correct facets' connectivity according to normals
    global_connectivity = new ByElementTypeUInt("global_connectivity", id);
    ByElementTypeUInt & g_connectivity = *global_connectivity;

    mesh_facets.initByElementTypeArray(g_connectivity, 1,
				       spatial_dimension - 1);

    /// copy facet normals for not ghost elements
    fillGlobalConnectivity(g_connectivity, _not_ghost, false);
    fillGlobalConnectivity(g_connectivity, _ghost, true);

    /// communicate
    synch_registry->synchronize(_gst_smmc_facets_conn);

    /// flip facets
    MeshUtils::flipFacets(mesh_facets, *global_connectivity, _ghost);

    delete global_connectivity;
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::updateFacetSynchronizer() {
  /// update synchronizer if needed
  if (facet_synchronizer != NULL) {
    DataAccessor * data_accessor = this;

    facet_synchronizer->updateDistributedSynchronizer(*cohesive_distributed_synchronizer,
						      *data_accessor,
                                                      cohesive_el_to_facet);

    for (ghost_type_t::iterator gt = ghost_type_t::begin();
         gt != ghost_type_t::end(); ++gt) {

      GhostType gt_facet = *gt;

      Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1, gt_facet);
      Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1, gt_facet);

      for(; it != end; ++it) {
        ElementType type_facet = *it;
        Array<UInt> & cohesive_el_to_f = cohesive_el_to_facet(type_facet, gt_facet);

        for (UInt f = 0; f < cohesive_el_to_f.getSize(); ++f)
          cohesive_el_to_f(f) = std::numeric_limits<UInt>::max();
      }
    }
  }
}

__END_AKANTU__
