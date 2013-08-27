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
    synch_registry->registerSynchronizer(*facet_synchronizer, _gst_smmc_normals);

    facet_stress_synchronizer =
      FacetStressSynchronizer::createFacetStressSynchronizer(*facet_synchronizer,
							     mesh_facets);

    synch_registry->registerSynchronizer(*facet_stress_synchronizer,
					 _gst_smmc_facets_stress);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::synchronizeCohesiveElements() {
  AKANTU_DEBUG_IN();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();

  if (psize > 1) {
    if (spatial_dimension == 3) AKANTU_DEBUG_TO_IMPLEMENT();

    /// correct facets' connectivity according to normals
    facet_normals = new ByElementTypeReal("facet_normals", id);
    ByElementTypeReal & f_normals = *facet_normals;

    mesh_facets.initByElementTypeArray(f_normals,
				       spatial_dimension,
				       spatial_dimension - 1,
				       false, _ek_regular, true);

    /// compute facet normals for not ghost elements
    MeshUtils::computeFacetNormals(mesh_facets, f_normals);

    /// communicate
    synch_registry->synchronize(_gst_smmc_normals);

    /// compute local normals for ghost facets
    ByElementTypeReal ghost_normals("ghost_normals", id);

    mesh_facets.initByElementTypeArray(ghost_normals,
				       spatial_dimension,
				       spatial_dimension - 1);

    MeshUtils::computeFacetNormals(mesh_facets, ghost_normals, _ghost);

    GhostType gt_facet = _ghost;

    Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1, gt_facet);
    Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1, gt_facet);

    /// loop on every ghost facet
    for(; it != end; ++it) {
      ElementType type_facet = *it;

      Array<UInt> & connectivity = mesh_facets.getConnectivity(type_facet, gt_facet);
      Array<std::vector<Element> > & el_to_f =
	mesh_facets.getElementToSubelement(type_facet, gt_facet);
      Array<Element> & subfacet_to_facet =
	mesh_facets.getSubelementToElement(type_facet, gt_facet);

      UInt nb_nodes_per_facet = connectivity.getNbComponent();
      UInt nb_subfacet_per_facet = subfacet_to_facet.getNbComponent();
      UInt nb_facet = connectivity.getSize();
      Array<Real> & recv_normals = f_normals(type_facet, gt_facet);
      Array<Real> & computed_normals = ghost_normals(type_facet, gt_facet);

      Array<Real>::iterator<Vector<Real> > recv_it =
	recv_normals.begin(spatial_dimension);
      Array<Real>::iterator<Vector<Real> > computed_it =
	computed_normals.begin(spatial_dimension);
      Array<UInt>::iterator<Vector<UInt> > conn_it =
	connectivity.begin(nb_nodes_per_facet);
      Array<Element>::iterator<Vector<Element> > subf_to_f =
	subfacet_to_facet.begin(nb_subfacet_per_facet);

      Vector<UInt> conn_tmp(nb_nodes_per_facet);
      std::vector<Element> el_tmp(2);
      Vector<Element> subf_tmp(nb_subfacet_per_facet);

      for (UInt f = 0; f < nb_facet; ++f, ++recv_it, ++computed_it, ++conn_it,
	     ++subf_to_f) {
	Real product = recv_it->dot( (*computed_it) );

	/// if product is negative, facets must be flipped
	if (product < 0) {
	  conn_tmp = (*conn_it);
	  (*conn_it)(0) = conn_tmp(1);
	  (*conn_it)(1) = conn_tmp(0);

	  el_tmp = el_to_f(f);
	  el_to_f(f)[0] = el_tmp[1];
	  el_to_f(f)[1] = el_tmp[0];

	  subf_tmp = (*subf_to_f);
	  (*subf_to_f)(0) = subf_tmp(1);
	  (*subf_to_f)(1) = subf_tmp(0);
	}
      }
    }

    delete facet_normals;
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
