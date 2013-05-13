/**
 * @file   solid_mechanics_model_cohesive.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Tue May 08 13:01:18 2012
 *
 * @brief  Solid mechanics model for cohesive elements
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

#include "shape_cohesive.hh"
#include "solid_mechanics_model_cohesive.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */

SolidMechanicsModelCohesive::SolidMechanicsModelCohesive(Mesh & mesh,
                                                         UInt dim,
                                                         const ID & id,
                                                         const MemoryID & memory_id)
  : SolidMechanicsModel(mesh, dim, id, memory_id),
    mesh_facets(mesh.getSpatialDimension(),
                mesh.getNodes().getID(),
                "mesh_facets", memory_id),
    tangents("tangents", id),
    stress_on_facet("stress_on_facet", id),
    facet_stress("facet_stress", id),
    fragment_to_element("fragment_to_element", id),
    fragment_velocity(0, spatial_dimension, "fragment_velocity"),
    fragment_center(0, spatial_dimension, "fragment_center"),
    facet_insertion("facet_insertion", id),
    doubled_facets("doubled_facets", id),
    facets_to_cohesive_el("facets_to_cohesive_el", id),
    cohesive_el_to_facet("cohesive_el_to_facet", id),
    facet_material("facet_material", id) {
  AKANTU_DEBUG_IN();

  facet_generated = false;
  facet_synchronizer = NULL;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::updateResidual(bool need_initialize) {
  AKANTU_DEBUG_IN();

  if (need_initialize) initializeUpdateResidualData();

  // f -= fint
  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    try {
      MaterialCohesive & mat = dynamic_cast<MaterialCohesive &>(**mat_it);
      mat.computeTraction(_not_ghost);
    } catch (std::bad_cast & bce) { }
  }

  SolidMechanicsModel::updateResidual(false);

  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    try {
      MaterialCohesive & mat = dynamic_cast<MaterialCohesive &>(**mat_it);
      mat.computeEnergies();
    } catch (std::bad_cast & bce) { }
  }


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::initFull(std::string material_file,
                                           AnalysisMethod method,
                                           bool extrinsic,
					   bool init_facet_filter) {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel::initFull(material_file, method);

  if(material_file != "")
    initCohesiveMaterial();

  if (extrinsic)
    initExtrinsic(material_file, init_facet_filter);

  mesh_facets.initByElementTypeArray(facet_material, 1, spatial_dimension - 1);
  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end();
       ++gt) {
    Mesh::type_iterator first = mesh_facets.firstType(spatial_dimension - 1, *gt);
    Mesh::type_iterator last  = mesh_facets.lastType(spatial_dimension - 1, *gt);
    for(;first != last; ++first) {
      Array<UInt> & f_material = facet_material(*first, *gt);
      f_material.resize(mesh_facets.getNbElement(*first, *gt));
      f_material.set(cohesive_index);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initParallel(MeshPartition * partition,
                                               DataAccessor * data_accessor,
                                               bool extrinsic) {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel::initParallel(partition, data_accessor);

  // /// create the distributed synchronizer for cohesive elements
  // cohesive_distributed_synchronizer =
  //   DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, NULL, 0,
  // 							       "cohesive_distributed_synchronizer");

  // synch_registry->registerSynchronizer(*cohesive_distributed_synchronizer,
  // 				       _gst_material_id);

  // synch_registry->registerSynchronizer(*cohesive_distributed_synchronizer,
  // 				       _gst_smm_stress);

  // synch_registry->registerSynchronizer(*cohesive_distributed_synchronizer,
  // 				       _gst_smm_boundary);

  /// create the facet synchronizer for extrinsic simulations
  if (extrinsic) {
    DistributedSynchronizer & distributed_synchronizer =
      dynamic_cast<DistributedSynchronizer &>(*synch_parallel);

    ByElementTypeUInt prank_to_element("prank_to_element", id);

    distributed_synchronizer.buildPrankToElement(prank_to_element);

    MeshUtils::buildAllFacetsParallel(mesh, mesh_facets, prank_to_element);
    facet_generated = true;

    facet_synchronizer =
      FacetSynchronizer::createFacetSynchronizer(distributed_synchronizer,
                                                 mesh_facets);

    synch_registry->registerSynchronizer(*facet_synchronizer, _gst_smmc_facets);
    synch_registry->registerSynchronizer(*facet_synchronizer, _gst_smmc_normals);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::initCohesiveMaterial() {

  AKANTU_DEBUG_IN();

  /// find the cohesive index
  cohesive_index = 0;

  while ((dynamic_cast<MaterialCohesive *>(materials[cohesive_index]) == NULL)
         && cohesive_index <= materials.size())
    ++cohesive_index;

  AKANTU_DEBUG_ASSERT(cohesive_index != materials.size(),
                      "No cohesive materials in the material input file");

  /// initialize element_index_by_material to to cohesive_index
  Material ** mat_val = &(materials.at(0));

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    /// fill the element filters of the materials using the element_material arrays
    Mesh::type_iterator it  = mesh.firstType(spatial_dimension, gt, _ek_cohesive);
    Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt, _ek_cohesive);

    for(; it != end; ++it) {
      UInt nb_element = mesh.getNbElement(*it, gt);

      Array<UInt> & el_id_by_mat = element_index_by_material(*it, gt);
      for (UInt el = 0; el < nb_element; ++el) {
        el_id_by_mat(el, 1) = cohesive_index;
        UInt index = mat_val[cohesive_index]->addElement(*it, el, gt);
        el_id_by_mat(el, 0) = index;
      }
    }
  }

  /// resize cohesive material vectors
  for (UInt m = 0; m < materials.size(); ++m) {
    try {
      MaterialCohesive & mat = dynamic_cast<MaterialCohesive &>(*materials[m]);
      mat.resizeCohesiveArrays();
    } catch(std::bad_cast&) { }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Initialize the model,basically it  pre-compute the shapes, shapes derivatives
 * and jacobian
 *
 */
void SolidMechanicsModelCohesive::initModel() {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel::initModel();

  registerFEMObject<MyFEMCohesiveType>("CohesiveFEM", mesh, spatial_dimension);

  /// add cohesive type connectivity
  ElementType type = _not_defined;

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType type_ghost = *gt;

    Mesh::type_iterator it   = mesh.firstType(spatial_dimension, type_ghost);
    Mesh::type_iterator last = mesh.lastType(spatial_dimension, type_ghost);

    for (; it != last; ++it) {
      const Array<UInt> & connectivity = mesh.getConnectivity(*it, type_ghost);
      if (connectivity.getSize() != 0) {
        type = *it;
        ElementType type_facet = Mesh::getFacetType(type);
        ElementType type_cohesive = FEM::getCohesiveElementType(type_facet);
        mesh.addConnectivityType(type_cohesive, type_ghost);
      }
    }
  }

  AKANTU_DEBUG_ASSERT(type != _not_defined, "No elements in the mesh");

  getFEM("CohesiveFEM").initShapeFunctions(_not_ghost);
  getFEM("CohesiveFEM").initShapeFunctions(_ghost);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::initExtrinsic(std::string material_file,
						bool init_facet_filter) {
  AKANTU_DEBUG_IN();

  if (!facet_generated)
    MeshUtils::buildAllFacets(mesh, mesh_facets);

  /// find element type for _not_ghost facet
  Mesh::type_iterator it_f   = mesh_facets.firstType(spatial_dimension - 1);
  Mesh::type_iterator last_f = mesh_facets.lastType(spatial_dimension - 1);

  internal_facet_type = _not_defined;

  for (; it_f != last_f; ++it_f) {
    UInt nb_facet = mesh_facets.getNbElement(*it_f);
    if (nb_facet != 0) {
      internal_facet_type = *it_f;
      break;
    }
  }

  AKANTU_DEBUG_ASSERT(internal_facet_type != _not_defined, "No facets in the mesh");

  for (UInt m = 0; m < materials.size(); ++m) {
    try {
      MaterialCohesive & mat = dynamic_cast<MaterialCohesive &>(*materials[m]);
      mat.initInsertionArrays(mesh_facets);
    } catch(std::bad_cast&) { }
  }


  /// init facet_filter
  if (init_facet_filter) {
    MaterialCohesive & mat_cohesive =
      dynamic_cast<MaterialCohesive &>(*materials[cohesive_index]);

    Array<UInt> & facet_filter = mat_cohesive.getFacetFilter(internal_facet_type);

    UInt nb_facet = mesh_facets.getNbElement(internal_facet_type);
    facet_filter.resize(nb_facet);
    for (UInt f = 0; f < nb_facet; ++f)
      facet_filter(f) = f;

    mat_cohesive.generateRandomDistribution(mesh_facets);
  }

  /// initialize facet insertion array
  mesh_facets.initByElementTypeArray(facet_insertion, 1,
                                     spatial_dimension - 1,
                                     false,
                                     _ek_regular,
                                     true);

  mesh.initByElementTypeArray(facets_to_cohesive_el, 2,
                              spatial_dimension, false,
                              _ek_cohesive);

  mesh_facets.initByElementTypeArray(cohesive_el_to_facet, 1, spatial_dimension - 1);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType gt_facet = *gt;

    Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1, gt_facet);
    Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1, gt_facet);

    for(; it != end; ++it) {
      ElementType type_facet = *it;
      UInt nb_facet = mesh_facets.getNbElement(type_facet, gt_facet);
      Array<UInt> & cohesive_el_to_f = cohesive_el_to_facet(type_facet, gt_facet);
      cohesive_el_to_f.resize(nb_facet);

      for (UInt f = 0; f < nb_facet; ++f)
        cohesive_el_to_f(f) = std::numeric_limits<UInt>::max();
    }
  }

  mesh_facets.initByElementTypeArray(doubled_facets, 2, spatial_dimension - 1);

  // if(material_file != "")
  //   buildFragmentsList();

  UInt nb_internal_facet = mesh_facets.getNbElement(internal_facet_type);

  mesh_facets.initByElementTypeArray(facet_stress,
				     spatial_dimension * spatial_dimension,
				     spatial_dimension - 1);

  const Array< std::vector<Element> > & element_to_facet
    = mesh_facets.getElementToSubelement(internal_facet_type);
  facets_check.resize(nb_internal_facet);

  for (UInt f = 0; f < nb_internal_facet; ++f) {
    if (element_to_facet(f)[1] != ElementNull) {
      facets_check(f) = true;
    }
    else facets_check(f) = false;
  }

  registerFEMObject<MyFEMType>("FacetsFEM", mesh_facets, spatial_dimension - 1);
  getFEM("FacetsFEM").initShapeFunctions(_not_ghost);
  getFEM("FacetsFEM").initShapeFunctions(_ghost);

  /// compute normals on facets
  computeSynchronizeNormals();

  /// compute quadrature points coordinates on facets
  Array<Real> & position = mesh.getNodes();

  ByElementTypeReal quad_facets("quad_facets", id);
  mesh_facets.initByElementTypeArray(quad_facets,
                                     spatial_dimension, spatial_dimension - 1);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType facet_gt = *gt;
    Mesh::type_iterator it   = mesh_facets.firstType(spatial_dimension - 1, facet_gt);
    Mesh::type_iterator last = mesh_facets.lastType(spatial_dimension - 1, facet_gt);

    for (; it != last; ++it) {
      ElementType facet_type = *it;

      UInt nb_quad_per_facet =
        getFEM("FacetsFEM").getNbQuadraturePoints(facet_type);

      UInt nb_facet = mesh_facets.getNbElement(facet_type, facet_gt);
      UInt nb_tot_quad = nb_quad_per_facet * nb_facet;

      Array<Real> & quad_f = quad_facets(facet_type, facet_gt);
      quad_f.resize(nb_tot_quad);

      getFEM("FacetsFEM").interpolateOnQuadraturePoints(position,
                                                        quad_f,
                                                        spatial_dimension,
                                                        facet_type,
                                                        facet_gt);
    }
  }

  /// compute elements quadrature point positions and build
  /// element-facet quadrature points data structure
  ByElementTypeReal elements_quad_facets("elements_quad_facets", id);

  mesh.initByElementTypeArray(elements_quad_facets,
                              spatial_dimension,
                              spatial_dimension);

  mesh.initByElementTypeArray(stress_on_facet,
                              spatial_dimension * spatial_dimension,
                              spatial_dimension);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType facet_gt = *gt;

    Mesh::type_iterator it   = mesh.firstType(spatial_dimension, facet_gt);
    Mesh::type_iterator last = mesh.lastType(spatial_dimension, facet_gt);

    for (; it != last; ++it) {
      ElementType type = *it;
      UInt nb_element = mesh.getNbElement(type, facet_gt);
      if (nb_element == 0) continue;

      /// compute elements' quadrature points and list of facet
      /// quadrature points positions by element
      Array<Element> & facet_to_element = mesh_facets.getSubelementToElement(type,
                                                                             facet_gt);
      UInt nb_facet_per_elem = facet_to_element.getNbComponent();

      Array<Real> & el_q_facet = elements_quad_facets(type, facet_gt);

      ElementType facet_type = Mesh::getFacetType(type);

      UInt nb_quad_per_facet =
        getFEM("FacetsFEM").getNbQuadraturePoints(facet_type);

      el_q_facet.resize(nb_element * nb_facet_per_elem * nb_quad_per_facet);

      stress_on_facet(type, facet_gt).resize(el_q_facet.getSize());

      for (UInt el = 0; el < nb_element; ++el) {
        for (UInt f = 0; f < nb_facet_per_elem; ++f) {
          Element global_facet_elem = facet_to_element(el, f);
          UInt global_facet = global_facet_elem.element;
          GhostType facet_gt = global_facet_elem.ghost_type;
          const Array<Real> & quad_f = quad_facets(facet_type, facet_gt);

          for (UInt q = 0; q < nb_quad_per_facet; ++q) {
            for (UInt s = 0; s < spatial_dimension; ++s) {
              el_q_facet(el * nb_facet_per_elem * nb_quad_per_facet
                         + f * nb_quad_per_facet + q, s)
                = quad_f(global_facet * nb_quad_per_facet + q, s);
            }
          }
        }
      }
    }
  }

  /// loop over materials
  for (UInt m = 0; m < materials.size(); ++m) {
    try {
      MaterialCohesive & mat __attribute__((unused)) =
	dynamic_cast<MaterialCohesive &>(*materials[m]);
    } catch(std::bad_cast&) {
      /// initialize the interpolation function
      materials[m]->initElementalFieldInterpolation(elements_quad_facets);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::computeSynchronizeNormals() {
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

  getFEM("FacetsFEM").computeNormalsOnControlPoints();

  /**
   *  @todo store tangents while computing normals instead of
   *  recomputing them as follows:
   */
  /* ------------------------------------------------------------------------ */
  const Array<Real> & normals =
    getFEM("FacetsFEM").getNormalsOnQuadPoints(internal_facet_type);

  Array<Real>::const_iterator< Vector<Real> > normal_it =
    normals.begin(spatial_dimension);

  mesh_facets.initByElementTypeArray(tangents, spatial_dimension, spatial_dimension - 1);
  Array<Real> & tang = tangents(internal_facet_type);
  tang.resize(normals.getSize());

  Array<Real>::iterator< Vector<Real> > tangent_it =
    tang.begin(spatial_dimension);

  for (UInt i = 0; i < normals.getSize(); ++i, ++normal_it, ++tangent_it) {
    Math::normal2( (*normal_it).storage(), (*tangent_it).storage() );
  }
  /* ------------------------------------------------------------------------ */

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::initFacetFilter() {
  AKANTU_DEBUG_IN();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end();
       ++gt) {
    GhostType gt_facet = *gt;

    Mesh::type_iterator first = mesh_facets.firstType(spatial_dimension - 1, gt_facet);
    Mesh::type_iterator last  = mesh_facets.lastType(spatial_dimension - 1, gt_facet);
    for(;first != last; ++first) {
      ElementType type_facet = *first;
      Array<UInt> & f_material = facet_material(type_facet, gt_facet);
      UInt nb_facet = f_material.getSize();

      for (UInt f = 0; f < nb_facet; ++f) {
	try {
	  MaterialCohesive & mat =
	    dynamic_cast<MaterialCohesive &>(*materials[f_material(f)]);

	  Array<UInt> & facet_filter = mat.getFacetFilter(type_facet, gt_facet);
	  facet_filter.push_back(f);

	} catch(std::bad_cast&) { }
      }
    }
  }

  for (UInt m = 0; m < materials.size(); ++m) {
    try {
      MaterialCohesive & mat_cohesive = dynamic_cast<MaterialCohesive &>(*materials[m]);
      mat_cohesive.generateRandomDistribution(mesh_facets);
    } catch(std::bad_cast&) { }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::checkCohesiveStress() {
  AKANTU_DEBUG_IN();

  UInt nb_materials = materials.size();

  UInt nb_quad_per_facet = getFEM("FacetsFEM").getNbQuadraturePoints(internal_facet_type);
  UInt nb_facet = mesh_facets.getNbElement(internal_facet_type);
  UInt sp2 = spatial_dimension * spatial_dimension;
  UInt nb_quad_f_two = nb_quad_per_facet * 2;

  /// vector containing stresses coming from the two elements of each facets
  Array<Real> & f_stress = facet_stress(internal_facet_type);
  f_stress.resize(2 * nb_facet * nb_quad_per_facet);

  facet_stress_count.resize(nb_facet);
  facet_stress_count.clear();

  /// loop over materials
  for (UInt m = 0; m < nb_materials; ++m) {
    try {
      MaterialCohesive & mat __attribute__((unused)) =
	dynamic_cast<MaterialCohesive &>(*materials[m]);
    } catch(std::bad_cast&) {
      /// interpolate stress on facet quadrature points positions
      materials[m]->interpolateStress(stress_on_facet);

      for (ghost_type_t::iterator gt = ghost_type_t::begin();
           gt != ghost_type_t::end(); ++gt) {

        GhostType ghost_type = *gt;

        Mesh::type_iterator it   = mesh.firstType(spatial_dimension, ghost_type);
        Mesh::type_iterator last = mesh.lastType(spatial_dimension, ghost_type);

        /// loop over element type
        for (; it != last; ++it) {
          ElementType type = *it;
          UInt nb_element = mesh.getNbElement(type, ghost_type);
          if (nb_element == 0) continue;

          Array<Real> & stress_on_f = stress_on_facet(type, ghost_type);

          /// store the interpolated stresses on the facet_stress vector
          Array<Element> & facet_to_element =
            mesh_facets.getSubelementToElement(type, ghost_type);

          UInt nb_facet_per_elem = facet_to_element.getNbComponent();

          Array<Element>::iterator<Vector<Element> > facet_to_el_it =
            facet_to_element.begin(nb_facet_per_elem);

          Array<Real>::iterator< Matrix<Real> > stress_on_f_it =
            stress_on_f.begin(spatial_dimension, spatial_dimension);

          for (UInt el = 0; el < nb_element; ++el, ++facet_to_el_it) {
            for (UInt f = 0; f < nb_facet_per_elem; ++f) {
              Element global_facet_elem = (*facet_to_el_it)(f);
              UInt global_facet = global_facet_elem.element;
              if (global_facet_elem.ghost_type == _ghost) {
                stress_on_f_it += nb_quad_per_facet;
                continue;
              }

              for (UInt q = 0; q < nb_quad_per_facet; ++q, ++stress_on_f_it) {

                if (facets_check(global_facet) == true) {
                  Matrix<Real> facet_stress_local(f_stress.storage()
                                                  + (global_facet * nb_quad_f_two
                                                     + q * 2
                                                     + facet_stress_count(global_facet))
                                                  * sp2,
                                                  spatial_dimension,
                                                  spatial_dimension);

                  facet_stress_local = *stress_on_f_it;
                }

              }
              facet_stress_count(global_facet) = true;
            }
          }
        }
      }
    }
  }


  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    GhostType facet_gt = *gt;
    Mesh::type_iterator it   = mesh_facets.firstType(spatial_dimension - 1, facet_gt);
    Mesh::type_iterator last = mesh_facets.lastType(spatial_dimension - 1, facet_gt);
    for (; it != last; ++it)
      facet_insertion(*it, facet_gt).clear();
  }

  for (UInt m = 0; m < nb_materials; ++m) {
    try {
      MaterialCohesive & mat_cohesive = dynamic_cast<MaterialCohesive &>(*materials[m]);
      /// check which not ghost cohesive elements are to be created
      mat_cohesive.checkInsertion(facet_stress, mesh_facets, facet_insertion);
    } catch(std::bad_cast&) { }
  }

  /// communicate data among processors
  synch_registry->synchronize(_gst_smmc_facets);

  /// insert cohesive elements
  MeshUtils::insertCohesiveElements(mesh,
                                    mesh_facets,
                                    facet_insertion,
                                    doubled_facets,
                                    true);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::onElementsAdded(__attribute__((unused)) const Array<Element> & element_list,
                                                  __attribute__((unused)) const NewElementsEvent & event) {
  AKANTU_DEBUG_IN();

  /// update model data for new cohesive elements
  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType gt_facet = *gt;

    Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1, gt_facet);
    Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1, gt_facet);

    for(; it != end; ++it) {

      ElementType type_facet = *it;
      Array<UInt> & doubled_f = doubled_facets(type_facet, gt_facet);

      if(doubled_f.getSize() == 0) continue;

      /// define facet material vector
      ElementType type_cohesive = FEM::getCohesiveElementType(type_facet);

      Array<UInt> & el_id_by_mat = element_index_by_material(type_cohesive,
                                                             gt_facet);

      UInt old_nb_cohesive_elements = el_id_by_mat.getSize();
      UInt new_cohesive_elements = doubled_f.getSize();
      UInt total_nb_cohesive_elements = old_nb_cohesive_elements + new_cohesive_elements;

      el_id_by_mat.resize(total_nb_cohesive_elements);

      if (facets_check.getSize() > 0 && gt_facet == _not_ghost)
        facets_check.resize(facets_check.getSize() + new_cohesive_elements);

      Array<UInt> & f_to_cohesive_el = facets_to_cohesive_el(type_cohesive, gt_facet);
      f_to_cohesive_el.resize(total_nb_cohesive_elements);

      Array<UInt> & cohesive_el_to_f = cohesive_el_to_facet(type_facet, gt_facet);

      const Array<UInt> & facet_material_by_type = facet_material(type_facet, gt_facet);

      Array<UInt>::iterator<Vector<UInt> > d_f_it = doubled_f.begin(2);

      for (UInt f = 0; f < doubled_f.getSize(); ++f, ++d_f_it) {
        UInt old_facet = (*d_f_it)(0);
        UInt new_facet = (*d_f_it)(1);

        /// assign the cohesive material to each new cohesive element
        UInt nb_cohesive_elements = old_nb_cohesive_elements + f;

	UInt material_id = facet_material_by_type(old_facet);

        el_id_by_mat(nb_cohesive_elements, 0) =
          materials[material_id]->addElement(type_cohesive,
					     nb_cohesive_elements,
					     gt_facet);

        el_id_by_mat(nb_cohesive_elements, 1) = material_id;

        /// update facets_check vector
        if (facets_check.getSize() > 0 && gt_facet == _not_ghost) {
          facets_check(old_facet) = false;
          facets_check(new_facet) = false;
        }

        /// update facets_to_cohesive_el vector
        f_to_cohesive_el(nb_cohesive_elements, 0) = old_facet;
        f_to_cohesive_el(nb_cohesive_elements, 1) = new_facet;

        if (facet_synchronizer != NULL)
          cohesive_el_to_f(old_facet) = nb_cohesive_elements;
      }
    }

    /// update shape functions
    getFEM("CohesiveFEM").initShapeFunctions(gt_facet);
  }

  /// resize cohesive material vectors
  for (UInt m = 0; m < materials.size(); ++m) {
    try {
      MaterialCohesive & mat = dynamic_cast<MaterialCohesive &>(*materials[m]);
      mat.resizeCohesiveArrays();
    } catch(std::bad_cast&) { }
  }

  assembleMassLumped();

  /// update synchronizer if needed
  if (facet_synchronizer != NULL) {
    DataAccessor * data_accessor = this;

    DistributedSynchronizer & distributed_synchronizer =
      dynamic_cast<DistributedSynchronizer &>(*synch_parallel);

    facet_synchronizer->updateDistributedSynchronizer(distributed_synchronizer,
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

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::onNodesAdded(const Array<UInt> & doubled_nodes,
                                               __attribute__((unused)) const NewNodesEvent & event) {
  AKANTU_DEBUG_IN();

  UInt nb_new_nodes = doubled_nodes.getSize();
  Array<UInt> nodes_list(nb_new_nodes);

  for (UInt n = 0; n < nb_new_nodes; ++n)
    nodes_list(n) = doubled_nodes(n, 1);

  SolidMechanicsModel::onNodesAdded(nodes_list, event);

  for (UInt n = 0; n < nb_new_nodes; ++n) {

    UInt old_node = doubled_nodes(n, 0);
    UInt new_node = doubled_nodes(n, 1);

    for (UInt dim = 0; dim < spatial_dimension; ++dim) {
      (*displacement)(new_node, dim) = (*displacement)(old_node, dim);
      (*velocity)(new_node, dim) = (*velocity)(old_node, dim);
      (*acceleration)(new_node, dim) = (*acceleration)(old_node, dim);
      (*boundary)(new_node, dim) = (*boundary)(old_node, dim);

      if (current_position)
        (*current_position)(new_node, dim) = (*current_position)(old_node, dim);

      if (increment_acceleration)
        (*increment_acceleration)(new_node, dim)
          = (*increment_acceleration)(old_node, dim);

      if (increment)
        (*increment)(new_node, dim) = (*increment)(old_node, dim);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::buildFragmentsList() {
  AKANTU_DEBUG_IN();

  Mesh::type_iterator it   = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last = mesh.lastType(spatial_dimension);

  const UInt max = std::numeric_limits<UInt>::max();

  for (; it != last; ++it) {
    UInt nb_element = mesh.getNbElement(*it);
    if (nb_element != 0) {
      /// initialize the list of checked elements and the list of
      /// elements to be checked
      fragment_to_element.alloc(nb_element, 1, *it);
      Array<UInt> & frag_to_el = fragment_to_element(*it);

      for (UInt el = 0; el < nb_element; ++el)
        frag_to_el(el) = max;
    }
  }

  Array< std::vector<Element> > & element_to_facet
    = mesh_facets.getElementToSubelement(internal_facet_type);

  ElementType type_cohesive = FEM::getCohesiveElementType(internal_facet_type);
  const Array<UInt> & f_to_cohesive_el = facets_to_cohesive_el(type_cohesive);

  nb_fragment = 0;
  it = mesh.firstType(spatial_dimension);

  MaterialCohesive * mat_cohesive
    = dynamic_cast<MaterialCohesive*>(materials[cohesive_index]);
  const Array<Real> & damage = mat_cohesive->getDamage(type_cohesive);
  UInt nb_quad_cohesive = getFEM("CohesiveFEM").getNbQuadraturePoints(type_cohesive);

  Real epsilon = std::numeric_limits<Real>::epsilon();

  for (; it != last; ++it) {
    Array<UInt> & checked_el = fragment_to_element(*it);
    UInt nb_element = checked_el.getSize();
    if (nb_element != 0) {

      Array<Element> & facet_to_element = mesh_facets.getSubelementToElement(*it);
      UInt nb_facet_per_elem = facet_to_element.getNbComponent();

      /// loop on elements
      for (UInt el = 0; el < nb_element; ++el) {
        if (checked_el(el) == max) {
          /// build fragment
          ++nb_fragment;
          checked_el(el) = nb_fragment - 1;
          Array<Element> elem_to_check;

          Element first_el(*it, el);
          elem_to_check.push_back(first_el);

          /// keep looping while there are elements to check
          while (elem_to_check.getSize() != 0) {
            UInt nb_elem_check = elem_to_check.getSize();

            for (UInt el_check = 0; el_check < nb_elem_check; ++el_check) {
              Element current_el = elem_to_check(el_check);

              for (UInt f = 0; f < nb_facet_per_elem; ++f) {
                /// find adjacent element on current facet
                UInt global_facet = facet_to_element(current_el.element, f).element;
                Element next_el;
                for (UInt i = 0; i < 2; ++i) {
                  next_el = element_to_facet(global_facet)[i];
                  if (next_el != current_el) break;
                }

                if (next_el.kind == _ek_cohesive) {
                  /// fragmention occurs when the cohesive element has
                  /// reached damage = 1 on every quadrature point
                  UInt q = 0;
                  while (q < nb_quad_cohesive &&
                         std::abs(damage(next_el.element * nb_quad_cohesive + q) - 1)
                         <= epsilon) ++q;

                  if (q == nb_quad_cohesive)
                    next_el = ElementNull;
                  else {
                    /// check which facet is the correct one
                    UInt other_facet_index
                      = f_to_cohesive_el(next_el.element, 0) == global_facet;
                    UInt other_facet
                      = f_to_cohesive_el(next_el.element, other_facet_index);

                    /// get the other regualar element
                    next_el = element_to_facet(other_facet)[0];
                  }
                }

                /// if it exists, add it to the fragment list
                if (next_el != ElementNull) {
                  Array<UInt> & checked_next_el = fragment_to_element(next_el.type);
                  /// check if the element isn't already part of a fragment
                  if (checked_next_el(next_el.element) == max) {
                    checked_next_el(next_el.element) = nb_fragment - 1;
                    elem_to_check.push_back(next_el);
                  }
                }
              }
            }

            /// erase elements that have already been checked
            for (UInt el_check = nb_elem_check; el_check > 0; --el_check)
              elem_to_check.erase(el_check - 1);
          }
        }
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::computeFragmentsMV() {
  AKANTU_DEBUG_IN();

  fragment_mass.resize(nb_fragment);
  fragment_mass.clear();

  fragment_velocity.resize(nb_fragment);
  fragment_velocity.clear();

  fragment_center.resize(nb_fragment);
  fragment_center.clear();

  UInt nb_nodes = mesh.getNbNodes();
  Array<bool> checked_node(nb_nodes);
  checked_node.clear();

  const Array<Real> & position = mesh.getNodes();

  Mesh::type_iterator it   = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last = mesh.lastType(spatial_dimension);

  for (; it != last; ++it) {
    UInt nb_element = mesh.getNbElement(*it);
    if (nb_element != 0) {
      const Array<UInt> & frag_to_el = fragment_to_element(*it);
      const Array<UInt> & connectivity = mesh.getConnectivity(*it);
      UInt nb_nodes_per_elem = connectivity.getNbComponent();

      /// loop over each node of every element
      for (UInt el = 0; el < nb_element; ++el) {
        for (UInt n = 0; n < nb_nodes_per_elem; ++n) {
          UInt node = connectivity(el, n);
          /// if the node hasn't been checked, store its data
          if (!checked_node(node)) {
            fragment_mass(frag_to_el(el)) += (*mass)(node);

            for (UInt s = 0; s < spatial_dimension; ++s) {
              fragment_velocity(frag_to_el(el), s)
                += (*mass)(node) * (*velocity)(node, s);

              fragment_center(frag_to_el(el), s)
                += (*mass)(node) * position(node, s);
            }

            checked_node(node) = true;
          }
        }
      }
    }
  }

  for (UInt frag = 0; frag < nb_fragment; ++frag) {
    for (UInt s = 0; s < spatial_dimension; ++s) {
      fragment_velocity(frag, s) /= fragment_mass(frag);
      fragment_center(frag, s) /= fragment_mass(frag);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::computeFragmentsData() {
  AKANTU_DEBUG_IN();

  buildFragmentsList();
  computeFragmentsMV();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "SolidMechanicsModelCohesive [" << std::endl;

  SolidMechanicsModel::printself(stream, indent + 1);

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */


__END_AKANTU__
