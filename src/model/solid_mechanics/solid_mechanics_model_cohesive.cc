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
#include <algorithm>
#include "shape_cohesive.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "material_cohesive.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

const SolidMechanicsModelCohesiveOptions default_solid_mechanics_model_cohesive_options(_explicit_lumped_mass,
											false,
											false,
											true,
											true);

/* -------------------------------------------------------------------------- */

SolidMechanicsModelCohesive::SolidMechanicsModelCohesive(Mesh & mesh,
                                                         UInt dim,
                                                         const ID & id,
                                                         const MemoryID & memory_id) :
  SolidMechanicsModel(mesh, dim, id, memory_id),
  mesh_facets(mesh.initMeshFacets("mesh_facets")),
  tangents("tangents", id),
  stress_on_facet("stress_on_facet", id),
  facet_stress("facet_stress", id),
  fragment_mass(0, spatial_dimension, "fragment_mass"),
  fragment_velocity(0, spatial_dimension, "fragment_velocity"),
  fragment_center(0, spatial_dimension, "fragment_center"),
  facet_material("facet_material", id),
  inserter(mesh, mesh_facets) {
  AKANTU_DEBUG_IN();

  facet_generated = false;

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  facet_synchronizer = NULL;
  facet_stress_synchronizer = NULL;
  cohesive_distributed_synchronizer = NULL;
  global_connectivity = NULL;
#endif

  delete material_selector;
  material_selector = new DefaultMaterialCohesiveSelector(*this);

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
SolidMechanicsModelCohesive::~SolidMechanicsModelCohesive() {
  AKANTU_DEBUG_IN();

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  delete cohesive_distributed_synchronizer;
  delete facet_synchronizer;
  delete facet_stress_synchronizer;
#endif

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initFull(std::string material_file,
					   const ModelOptions & options) {
  AKANTU_DEBUG_IN();

  const SolidMechanicsModelCohesiveOptions & smmc_options =
    dynamic_cast<const SolidMechanicsModelCohesiveOptions &>(options);

  this->stress_interpolation = smmc_options.stress_interpolation;
  this->is_extrinsic = smmc_options.extrinsic;

  if (is_extrinsic) {
    if (!facet_generated)
      MeshUtils::buildAllFacets(mesh, mesh_facets);
  }

  SolidMechanicsModel::initFull(material_file, options);

  if (is_extrinsic)
    initAutomaticInsertion();

  // if(material_file != "")
  //   initCohesiveMaterials(smmc_options.extrinsic, smmc_options.init_facet_filter);

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initMaterials() {
  AKANTU_DEBUG_IN();

  /// init element inserter
  inserter.init();

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  if (facet_synchronizer != NULL)
    inserter.initParallel(facet_synchronizer);

  if (facet_stress_synchronizer != NULL) {
    DataAccessor * data_accessor = this;
    const ByElementTypeUInt & rank_to_element = synch_parallel->getPrankToElement();

    facet_stress_synchronizer->updateFacetStressSynchronizer(inserter,
							     rank_to_element,
							     *data_accessor);
  }
#endif

  // make sure the material are instantiated
  if(!are_materials_instantiated) instantiateMaterials();

  /// find the first cohesive material
  UInt cohesive_index = 0;

  while ((dynamic_cast<MaterialCohesive *>(materials[cohesive_index]) == NULL)
         && cohesive_index <= materials.size())
    ++cohesive_index;

  AKANTU_DEBUG_ASSERT(cohesive_index != materials.size(),
                      "No cohesive materials in the material input file");

  material_selector->setFallback(cohesive_index);

  // set the facet information in the material in case of dynamic insertion
  mesh_facets.initByElementTypeArray(facet_material, 1, spatial_dimension - 1);

  if (is_extrinsic) {
    Element element;
    for (ghost_type_t::iterator gt = ghost_type_t::begin(); gt != ghost_type_t::end(); ++gt) {
      element.ghost_type = *gt;
      Mesh::type_iterator first = mesh_facets.firstType(spatial_dimension - 1, *gt);
      Mesh::type_iterator last  = mesh_facets.lastType(spatial_dimension - 1, *gt);
      for(;first != last; ++first) {
	element.type = *first;
	Array<UInt> & f_material = facet_material(*first, *gt);
	UInt nb_element = mesh_facets.getNbElement(*first, *gt);
	f_material.resize(nb_element);
	f_material.set(cohesive_index);
	for (UInt el = 0; el < nb_element; ++el) {
	  element.element = el;
	  UInt mat_index = (*material_selector)(element);
	  f_material(el) = mat_index;
	  MaterialCohesive & mat = dynamic_cast<MaterialCohesive &>(*materials[mat_index]);
	  mat.addFacet(element);
	}
      }
    }
  } else {
    for (ghost_type_t::iterator gt = ghost_type_t::begin(); gt != ghost_type_t::end(); ++gt) {
      Mesh::type_iterator first = mesh.firstType(spatial_dimension, *gt, _ek_cohesive);
      Mesh::type_iterator last  = mesh.lastType(spatial_dimension, *gt, _ek_cohesive);

      for(;first != last; ++first) {
	Array<UInt> & el_id_by_mat = element_index_by_material(*first, *gt);
	Vector<UInt> el_mat(2); el_mat(0) = cohesive_index; el_mat(1) = 0;
	el_id_by_mat.set(el_mat);
      }
    }
  }

  SolidMechanicsModel::initMaterials();

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

  registerFEMObject<MyFEMType>("FacetsFEM", mesh_facets, spatial_dimension - 1);
  if (is_extrinsic) {
    getFEM("FacetsFEM").initShapeFunctions(_not_ghost);
    getFEM("FacetsFEM").initShapeFunctions(_ghost);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initAutomaticInsertion() {
  AKANTU_DEBUG_IN();

  mesh_facets.initByElementTypeArray(facet_stress,
				     2 * spatial_dimension * spatial_dimension,
				     spatial_dimension - 1);

  resizeFacetStress();

  /// compute normals on facets
  computeNormals();

  /// allocate stress_on_facet to store element stress on facets
  mesh.initByElementTypeArray(stress_on_facet,
                              spatial_dimension * spatial_dimension,
                              spatial_dimension);

  Mesh::type_iterator it   = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last = mesh.lastType(spatial_dimension);

  for (; it != last; ++it) {
    ElementType type = *it;
    UInt nb_element = mesh.getNbElement(type);

    UInt nb_facet_per_elem = Mesh::getNbFacetsPerElement(type);
    ElementType type_facet = Mesh::getFacetType(type);
    UInt nb_quad_per_facet = getFEM("FacetsFEM").getNbQuadraturePoints(type_facet);

    stress_on_facet(type).resize(nb_quad_per_facet * nb_facet_per_elem * nb_element);
  }

  if (stress_interpolation)
    initStressInterpolation();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::updateAutomaticInsertion() {
  AKANTU_DEBUG_IN();

  inserter.limitCheckFacets();

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
    if (facet_stress_synchronizer != NULL) {
      DataAccessor * data_accessor = this;
      const ByElementTypeUInt & rank_to_element = synch_parallel->getPrankToElement();

      facet_stress_synchronizer->updateFacetStressSynchronizer(inserter,
							       rank_to_element,
							       *data_accessor);
    }
#endif

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initStressInterpolation() {
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

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType elem_gt = *gt;

    Mesh::type_iterator it   = mesh.firstType(spatial_dimension, elem_gt);
    Mesh::type_iterator last = mesh.lastType(spatial_dimension, elem_gt);

    for (; it != last; ++it) {
      ElementType type = *it;
      UInt nb_element = mesh.getNbElement(type, elem_gt);
      if (nb_element == 0) continue;

      /// compute elements' quadrature points and list of facet
      /// quadrature points positions by element
      Array<Element> & facet_to_element = mesh_facets.getSubelementToElement(type,
                                                                             elem_gt);
      UInt nb_facet_per_elem = facet_to_element.getNbComponent();

      Array<Real> & el_q_facet = elements_quad_facets(type, elem_gt);

      ElementType facet_type = Mesh::getFacetType(type);

      UInt nb_quad_per_facet =
        getFEM("FacetsFEM").getNbQuadraturePoints(facet_type);

      el_q_facet.resize(nb_element * nb_facet_per_elem * nb_quad_per_facet);

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

  /// loop over non cohesive materials
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
void SolidMechanicsModelCohesive::computeNormals() {
  AKANTU_DEBUG_IN();

  getFEM("FacetsFEM").computeNormalsOnControlPoints(_not_ghost);

  /**
   *  @todo store tangents while computing normals instead of
   *  recomputing them as follows:
   */
  /* ------------------------------------------------------------------------ */
  UInt tangent_components = spatial_dimension * (spatial_dimension - 1);

  mesh_facets.initByElementTypeArray(tangents,
				     tangent_components,
				     spatial_dimension - 1);

  Mesh::type_iterator it   = mesh_facets.firstType(spatial_dimension - 1);
  Mesh::type_iterator last = mesh_facets.lastType(spatial_dimension - 1);

  for (; it != last; ++it) {
    ElementType facet_type = *it;

    const Array<Real> & normals =
      getFEM("FacetsFEM").getNormalsOnQuadPoints(facet_type);

    UInt nb_quad = normals.getSize();

    Array<Real> & tang = tangents(facet_type);
    tang.resize(nb_quad);

    Real * normal_it = normals.storage();
    Real * tangent_it = tang.storage();

    /// compute first tangent
    for (UInt q = 0; q < nb_quad; ++q) {

      /// if normal is orthogonal to xy plane, arbitrarly define tangent
      if ( Math::are_float_equal(Math::norm2(normal_it), 0) )
	tangent_it[0] = 1;
      else
	Math::normal2(normal_it, tangent_it);

      normal_it += spatial_dimension;
      tangent_it += tangent_components;
    }

    /// compute second tangent (3D case)
    if (spatial_dimension == 3) {
      normal_it = normals.storage();
      tangent_it = tang.storage();

      for (UInt q = 0; q < nb_quad; ++q) {
	Math::normal3(normal_it, tangent_it, tangent_it + spatial_dimension);
	normal_it += spatial_dimension;
	tangent_it += tangent_components;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::averageStressOnFacets(UInt material_index) {
  AKANTU_DEBUG_IN();

  Mesh::type_iterator it   = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last = mesh.lastType(spatial_dimension);

  /// loop over element type
  for (; it != last; ++it) {
    ElementType type = *it;
    UInt nb_element = mesh.getNbElement(type);
    UInt nb_quad_per_element = getFEM().getNbQuadraturePoints(type);
    const Array<Real> & stress = materials[material_index]->getStress(type);
    Array<Real> & s_on_facet = stress_on_facet(type);

    UInt nb_facet_per_elem = Mesh::getNbFacetsPerElement(type);
    ElementType type_facet = Mesh::getFacetType(type);
    UInt nb_quad_per_facet = getFEM("FacetsFEM").getNbQuadraturePoints(type_facet);
    UInt nb_facet_quad_per_elem = nb_quad_per_facet * nb_facet_per_elem;

    Array<Real>::const_iterator<Matrix<Real> > stress_it
      = stress.begin(spatial_dimension, spatial_dimension);
    Array<Real>::iterator<Matrix<Real> > s_on_facet_it
      = s_on_facet.begin(spatial_dimension, spatial_dimension);

    Matrix<Real> average_stress(spatial_dimension, spatial_dimension);

    for (UInt el = 0; el < nb_element; ++el) {

      /// compute average stress
      average_stress.clear();

      for (UInt q = 0; q < nb_quad_per_element; ++q, ++stress_it)
	average_stress += *stress_it;

      average_stress /= nb_quad_per_element;

      /// store average stress
      for (UInt q = 0; q < nb_facet_quad_per_elem; ++q, ++s_on_facet_it)
	*s_on_facet_it = average_stress;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::fillStressOnFacet() {
  AKANTU_DEBUG_IN();

  UInt sp2 = spatial_dimension * spatial_dimension;
  UInt sp4 = sp2 * 2;

  /// loop over materials
  for (UInt m = 0; m < materials.size(); ++m) {
    try {
      MaterialCohesive & mat __attribute__((unused)) =
        dynamic_cast<MaterialCohesive &>(*materials[m]);
    } catch(std::bad_cast&) {

      if (stress_interpolation)
	/// interpolate stress on facet quadrature points positions
	materials[m]->interpolateStress(stress_on_facet);
      else
	averageStressOnFacets(m);

      GhostType ghost_type = _not_ghost;

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

        ElementType facet_type = _not_defined;
        GhostType facet_gt = _casper;
        Array<std::vector<Element> > * element_to_facet = NULL;
        Array<Real> * f_stress = NULL;
        Array<bool> * f_check = NULL;
        UInt nb_quad_per_facet = 0;
        UInt element_rank = 0;

        Element current_el(type, 0, ghost_type);

        for (UInt el = 0; el < nb_element; ++el, ++facet_to_el_it) {
          current_el.element = el;
          for (UInt f = 0; f < nb_facet_per_elem; ++f) {
            Element global_facet_elem = (*facet_to_el_it)(f);
            UInt global_facet = global_facet_elem.element;

            if (facet_type != global_facet_elem.type ||
                facet_gt != global_facet_elem.ghost_type) {
              facet_type = global_facet_elem.type;
              facet_gt = global_facet_elem.ghost_type;

              element_to_facet =
                &( mesh_facets.getElementToSubelement(facet_type, facet_gt) );
              f_stress = &( facet_stress(facet_type, facet_gt) );

              nb_quad_per_facet =
                getFEM("FacetsFEM").getNbQuadraturePoints(facet_type, facet_gt);

              f_check = &( inserter.getCheckFacets(facet_type, facet_gt) );
            }

	    if (!(*f_check)(global_facet)) {
	      stress_on_f_it += nb_quad_per_facet;
	      continue;
	    }

            for (UInt q = 0; q < nb_quad_per_facet; ++q, ++stress_on_f_it) {

	      element_rank = (*element_to_facet)(global_facet)[0] != current_el;

	      Matrix<Real> facet_stress_local(f_stress->storage()
					      + (global_facet * nb_quad_per_facet + q) * sp4
					      + element_rank * sp2,
					      spatial_dimension,
					      spatial_dimension);

	      facet_stress_local = *stress_on_f_it;
            }
          }
        }
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::checkCohesiveStress() {
  AKANTU_DEBUG_IN();

  fillStressOnFacet();

#if defined(AKANTU_DEBUG_TOOLS)
  debug::element_manager.printData(debug::_dm_model_cohesive,
				   "Interpolated stresses before",
				   facet_stress);
#endif

  synch_registry->synchronize(_gst_smmc_facets_stress);

#if defined(AKANTU_DEBUG_TOOLS)
  debug::element_manager.printData(debug::_dm_model_cohesive,
				   "Interpolated stresses",
				   facet_stress);
#endif

  for (UInt m = 0; m < materials.size(); ++m) {
    try {
      MaterialCohesive & mat_cohesive = dynamic_cast<MaterialCohesive &>(*materials[m]);
      /// check which not ghost cohesive elements are to be created
      mat_cohesive.checkInsertion();
    } catch(std::bad_cast&) { }
  }

  /// communicate data among processors
  synch_registry->synchronize(_gst_smmc_facets);

  /// insert cohesive elements
  inserter.insertExtrinsicElements();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::onElementsAdded(const Array<Element> & element_list,
                                                  const NewElementsEvent & event) {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel::onElementsAdded(element_list, event);

  /// update shape functions
  getFEM("CohesiveFEM").initShapeFunctions(_not_ghost);
  getFEM("CohesiveFEM").initShapeFunctions(_ghost);
  if(is_extrinsic) {
    getFEM("FacetsFEM").initShapeFunctions(_not_ghost);
    getFEM("FacetsFEM").initShapeFunctions(_ghost);
  }

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  updateFacetSynchronizers();
#endif

  resizeFacetStress();

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
      (*velocity)    (new_node, dim) = (*velocity)    (old_node, dim);
      (*acceleration)(new_node, dim) = (*acceleration)(old_node, dim);
      (*boundary)    (new_node, dim) = (*boundary)    (old_node, dim);

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
class CohesiveElementFilter : public GroupManager::ClusteringFilter {
public:
  CohesiveElementFilter(const SolidMechanicsModelCohesive & model,
			const Real max_damage = 1.) :
    model(model), is_unbroken(max_damage) {}

  bool operator() (const Element & el) const {
    if (el.kind == _ek_regular)
      return true;

    const Array<UInt> & el_id_by_mat = model.getElementIndexByMaterial(el.type,
								       el.ghost_type);

    const MaterialCohesive & mat
      = static_cast<const MaterialCohesive &>
      (model.getMaterial(el_id_by_mat(el.element, 0)));

    UInt el_index = el_id_by_mat(el.element, 1);
    UInt nb_quad_per_element
      = model.getFEM("CohesiveFEM").getNbQuadraturePoints(el.type, el.ghost_type);

    Vector<Real> damage
      = mat.getDamage(el.type, el.ghost_type).begin(nb_quad_per_element)[el_index];

    UInt unbroken_quads = std::count_if(damage.storage(),
					damage.storage() + nb_quad_per_element,
					is_unbroken);

    if (unbroken_quads > 0)
      return true;

    return false;
  }

private:

  struct IsUnbrokenFunctor {
    IsUnbrokenFunctor(const Real & max_damage) : max_damage(max_damage) {}
    bool operator() (const Real & x) {return x < max_damage;}
    const Real max_damage;
  };

  const SolidMechanicsModelCohesive & model;
  const IsUnbrokenFunctor is_unbroken;
};

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::buildFragments() {
  AKANTU_DEBUG_IN();

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  if (cohesive_distributed_synchronizer) {
    cohesive_distributed_synchronizer->computeBufferSize(*this, _gst_smmc_damage);
    cohesive_distributed_synchronizer->asynchronousSynchronize(*this, _gst_smmc_damage);
    cohesive_distributed_synchronizer->waitEndSynchronize(*this, _gst_smmc_damage);
  }
#endif

  nb_fragment = mesh.createClusters(spatial_dimension, "fragment",
				    CohesiveElementFilter(*this),
				    synch_parallel,
				    &mesh_facets);

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

  Array<Real>::const_iterator< Vector<Real> > it_mass =
    mass->begin(spatial_dimension);
  Array<Real>::const_iterator< Vector<Real> > it_velocity =
    velocity->begin(spatial_dimension);
  Array<Real>::const_iterator< Vector<Real> > it_position =
    mesh.getNodes().begin(spatial_dimension);

  Array<Real>::iterator< Vector<Real> > it_frag_mass =
    fragment_mass.begin(spatial_dimension);
  Array<Real>::iterator< Vector<Real> > it_frag_velocity =
    fragment_velocity.begin(spatial_dimension);
  Array<Real>::iterator< Vector<Real> > it_frag_center =
    fragment_center.begin(spatial_dimension);


  /// compute data for local fragments
  for (UInt frag = 0; frag < nb_fragment; ++frag) {
    std::stringstream sstr;
    sstr << "fragment_" << frag;

    GroupManager::const_element_group_iterator eg_it
      = mesh.element_group_find(sstr.str());

    if (eg_it == mesh.element_group_end()) continue;

    Vector<Real> & frag_mass = it_frag_mass[frag];
    Vector<Real> & frag_velocity = it_frag_velocity[frag];
    Vector<Real> & frag_center = it_frag_center[frag];

    const NodeGroup & node_group = eg_it->second->getNodeGroup();
    NodeGroup::const_node_iterator it = node_group.begin();
    NodeGroup::const_node_iterator end = node_group.end();

    for (; it != end; ++it) {
      UInt node = *it;
      if (mesh.isLocalOrMasterNode(node)) {
	const Vector<Real> & node_mass = it_mass[node];
	const Vector<Real> & node_velocity = it_velocity[node];
	const Vector<Real> & node_position = it_position[node];

	frag_mass += node_mass;

	for (UInt s = 0; s < spatial_dimension; ++s) {
	  frag_velocity(s)
	    += node_mass(s) * node_velocity(s);

	  frag_center(s)
	    += node_mass(s) * node_position(s);
	}
      }
    }
  }

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  /// reduce all data
  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  comm.allReduce(fragment_mass.storage()    , nb_fragment * spatial_dimension, _so_sum);
  comm.allReduce(fragment_velocity.storage(), nb_fragment * spatial_dimension, _so_sum);
  comm.allReduce(fragment_center.storage()  , nb_fragment * spatial_dimension, _so_sum);
#endif

  /// compute averages
  it_frag_mass = fragment_mass.begin(spatial_dimension);
  it_frag_velocity = fragment_velocity.begin(spatial_dimension);
  it_frag_center = fragment_center.begin(spatial_dimension);

  for (UInt frag = 0;
       frag < nb_fragment;
       ++frag, ++it_frag_mass, ++it_frag_velocity, ++it_frag_center) {

    for (UInt s = 0; s < spatial_dimension; ++s) {
      (*it_frag_velocity)(s) /= (*it_frag_mass)(s);
      (*it_frag_center)(s) /= (*it_frag_mass)(s);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::computeFragmentsData() {
  AKANTU_DEBUG_IN();

  buildFragments();
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
void SolidMechanicsModelCohesive::resizeFacetStress() {
  AKANTU_DEBUG_IN();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end();
       ++gt) {
    GhostType ghost_type = *gt;

    Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1, ghost_type);
    Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1, ghost_type);
    for(; it != end; ++it) {
      ElementType type = *it;

      UInt nb_facet = mesh_facets.getNbElement(type, ghost_type);

      UInt nb_quadrature_points =
	getFEM("FacetsFEM").getNbQuadraturePoints(type, ghost_type);

      UInt new_size = nb_facet * nb_quadrature_points;

      facet_stress(type, ghost_type).resize(new_size);
    }
  }

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
