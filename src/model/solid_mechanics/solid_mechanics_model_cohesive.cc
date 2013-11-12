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
#include "material_cohesive.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */

SolidMechanicsModelCohesive::SolidMechanicsModelCohesive(Mesh & mesh,
                                                         UInt dim,
                                                         const ID & id,
                                                         const MemoryID & memory_id)
  : SolidMechanicsModel(mesh, dim, id, memory_id),
    mesh_facets(mesh.initMeshFacets()),
    facets_check("facets_check", id),
    tangents("tangents", id),
    stress_on_facet("stress_on_facet", id),
    facet_stress("facet_stress", id),
    fragment_to_element("fragment_to_element", id),
    fragment_velocity(0, spatial_dimension, "fragment_velocity"),
    fragment_center(0, spatial_dimension, "fragment_center"),
    facet_insertion("facet_insertion", id),
    cohesive_el_to_facet("cohesive_el_to_facet", id),
    facet_material("facet_material", id) {
  AKANTU_DEBUG_IN();

  facet_generated = false;
#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  facet_synchronizer = NULL;
  facet_stress_synchronizer = NULL;
  cohesive_distributed_synchronizer = NULL;
  global_connectivity = NULL;
  rank_to_element = NULL;
#endif

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
SolidMechanicsModelCohesive::~SolidMechanicsModelCohesive() {
  AKANTU_DEBUG_IN();

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  delete cohesive_distributed_synchronizer;
  delete facet_synchronizer;
  delete facet_stress_synchronizer;
  delete rank_to_element;
#endif

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

  if (extrinsic) {
    if (!facet_generated)
      MeshUtils::buildAllFacets(mesh, mesh_facets);
    initAutomaticInsertion();
  }

  mesh_facets.initByElementTypeArray(facets_check, 1, spatial_dimension - 1);
  if (extrinsic) {
    initFacetsCheck();

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
    if (facet_stress_synchronizer != NULL) {
      DataAccessor * data_accessor = this;
      facet_stress_synchronizer->updateFacetStressSynchronizer(facets_check,
							       *rank_to_element,
							       *data_accessor);
    }
#endif

  }

  if(material_file != "")
    initCohesiveMaterials(extrinsic, init_facet_filter);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::initCohesiveMaterials(bool extrinsic,
							bool init_facet_filter) {

  AKANTU_DEBUG_IN();

  /// find the first cohesive material
  UInt cohesive_index = 0;

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
      if (extrinsic) mat.initInsertionArrays(mesh_facets);
    } catch(std::bad_cast&) { }
  }

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

  if (init_facet_filter) initFacetFilter();

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
void SolidMechanicsModelCohesive::initAutomaticInsertion() {
  AKANTU_DEBUG_IN();

  MeshUtils::resetFacetToDouble(mesh_facets);

  /// initialize facet insertion array
  mesh_facets.initByElementTypeArray(facet_insertion, 1,
                                     spatial_dimension - 1,
                                     false,
                                     _ek_regular,
                                     true);

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

  registerFEMObject<MyFEMType>("FacetsFEM", mesh_facets, spatial_dimension - 1);
  getFEM("FacetsFEM").initShapeFunctions(_not_ghost);
  getFEM("FacetsFEM").initShapeFunctions(_ghost);

  mesh_facets.initByElementTypeArray(facet_stress,
				     2 * spatial_dimension * spatial_dimension,
				     spatial_dimension - 1);

  /// compute normals on facets
  computeNormals();

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

      stress_on_facet(type, elem_gt).resize(el_q_facet.getSize());

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
  UInt sp2 = spatial_dimension * spatial_dimension;
  UInt sp4 = sp2 * 2;

  resizeFacetArray(facet_stress);

  /// loop over materials
  for (UInt m = 0; m < nb_materials; ++m) {
    try {
      MaterialCohesive & mat __attribute__((unused)) =
        dynamic_cast<MaterialCohesive &>(*materials[m]);
    } catch(std::bad_cast&) {
      /// interpolate stress on facet quadrature points positions
      materials[m]->interpolateStress(stress_on_facet);

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

              f_check = &( facets_check(facet_type, facet_gt) );
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

  for (UInt m = 0; m < nb_materials; ++m) {
    try {
      MaterialCohesive & mat_cohesive = dynamic_cast<MaterialCohesive &>(*materials[m]);
      /// check which not ghost cohesive elements are to be created
      mat_cohesive.checkInsertion(facet_stress, facet_insertion);
    } catch(std::bad_cast&) { }
  }

  /// communicate data among processors
  synch_registry->synchronize(_gst_smmc_facets);

  /// insert cohesive elements
  MeshUtils::insertCohesiveElements(mesh,
                                    mesh_facets,
                                    facet_insertion,
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

      /// define facet material vector
      ElementType type_cohesive = FEM::getCohesiveElementType(type_facet);

      Array<UInt> & el_id_by_mat = element_index_by_material(type_cohesive,
                                                             gt_facet);

      UInt old_nb_cohesive_elements = el_id_by_mat.getSize();
      UInt total_nb_cohesive_elements = mesh.getNbElement(type_cohesive, gt_facet);
      UInt new_cohesive_elements
	= total_nb_cohesive_elements - old_nb_cohesive_elements;

      if (new_cohesive_elements == 0) continue;

      const Array<Element> & facet_to_coh_element
	= mesh.getSubelementToElement(type_cohesive, gt_facet);

      el_id_by_mat.resize(total_nb_cohesive_elements);

      Array<bool> & f_check = facets_check(type_facet, gt_facet);

      if (f_check.getSize() > 0)
        f_check.resize(f_check.getSize() + new_cohesive_elements);

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
      Array<UInt> & cohesive_el_to_f = cohesive_el_to_facet(type_facet, gt_facet);
#endif

      const Array<UInt> & facet_material_by_type = facet_material(type_facet, gt_facet);

      bool third_dimension = spatial_dimension == 3;

      for (UInt cel = old_nb_cohesive_elements;
	   cel < total_nb_cohesive_elements;
	   ++cel) {
        UInt old_facet = facet_to_coh_element(cel,  third_dimension).element;
        UInt new_facet = facet_to_coh_element(cel, !third_dimension).element;

        /// assign the cohesive material to each new cohesive element
	UInt material_id = facet_material_by_type(old_facet);

        el_id_by_mat(cel, 0) =
          materials[material_id]->addElement(type_cohesive,
					     cel,
					     gt_facet);

        el_id_by_mat(cel, 1) = material_id;

        /// update facets_check vector
        if (f_check.getSize() > 0) {
          f_check(old_facet) = false;
          f_check(new_facet) = false;
        }

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
        if (facet_synchronizer != NULL)
          cohesive_el_to_f(old_facet) = cel;
#endif
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

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  updateFacetSynchronizer();
    if (facet_stress_synchronizer != NULL) {
      DataAccessor * data_accessor = this;
      facet_stress_synchronizer->updateFacetStressSynchronizer(facets_check,
							       *rank_to_element,
							       *data_accessor);
    }
#endif
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

void SolidMechanicsModelCohesive::buildFragmentsList() {
  AKANTU_DEBUG_IN();

  UInt cohesive_index = 0;

  while ((dynamic_cast<MaterialCohesive *>(materials[cohesive_index]) == NULL)
         && cohesive_index <= materials.size())
    ++cohesive_index;

  /// find element type for _not_ghost facet
  Mesh::type_iterator it_f   = mesh_facets.firstType(spatial_dimension - 1);
  Mesh::type_iterator last_f = mesh_facets.lastType(spatial_dimension - 1);

  ElementType internal_facet_type = _not_defined;

  for (; it_f != last_f; ++it_f) {
    UInt nb_facet = mesh_facets.getNbElement(*it_f);
    if (nb_facet != 0) {
      internal_facet_type = *it_f;
      break;
    }
  }

  AKANTU_DEBUG_ASSERT(internal_facet_type != _not_defined, "No facets in the mesh");


  Mesh::type_iterator it   = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last = mesh.lastType(spatial_dimension);

  UInt max = std::numeric_limits<UInt>::max();

  for (; it != last; ++it) {
    UInt nb_element = mesh.getNbElement(*it);
    if (nb_element != 0) {
      /// initialize the list of checked elements and the list of
      /// elements to be checked
      fragment_to_element.alloc(nb_element, 1, *it);
      fragment_to_element(*it).set(max);
    }
  }

  Array< std::vector<Element> > & element_to_facet
    = mesh_facets.getElementToSubelement(internal_facet_type);

  ElementType type_cohesive = FEM::getCohesiveElementType(internal_facet_type);
  if (mesh.getNbElement(type_cohesive) == 0) return;
  const Array<Element> & f_to_cohesive_el
    = mesh.getSubelementToElement(type_cohesive);

  nb_fragment = 0;
  it = mesh.firstType(spatial_dimension);

  MaterialCohesive * mat_cohesive
    = dynamic_cast<MaterialCohesive*>(materials[cohesive_index]);
  const Array<Real> & damage = mat_cohesive->getDamage(type_cohesive);
  UInt nb_quad_cohesive = getFEM("CohesiveFEM").getNbQuadraturePoints(type_cohesive);

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
			 Math::are_float_equal(damage(next_el.element
						      * nb_quad_cohesive + q), 1)) ++q;

                  if (q == nb_quad_cohesive)
                    next_el = ElementNull;
                  else {
                    /// check which facet is the correct one
                    UInt other_facet_index
                      = f_to_cohesive_el(next_el.element, 0).element == global_facet;
                    UInt other_facet
                      = f_to_cohesive_el(next_el.element, other_facet_index).element;

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

void SolidMechanicsModelCohesive::initFacetsCheck() {
  AKANTU_DEBUG_IN();

  resizeFacetArray(facets_check, false);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType facet_gt = *gt;
    Mesh::type_iterator it   = mesh_facets.firstType(spatial_dimension - 1, facet_gt);
    Mesh::type_iterator last = mesh_facets.lastType(spatial_dimension - 1, facet_gt);

    for (; it != last; ++it) {
      ElementType facet_type = *it;

      Array<bool> & f_check = facets_check(facet_type, facet_gt);

      const Array< std::vector<Element> > & element_to_facet
	= mesh_facets.getElementToSubelement(facet_type, facet_gt);

      for (UInt f = 0; f < f_check.getSize(); ++f) {
	if (element_to_facet(f)[1] == ElementNull ||
	    (element_to_facet(f)[0].ghost_type == _ghost &&
	     element_to_facet(f)[1].ghost_type == _ghost)) {
	  f_check(f) = false;
	}
	else f_check(f) = true;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::enableFacetsCheckOnArea(const Array<Real> & limits) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(limits.getNbComponent() == 2,
		      "Number of components for limits array must be 2");

  AKANTU_DEBUG_ASSERT(limits.getSize() == spatial_dimension,
		      "Limits array size must be equal to spatial dimension");

  Real * bary_facet = new Real[spatial_dimension];

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end();
       ++gt) {
    GhostType ghost_type = *gt;

    Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1, ghost_type);
    Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1, ghost_type);
    for(; it != end; ++it) {
      ElementType type = *it;
      Array<bool> & f_check = facets_check(type, ghost_type);
      UInt nb_facet = mesh_facets.getNbElement(type, ghost_type);

      for (UInt f = 0; f < nb_facet; ++f) {
	if (f_check(f)) {

	  mesh_facets.getBarycenter(f, type, bary_facet, ghost_type);

	  UInt coord_in_limit = 0;

	  while (coord_in_limit < spatial_dimension &&
		 bary_facet[coord_in_limit] > limits(coord_in_limit, 0) &&
		 bary_facet[coord_in_limit] < limits(coord_in_limit, 1))
	    ++coord_in_limit;

	  if (coord_in_limit != spatial_dimension)
	    f_check(f) = false;

	}
      }
    }
  }

  delete [] bary_facet;

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
    if (facet_stress_synchronizer != NULL) {
      DataAccessor * data_accessor = this;
      facet_stress_synchronizer->updateFacetStressSynchronizer(facets_check,
							       *rank_to_element,
							       *data_accessor);
    }
#endif

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template<typename T>
void SolidMechanicsModelCohesive::resizeFacetArray(ByElementTypeArray<T> & by_el_type_array,
						   bool per_quadrature_point_data) const {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = 0;

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end();
       ++gt) {
    GhostType ghost_type = *gt;

    Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1, ghost_type);
    Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1, ghost_type);
    for(; it != end; ++it) {
      ElementType type = *it;

      UInt nb_facet = mesh_facets.getNbElement(type, ghost_type);

      if (per_quadrature_point_data)
	nb_quadrature_points =
	  getFEM("FacetsFEM").getNbQuadraturePoints(type, ghost_type);
      else nb_quadrature_points = 1;

      UInt new_size = nb_facet * nb_quadrature_points;

      Array<T> & vect = by_el_type_array(type, ghost_type);
      vect.resize(new_size);
    }
  }

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
