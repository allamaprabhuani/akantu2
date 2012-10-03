/**
 * @file   solid_mechanics_model_cohesive.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Thu Apr 19 10:42:11 2012
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
		id, memory_id),
    stress_on_facet("stress_on_facet", id),
    facet_stress(0, spatial_dimension * spatial_dimension, "facet_stress"),
    facets_to_cohesive_el(0, 2, "facets_to_cohesive_el"),
    fragment_to_element("fragment_to_element", id),
    fragment_velocity(0, spatial_dimension, "fragment_velocity"),
    fragment_center(0, spatial_dimension, "fragment_center") {
  AKANTU_DEBUG_IN();

  registerFEMObject<MyFEMCohesiveType>("CohesiveFEM", mesh, spatial_dimension);

  MeshUtils::buildAllFacets(mesh, mesh_facets);

  /// assign cohesive type
  Mesh::type_iterator it   = mesh_facets.firstType(spatial_dimension - 1);
  Mesh::type_iterator last = mesh_facets.lastType(spatial_dimension - 1);

  for (; it != last; ++it) {
    const Vector<UInt> & connectivity = mesh_facets.getConnectivity(*it);
    if (connectivity.getSize() != 0) {
      type_facet = *it;
      break;
    }
  }

  if (type_facet == _segment_2) type_cohesive = _cohesive_2d_4;
  else if (type_facet == _segment_3) type_cohesive = _cohesive_2d_6;

  mesh.addConnectivityType(type_cohesive);

  nb_fragment = 1;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::updateResidual(bool need_initialize) {
  AKANTU_DEBUG_IN();

  // f = f_ext
  if (need_initialize) initializeUpdateResidualData();

  // f -= fint
  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    Material & mat = **mat_it;
    mat.updateResidual(_not_ghost);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::initFull(std::string material_file,
					   AnalysisMethod method) {

  SolidMechanicsModel::initFull(material_file, method);

  if(material_file != "") {
    initCohesiveMaterial();
  }
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

  AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
/**
 * Initialize the model,basically it  pre-compute the shapes, shapes derivatives
 * and jacobian
 *
 */
void SolidMechanicsModelCohesive::initModel() {
  SolidMechanicsModel::initModel();
  getFEM("CohesiveFEM").initShapeFunctions(_not_ghost);
  getFEM("CohesiveFEM").initShapeFunctions(_ghost);
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::initExtrinsic() {
  AKANTU_DEBUG_IN();

  /// assign sigma_c to each facet
  MaterialCohesive * mat_cohesive
    = dynamic_cast<MaterialCohesive*>(materials[cohesive_index]);

  UInt nb_facet = mesh_facets.getNbElement(type_facet);
  sigma_lim.resize(nb_facet);

  mat_cohesive->generateRandomDistribution(sigma_lim);

  if (facets_check.getSize() < 1) {
    const Vector<Vector<Element> > & element_to_facet
      = mesh_facets.getElementToSubelement(type_facet);
    facets_check.resize(nb_facet);

    for (UInt f = 0; f < nb_facet; ++f) {
      if (element_to_facet(f)(1) != ElementNull) facets_check(f) = true;
      else facets_check(f) = false;
    }
  }

  /// compute normals on facets
  registerFEMObject<MyFEMType>("FacetsFEM", mesh_facets, spatial_dimension-1);
  getFEM("FacetsFEM").initShapeFunctions();
  getFEM("FacetsFEM").computeNormalsOnControlPoints();

  /**
   *  @todo store tangents while computing normals instead of
   *  recomputing them as follows:
   */
  /* ------------------------------------------------------------------------ */
  const Vector<Real> & normals = getFEM("FacetsFEM").getNormalsOnQuadPoints(type_facet);

  Vector<Real>::const_iterator<types::RVector> normal_it =
    normals.begin(spatial_dimension);

  tangents.resize(normals.getSize());
  tangents.extendComponentsInterlaced(spatial_dimension, tangents.getNbComponent());

  Vector<Real>::iterator<types::RVector> tangent_it =
    tangents.begin(spatial_dimension);

  for (UInt i = 0; i < normals.getSize(); ++i, ++normal_it, ++tangent_it) {
    Math::normal2( (*normal_it).storage(), (*tangent_it).storage() );
  }
  /* ------------------------------------------------------------------------ */

  /// compute quadrature points coordinates on facets
  Vector<Real> & position = mesh.getNodes();

  UInt nb_quad_per_facet = getFEM("FacetsFEM").getNbQuadraturePoints(type_facet);
  UInt nb_tot_quad = nb_quad_per_facet * nb_facet;

  Vector<Real> quad_facets(nb_tot_quad, spatial_dimension);

  getFEM("FacetsFEM").interpolateOnQuadraturePoints(position,
						    quad_facets,
						    spatial_dimension,
						    type_facet);


  /// compute elements quadrature point positions and build
  /// element-facet quadrature points data structure
  ByElementTypeReal elements_quad_facets;

  Mesh::type_iterator it   = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last = mesh.lastType(spatial_dimension);

  for (; it != last; ++it) {
    UInt nb_element = mesh.getNbElement(*it);
    if (nb_element == 0) continue;

    /// compute elements' quadrature points and list of facet
    /// quadrature points positions by element
    Vector<Element> & facet_to_element = mesh_facets.getSubelementToElement(*it);
    UInt nb_facet_per_elem = facet_to_element.getNbComponent();

    elements_quad_facets.alloc(nb_element * nb_facet_per_elem * nb_quad_per_facet,
			       spatial_dimension,
			       *it);

    Vector<Real> & el_q_facet = elements_quad_facets(*it);

    stress_on_facet.alloc(el_q_facet.getSize(),
			  spatial_dimension * spatial_dimension,
			  *it);

    for (UInt el = 0; el < nb_element; ++el) {
      for (UInt f = 0; f < nb_facet_per_elem; ++f) {
	UInt global_facet = facet_to_element(el, f).element;

	for (UInt q = 0; q < nb_quad_per_facet; ++q) {
	  for (UInt s = 0; s < spatial_dimension; ++s) {
	    el_q_facet(el * nb_facet_per_elem * nb_quad_per_facet
		       + f * nb_quad_per_facet + q, s)
	      = quad_facets(global_facet * nb_quad_per_facet + q, s);
	  }
	}
      }
    }
  }

  /// initialize the interpolation function
  materials[0]->initElementalFieldInterpolation(elements_quad_facets);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::checkCohesiveStress() {
  AKANTU_DEBUG_IN();

  Vector<UInt> facet_insertion;

  UInt nb_materials = materials.size();

  UInt nb_quad_per_facet = getFEM("FacetsFEM").getNbQuadraturePoints(type_facet);
  UInt nb_facet = mesh_facets.getNbElement(type_facet);

  /// vector containing stresses coming from the two elements of each facets
  facet_stress.resize(2 * nb_facet * nb_quad_per_facet);

  facet_stress_count.resize(nb_facet);
  facet_stress_count.clear();

  /// loop over materials
  for (UInt m = 0; m < nb_materials; ++m) {
    if (m == cohesive_index) continue;

    Mesh::type_iterator it   = mesh.firstType(spatial_dimension);
    Mesh::type_iterator last = mesh.lastType(spatial_dimension);

    /// loop over element type
    for (; it != last; ++it) {
      UInt nb_element = mesh.getNbElement(*it);
      if (nb_element == 0) continue;

      Vector<Real> & stress_on_f = stress_on_facet(*it);

      /// interpolate stress on facet quadrature points positions
      materials[m]->interpolateStress(*it,
				      stress_on_f);

      /// store the interpolated stresses on the facet_stress vector
      Vector<Element> & facet_to_element = mesh_facets.getSubelementToElement(*it);
      UInt nb_facet_per_elem = facet_to_element.getNbComponent();

      Vector<Element>::iterator<types::Vector<Element> > facet_to_el_it =
	facet_to_element.begin(nb_facet_per_elem);

      Vector<Real>::iterator<types::Matrix> stress_on_f_it =
	stress_on_f.begin(spatial_dimension, spatial_dimension);

      UInt sp2 = spatial_dimension * spatial_dimension;
      UInt nb_quad_f_two = nb_quad_per_facet * 2;

      for (UInt el = 0; el < nb_element; ++el, ++facet_to_el_it) {
	for (UInt f = 0; f < nb_facet_per_elem; ++f) {
	  UInt global_facet = (*facet_to_el_it)(f).element;

	  for (UInt q = 0; q < nb_quad_per_facet; ++q, ++stress_on_f_it) {

	    if (facets_check(global_facet) == true) {
	      types::Matrix facet_stress_local(facet_stress.storage()
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

  MaterialCohesive * mat_cohesive
    = dynamic_cast<MaterialCohesive*>(materials[cohesive_index]);

  mat_cohesive->checkInsertion(facet_stress, facet_insertion);

  if (facet_insertion.getSize() != 0)
    insertCohesiveElements(facet_insertion);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::insertCohesiveElements(const Vector<UInt> & facet_insertion) {
  AKANTU_DEBUG_IN();
  if(facet_insertion.getSize() == 0) return;

  Vector<UInt> facet_material(facet_insertion.getSize(), 1, cohesive_index);
  insertCohesiveElements(facet_insertion, facet_material);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::insertCohesiveElements(const Vector<UInt> & facet_insertion,
							 const Vector<UInt> & facet_material) {
  AKANTU_DEBUG_IN();

  if(facet_insertion.getSize() == 0) return;

  Vector<UInt> doubled_nodes(0, 2);
  Vector<UInt> doubled_facets(0, 2);

  /// update mesh
  for (UInt f = 0; f < facet_insertion.getSize(); ++f) {
    Element facet(type_facet, facet_insertion(f));
    doubleFacet(facet, doubled_nodes, doubled_facets);
  }

  /// double middle nodes if it's the case
  if (type_facet == _segment_3)
    doubleMiddleNode(doubled_nodes, doubled_facets);

  /// loop over doubled facets to insert cohesive elements
  Vector<UInt> & conn_cohesive = mesh.getConnectivity(type_cohesive);
  const Vector<UInt> & conn_facet = mesh_facets.getConnectivity(type_facet);
  UInt nb_nodes_per_facet = conn_facet.getNbComponent();
  const Vector<Real> & position = mesh.getNodes();
  Vector<UInt> & element_mat = getElementMaterial(type_cohesive);
  Vector<Vector<Element> > & element_to_facet
    = mesh_facets.getElementToSubelement(type_facet);

  const Real epsilon = std::numeric_limits<Real>::epsilon();

  for (UInt f = 0; f < doubled_facets.getSize(); ++f) {
    UInt nb_cohesive_elements = conn_cohesive.getSize();
    conn_cohesive.resize(nb_cohesive_elements + 1);

    UInt first_facet = doubled_facets(f, 0);
    UInt second_facet = doubled_facets(f, 1);

    /// copy first facet's connectivity
    for (UInt n = 0; n < nb_nodes_per_facet; ++n)
      conn_cohesive(nb_cohesive_elements, n) = conn_facet(first_facet, n);

    /// check if first nodes of the two facets are coincident or not
    UInt first_facet_node = conn_facet(first_facet, 0);
    UInt second_facet_node = conn_facet(second_facet, 0);
    bool second_facet_inversion = false;

    for (UInt dim = 0; dim < mesh.getSpatialDimension(); ++dim) {
      if (std::abs( (position(first_facet_node, dim) - position(second_facet_node, dim))
		 / position(second_facet_node, dim)) >= epsilon) {
	second_facet_inversion = true;
	break;
      }
    }

    /// if the two nodes are coincident second facet connectivity is
    /// normally copied, otherwise the copy is reverted
    if (!second_facet_inversion) {
      for (UInt n = 0; n < nb_nodes_per_facet; ++n)
	conn_cohesive(nb_cohesive_elements, n + nb_nodes_per_facet)
	  = conn_facet(second_facet, n);
    }
    else {
      for (UInt n = 0; n < nb_nodes_per_facet; ++n)
	conn_cohesive(nb_cohesive_elements, n + nb_nodes_per_facet)
	  = conn_facet(second_facet, nb_nodes_per_facet - n - 1);
    }

    /// assign the cohesive material to the new element
    element_mat.resize(nb_cohesive_elements + 1);
    element_mat(nb_cohesive_elements) = facet_material(f);
    materials[cohesive_index]->addElement(type_cohesive, nb_cohesive_elements, _not_ghost);

    /// update element_to_facet vectors
    Element cohesive_element(type_cohesive, nb_cohesive_elements,
			     _not_ghost, _ek_cohesive);
    element_to_facet(first_facet)(1) = cohesive_element;
    element_to_facet(second_facet)(1) = cohesive_element;

    /// update facets_to_cohesive_el vector
    facets_to_cohesive_el.resize(nb_cohesive_elements + 1);
    facets_to_cohesive_el(nb_cohesive_elements, 0) = first_facet;
    facets_to_cohesive_el(nb_cohesive_elements, 1) = second_facet;
  }

  /// update shape functions
  getFEM("CohesiveFEM").initShapeFunctions(_not_ghost);

  /// resize cohesive material vectors
  MaterialCohesive * mat_cohesive
    = dynamic_cast<MaterialCohesive*>(materials[cohesive_index]);
  AKANTU_DEBUG_ASSERT(mat_cohesive, "No cohesive materials in the materials vector");
  mat_cohesive->resizeCohesiveVectors();

  /// update nodal values
  updateDoubledNodes(doubled_nodes);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::doubleMiddleNode(Vector<UInt> & doubled_nodes,
						   const Vector<UInt> & doubled_facets) {

  AKANTU_DEBUG_IN();

  Vector<UInt> & conn_facet = mesh_facets.getConnectivity(type_facet);
  Vector<Real> & position = mesh.getNodes();

  Vector<Vector<Element> > & elem_to_facet
    = mesh_facets.getElementToSubelement(type_facet);

  for (UInt f = 0; f < doubled_facets.getSize(); ++f) {

    UInt facet_first = doubled_facets(f, 0);
    UInt facet_second = doubled_facets(f, 1);

    UInt new_node = position.getSize();

    /// store doubled nodes
    UInt nb_doubled_nodes = doubled_nodes.getSize();
    doubled_nodes.resize(nb_doubled_nodes + 1);

    UInt old_node = conn_facet(facet_first, 2);

    doubled_nodes(nb_doubled_nodes, 0) = old_node;
    doubled_nodes(nb_doubled_nodes, 1) = new_node;

    /// update position
    position.resize(new_node + 1);
    for (UInt dim = 0; dim < spatial_dimension; ++dim)
      position(new_node, dim) = position(old_node, dim);

    /// update facet connectivity
    conn_facet(facet_second, 2) = new_node;

    /// update element connectivity
    for (UInt el = 0; el < elem_to_facet(facet_second).getSize(); ++el) {
      const ElementType type_elem = elem_to_facet(facet_second)(el).type;
      if (type_elem != _not_defined) {
	UInt elem_global = elem_to_facet(facet_second)(el).element;
	Vector<UInt> & conn_elem = mesh.getConnectivity(type_elem);

	for (UInt n = 0; n < conn_elem.getNbComponent(); ++n) {
	  if (conn_elem(elem_global, n) == old_node)
	    conn_elem(elem_global, n) = new_node;
	}
      }
    }

  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::updateDoubledNodes(const Vector<UInt> & doubled_nodes) {

  AKANTU_DEBUG_IN();

  UInt nb_old_nodes = displacement->getSize();
  UInt nb_new_nodes = nb_old_nodes + doubled_nodes.getSize();
  displacement          ->resize(nb_new_nodes);
  velocity              ->resize(nb_new_nodes);
  acceleration          ->resize(nb_new_nodes);
  if (increment_acceleration)
    increment_acceleration->resize(nb_new_nodes);
  force                 ->resize(nb_new_nodes);
  residual              ->resize(nb_new_nodes);
  boundary              ->resize(nb_new_nodes);
  mass                  ->resize(nb_new_nodes);
  if (increment)
    increment->resize(nb_new_nodes);


  //initsolver();
  /**
   * @todo temporary patch, should be done in a better way that works
   * always (pbc, parallel, ...)
   **/
  delete dof_synchronizer;
  dof_synchronizer = new DOFSynchronizer(mesh, spatial_dimension);
  dof_synchronizer->initLocalDOFEquationNumbers();
  dof_synchronizer->initGlobalDOFEquationNumbers();

  if (method != _explicit_dynamic) {
    if(method == _implicit_dynamic)
      delete stiffness_matrix;

    delete jacobian_matrix;
    delete solver;

    SolverOptions  solver_options;
    initImplicit((method == _implicit_dynamic), solver_options);
  }

  for (UInt n = 0; n < doubled_nodes.getSize(); ++n) {

    UInt old_node = doubled_nodes(n, 0);
    UInt new_node = doubled_nodes(n, 1);

    for (UInt dim = 0; dim < displacement->getNbComponent(); ++dim) {
      (*displacement)(new_node, dim) = (*displacement)(old_node, dim);
    }

    for (UInt dim = 0; dim < velocity->getNbComponent(); ++dim) {
      (*velocity)(new_node, dim) = (*velocity)(old_node, dim);
    }

    for (UInt dim = 0; dim < acceleration->getNbComponent(); ++dim) {
      (*acceleration)(new_node, dim) = (*acceleration)(old_node, dim);
    }

    if (increment_acceleration) {
      for (UInt dim = 0; dim < increment_acceleration->getNbComponent(); ++dim) {
	(*increment_acceleration)(new_node, dim)
	  = (*increment_acceleration)(old_node, dim);
      }
    }
  }

  if (increment) {
    for (UInt n = 0; n < doubled_nodes.getSize(); ++n) {

      UInt old_node = doubled_nodes(n, 0);
      UInt new_node = doubled_nodes(n, 1);

      for (UInt dim = 0; dim < increment->getNbComponent(); ++dim)
	(*increment)(new_node, dim) = (*increment)(old_node, dim);
    }
  }

  assembleMassLumped();
  updateCurrentPosition();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::doubleFacet(Element & facet,
					      Vector<UInt> & doubled_nodes,
					      Vector<UInt> & doubled_facets) {
  AKANTU_DEBUG_IN();

  const UInt f_index = facet.element;
  const ElementType type_facet = facet.type;

  const ElementType type_subfacet = mesh.getFacetElementType(type_facet);
  const UInt nb_subfacet = mesh.getNbFacetsPerElement(type_facet);

  Vector<Vector<Element> > & facet_to_subfacet
    = mesh_facets.getElementToSubelement(type_subfacet);

  Vector<Vector<Element> > & element_to_facet
    = mesh_facets.getElementToSubelement(type_facet);

  Vector<Element> & subfacet_to_facet
    = mesh_facets.getSubelementToElement(type_facet);

  /// adding a new facet by copying original one

  /// create new connectivity
  Vector<UInt> & conn_facet = mesh_facets.getConnectivity(type_facet);
  UInt nb_facet = conn_facet.getSize();
  conn_facet.resize(nb_facet + 1);
  for (UInt n = 0; n < conn_facet.getNbComponent(); ++n)
    conn_facet(nb_facet, n) = conn_facet(f_index, n);

  /// update facets_check vector
  if (facets_check.getSize() > 0) {
    facets_check.resize(nb_facet + 1);
    facets_check(nb_facet) = false;
    facets_check(f_index) = false;
  }

  /// store doubled facets
  UInt nb_doubled_facets = doubled_facets.getSize();
  doubled_facets.resize(nb_doubled_facets + 1);
  doubled_facets(nb_doubled_facets, 0) = f_index;
  doubled_facets(nb_doubled_facets, 1) = nb_facet;

  /// update elements connected to facet
  Vector<Element> first_facet_list = element_to_facet(f_index);
  element_to_facet.push_back(first_facet_list);

  /// set new and original facets as boundary facets
  AKANTU_DEBUG_ASSERT(element_to_facet(f_index)(1) != ElementNull,
		      "can't double a facet on the boundary!");

  element_to_facet(f_index)(1) = ElementNull;

  element_to_facet(nb_facet)(0) = element_to_facet(nb_facet)(1);
  element_to_facet(nb_facet)(1) = ElementNull;

  /// update facet_to_element vector
  ElementType type = element_to_facet(nb_facet)(0).type;
  UInt el = element_to_facet(nb_facet)(0).element;
  Vector<Element> & facet_to_element = mesh_facets.getSubelementToElement(type);

  UInt i;
  for (i = 0; facet_to_element(el, i).element != f_index
	 && i <= facet_to_element.getNbComponent(); ++i);

  facet_to_element(el, i).element = nb_facet;

  /// create new links to subfacets and update list of facets
  /// connected to subfacets
  subfacet_to_facet.resize(nb_facet + 1);
  for (UInt sf = 0; sf < nb_subfacet; ++sf) {
    subfacet_to_facet(nb_facet, sf) = subfacet_to_facet(f_index, sf);

    UInt sf_index = subfacet_to_facet(f_index, sf).element;

    /// find index to start looping around facets connected to current
    /// subfacet
    UInt start = 0;
    UInt nb_connected_facets = facet_to_subfacet(sf_index).getSize();

    while (facet_to_subfacet(sf_index)(start).element != f_index
    	   && start <= facet_to_subfacet(sf_index).getSize()) ++start;

    /// add the new facet to the list next to the original one
    ++nb_connected_facets;
    facet_to_subfacet(sf_index).resize(nb_connected_facets);

    for (UInt f = nb_connected_facets - 1; f > start; --f) {
      facet_to_subfacet(sf_index)(f) = facet_to_subfacet(sf_index)(f - 1);
    }


    /// check if the new facet should be inserted before or after the
    /// original one: the second element connected to both original
    /// and new facet will be _not_defined, so I check if the first
    /// one is equal to one of the elements connected to the following
    /// facet in the facet_to_subfacet vector
    UInt f_start = facet_to_subfacet(sf_index)(start).element;
    UInt f_next;
    if (start + 2 == nb_connected_facets)
      f_next = facet_to_subfacet(sf_index)(0).element;
    else
      f_next = facet_to_subfacet(sf_index)(start + 2).element;

    if ((element_to_facet(f_start)(0) == element_to_facet(f_next)(0))
	|| ( element_to_facet(f_start)(0) == element_to_facet(f_next)(1)))
      facet_to_subfacet(sf_index)(start).element = nb_facet;
    else
      facet_to_subfacet(sf_index)(start + 1).element = nb_facet;

    /// loop on every facet connected to the current subfacet
    for (UInt f = start + 2; ; ++f) {

      /// reset f in order to continue looping from the beginning
      if (f == nb_connected_facets) f = 0;
      /// exit loop if it reaches the end
      if (f == start) break;

      /// if current loop facet is on the boundary, double subfacet
      UInt f_global = facet_to_subfacet(sf_index)(f).element;
      if (element_to_facet(f_global)(1).type == _not_defined ||
	  element_to_facet(f_global)(1).kind == _ek_cohesive) {
	doubleSubfacet(subfacet_to_facet(f_index, sf), start, f, doubled_nodes);
	break;
      }

    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::doubleSubfacet(const Element & subfacet,
						 UInt start,
						 UInt end,
						 Vector<UInt> & doubled_nodes) {
  AKANTU_DEBUG_IN();

  const UInt sf_index = subfacet.element;
  const ElementType type_subfacet = subfacet.type;

  Vector<Vector<Element> > & facet_to_subfacet
    = mesh_facets.getElementToSubelement(type_subfacet);
  UInt nb_subfacet = facet_to_subfacet.getSize();

  Vector<UInt> & conn_point = mesh_facets.getConnectivity(_point);
  Vector<Real> & position = mesh.getNodes();

  /// add the new subfacet
  if (spatial_dimension == 2) {

    UInt new_node = position.getSize();

    /// add new node in connectivity
    UInt new_subfacet = conn_point.getSize();
    conn_point.resize(new_subfacet + 1);
    conn_point(new_subfacet) = new_node;

    /// store doubled nodes
    UInt nb_doubled_nodes = doubled_nodes.getSize();
    doubled_nodes.resize(nb_doubled_nodes + 1);

    UInt old_node = doubled_nodes(nb_doubled_nodes, 0)
      = conn_point(sf_index);
    doubled_nodes(nb_doubled_nodes, 1) = new_node;

    /// update position
    position.resize(new_node + 1);
    for (UInt dim = 0; dim < spatial_dimension; ++dim)
      position(new_node, dim) = position(old_node, dim);
  }

  /// create a vector for the new subfacet in facet_to_subfacet
  facet_to_subfacet.resize(nb_subfacet + 1);

  UInt nb_connected_facets = facet_to_subfacet(sf_index).getSize();
  /// loop over facets from start to end
  for (UInt f = start + 1; ; ++f) {

    /// reset f in order to continue looping from the beginning
    if (f == nb_connected_facets) f = 0;

    UInt f_global = facet_to_subfacet(sf_index)(f).element;
    ElementType type_facet = facet_to_subfacet(sf_index)(f).type;
    Vector<UInt> & conn_facet = mesh_facets.getConnectivity(type_facet);

    UInt old_node = conn_point(sf_index);
    UInt new_node = conn_point(conn_point.getSize() - 1);

    /// update facet connectivity
    UInt i;
    for (i = 0; conn_facet(f_global, i) != old_node
	   && i <= conn_facet.getNbComponent(); ++i);
    conn_facet(f_global, i) = new_node;

    /// update element connectivity
    Vector<Vector<Element> > & elem_to_facet
      = mesh_facets.getElementToSubelement(type_facet);

    for (UInt el = 0; el < elem_to_facet(f_global).getSize(); ++el) {
      const ElementType type_elem = elem_to_facet(f_global)(el).type;
      if (type_elem != _not_defined) {
	UInt elem_global = elem_to_facet(f_global)(el).element;
	Vector<UInt> & conn_elem = mesh.getConnectivity(type_elem);

	for (UInt n = 0; n < conn_elem.getNbComponent(); ++n) {
	  if (conn_elem(elem_global, n) == old_node)
	    conn_elem(elem_global, n) = new_node;
	}
      }
    }

    /// update subfacet_to_facet vector
    Vector<Element> & subfacet_to_facet
      = mesh_facets.getSubelementToElement(type_facet);

    for (i = 0; subfacet_to_facet(f_global, i).element != sf_index
	   && i <= subfacet_to_facet.getNbComponent(); ++i);
    subfacet_to_facet(f_global, i).element = nb_subfacet;

    /// add current facet to facet_to_subfacet last position
    Element current_facet(type_facet, f_global);
    facet_to_subfacet(nb_subfacet).push_back(current_facet);

    /// exit loop if it reaches the end
    if (f == end) break;
  }

  /// rearrange the facets connected to the original subfacet and
  /// compute the new number of facets connected to it
  if (end < start) {
    for (UInt f = 0; f < start - end; ++f)
      facet_to_subfacet(sf_index)(f) = facet_to_subfacet(sf_index)(f + end + 1);

    nb_connected_facets = start - end;
  }
  else {
    for (UInt f = 1; f < nb_connected_facets - end; ++f)
      facet_to_subfacet(sf_index)(start + f) = facet_to_subfacet(sf_index)(end + f);

    nb_connected_facets -= end - start;
  }

  /// resize list of facets of the original subfacet
  facet_to_subfacet(sf_index).resize(nb_connected_facets);

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
      Vector<UInt> & frag_to_el = fragment_to_element(*it);

      for (UInt el = 0; el < nb_element; ++el)
	frag_to_el(el) = max;
    }
  }

  Vector<Vector<Element> > & element_to_facet
    = mesh_facets.getElementToSubelement(type_facet);

  nb_fragment = 0;
  it = mesh.firstType(spatial_dimension);

  MaterialCohesive * mat_cohesive
    = dynamic_cast<MaterialCohesive*>(materials[cohesive_index]);
  const Vector<Real> & damage = mat_cohesive->getDamage(type_cohesive);
  UInt nb_quad_cohesive = getFEM("CohesiveFEM").getNbQuadraturePoints(type_cohesive);

  Real epsilon = std::numeric_limits<Real>::epsilon();

  for (; it != last; ++it) {
    Vector<UInt> & checked_el = fragment_to_element(*it);
    UInt nb_element = checked_el.getSize();
    if (nb_element != 0) {

      Vector<Element> & facet_to_element = mesh_facets.getSubelementToElement(*it);
      UInt nb_facet_per_elem = facet_to_element.getNbComponent();

      /// loop on elements
      for (UInt el = 0; el < nb_element; ++el) {
	if (checked_el(el) == max) {
	  /// build fragment
	  ++nb_fragment;
	  checked_el(el) = nb_fragment - 1;
	  Vector<Element> elem_to_check;

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
		  next_el = element_to_facet(global_facet)(i);
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
		      = facets_to_cohesive_el(next_el.element, 0) == global_facet;
		    UInt other_facet
		      = facets_to_cohesive_el(next_el.element, other_facet_index);

		    /// get the other regualar element
		    next_el = element_to_facet(other_facet)(0);
		  }
		}

		/// if it exists, add it to the fragment list
		if (next_el != ElementNull) {
		  Vector<UInt> & checked_next_el = fragment_to_element(next_el.type);
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
  Vector<bool> checked_node(nb_nodes);
  checked_node.clear();

  const Vector<Real> & position = mesh.getNodes();

  Mesh::type_iterator it   = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last = mesh.lastType(spatial_dimension);

  for (; it != last; ++it) {
    UInt nb_element = mesh.getNbElement(*it);
    if (nb_element != 0) {
      const Vector<UInt> & frag_to_el = fragment_to_element(*it);
      const Vector<UInt> & connectivity = mesh.getConnectivity(*it);
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
