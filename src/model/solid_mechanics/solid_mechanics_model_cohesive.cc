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
#include <cstdlib>
#include <algorithm>
#include <ctime>

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
		id, memory_id) {
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

  MaterialCohesive * mat_cohesive
    = dynamic_cast<MaterialCohesive*>(materials[cohesive_index]);

  UInt nb_facet = mesh_facets.getNbElement(type_facet);
  sigma_lim.resize(nb_facet);
  const Real sigma_c = mat_cohesive->getSigmaC();
  const Real rand = mat_cohesive->getRandFactor();
  std::srand(time(NULL));

  const Vector<Vector<Element> > & element_to_facet
    = mesh_facets.getElementToSubelement(type_facet);
  facets_check.resize(nb_facet);

  for (UInt f = 0; f < nb_facet; ++f) {
    if (element_to_facet(f)(1) != ElementNull) {
      facets_check(f) = true;
      sigma_lim(f) = sigma_c * (1 + std::rand()/(Real)RAND_MAX * rand);
    }
    else {
      facets_check(f) = false;
      sigma_lim(f) = std::numeric_limits<Real>::max();
    }
  }

  registerFEMObject<MyFEMType>("FacetsFEM", mesh_facets, spatial_dimension-1);
  getFEM("FacetsFEM").initShapeFunctions();
  getFEM("FacetsFEM").computeNormalsOnControlPoints();

  /// THIS HAS TO BE CHANGED:
  const Vector<Real> & normals = getFEM("FacetsFEM").getNormalsOnQuadPoints(type_facet);

  Vector<Real>::const_iterator<types::RVector> normal_it =
    normals.begin(spatial_dimension);

  tangents.resize(normals.getSize());
  tangents.extendComponents(spatial_dimension);

  Vector<Real>::iterator<types::RVector> tangent_it =
    tangents.begin(spatial_dimension);

  for (UInt i = 0; i < normals.getSize(); ++i, ++normal_it, ++tangent_it) {
    Math::normal2( (*normal_it).storage(), (*tangent_it).storage() );
  }


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::checkCohesiveStress() {
  AKANTU_DEBUG_IN();

  Vector<UInt> facet_insertion;

  MaterialCohesive * mat_cohesive
    = dynamic_cast<MaterialCohesive*>(materials[cohesive_index]);

  mat_cohesive->checkInsertion(facet_insertion);

  if (facet_insertion.getSize() != 0)
    insertCohesiveElements(facet_insertion);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::insertCohesiveElements(const Vector<UInt> & facet_insertion) {

  AKANTU_DEBUG_IN();

  if(facet_insertion.getSize() == 0) return;

  Vector<UInt> doubled_nodes(0, 2);
  Vector<UInt> doubled_facets(0, 2);

  /// update mesh
  for (UInt f = 0; f < facet_insertion.getSize(); ++f) {
      Element facet;
      facet.element = facet_insertion(f);
      facet.type = type_facet;
      doubleFacet(facet, doubled_nodes, doubled_facets);
  }

  /// double middle nodes if it's the case
  if (type_facet == _segment_3) {
    doubleMiddleNode(doubled_nodes, doubled_facets);
  }

  /// update nodal values
  updateDoubledNodes(doubled_nodes);

  /// loop over doubled facets to insert cohesive elements
  Vector<UInt> & conn_cohesive = mesh.getConnectivity(type_cohesive);
  const Vector<UInt> & conn_facet = mesh_facets.getConnectivity(type_facet);
  UInt nb_nodes_per_facet = conn_facet.getNbComponent();
  const Vector<Real> & position = mesh.getNodes();
  Vector<UInt> & element_mat = getElementMaterial(type_cohesive);

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
      if (position(first_facet_node, dim) != position(second_facet_node, dim)) {
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
    element_mat(nb_cohesive_elements) = cohesive_index;
    materials[cohesive_index]->addElement(type_cohesive, nb_cohesive_elements, _not_ghost);
  }

  /// update shape functions
  getFEM("CohesiveFEM").initShapeFunctions(_not_ghost);

  /// resize cohesive material vectors
  MaterialCohesive * mat_cohesive
    = dynamic_cast<MaterialCohesive*>(materials[cohesive_index]);
  AKANTU_DEBUG_ASSERT(mat_cohesive, "No cohesive materials in the materials vector");
  mat_cohesive->resizeCohesiveVectors();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModelCohesive::doubleMiddleNode(Vector<UInt> & doubled_nodes, const Vector<UInt> & doubled_facets) {

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
  increment_acceleration->resize(nb_new_nodes);
  force                 ->resize(nb_new_nodes);
  residual              ->resize(nb_new_nodes);
  boundary              ->resize(nb_new_nodes);
  mass                  ->resize(nb_new_nodes);

  /**
   * @todo temporary patch, should be done in a better way that works
   * always (pbc, parallel, ...)
   **/
  delete dof_synchronizer;
  dof_synchronizer = new DOFSynchronizer(mesh, spatial_dimension);
  dof_synchronizer->initLocalDOFEquationNumbers();
  dof_synchronizer->initGlobalDOFEquationNumbers();

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

    for (UInt dim = 0; dim < increment_acceleration->getNbComponent(); ++dim) {
      (*increment_acceleration)(new_node, dim)
	= (*increment_acceleration)(old_node, dim);
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

  /// store doubled facets
  UInt nb_doubled_facets = doubled_facets.getSize();
  doubled_facets.resize(nb_doubled_facets + 1);
  doubled_facets(nb_doubled_facets, 0) = f_index;
  doubled_facets(nb_doubled_facets, 1) = nb_facet;

  /// update elements connected to facet
  element_to_facet.push_back(element_to_facet(f_index));

  /// set new and original facets as boundary facets
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
      if (element_to_facet(f_global)(1).type == _not_defined) {
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
    Element current_facet;
    current_facet.element = f_global;
    current_facet.type = type_facet;
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

Real SolidMechanicsModelCohesive::getReversibleEnergy() {
  AKANTU_DEBUG_IN();

  MaterialCohesive * mat_cohesive
    = dynamic_cast<MaterialCohesive*>(materials[cohesive_index]);

  return mat_cohesive->getReversibleEnergy();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

Real SolidMechanicsModelCohesive::getDissipatedEnergy() {
  AKANTU_DEBUG_IN();

  MaterialCohesive * mat_cohesive
    = dynamic_cast<MaterialCohesive*>(materials[cohesive_index]);

  return mat_cohesive->getDissipatedEnergy();

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
