/**
 * @file   solid_mechanics_model_mass.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Oct  4 15:21:23 2010
 *
 * @brief  function handling mass computation
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
#include "material.hh"

//#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMassDiagonal() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = fem->getMesh().getNbNodes();
  memset(mass->values, 0, nb_nodes*sizeof(Real));

  assembleMassDiagonal(_not_ghost);
  assembleMassDiagonal(_ghost);

  /// for not connected nodes put mass to one in order to avoid
  /// wrong range in paraview
  Real * mass_values = mass->values;
  for (UInt i = 0; i < nb_nodes; ++i) {
    if (!mass_values[i] || isnan(mass_values[i]))
      mass_values[i] = 1;
  }


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMassDiagonal(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const Mesh::ConnectivityTypeList & type_list = fem->getMesh().getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;
    switch(*it) {
    case _triangle_2:
      assembleMassDiagonalTriangle2(ghost_type);
      break;
    default: assembleMassDiagonalGeneric(ghost_type, *it);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMassDiagonalGeneric(GhostType ghost_type, ElementType type) {
  AKANTU_DEBUG_IN();
  Material ** mat_val = &(materials.at(0));

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points = FEM::getNbQuadraturePoints(type);
  UInt shape_size           = FEM::getShapeSize(type);

  UInt nb_element;
  const Vector<Real> * shapes;
  UInt * elem_mat_val;

  if(ghost_type == _not_ghost) {
    nb_element   = fem->getMesh().getNbElement(type);
    shapes       = &(fem->getShapes(type));
    elem_mat_val = element_material[type]->values;
  } else {
    nb_element   = fem->getMesh().getNbGhostElement(type);
    shapes       = &(fem->getGhostShapes(type));
    elem_mat_val = ghost_element_material[type]->values;
  }

  if (nb_element == 0) {
    AKANTU_DEBUG_OUT();
    return;
  }

  Vector<Real> * rho_phi_i = new Vector<Real>(nb_element, shape_size, "rho_x_shapes");

  Real * rho_phi_i_val = rho_phi_i->values;
  Real * shapes_val = shapes->values;

  /// compute @f$ rho * \phi_i @f$ for each nodes of each element
  for (UInt el = 0; el < nb_element; ++el) {
    Real rho = mat_val[elem_mat_val[el]]->getRho();
    for (UInt n = 0; n < shape_size; ++n) {
      *rho_phi_i_val++ = rho * *shapes_val++;
    }
  }

  Vector<Real> * int_rho_phi_i = new Vector<Real>(nb_element, shape_size / nb_quadrature_points,
						  "inte_rho_x_shapes");
  fem->integrate(*rho_phi_i, *int_rho_phi_i, nb_nodes_per_element, type, ghost_type);
  delete rho_phi_i;

  fem->assembleVector(*int_rho_phi_i, *mass, 1, type, ghost_type);
  delete int_rho_phi_i;
 
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMassDiagonalTriangle2(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  Material ** mat_val = &(materials.at(0));

  ElementType type = _triangle_2;
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points = FEM::getNbQuadraturePoints(type);
  UInt shape_size           = FEM::getShapeSize(type);

  UInt nb_element;
  const Vector<Real> * shapes;
  UInt * elem_mat_val;

  if(ghost_type == _not_ghost) {
    nb_element   = fem->getMesh().getNbElement(type);
    elem_mat_val = element_material[type]->values;
  } else {
    nb_element   = fem->getMesh().getNbGhostElement(type);
    elem_mat_val = ghost_element_material[type]->values;
  }

  if (nb_element == 0) {
    AKANTU_DEBUG_OUT();
    return;
  }

  Vector<Real> * rho_1 = new Vector<Real>(nb_element, nb_quadrature_points, "rho_x_1");

  Real * rho_1_val = rho_1->values;

  /// compute @f$ rho @f$ for each nodes of each element
  for (UInt el = 0; el < nb_element; ++el) {
    Real rho = mat_val[elem_mat_val[el]]->getRho();
    for (UInt n = 0; n < nb_quadrature_points; ++n) {
      *rho_1_val++ = rho;
    }
  }

  /// compute @f$ \int \rho dV @f$ for each element
  Vector<Real> * int_rho_1 = new Vector<Real>(nb_element, 1,
					      "inte_rho_x_1");
  fem->integrate(*rho_1, *int_rho_1, 1, type, ghost_type);
  delete rho_1;

  /// distribute the mass of the element to the nodes
  Vector<Real> * mass_per_node = new Vector<Real>(nb_element, nb_nodes_per_element, "mass_per_node");
  Real * int_rho_1_val = int_rho_1->values;
  Real * mass_per_node_val = mass_per_node->values;
  for (UInt e = 0; e < nb_element; ++e) {
    Real lmass = *int_rho_1_val / 12;
    for (UInt n = 0; n < 3; ++n) *mass_per_node_val++ = lmass; /// corner points 1/12 of the element mass

    lmass = *int_rho_1_val / 4;
    for (UInt n = 3; n < 6; ++n) *mass_per_node_val++ = lmass; /// mid points 1/4 of the element mass

    int_rho_1_val++;
  }
  delete int_rho_1;

  fem->assembleVector(*mass_per_node, *mass, 1, type, ghost_type);
  delete mass_per_node;

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
