/**
 * @file   solid_mechanics_model_mass.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Oct  4 15:21:23 2010
 *
 * @brief  function handling mass computation
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
#include "material.hh"

//#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMassLumped() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = fem->getMesh().getNbNodes();
  memset(mass->values, 0, nb_nodes*sizeof(Real));

  assembleMassLumped(_not_ghost);
  assembleMassLumped(_ghost);

  /// for not connected nodes put mass to one in order to avoid
  /// wrong range in paraview
  Real * mass_values = mass->values;
  for (UInt i = 0; i < nb_nodes; ++i) {
    if (fabs(mass_values[i]) < std::numeric_limits<Real>::epsilon() || isnan(mass_values[i]))
      mass_values[i] = 1;
  }


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMassLumped(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const Mesh::ConnectivityTypeList & type_list = fem->getMesh().getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;
    switch(*it) {
    case _triangle_6:
    case _tetrahedron_10:
      assembleMassLumpedDiagonalScaling(ghost_type, *it);
      break;
    default: assembleMassLumpedRowSum(ghost_type, *it);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * @f$ \tilde{M}_{i} = \sum_j M_{ij} = \sum_j \int \rho \varphi_i \varphi_j dV = \int \rho \varphi_i dV @f$
 */
void SolidMechanicsModel::assembleMassLumpedRowSum(GhostType ghost_type, ElementType type) {
  AKANTU_DEBUG_IN();
  Material ** mat_val = &(materials.at(0));

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

  Vector<Real> * rho_phi_i = new Vector<Real>(nb_element * nb_quadrature_points, shape_size, "rho_x_shapes");

  Real * rho_phi_i_val = rho_phi_i->values;
  Real * shapes_val = shapes->values;

  /// compute @f$ rho * \varphi_i @f$ for each nodes of each element
  for (UInt el = 0; el < nb_element; ++el) {
    Real rho = mat_val[elem_mat_val[el]]->getRho();
    for (UInt n = 0; n < shape_size * nb_quadrature_points; ++n) {
      *rho_phi_i_val++ = rho * *shapes_val++;
    }
  }

  Vector<Real> * int_rho_phi_i = new Vector<Real>(0, shape_size,
						  "inte_rho_x_shapes");
  fem->integrate(*rho_phi_i, *int_rho_phi_i, shape_size, type, ghost_type);
  delete rho_phi_i;

  fem->assembleVector(*int_rho_phi_i, *mass, 1, type, ghost_type);
  delete int_rho_phi_i;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
/**
 * @f$ \tilde{M}_{i} = c * M_{ii} = \int_{V_e} \rho dV @f$
 */
void SolidMechanicsModel::assembleMassLumpedDiagonalScaling(GhostType ghost_type, ElementType type) {
  AKANTU_DEBUG_IN();
  Material ** mat_val = &(materials.at(0));

  UInt nb_nodes_per_element_p1 = Mesh::getNbNodesPerElement(Mesh::getP1ElementType(type));
  UInt nb_nodes_per_element    = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points    = FEM::getNbQuadraturePoints(type);

  UInt nb_element;
  UInt * elem_mat_val;

  Real corner_factor = 0;
  Real mid_factor    = 0;
  switch(type){
  case _triangle_6 :
    corner_factor = 1./12.;
    mid_factor    = 1./4.;
    break;
  case _tetrahedron_10:
    corner_factor = 1./32.;
    mid_factor    = 7./48.;
    break;
  default:
    corner_factor = 0;
    mid_factor    = 0;
  }

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

  Vector<Real> * rho_1 = new Vector<Real>(nb_element * nb_quadrature_points, 1, "rho_x_1");

  Real * rho_1_val = rho_1->values;

  /// compute @f$ rho @f$ for each nodes of each element
  for (UInt el = 0; el < nb_element; ++el) {
    Real rho = mat_val[elem_mat_val[el]]->getRho(); /// here rho is constant in an element
    for (UInt n = 0; n < nb_quadrature_points; ++n) {
      *rho_1_val++ = rho;
    }
  }

  /// compute @f$ \int \rho dV = \rho V @f$ for each element
  Vector<Real> * int_rho_1 = new Vector<Real>(nb_element * nb_quadrature_points, 1,
					      "inte_rho_x_1");
  fem->integrate(*rho_1, *int_rho_1, 1, type, ghost_type);
  delete rho_1;

  /// distribute the mass of the element to the nodes
  Vector<Real> * mass_per_node = new Vector<Real>(nb_element, nb_nodes_per_element, "mass_per_node");
  Real * int_rho_1_val = int_rho_1->values;
  Real * mass_per_node_val = mass_per_node->values;

  for (UInt e = 0; e < nb_element; ++e) {
    Real lmass = *int_rho_1_val * corner_factor;
    for (UInt n = 0; n < nb_nodes_per_element_p1; ++n)
      *mass_per_node_val++ = lmass; /// corner points

    lmass = *int_rho_1_val * mid_factor;
    for (UInt n = nb_nodes_per_element_p1; n < nb_nodes_per_element; ++n)
      *mass_per_node_val++ = lmass; /// mid points

    int_rho_1_val++;
  }
  delete int_rho_1;

  fem->assembleVector(*mass_per_node, *mass, 1, type, ghost_type);
  delete mass_per_node;

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
