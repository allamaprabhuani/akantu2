/**
 * @file   solid_mechanics_model_mass.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Oct  4 15:21:23 2010
 *
 * @brief  function handling mass computation
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
#include "solid_mechanics_model.hh"
#include "material.hh"

//#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMassLumped() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = mesh.getNbNodes();
  memset(mass->values, 0, nb_nodes*sizeof(Real));

  assembleMassLumped(_not_ghost);
  assembleMassLumped(_ghost);

  /// for not connected nodes put mass to one in order to avoid
  /// wrong range in paraview
  Real * mass_values = mass->values;
  for (UInt i = 0; i < nb_nodes; ++i) {
    if (fabs(mass_values[i]) < std::numeric_limits<Real>::epsilon() || Math::isnan(mass_values[i]))
      mass_values[i] = 1.;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMass() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(equation_number != NULL, "Call first initImplicitSolver()!");

  assembleMass(_not_ghost);
  //  assembleMass(_ghost);

  UInt nb_global_node = mesh.getNbGlobalNodes();

  std::stringstream sstr; sstr << id << ":mass_matrix";
  mass_matrix = new SparseMatrix(nb_global_node * spatial_dimension, _symmetric,
				 spatial_dimension, sstr.str(), memory_id);


  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMassLumped(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  FEM & fem = getFEM();

  Vector<Real> rho_1(0,1);

  const Mesh::ConnectivityTypeList & type_list = fem.getMesh().getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;
    ElementType type = *it;

    computeRho(rho_1, type, ghost_type);

    fem.assembleFieldLumped(rho_1, *mass, type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMass(GhostType ghost_type) {
  AKANTU_DEBUG_IN();


  MyFEMType & fem = getFEMClass<MyFEMType>();

  Vector<Real> rho_1(0,1);
  //UInt nb_element;

  const Mesh::ConnectivityTypeList & type_list = fem.getMesh().getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;
    ElementType type = *it;

    // const Vector<Real> * shapes;
    // if(ghost_type == _not_ghost) {
    //   nb_element   = fem.getMesh().getNbElement(type);
    //   shapes       = &(fem.getShapes(type));
    // } else {
    //   nb_element   = fem.getMesh().getNbGhostElement(type);
    //   shapes       = &(fem.getGhostShapes(type));
    // }

    computeRho(rho_1, type, ghost_type);
    fem.assembleFieldMatrix(rho_1, spatial_dimension, *equation_number, *mass_matrix, type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::computeRho(Vector<Real> & rho,
				     ElementType type,
				     GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Material ** mat_val = &(materials.at(0));
  UInt * elem_mat_val;

  UInt nb_element;

  FEM & fem = getFEM();

  if(ghost_type == _not_ghost) {
    nb_element   = fem.getMesh().getNbElement(type);
    elem_mat_val = element_material[type]->values;
  } else {
    nb_element   = fem.getMesh().getNbGhostElement(type);
    elem_mat_val = ghost_element_material[type]->values;
  }

  UInt nb_quadrature_points = fem.getNbQuadraturePoints(type);

  rho.resize(nb_element * nb_quadrature_points);
  Real * rho_1_val = rho.values;

  /// compute @f$ rho @f$ for each nodes of each element
  for (UInt el = 0; el < nb_element; ++el) {
    Real mat_rho = mat_val[elem_mat_val[el]]->getRho(); /// here rho is constant in an element
    for (UInt n = 0; n < nb_quadrature_points; ++n) {
      *rho_1_val++ = mat_rho;
    }
  }

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
