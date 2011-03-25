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

  UInt nb_nodes = getFEM().getMesh().getNbNodes();
  memset(mass->values, 0, nb_nodes*sizeof(Real));

  assembleMassLumped(_not_ghost);
  assembleMassLumped(_ghost);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMassLumped(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Material ** mat_val = &(materials.at(0));

  UInt nb_element;
  UInt * elem_mat_val;

  ByElementTypeReal rho_1;

  const Mesh::ConnectivityTypeList & type_list = getFEM().getMesh().getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;

    if(ghost_type == _not_ghost) {
      nb_element   = getFEM().getMesh().getNbElement(type);
      elem_mat_val = element_material[type]->values;
    } else {
      nb_element   = getFEM().getMesh().getNbGhostElement(type);
      elem_mat_val = ghost_element_material[type]->values;
    }

    UInt nb_quadrature_points = getFEM().getNbQuadraturePoints(type);
    rho_1[type] = new Vector<Real>(nb_element * nb_quadrature_points, 1, "rho_x_1");
    Real * rho_1_val = rho_1[type]->values;

    /// compute @f$ rho @f$ for each nodes of each element
    for (UInt el = 0; el < nb_element; ++el) {
      Real rho = mat_val[elem_mat_val[el]]->getRho(); /// here rho is constant in an element
      for (UInt n = 0; n < nb_quadrature_points; ++n) {
	*rho_1_val++ = rho;
      }
    }
  }

  getFEM().assembleFieldLumped(rho_1, *mass, ghost_type);

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    delete rho_1[*it];
  }

  AKANTU_DEBUG_OUT();
}
__END_AKANTU__
