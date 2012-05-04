/**
 * @file   material_cohesive_bilinear.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Thu Feb 16 16:44:28 2012
 *
 * @brief  Bilinear cohesive constitutive law
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
#include "material_cohesive_bilinear.hh"
#include "solid_mechanics_model.hh"
#include "sparse_matrix.hh"
#include "dof_synchronizer.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MaterialCohesiveBilinear::MaterialCohesiveBilinear(Model & model, const ID & id) :
  MaterialCohesiveLinear(model,id) {
  AKANTU_DEBUG_IN();

  delta_0   = 0;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MaterialCohesiveBilinear::~MaterialCohesiveBilinear() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesiveBilinear::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialCohesiveLinear::initMaterial();

  /**
   * Recompute sigma_c as
   * @f$ {\sigma_c}_\textup{new} =
   * \frac{{\sigma_c}_\textup{old} \delta_c} {\delta_c - \delta_0} @f$
   */

  AKANTU_DEBUG_ASSERT(delta_c != delta_0, "Check your material.dat");

  sigma_c *= delta_c / (delta_c - delta_0);

  updateDeltaMax(_ghost);
  updateDeltaMax(_not_ghost);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesiveBilinear::resizeCohesiveVectors() {
  MaterialCohesiveLinear::resizeCohesiveVectors();
  updateDeltaMax(_ghost);
  updateDeltaMax(_not_ghost);
}

/* -------------------------------------------------------------------------- */
void MaterialCohesiveBilinear::updateDeltaMax(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & mesh = model->getFEM().getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type, _ek_cohesive);

  for(; it != last_type; ++it) {
    Vector<Real>::iterator<Real>delta_max_it =
      delta_max(*it, ghost_type).begin();

    Vector<Real>::iterator<Real>delta_max_end =
      delta_max(*it, ghost_type).end();

    for (; delta_max_it != delta_max_end; ++delta_max_it) {
      *delta_max_it = delta_0;
    }
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
bool MaterialCohesiveBilinear::setParam(const std::string & key, 
				       const std::string & value,
				       const ID & id) {
  std::stringstream sstr(value);
  if(key == "delta_0") { sstr >> delta_0; }
  else { return MaterialCohesiveLinear::setParam(key, value, id); }
  return true;
}

/* -------------------------------------------------------------------------- */
void MaterialCohesiveBilinear::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material<_cohesive_bilinear> [" << std::endl;
  stream << space << " + delta_0      : " << delta_0 << std::endl;
  MaterialCohesiveLinear::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */

__END_AKANTU__
