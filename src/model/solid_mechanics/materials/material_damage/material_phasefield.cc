/**
 * @file   material_phasefield.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Fri Apr 02 2021
 *
 * @brief  Specialization of the material class for the phasefield material
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_phasefield.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialPhaseField<dim>::MaterialPhaseField(SolidMechanicsModel & model,
                                            const ID & id)
    : Parent(model, id) {

  AKANTU_DEBUG_IN();

  this->registerParam("eta", eta, Real(0.), _pat_parsable, "eta");
  this->damage.initialize(0);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialPhaseField<dim>::computeStress(ElementType el_type,
                                            GhostType ghost_type) {
  for (auto && args : Parent::getArguments(el_type, ghost_type)) {
    computeStressOnQuad(args);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialPhaseField<dim>::computeTangentModuli(ElementType el_type,
                                                   Array<Real> & tangent_matrix,
                                                   GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  for (auto && args :
       Parent::getArgumentsTangent(tangent_matrix, el_type, ghost_type)) {
    computeTangentModuliOnQuad(args);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(phasefield, MaterialPhaseField);

} // namespace akantu
