/**
 * @file   phasefield_exponential.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Fri Jun 19 2020
 * @date last modification: Wed Jun 23 2021
 *
 * @brief  Specialization of the phasefield law class for exponential type
 * law
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "phasefield_exponential.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
PhaseFieldExponential::PhaseFieldExponential(PhaseFieldModel & model,
                                             const ID & id)
    : PhaseField(model, id) {}

/* -------------------------------------------------------------------------- */
void PhaseFieldExponential::updateInternalParameters() {
  PhaseField::updateInternalParameters();

  Matrix<Real> d(spatial_dimension, spatial_dimension);
  d.eye(this->g_c * this->l0);
  damage_energy.set(d);
}

/* -------------------------------------------------------------------------- */
void PhaseFieldExponential::computeDrivingForce(const ElementType & el_type,
                                                GhostType ghost_type) {
  for (auto && tuple : zip(this->phi(el_type, ghost_type),
                           this->phi.previous(el_type, ghost_type),
                           this->driving_force(el_type, ghost_type),
                           this->damage_energy_density(el_type, ghost_type),
                           make_view(this->strain(el_type, ghost_type),
                                     spatial_dimension, spatial_dimension))) {
    computePhiOnQuad(std::get<4>(tuple), std::get<0>(tuple),
                     std::get<1>(tuple));
    computeDamageEnergyDensityOnQuad(std::get<0>(tuple), std::get<3>(tuple));
    computeDrivingForceOnQuad(std::get<0>(tuple), std::get<2>(tuple));
  }
}

  
/* -------------------------------------------------------------------------- */
void PhaseFieldExponential::computeDissipatedEnergy(ElementType el_type) {
  AKANTU_DEBUG_IN();

  for (auto && tuple : zip(this->dissipated_energy(el_type, _not_ghost),
			   this->damage_on_qpoints(el_type, _not_ghost),
			   make_view(this->gradd(el_type, _not_ghost),
				     spatial_dimension))) {

    this->computeDissipatedEnergyOnQuad(std::get<1>(tuple),
					std::get<2>(tuple),
					std::get<0>(tuple));
  }
  
  AKANTU_DEBUG_OUT();
}



/* -------------------------------------------------------------------------- */
void PhaseFieldExponential::computeDissipatedEnergyByElement(
      ElementType type, UInt index, Vector<Real> & edis_on_quad_points) {
  auto gradd_it = this->gradd(type).begin(spatial_dimension);
  auto gradd_end = this->gradd(type).begin(spatial_dimension);
  auto damage_it = this->damage_on_qpoints(type).begin();


  UInt nb_quadrature_points = fem.getNbIntegrationPoints(type);

  gradd_it += index * nb_quadrature_points;
  gradd_end += (index + 1) * nb_quadrature_points;
  damage_it += index * nb_quadrature_points;

  Real * edis_quad = edis_on_quad_points.storage();

  for (; gradd_it != gradd_end; ++gradd_it, ++damage_it, ++edis_quad) {
    this->computeDissipatedEnergyOnQuad(*damage_it,*gradd_it, *edis_quad);
  }
}
  


INSTANTIATE_PHASEFIELD(exponential, PhaseFieldExponential);

} // namespace akantu
