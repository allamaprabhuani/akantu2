/**
 * @file   material_iterative_stiffness_reduction.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Thu Feb 18 16:03:56 2016
 *
 * @brief  Implementation of material iterative stiffness reduction
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
#include "communicator.hh"
#include "material_iterative_stiffness_reduction.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, class IsotropicParent>
MaterialIterativeStiffnessReduction<spatial_dimension, IsotropicParent>::
    MaterialIterativeStiffnessReduction(SolidMechanicsModel & model,
                                        const ID & id)
    : IsotropicParent(model, id), eps_u("eps_u", *this),
      Sc_init("Sc_init", *this), D("D", *this), Gf(0.),
      crack_band_width(0.), reduction_constant(0.) {
  AKANTU_DEBUG_IN();

  this->registerParam("Gf", Gf, _pat_parsable | _pat_modifiable,
                      "fracture energy");
  this->registerParam("crack_band_width", crack_band_width,
                      _pat_parsable | _pat_modifiable, "crack_band_width");
  this->registerParam("reduction_constant", reduction_constant, 2.,
                      _pat_parsable | _pat_modifiable, "reduction constant");

  this->eps_u.initialize(1);
  this->Sc_init.initialize(1);
  this->D.initialize(1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, class IsotropicParent>
void MaterialIterativeStiffnessReduction<spatial_dimension,
                                         IsotropicParent>::initMaterial() {
  AKANTU_DEBUG_IN();
  IsotropicParent::initMaterial();
  this->Sc_init.copy(this->Sc);

  for (auto ghost_type : ghost_types) {
    /// loop over all types in the filter
    for (auto & el_type :
         this->element_filter.elementTypes(_ghost_type = ghost_type)) {
      /// get the stiffness on each quad point
      auto Sc_it = this->Sc(el_type, ghost_type).begin();
      /// get the tangent of the tensile softening on each quad point
      auto D_it = this->D(el_type, ghost_type).begin();
      auto D_end = this->D(el_type, ghost_type).end();
      /// get the ultimate strain on each quad
      auto eps_u_it = this->eps_u(el_type, ghost_type).begin();
      // compute the tangent and the ultimate strain for each quad
      for (; D_it != D_end; ++Sc_it, ++D_it, ++eps_u_it) {
        *eps_u_it = ((2. * this->Gf) / (*Sc_it * this->crack_band_width));
        *D_it = *(Sc_it) / ((*eps_u_it) - ((*Sc_it) / this->E));
      }
    }
  }
  AKANTU_DEBUG_OUT();
} // namespace akantu

/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension, class IsotropicParent>
UInt MaterialIterativeStiffnessReduction<spatial_dimension,
                                         IsotropicParent>::updateDamage() {
  UInt nb_damaged_elements = 0;

  if (this->norm_max_equivalent_stress >= 1.) {

    AKANTU_DEBUG_IN();

    /// update the damage only on non-ghosts elements! Doesn't make sense to
    /// update on ghost.
    GhostType ghost_type = _not_ghost;

    for (auto && el_type :
         this->element_filter.elementTypes(spatial_dimension, ghost_type)) {
      /// loop over all the elements

      /// get iterators on the needed internal fields
      auto equivalent_stress_it =
          this->equivalent_stress(el_type, ghost_type).begin();
      auto equivalent_stress_end =
          this->equivalent_stress(el_type, ghost_type).end();
      auto dam_it = this->damage(el_type, ghost_type).begin();
      auto reduction_it = this->reduction_step(el_type, ghost_type).begin();
      auto eps_u_it = this->eps_u(el_type, ghost_type).begin();
      auto Sc_it = this->Sc(el_type, ghost_type).begin();
      auto Sc_init_it = this->Sc_init(el_type, ghost_type).begin();
      auto D_it = this->D(el_type, ghost_type).begin();

      /// EXPERIMENTAL
      Real E_ef = this->E;
      Real dt = this->model.getTimeStep();
      Real E_sum = E_ef;

      if (IsotropicParent::getID().find("viscoelastic_maxwell") != std::string::npos) {
        E_ef = this->getParam("Einf");
        const Vector<Real> & Eta = this->getParam("Eta");
        const Vector<Real> & Ev = this->getParam("Ev");

        for (UInt k = 0; k < Eta.size(); ++k) {
          Real lambda = Eta(k) / Ev(k);
          Real exp_dt_lambda = exp(-dt / lambda);
          Real gamma = lambda * (1 - exp_dt_lambda) / dt;
          E_sum += Ev(k);

          if (Math::are_float_equal(exp_dt_lambda, 1)) {
            E_ef += Ev(k);
          } else {
            E_ef += std::min(gamma * Ev(k), Ev(k));
          }
        }
      }
      /// loop over all the quads of the given element type
      for (; equivalent_stress_it != equivalent_stress_end;
           ++equivalent_stress_it, ++dam_it, ++reduction_it, ++eps_u_it,
           ++Sc_it, ++Sc_init_it, ++D_it) {

        if (Material::getID().find("viscoelastic_maxwell") != std::string::npos) {
          *eps_u_it = 2. * this->Gf / (*Sc_init_it * this->crack_band_width) +
                      *Sc_init_it * (1 / E_ef - 1 / E_sum);
          *D_it = *(Sc_init_it) / ((*eps_u_it) - ((*Sc_init_it) / E_ef));
        }

        /// check if damage occurs
        if (*equivalent_stress_it >=
            (1 - this->dam_tolerance) * this->norm_max_equivalent_stress) {

          /// check if this element can still be damaged
          if (*reduction_it == this->max_reductions)
            continue;

          /// increment the counter of stiffness reduction steps
          *reduction_it += 1;

          if (*reduction_it == this->max_reductions)
            *dam_it = this->max_damage;
          else {
            /// update the damage on this quad
            *dam_it =
                1. - (1. / std::pow(this->reduction_constant, *reduction_it));
            /// update the stiffness on this quad
            *Sc_it = (*eps_u_it) * (1. - (*dam_it)) * E_ef * (*D_it) /
                     ((1. - (*dam_it)) * E_ef + (*D_it));
          }
          nb_damaged_elements += 1;
        }
      }
    }
  }

  // auto rve_model = dynamic_cast<SolidMechanicsModelRVE *>(&this->model);
  // if (rve_model == NULL) {
  const auto & comm = this->model.getMesh().getCommunicator();
  comm.allReduce(nb_damaged_elements, SynchronizerOperation::_sum);
  //}

  AKANTU_DEBUG_OUT();
  return nb_damaged_elements;
} // namespace akantu

/* --------------------------------------------------------------------------
 */

} // namespace akantu
