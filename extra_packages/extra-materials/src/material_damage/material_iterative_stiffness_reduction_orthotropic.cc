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
#include "material_iterative_stiffness_reduction_orthotropic.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

template <UInt spatial_dimension>
MaterialIterativeStiffnessReductionOrthotropic<spatial_dimension>::
    MaterialIterativeStiffnessReductionOrthotropic(SolidMechanicsModel & model,
                                                   const ID & id)
    : parent(model, id) {}

/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension>
UInt MaterialIterativeStiffnessReductionOrthotropic<
    spatial_dimension>::updateDamage() {
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
          this->template getInternal<Real>("equivalent_stress")(el_type,
                                                                ghost_type)
              .begin();
      auto equivalent_stress_end =
          this->template getInternal<Real>("equivalent_stress")(el_type,
                                                                ghost_type)
              .end();
      auto dam_it =
          this->template getInternal<Real>("damage")(el_type, ghost_type)
              .begin();
      auto reduction_it = this->template getInternal<UInt>("reduction_step")(
                                  el_type, ghost_type)
                              .begin();
      auto eps_u_it =
          this->template getInternal<Real>("eps_u")(el_type, ghost_type)
              .begin();

      auto Sc_it =
          this->template getInternal<Real>("Sc")(el_type, ghost_type).begin();
      auto D_it =
          this->template getInternal<Real>("D")(el_type, ghost_type).begin();
      auto crack_norm_it = this->crack_normals(el_type, ghost_type)
                               .begin(spatial_dimension, spatial_dimension);
      auto dir_vecs_it = this->template getInternal<Real>("dir_vecs_field")(
                                 el_type, ghost_type)
                             .begin(spatial_dimension, spatial_dimension);

      /// loop over all the quads of the given element type
      for (; equivalent_stress_it != equivalent_stress_end;
           ++equivalent_stress_it, ++dam_it, ++reduction_it, ++eps_u_it,
           ++Sc_it, ++D_it, ++crack_norm_it, ++dir_vecs_it) {

        /// check if damage occurs
        if (*equivalent_stress_it >=
            (1 - this->dam_tolerance) * this->norm_max_equivalent_stress) {

          /// if element is damaged first time -> note down the directions
          if (*reduction_it == 0) {
            *dir_vecs_it = *crack_norm_it;
          }
          /// check if this element can still be damaged
          if (*reduction_it == this->max_reductions)
            continue;

          /// increment the counter of stiffness reduction steps
          *reduction_it += 1;

          if (*reduction_it == this->max_reductions)
            *dam_it = this->max_damage;
          else {
            /// update the damage on this quad
            Real red_cons = this->getParam("reduction_constant");
            *dam_it = 1. - (1. / std::pow(red_cons, *reduction_it));
            /// update the stiffness on this quad
            *Sc_it = (*eps_u_it) * (1. - (*dam_it)) * this->E * (*D_it) /
                     ((1. - (*dam_it)) * this->E + (*D_it));
          }
          nb_damaged_elements += 1;
        }
      }
    }
  }

  const auto & comm = this->model.getMesh().getCommunicator();
  comm.allReduce(nb_damaged_elements, SynchronizerOperation::_sum);

  AKANTU_DEBUG_OUT();
  return nb_damaged_elements;
}
/* --------------------------------------------------------------------------
 */

INSTANTIATE_MATERIAL(iterative_stiffness_reduction_orthotropic,
                     MaterialIterativeStiffnessReductionOrthotropic);

} // namespace akantu
