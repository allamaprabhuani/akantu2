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
#include "material_damage_iterative_orthotropic.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialDamageIterativeOrthotropic<spatial_dimension>::
    MaterialDamageIterativeOrthotropic(SolidMechanicsModel & model,
                                       const ID & id)
    : parent(model, id), nb_state_changes("nb_state_changes", *this), E(0.) {
  this->registerParam("max_state_changes_allowed", max_state_changes_allowed,
                      UInt(5), _pat_parsmod,
                      "How many times an element can change between tension "
                      "and compression stiffness");
  this->registerParam("contact", contact, true, _pat_parsmod,
                      "Account for contact by recovering initial stiffness in "
                      "direction orthogonal to crack");
  this->registerParam("iso_damage", iso_damage, false, _pat_parsmod,
                      "Reduce stiffness in all directions");
}
/* -------------------------------------------------------------------------- */

template <UInt spatial_dimension>
void MaterialDamageIterativeOrthotropic<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  parent::initMaterial();
  this->nb_state_changes.initialize(1);
  this->E = this->E1;
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeOrthotropic<spatial_dimension>::computeStress(
    const ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  parent::computeStress(el_type, ghost_type);

  Real * dam = this->damage(el_type, ghost_type).storage();
  auto E1_it =
      this->template getInternal<Real>("E1_field")(el_type, ghost_type).begin();
  auto E2_it =
      this->template getInternal<Real>("E2_field")(el_type, ghost_type).begin();
  auto E3_it =
      this->template getInternal<Real>("E3_field")(el_type, ghost_type).begin();
  auto nu12_it =
      this->template getInternal<Real>("nu12_field")(el_type, ghost_type)
          .begin();
  auto nu13_it =
      this->template getInternal<Real>("nu13_field")(el_type, ghost_type)
          .begin();
  auto nu23_it =
      this->template getInternal<Real>("nu23_field")(el_type, ghost_type)
          .begin();
  auto G12_it =
      this->template getInternal<Real>("G12_field")(el_type, ghost_type)
          .begin();
  auto G13_it =
      this->template getInternal<Real>("G13_field")(el_type, ghost_type)
          .begin();
  auto G23_it =
      this->template getInternal<Real>("G23_field")(el_type, ghost_type)
          .begin();
  auto Cprime_it =
      this->template getInternal<Real>("Cprime_field")(el_type, ghost_type)
          .begin(spatial_dimension * spatial_dimension,
                 spatial_dimension * spatial_dimension);
  auto C_it = this->template getInternal<Real>("C_field")(el_type, ghost_type)
                  .begin(voigt_h::size, voigt_h::size);
  auto eigC_it =
      this->template getInternal<Real>("eigC_field")(el_type, ghost_type)
          .begin(voigt_h::size);
  auto dir_vecs_it =
      this->template getInternal<Real>("dir_vecs_field")(el_type, ghost_type)
          .begin(spatial_dimension, spatial_dimension);
  auto flick_it = this->nb_state_changes(el_type, ghost_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  computeOrthotropicStress(sigma, grad_u, *dam, *E1_it, *E2_it, *E3_it,
                           *nu12_it, *nu13_it, *nu23_it, *G12_it, *G13_it,
                           *G23_it, *Cprime_it, *C_it, *eigC_it, *dir_vecs_it,
                           *flick_it);

  ++dam;
  ++E1_it;
  ++E2_it;
  ++E3_it;
  ++nu12_it;
  ++nu13_it;
  ++nu23_it;
  ++G12_it;
  ++G13_it;
  ++G23_it;
  ++Cprime_it;
  ++C_it;
  ++eigC_it;
  ++dir_vecs_it;
  ++flick_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  this->computeNormalizedEquivalentStress(this->gradu(el_type, ghost_type),
                                          el_type, ghost_type);
  this->norm_max_equivalent_stress = 0;
  this->findMaxNormalizedEquivalentStress(el_type, ghost_type);

  AKANTU_DEBUG_OUT();
}
/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension>
inline void
MaterialDamageIterativeOrthotropic<spatial_dimension>::computeOrthotropicStress(
    Matrix<Real> & sigma, Matrix<Real> & grad_u, Real & dam, Real & _E1,
    Real & _E2, Real & _E3, Real & _nu12, Real & _nu13, Real & _nu23,
    Real & _G12, Real & _G13, Real & _G23, Matrix<Real> & _Cprime,
    Matrix<Real> & _C, Vector<Real> & _eigC, Matrix<Real> & _dir_vecs,
    UInt & nb_flicks) {

  /// parameters reduction & update of C and Cprime only in case of damage
  if (dam > 0) {

    /// detect compression in the normal to crack plane direction
    Vector<Real> normal_to_crack(spatial_dimension);
    /// normals are stored row-wise
    for (auto i : arange(spatial_dimension))
      normal_to_crack(i) = _dir_vecs(0, i);
    auto traction_on_crack = sigma * normal_to_crack;
    auto stress_normal_to_crack = normal_to_crack.dot(traction_on_crack);

    /// elements who exceed max allowed number of state changes get the average
    /// elastic properties
    if (nb_flicks == this->max_state_changes_allowed) {
      _E1 = std::sqrt(this->E1 * this->E1 * (1 - dam));
      _nu12 = std::sqrt(this->nu12 * this->nu12 * (1 - dam));
      _nu13 = std::sqrt(this->nu13 * this->nu13 * (1 - dam));
      // _E1 = this->E1 * (1 - dam);
      // _nu12 = this->nu12 * (1 - dam);
      // _nu13 = this->nu13 * (1 - dam);
    } else {
      /// recover stiffness only when compressive stress is considerable
      if (this->contact && std::abs(stress_normal_to_crack) > this->E / 1e9 &&
          stress_normal_to_crack < 0) {
        _E1 = this->E1;
        _nu12 = this->nu12;
        _nu13 = this->nu13;
        // _G12 = this->G12;
        // _G13 = this->G13;
        ++nb_flicks;
      } else {
        _E1 = this->E1 * (1 - dam);
        _nu12 = this->nu12 * (1 - dam);
        _nu13 = this->nu13 * (1 - dam);

        if (this->iso_damage) {
          _E2 = this->E2 * (1 - dam);
          _nu23 = this->nu23 * (1 - dam);
        }
        /// update shear moduli (Stable orthotropic materials (Li & Barbic))
        // _G12 = std::sqrt(_E1 * _E2) / 2 / (1 + std::sqrt(_nu12 *
        // this->nu12)); _G13 = std::sqrt(_E1 * _E3) / 2 / (1 + std::sqrt(_nu13
        // * this->nu13));
        _G12 = this->G12 * (1 - dam);
        _G13 = this->G13 * (1 - dam);
        // _G12 = 0.;
        // _G13 = 0.;
      }
    }
    // calling on quad function of mat_orthotropic_heterogeneous
    this->updateInternalParametersOnQuad(_E1, _E2, _E3, _nu12, _nu13, _nu23,
                                         _G12, _G13, _G23, _Cprime, _C, _eigC,
                                         _dir_vecs);
  }
  // calling on quad function of mat_anisotropic_heterogeneous
  this->computeStressOnQuad(grad_u, sigma, _C);
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeOrthotropic<
    spatial_dimension>::computeTangentModuli(const ElementType & el_type,
                                             Array<Real> & tangent_matrix,
                                             GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  OrthotropicParent::computeTangentModuli(el_type, tangent_matrix, ghost_type);

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeOrthotropic<spatial_dimension>::beforeSolveStep() {

  parent::beforeSolveStep();

  for (auto & el_type : this->element_filter.elementTypes(
           _all_dimensions, _not_ghost, _ek_not_defined)) {
    this->nb_state_changes(el_type, _not_ghost).clear();
  }
}
/* --------------------------------------------------------------------------
 */

} // namespace akantu
