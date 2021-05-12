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
#include "non_linear_solver.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialDamageIterativeOrthotropic<spatial_dimension>::
    MaterialDamageIterativeOrthotropic(SolidMechanicsModel & model,
                                       const ID & id)
    : parent(model, id), nb_state_changes("nb_state_changes", *this),
      damage_prev_iteration("damage_prev_iteration", *this),
      in_tension("in_tension", *this) {
  this->registerParam("max_state_changes_allowed", max_state_changes_allowed,
                      Real(5), _pat_parsmod,
                      "How many times an element can change between tension "
                      "and compression stiffness");
  this->registerParam("contact", contact, true, _pat_parsmod,
                      "Account for contact by recovering initial stiffness in "
                      "direction orthogonal to crack");
  this->registerParam("iso_damage", iso_damage, false, _pat_parsmod,
                      "Reduce stiffness in all directions");
  this->registerParam("loading_test", loading_test, false, _pat_modifiable,
                      "Indicate to material that it's in the loading test");
}
/* -------------------------------------------------------------------------- */

template <UInt spatial_dimension>
void MaterialDamageIterativeOrthotropic<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  parent::initMaterial();
  this->nb_state_changes.setDefaultValue(0.);
  this->nb_state_changes.initialize(1);
  this->damage_prev_iteration.initialize(1);
  this->in_tension.setDefaultValue(true);
  this->in_tension.initialize(1);
  this->E = this->E1;
  this->nu = this->nu12;

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeOrthotropic<spatial_dimension>::computeStress(
    const ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  PlaneStressToolbox<spatial_dimension, MaterialThermal<spatial_dimension>>::
      computeStress(el_type, ghost_type);
  /// if in loading test stop updating stiffness after 1st iteration
  Int nb_iter = this->model.getDOFManager().getNonLinearSolver("static").get(
      "nb_iterations");

  if (not(nb_iter > 0 and this->loading_test)) {
    computeC(el_type, ghost_type);
  }

  auto C_it =
      this->C_field(el_type, ghost_type).begin(voigt_h::size, voigt_h::size);
  auto sigma_th_it = make_view(this->sigma_th(el_type, ghost_type)).begin();
  const auto & elem_filter = this->element_filter(el_type);
  auto & extra_vol = this->extra_volume(el_type);
  Array<Real> volumetric_strain(
      elem_filter.size() * this->fem.getNbIntegrationPoints(el_type), 1, 0.);
  auto vol_strain_it = volumetric_strain.begin();
  auto dam_it = this->damage(el_type, ghost_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
  /// compute stress according to anisotropic material law
  this->computeStressOnQuad(grad_u, sigma, *C_it, *sigma_th_it);
  /// compute volumetric strain in record it into extra_volume array
  if (*dam_it) {
    *vol_strain_it = grad_u.trace();
    *vol_strain_it *= int(*vol_strain_it > 0);
  }

  ++C_it;
  ++sigma_th_it;
  ++vol_strain_it;
  ++dam_it;
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  /// integrate volumetric strain across elements, subtract initial volume
  this->fem.integrate(volumetric_strain, extra_vol, 1, el_type, ghost_type,
                      elem_filter);

  this->computeNormalizedEquivalentStress(el_type, ghost_type);
  this->norm_max_equivalent_stress = 0;
  this->findMaxNormalizedEquivalentStress(el_type, ghost_type);

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
inline void MaterialDamageIterativeOrthotropic<spatial_dimension>::computeC(
    const ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam = this->damage(el_type, ghost_type).storage();
  Real * dam_prev_iter =
      this->damage_prev_iteration(el_type, ghost_type).storage();
  auto E1_it = this->E1_field(el_type, ghost_type).begin();
  auto E2_it = this->E2_field(el_type, ghost_type).begin();
  auto E3_it = this->E3_field(el_type, ghost_type).begin();
  auto nu12_it = this->nu12_field(el_type, ghost_type).begin();
  auto nu13_it = this->nu13_field(el_type, ghost_type).begin();
  auto nu23_it = this->nu23_field(el_type, ghost_type).begin();
  auto G12_it = this->G12_field(el_type, ghost_type).begin();
  auto G13_it = this->G13_field(el_type, ghost_type).begin();
  auto G23_it = this->G23_field(el_type, ghost_type).begin();
  auto Cprime_it = this->Cprime_field(el_type, ghost_type)
                       .begin(spatial_dimension * spatial_dimension,
                              spatial_dimension * spatial_dimension);
  auto C_it =
      this->C_field(el_type, ghost_type).begin(voigt_h::size, voigt_h::size);
  auto eigC_it = this->eigC_field(el_type, ghost_type).begin(voigt_h::size);
  auto dir_vecs_it = this->dir_vecs_field(el_type, ghost_type)
                         .begin(spatial_dimension, spatial_dimension);
  auto flick_it = this->nb_state_changes(el_type, ghost_type).begin();
  auto tension_it = this->in_tension(el_type, ghost_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  /// parameters reduction & update of C and Cprime only if damage changed
  // if (*dam != *dam_prev_iter) {
  /// reduce or recover elastic moduli due to damage
  updateElasticModuli(sigma, grad_u, *dam, *E1_it, *E2_it, *E3_it, *nu12_it,
                      *nu13_it, *nu23_it, *G12_it, *G13_it, *G23_it,
                      *dir_vecs_it, *flick_it, *tension_it);
  /// construct the stiffness matrix with update parameters
  this->updateInternalParametersOnQuad(
      *E1_it, *E2_it, *E3_it, *nu12_it, *nu13_it, *nu23_it, *G12_it, *G13_it,
      *G23_it, *Cprime_it, *C_it, *eigC_it, *dir_vecs_it);
  // }

  /// update damage at previous iteration value
  *dam_prev_iter = *dam;

  ++dam;
  ++dam_prev_iter;
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
  ++tension_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  AKANTU_DEBUG_OUT();
}
/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension>
inline void
MaterialDamageIterativeOrthotropic<spatial_dimension>::updateElasticModuli(
    Matrix<Real> & /*sigma*/, Matrix<Real> & grad_u, Real & dam, Real & E1,
    Real & /*E2*/, Real & /*E3*/, Real & nu12, Real & nu13, Real & /*nu23*/,
    Real & G12, Real & G13, Real & /*G23*/, Matrix<Real> & dir_vecs,
    Real & nb_flicks, bool & in_tension) {

  if (dam == 0.) {
    E1 = this->E1;
    nu12 = this->nu12;
    nu13 = this->nu13;
  } else {
    Matrix<Real> strain(spatial_dimension, spatial_dimension);
    strain = 0.5 * (grad_u + grad_u.transpose());

    Vector<Real> normal_to_crack(spatial_dimension);
    /// normals are stored row-wise
    for (auto i : arange(spatial_dimension)) {
      normal_to_crack(i) = dir_vecs(0, i);
    }
    auto def_on_crack = strain * normal_to_crack;
    auto def_normal_to_crack = normal_to_crack.dot(def_on_crack);
    // auto traction_on_crack = sigma * normal_to_crack;
    // auto stress_normal_to_crack = normal_to_crack.dot(traction_on_crack);
    /// coefficients of the cubic spline function

    // if (nb_flicks == this->max_state_changes_allowed) {
    //   _E1 = std::sqrt(this->E1 * this->E1 * (1 - dam));
    //   _nu12 = std::sqrt(this->nu12 * this->nu12 * (1 - dam));
    //   _nu13 = std::sqrt(this->nu13 * this->nu13 * (1 - dam));
    //   // _G12 = std::sqrt(this->G12 * this->G12 * (1 - dam));
    //   // _G13 = std::sqrt(this->G13 * this->G13 * (1 - dam));
    //   std::cout << "Max number of state changes" << std::endl;
    // } else {
    if (def_normal_to_crack >= 0) {
      if (not in_tension) {
        ++nb_flicks;
      }
      in_tension = true;
      E1 = this->E1 * (1 - dam);
      nu12 = this->nu12 * (1 - dam);
      nu13 = this->nu13 * (1 - dam);
      G12 = this->G12 * (1 - dam);
      G13 = this->G13 * (1 - dam);
    } else {
      if (in_tension) {
        ++nb_flicks;
      }
      in_tension = false;
      // _E1 = std::min(-def_normal_to_crack * 1e13 + this->E1 * (1 - dam),
      //                this->E1);
      // _nu12 = std::min(this->nu12 * _E1 / this->E1, this->nu12);
      // _nu13 = std::min(this->nu13 * _E1 / this->E1, this->nu13);
      E1 = this->E1;
      nu12 = this->nu12;
      nu13 = this->nu13;
      // _G12 = this->G12;
      // _G13 = this->G13;
    }
  }
  // }
}
/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension>
inline void
MaterialDamageIterativeOrthotropic<spatial_dimension>::reduceInternalParameters(
    Matrix<Real> & sigma, Real & dam, Real & E1, Real & E2, Real & /*E3*/,
    Real & nu12, Real & nu13, Real & nu23, Real & /*G12*/, Real & /*G13*/,
    Real & /*G23*/, Matrix<Real> & dir_vecs, Real & nb_flicks) {

  /// detect compression in the normal to crack plane direction
  Vector<Real> normal_to_crack(spatial_dimension);
  /// normals are stored row-wise
  for (auto i : arange(spatial_dimension)) {
    normal_to_crack(i) = dir_vecs(0, i);
  }
  auto traction_on_crack = sigma * normal_to_crack;
  auto stress_normal_to_crack = normal_to_crack.dot(traction_on_crack);

  /// elements who exceed max allowed number of state changes get the average
  /// elastic properties
  if (nb_flicks == this->max_state_changes_allowed) {
    E1 = std::sqrt(this->E1 * this->E1 * (1 - dam));
    nu12 = std::sqrt(this->nu12 * this->nu12 * (1 - dam));
    nu13 = std::sqrt(this->nu13 * this->nu13 * (1 - dam));
    // G12 = std::sqrt(this->G12 * this->G12 * (1 - dam));
    // G13 = std::sqrt(this->G13 * this->G13 * (1 - dam));
  } else {
    /// recover stiffness only when compressive stress is considerable
    if (this->contact && std::abs(stress_normal_to_crack) > this->E / 1e9 &&
        stress_normal_to_crack < 0) {
      E1 = this->E1;
      nu12 = this->nu12;
      nu13 = this->nu13;
      // G12 = this->G12;
      // G13 = this->G13;
      ++nb_flicks;
    } else {
      E1 = this->E1 * (1 - dam);
      nu12 = this->nu12 * (1 - dam);
      nu13 = this->nu13 * (1 - dam);
      // _G12 = this->G12 * (1 - dam);
      // _G13 = this->G13 * (1 - dam);
      // _G12 = std::sqrt(this->G12 * this->G12 * (1 - dam));
      // _G13 = std::sqrt(this->G13 * this->G13 * (1 - dam));

      if (this->iso_damage) {
        E2 = this->E2 * (1 - dam);
        nu23 = this->nu23 * (1 - dam);
      }
      /// update shear moduli (Stable orthotropic materials (Li & Barbic))
      // _G12 = std::sqrt(_E1 * _E2) / 2 / (1 + std::sqrt(_nu12 *
      // this->nu12)); _G13 = std::sqrt(_E1 * _E3) / 2 / (1 + std::sqrt(_nu13
      // * this->nu13));
    }
  }
}
/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension>
void MaterialDamageIterativeOrthotropic<
    spatial_dimension>::computeTangentModuli(ElementType el_type,
                                             Array<Real> & tangent_matrix,
                                             GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto release = this->model.getDisplacementRelease();
  if (release != this->last_displacement_release) {
    /// if in loading test stop updating stiffness after 1st iteration
    Int nb_iter = this->model.getDOFManager().getNonLinearSolver("static").get(
        "nb_iterations");
    if (not(nb_iter > 0 and this->loading_test)) {
      computeC(el_type, ghost_type);
    }
  }
  OrthotropicParent::computeTangentModuli(el_type, tangent_matrix, ghost_type);
  AKANTU_DEBUG_OUT();
}
/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension>
void MaterialDamageIterativeOrthotropic<spatial_dimension>::beforeSolveStep() {

  parent::beforeSolveStep();

  for (auto & el_type : this->element_filter.elementTypes(
           _all_dimensions, _not_ghost, _ek_not_defined)) {
    this->nb_state_changes(el_type, _not_ghost).zero();
  }
}
/* --------------------------------------------------------------------------
 */

} // namespace akantu
