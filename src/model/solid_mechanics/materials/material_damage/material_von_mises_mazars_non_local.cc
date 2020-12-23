

/* -------------------------------------------------------------------------- */
#include "material_von_mises_mazars_non_local.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialVonMisesMazarsNonLocal<spatial_dimension>::MaterialVonMisesMazarsNonLocal(
    SolidMechanicsModel & model, const ID & id)
    : MaterialNonLocalParent(model, id), Ehat("epsilon_equ", *this),
      non_local_variable("mazars_non_local", *this) {
  AKANTU_DEBUG_IN();

  this->is_non_local = true;
  this->Ehat.initialize(1);
  this->non_local_variable.initialize(1);

  this->registerParam("average_on_damage", this->damage_in_compute_stress,
                      false, _pat_parsable | _pat_modifiable,
                      "Is D the non local variable");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialVonMisesMazarsNonLocal<spatial_dimension>::registerNonLocalVariables() {
  ID local;
  if (this->damage_in_compute_stress) {
    local = this->damage.getName();
  } else {
    local = this->Ehat.getName();
  }

  this->model.getNonLocalManager().registerNonLocalVariable(
      local, non_local_variable.getName(), 1);
  this->model.getNonLocalManager()
      .getNeighborhood(this->name)
      .registerNonLocalVariable(non_local_variable.getName());
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialVonMisesMazarsNonLocal<spatial_dimension>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * damage = this->damage(el_type, ghost_type).storage();
  Real * epsilon_equ = this->Ehat(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  MaterialVonMisesMazars<spatial_dimension>::computeStressOnQuad(grad_u, sigma, *damage,
                                                         *epsilon_equ);
  ++damage;
  ++epsilon_equ;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialVonMisesMazarsNonLocal<spatial_dimension>::computeNonLocalStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  auto & non_loc_var = non_local_variable(el_type, ghost_type);
  Real * damage;
  Real * epsilon_equ;
  if (this->damage_in_compute_stress) {
    damage = non_loc_var.storage();
    epsilon_equ = this->Ehat(el_type, ghost_type).storage();
  } else {
    damage = this->damage(el_type, ghost_type).storage();
    epsilon_equ = non_loc_var.storage();
  }

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
  this->computeDamageAndStressOnQuad(grad_u, sigma, *damage, *epsilon_equ);

  ++damage;
  ++epsilon_equ;
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

INSTANTIATE_MATERIAL(von_mises_mazars_non_local, MaterialVonMisesMazarsNonLocal);

} // namespace akantu
