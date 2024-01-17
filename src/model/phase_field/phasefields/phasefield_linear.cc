/* -------------------------------------------------------------------------- */
#include "phasefield_linear.hh"
#include "aka_common.hh"
#include <tuple>

namespace akantu {

/* -------------------------------------------------------------------------- */
PhaseFieldLinear::PhaseFieldLinear(PhaseFieldModel & model, const ID & id)
    : PhaseField(model, id) {
  registerParam("irreversibility_tol", tol_ir, Real(1e-2),
                _pat_parsable | _pat_readable, "Irreversibility tolerance");
}

/* -------------------------------------------------------------------------- */
void PhaseFieldLinear::initPhaseField() {
  PhaseField::initPhaseField();

  this->gamma = Real(this->g_c) / this->l0 * 27. / (64. * tol_ir * tol_ir);

  this->dim = spatial_dimension;
  if (spatial_dimension == 2 && !this->plane_stress) {
    this->dim = 3;
  }
}

/* -------------------------------------------------------------------------- */
void PhaseFieldLinear::updateInternalParameters() {
  PhaseField::updateInternalParameters();

  for (const auto & type :
       element_filter.elementTypes(spatial_dimension, _not_ghost)) {
    for (auto && tuple : zip(make_view(this->damage_energy(type, _not_ghost),
                                       spatial_dimension, spatial_dimension),
                             this->g_c(type, _not_ghost))) {
      Matrix<Real> d(spatial_dimension, spatial_dimension);
      // eye g_c * l0
      d.eye(3. / 4. * std::get<1>(tuple) * this->l0);
      std::get<0>(tuple) = d;
    }
  }
}

/* -------------------------------------------------------------------------- */
// void PhaseFieldLinear::computeDrivingForce(const ElementType & el_type,
//                                            GhostType ghost_type) {
//   for (auto && tuple : zip(this->phi(el_type, ghost_type),
//                            this->phi.previous(el_type, ghost_type),
//                            this->driving_force(el_type, ghost_type),
//                            this->damage_energy_density(el_type, ghost_type),
//                            make_view(this->strain(el_type, ghost_type),
//                                      spatial_dimension, spatial_dimension),
//                            this->g_c(el_type, ghost_type))) {
//     computePhiOnQuad(std::get<4>(tuple), std::get<0>(tuple),
//                      std::get<1>(tuple));
//     computeDamageEnergyDensityOnQuad(std::get<0>(tuple), std::get<3>(tuple),
//                                      std::get<5>(tuple));
//     computeDrivingForceOnQuad(std::get<0>(tuple), std::get<2>(tuple),
//                               std::get<5>(tuple));
//   }
// }

/* -------------------------------------------------------------------------- */
void PhaseFieldLinear::computeDrivingForce(ElementType el_type,
                                           GhostType ghost_type) {

  if (this->isotropic) {
    for (auto && tuple : zip(this->phi(el_type, ghost_type),
                             make_view(this->strain(el_type, ghost_type),
                                       spatial_dimension, spatial_dimension))) {
      auto & phi_quad = std::get<0>(tuple);
      auto & strain = std::get<1>(tuple);
      computePhiIsotropicOnQuad(strain, phi_quad);
    }
  } else {
    for (auto && tuple : zip(this->phi(el_type, ghost_type),
                             make_view(this->strain(el_type, ghost_type),
                                       spatial_dimension, spatial_dimension))) {
      auto & phi_quad = std::get<0>(tuple);
      auto & strain = std::get<1>(tuple);
      computePhiOnQuad(strain, phi_quad);
    }
  }

  for (auto && tuple :
       zip(this->phi(el_type, ghost_type),
           this->driving_force(el_type, ghost_type),
           this->damage_energy_density(el_type, ghost_type),
           this->damage_on_qpoints(el_type, _not_ghost),
           make_view(this->driving_energy(el_type, ghost_type),
                     spatial_dimension),
           make_view(this->damage_energy(el_type, ghost_type),
                     spatial_dimension, spatial_dimension),
           make_view(this->gradd(el_type, ghost_type), spatial_dimension),
           this->g_c(el_type, ghost_type),
           this->damage_on_qpoints.previous(el_type, ghost_type))) {
    auto & phi_quad = std::get<0>(tuple);
    auto & driving_force_quad = std::get<1>(tuple);
    auto & dam_energy_density_quad = std::get<2>(tuple);
    auto & dam_on_quad = std::get<3>(tuple);
    auto & driving_energy_quad = std::get<4>(tuple);
    auto & damage_energy_quad = std::get<5>(tuple);
    auto & gradd_quad = std::get<6>(tuple);
    auto & g_c_quad = std::get<7>(tuple);
    auto & dam_prev_quad = std::get<8>(tuple);

    computeDamageEnergyDensityOnQuad(phi_quad, dam_energy_density_quad,
                                     g_c_quad);
    Real penalization =
        this->gamma * std::min(Real(0.), dam_on_quad - dam_prev_quad);

    driving_force_quad = dam_on_quad * dam_energy_density_quad - 2 * phi_quad +
                         3 * g_c_quad / (8 * this->l0) + penalization;
    driving_energy_quad = damage_energy_quad * gradd_quad;

    dam_energy_density_quad += this->gamma * (dam_on_quad < dam_prev_quad);
  }
}

/* -------------------------------------------------------------------------- */
void PhaseFieldLinear::computeDissipatedEnergy(ElementType el_type) {
  AKANTU_DEBUG_IN();

  for (auto && tuple :
       zip(this->dissipated_energy(el_type, _not_ghost),
           this->damage_on_qpoints(el_type, _not_ghost),
           this->damage_on_qpoints.previous(el_type, _not_ghost),
           make_view(this->gradd(el_type, _not_ghost), spatial_dimension),
           this->g_c(el_type, _not_ghost))) {

    this->computeDissipatedEnergyOnQuad(std::get<1>(tuple), std::get<2>(tuple),
                                        std::get<3>(tuple), std::get<0>(tuple),
                                        std::get<4>(tuple));
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PhaseFieldLinear::computeDissipatedEnergyByElement(
    ElementType type, Idx index, Vector<Real> & edis_on_quad_points) {
  auto gradd_it = this->gradd(type).begin(spatial_dimension);
  auto gradd_end = this->gradd(type).begin(spatial_dimension);
  auto damage_it = this->damage_on_qpoints(type).begin();
  auto damage_prev_it = this->damage_on_qpoints.previous(type).begin();
  auto g_c_it = this->g_c(type).begin();

  UInt nb_quadrature_points = fem.getNbIntegrationPoints(type);

  gradd_it += index * nb_quadrature_points;
  gradd_end += (index + 1) * nb_quadrature_points;
  damage_it += index * nb_quadrature_points;
  damage_prev_it += index * nb_quadrature_points;
  g_c_it += index * nb_quadrature_points;

  Real * edis_quad = edis_on_quad_points.data();

  for (; gradd_it != gradd_end; ++gradd_it, ++damage_it, ++edis_quad) {
    this->computeDissipatedEnergyOnQuad(*damage_it, *damage_prev_it, *gradd_it,
                                        *edis_quad, *g_c_it);
  }
}

/* -------------------------------------------------------------------------- */
void PhaseFieldLinear::computeDissipatedEnergyByElement(
    const Element & element, Vector<Real> & edis_on_quad_points) {
  computeDissipatedEnergyByElement(element.type, element.element,
                                   edis_on_quad_points);
}

/* -------------------------------------------------------------------------- */
void PhaseFieldLinear::afterSolveStep() {
  // clamp negative damage to 0
  for (auto & dam : this->model.getDamage()) {
    dam = std::max(Real(0.), dam);
  }
}

INSTANTIATE_PHASEFIELD(linear, PhaseFieldLinear);

} // namespace akantu
