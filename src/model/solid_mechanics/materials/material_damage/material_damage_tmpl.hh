/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
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
 */

/* -------------------------------------------------------------------------- */
#include "material_damage.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class Parent>
MaterialDamage<dim, Parent>::MaterialDamage(SolidMechanicsModel & model,
                                            const ID & id)
    : Parent<dim>(model, id), damage("damage", *this),
      dissipated_energy("damage dissipated energy", *this),
      int_sigma("integral of sigma", *this) {
  AKANTU_DEBUG_IN();

  this->is_non_local = false;
  this->use_previous_stress = true;
  this->use_previous_gradu = true;

  this->damage.initialize(1);
  this->dissipated_energy.initialize(1);
  this->int_sigma.initialize(1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class Parent>
void MaterialDamage<dim, Parent>::initMaterial() {
  AKANTU_DEBUG_IN();
  Parent<dim>::initMaterial();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Compute the dissipated energy in  each element by a trapezoidal approximation
 * of
 * @f$ Ed = \int_0^{\epsilon}\sigma(\omega)d\omega -
 * \frac{1}{2}\sigma:\epsilon@f$
 */
template <Int dim, template <Int> class Parent>
void MaterialDamage<dim, Parent>::updateEnergies(ElementType el_type) {
  Parent<dim>::updateEnergies(el_type);

  this->computePotentialEnergy(el_type);

  auto epsilon_p = this->gradu.previous(el_type).begin(dim, dim);
  auto sigma_p = this->stress.previous(el_type).begin(dim, dim);

  auto epot = this->potential_energy(el_type).begin();
  auto ints = this->int_sigma(el_type).begin();
  auto ed = this->dissipated_energy(el_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

  Matrix<Real> delta_gradu(grad_u);
  delta_gradu -= *epsilon_p;

  Matrix<Real> sigma_h(sigma);
  sigma_h += *sigma_p;

  Real dint = .5 * sigma_h.doubleDot(delta_gradu);

  *ints += dint;
  *ed = *ints - *epot;

  ++epsilon_p;
  ++sigma_p;
  ++epot;
  ++ints;
  ++ed;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class Parent>
void MaterialDamage<dim, Parent>::computeTangentModuli(
    ElementType el_type, Array<Real> & tangent_matrix, GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  constexpr auto tangent_size = Material::getTangentStiffnessVoigtSize(dim);

  for (auto && data :
       zip_append(getArguments(el_type, ghost_type),
                  "tangent_moduli"_n =
                      make_view<tangent_size, tangent_size>(tangent_matrix))) {
    Parent<dim>::computeTangentModuliOnQuad(data);
    computeTangentModuliOnQuad(data);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class Parent>
template <class Args>
void MaterialDamage<dim, Parent>::computeTangentModuliOnQuad(
    Args && arguments) {
  arguments["tangent_moduli"_n].array() *= (1 - arguments["damage"_n]);
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class Parent>
auto MaterialDamage<dim, Parent>::getDissipatedEnergy() const -> Real {
  AKANTU_DEBUG_IN();

  Real de = 0.;

  /// integrate the dissipated energy for each type of elements
  for (auto && type : this->element_filter.elementTypes(dim, _not_ghost)) {
    de +=
        this->fem.integrate(dissipated_energy(type, _not_ghost), type,
                            _not_ghost, this->element_filter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return de;
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class Parent>
Real MaterialDamage<dim, Parent>::getEnergy(const std::string & type) {
  if (type == "dissipated") {
    return getDissipatedEnergy();
  }

  return Parent<dim>::getEnergy(type);
}

/* -------------------------------------------------------------------------- */
} // namespace akantu
