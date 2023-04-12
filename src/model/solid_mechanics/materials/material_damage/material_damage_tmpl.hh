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
    : Parent<dim>(model, id) {
  this->damage = this->registerInternal("damage", 1);
  this->dissipated_energy =
      this->registerInternal("damage dissipated energy", 1);
  this->int_sigma = this->registerInternal("integral of sigma", 1);
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class Parent>
void MaterialDamage<dim, Parent>::initMaterial() {
  Parent<dim>::initMaterial();
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

  for (auto && [args, epot, edis, ints] :
       zip(getArguments(el_type), (*this->potential_energy)(el_type),
           (*this->dissipated_energy)(el_type), (*this->int_sigma)(el_type))) {

    Matrix<Real, dim, dim> delta_gradu =
        args["grad_u"_n] - args["previous_grad_u"_n];
    Matrix<Real, dim, dim> sigma_h = args["sigma"_n] + args["previous_sigma"_n];
    auto dint = sigma_h.doubleDot(delta_gradu) / 2;

    ints += dint;
    edis = ints - epot;
  }
}

/* --------------------------------------------------------------------------
 */
template <Int dim, template <Int> class Parent>
void MaterialDamage<dim, Parent>::computeTangentModuli(
    ElementType el_type, Array<Real> & tangent_matrix, GhostType ghost_type) {
  constexpr auto tangent_size = Material::getTangentStiffnessVoigtSize(dim);

  for (auto && data :
       zip_append(getArguments(el_type, ghost_type),
                  "tangent_moduli"_n =
                      make_view<tangent_size, tangent_size>(tangent_matrix))) {
    Parent<dim>::computeTangentModuliOnQuad(data);
    computeTangentModuliOnQuad(data);
  }
}

/* --------------------------------------------------------------------------
 */
template <Int dim, template <Int> class Parent>
template <class Args>
void MaterialDamage<dim, Parent>::computeTangentModuliOnQuad(
    Args && arguments) {
  arguments["tangent_moduli"_n].array() *= (1 - arguments["damage"_n]);
}

/* --------------------------------------------------------------------------
 */
template <Int dim, template <Int> class Parent>
auto MaterialDamage<dim, Parent>::getDissipatedEnergy() const -> Real {
  Real de = 0.;

  /// integrate the dissipated energy for each type of elements
  for (auto && type : this->element_filter.elementTypes(dim, _not_ghost)) {
    de += this->getFEEngine().integrate((*dissipated_energy)(type, _not_ghost),
                                        type, _not_ghost,
                                        this->element_filter(type, _not_ghost));
  }
  return de;
}

/* --------------------------------------------------------------------------
 */
template <Int dim, template <Int> class Parent>
Real MaterialDamage<dim, Parent>::getEnergy(const std::string & type) {
  if (type == "dissipated") {
    return getDissipatedEnergy();
  }

  return Parent<dim>::getEnergy(type);
}

/* --------------------------------------------------------------------------
 */
} // namespace akantu
