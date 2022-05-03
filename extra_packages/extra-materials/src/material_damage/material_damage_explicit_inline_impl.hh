/**
 * @file   material_damage_explicit_inline_impl.hh
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @date   Tue Feb 12 2019
 *
 * @brief  Inline implementation of material damage explicit
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
#include "material_damage_explicit.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt dim, template <UInt> class ElasticParent>
MaterialDamageExplicit<dim, ElasticParent>::MaterialDamageExplicit(
    SolidMechanicsModel & model, const ID & id)
    : parent(model, id), Sc("Sc", *this), Gf(0.), crack_band_width(0.),
      max_damage(1.), Et("Et", *this) {

  AKANTU_DEBUG_IN();
  this->registerParam("Sc", Sc, _pat_parsable, "critical stress threshold");
  this->registerParam("Gf", Gf, _pat_parsable | _pat_modifiable,
                      "fracture energy");
  this->registerParam("crack_band_width", crack_band_width,
                      _pat_parsable | _pat_modifiable, "crack_band_width");
  this->registerParam("max_damage", max_damage, _pat_parsable | _pat_modifiable,
                      "maximum possible damage");

  this->Sc.initialize(1);
  this->Et.initialize(1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim, template <UInt> class ElasticParent>
void MaterialDamageExplicit<dim, ElasticParent>::initMaterial() {
  AKANTU_DEBUG_IN();

  parent::initMaterial();

  GhostType ghost_type = _not_ghost;

  for (auto && el_type : this->element_filter.elementTypes(dim, ghost_type)) {
    for (auto && data :
         zip(this->Sc(el_type, ghost_type), this->Et(el_type, ghost_type))) {
      auto & strength = std::get<0>(data);
      auto & soft_slope = std::get<1>(data);
      soft_slope = strength / (strength / this->E - 2 * this->Gf / strength /
                                                        this->crack_band_width);
      AKANTU_DEBUG_ASSERT(soft_slope < 0,
                          "The softening slope is positive or 0");
    }
  }
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt dim, template <UInt> class ElasticParent>
void MaterialDamageExplicit<dim, ElasticParent>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  parent::computeStress(el_type, ghost_type);

  auto dam = this->damage(el_type, ghost_type).begin();
  auto Et = this->Et(el_type, ghost_type).begin();
  auto Sc = this->Sc(el_type, ghost_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  computeDamageAndStressOnQuad(grad_u, sigma, *dam, *Et, *Sc);

  ++dam;
  ++Et;
  ++Sc;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim, template <UInt> class ElasticParent>
inline void
MaterialDamageExplicit<dim, ElasticParent>::computeDamageAndStressOnQuad(
    Matrix<Real> & grad_u, Matrix<Real> & sigma, Real & dam, const Real & Et,
    const Real & Sc) {
  Real w = 0;
  for (UInt i = 0; i < dim; ++i) {
    for (UInt j = 0; j < dim; ++j) {
      w += sigma(i, j) * (grad_u(i, j) + grad_u(j, i)) / 2.;
    }
  }
  w *= 0.5;

  Real wy = Sc * Sc / (2 * this->E);
  Real gamma = -Et / this->E;
  Real k = wy * pow(((1 + gamma) / (1 + gamma - dam)), 2);

  Real F = w - k;

  if (F > 0) {
    dam = (1 + gamma) * (1 - std::sqrt(wy / w));
    dam = std::min(dam, this->max_damage);
    dam = std::max(0., dam);
  }
  sigma *= 1 - dam;
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
