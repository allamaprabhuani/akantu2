/**
 * @file   material_damage_tmpl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Aug 18 2015
 *
 * @brief  Specialization of the material class for the damage material
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "material_damage.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class Parent>
MaterialDamage<spatial_dimension, Parent>::MaterialDamage(
    SolidMechanicsModel & model, const ID & id)
    : Parent<spatial_dimension>(model, id), damage("damage", *this),
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
template <UInt spatial_dimension, template <UInt> class Parent>
void MaterialDamage<spatial_dimension, Parent>::initMaterial() {
  AKANTU_DEBUG_IN();
  Parent<spatial_dimension>::initMaterial();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Compute the dissipated energy in  each element by a trapezoidal approximation
 * of
 * @f$ Ed = \int_0^{\epsilon}\sigma(\omega)d\omega -
 * \frac{1}{2}\sigma:\epsilon@f$
 */
template <UInt spatial_dimension, template <UInt> class Parent>
void MaterialDamage<spatial_dimension, Parent>::updateEnergies(
    ElementType el_type, GhostType ghost_type) {
  Parent<spatial_dimension>::updateEnergies(el_type, ghost_type);

  this->computePotentialEnergy(el_type, ghost_type);

  Array<Real>::matrix_iterator epsilon_p =
      this->gradu.previous(el_type, ghost_type)
          .begin(spatial_dimension, spatial_dimension);
  Array<Real>::matrix_iterator sigma_p =
      this->stress.previous(el_type, ghost_type)
          .begin(spatial_dimension, spatial_dimension);

  Array<Real>::const_scalar_iterator epot =
      this->potential_energy(el_type, ghost_type).begin();

  Array<Real>::scalar_iterator ints =
      this->int_sigma(el_type, ghost_type).begin();
  Array<Real>::scalar_iterator ed =
      this->dissipated_energy(el_type, ghost_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  Matrix<Real> delta_gradu_it(*gradu_it);
  delta_gradu_it -= *epsilon_p;

  Matrix<Real> sigma_h(sigma);
  sigma_h += *sigma_p;

  Real dint = .5 * sigma_h.doubleDot(delta_gradu_it);

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
template <UInt spatial_dimension, template <UInt> class Parent>
void MaterialDamage<spatial_dimension, Parent>::computeTangentModuli(
    const ElementType & el_type, Array<Real> & tangent_matrix,
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Parent<spatial_dimension>::computeTangentModuli(el_type, tangent_matrix,
                                                  ghost_type);

  Real * dam = this->damage(el_type, ghost_type).storage();

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
  computeTangentModuliOnQuad(tangent, *dam);

  ++dam;

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class Parent>
void MaterialDamage<spatial_dimension, Parent>::computeTangentModuliOnQuad(
    Matrix<Real> & tangent, Real & dam) {
  tangent *= (1 - dam);
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class Parent>
Real MaterialDamage<spatial_dimension, Parent>::getDissipatedEnergy() const {
  AKANTU_DEBUG_IN();

  Real de = 0.;

  /// integrate the dissipated energy for each type of elements
  for (auto & type :
       this->element_filter.elementTypes(spatial_dimension, _not_ghost)) {
    de +=
        this->fem.integrate(dissipated_energy(type, _not_ghost), type,
                            _not_ghost, this->element_filter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return de;
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class Parent>
Real MaterialDamage<spatial_dimension, Parent>::getEnergy(std::string type) {
  if (type == "dissipated")
    return getDissipatedEnergy();
  else
    return Parent<spatial_dimension>::getEnergy(type);
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
