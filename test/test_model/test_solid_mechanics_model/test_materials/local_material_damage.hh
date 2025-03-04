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
#include "aka_common.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_LOCAL_MATERIAL_DAMAGE_HH_
#define AKANTU_LOCAL_MATERIAL_DAMAGE_HH_

namespace akantu {

class LocalMaterialDamage : public Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  LocalMaterialDamage(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initMaterial() override;

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  /// constitutive law for a given quadrature point
  template <class D1, class D2>
  inline void computeStressOnQuad(Eigen::MatrixBase<D1> & grad_u,
                                  Eigen::MatrixBase<D2> & sigma, Real & dam);

  /// compute the potential energy for all elements
  void computePotentialEnergy(ElementType el_type) override;

  /// compute the potential energy for on element
  template <class D1, class D2>
  inline void computePotentialEnergyOnQuad(Eigen::MatrixBase<D1> & grad_u,
                                           Eigen::MatrixBase<D2> & sigma,
                                           Real & epot);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// compute the celerity of wave in the material
  [[nodiscard]] inline Real getCelerity(const Element & element) const override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Damage, damage, Real);

private:
  /// the young modulus
  Real E{};

  /// Poisson coefficient
  Real nu{};

  /// First Lamé coefficient
  Real lambda{};

  /// Second Lamé coefficient (shear modulus)
  Real mu{};

  /// resistance to damage
  Real Yd{};

  /// damage threshold
  Real Sd{};

  /// Bulk modulus
  Real kpa{};

  /// damage internal variable
  InternalField<Real> & damage;
};

} // namespace akantu

#include "local_material_damage_inline_impl.hh"

#endif /* AKANTU_LOCAL_MATERIAL_DAMAGE_HH_ */
