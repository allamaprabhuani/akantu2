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
#include "material.hh"
#include "material_damage.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_MAZARS_HH_
#define AKANTU_MATERIAL_MAZARS_HH_

namespace akantu {

/**
 * Material Mazars
 *
 * parameters in the material files :
 *   - rho  : density (default: 0)
 *   - E    : Young's modulus (default: 0)
 *   - nu   : Poisson's ratio (default: 1/2)
 *   - K0   : Damage threshold
 *   - At   : Parameter damage traction 1
 *   - Bt   : Parameter damage traction 2
 *   - Ac   : Parameter damage compression 1
 *   - Bc   : Parameter damage compression 2
 *   - beta : Parameter for shear
 */
template <Int dim, template <Int> class Parent = MaterialElastic>
class MaterialMazars : public MaterialDamage<dim, Parent> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  using parent_damage = MaterialDamage<dim, Parent>;

public:
  MaterialMazars(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

protected:
  /// constitutive law for a given quadrature point
  template <typename Args> inline void computeStressOnQuad(Args && args);

  template <typename Args>
  inline void computeDamageAndStressOnQuad(Args && args);

  template <typename Args, typename Derived>
  inline void
  computeDamageOnQuad(Args && args,
                      const Eigen::MatrixBase<Derived> & epsilon_princ);

public:
  decltype(auto) getArguments(ElementType el_type, GhostType ghost_type) {
    return zip_append(
        parent_damage::getArguments(el_type, ghost_type),
        "K0"_n = this->K0(el_type, ghost_type),
        "Ehat"_n =
            broadcast(this->Ehat, this->damage(el_type, ghost_type).size()));
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// damage threshold
  DefaultRandomInternalField<Real> & K0;
  /// parameter damage traction 1
  Real At{0.};
  /// parameter damage traction 2
  Real Bt{0.};
  /// parameter damage compression 1
  Real Ac{0.};
  /// parameter damage compression 2
  Real Bc{0.};
  /// parameter for shear
  Real beta{0.};

  /// specify the variable to average false = ehat, true = damage (only valid
  /// for non local version)
  bool damage_in_compute_stress{true};

  Real Ehat{0};
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "material_mazars_inline_impl.hh"

#endif /* __AKANTU_MATERIAL_MAZARS_HH__ */
