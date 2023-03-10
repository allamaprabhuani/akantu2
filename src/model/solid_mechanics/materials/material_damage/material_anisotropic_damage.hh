/**
 * @file   material_anisotropic_damage.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sat Feb 03 2018
 * @date last modification: Fri Jul 24 2020
 *
 * @brief  Base class for anisotropic damage materials
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_ANISOTROPIC_DAMAGE_HH_
#define AKANTU_MATERIAL_ANISOTROPIC_DAMAGE_HH_

namespace akantu {

template <Int dim, template <Int> class EquivalentStrain,
          template <Int> class DamageThreshold,
          template <Int> class Parent = MaterialElastic>
class MaterialAnisotropicDamage : public Parent<dim> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialAnisotropicDamage(SolidMechanicsModel & model, const ID & id = "");
  ~MaterialAnisotropicDamage() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void computeStress(ElementType el_type, GhostType ghost_type) override;

  decltype(auto) getArguments(ElementType el_type,
                              GhostType ghost_type = _not_ghost) {
    return zip_append(
        Parent<dim>::getArguments(el_type, ghost_type),
        "damage"_n = make_view<dim, dim>(this->damage(el_type, ghost_type)),
        "sigma_el"_n =
            make_view<dim, dim>(this->elastic_stress(el_type, ghost_type)),
        "epsilon_hat"_n = this->equivalent_strain(el_type, ghost_type),
        "TrD"_n = this->trace_damage(el_type, ghost_type),
        "TrD_n_1"_n = this->trace_damage.previous(el_type, ghost_type),
        "equivalent_strain_data"_n = equivalent_strain_function,
        "damage_threshold_data"_n = damage_threshold_function);
  }

  template <class Args> void computeStressOnQuad(Args && args);

private:
  template <class D1, class D2, class D3>
  void damageStress(Eigen::MatrixBase<D1> & sigma,
                    const Eigen::MatrixBase<D2> & sigma_el,
                    const Eigen::MatrixBase<D3> & D, Real TrD);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  Real Dc{0.99};

  /// damage internal variable
  InternalField<Real> damage;

  /// elastic stress
  InternalField<Real> elastic_stress;

  /// equivalent strain
  InternalField<Real> equivalent_strain;

  /// trace of the damageThreshold
  InternalField<Real> trace_damage;

  /// damage criteria
  EquivalentStrain<dim> equivalent_strain_function;

  /// damage evolution
  DamageThreshold<dim> damage_threshold_function;
};

} // namespace akantu

#include "material_anisotropic_damage_tmpl.hh"

#endif /* AKANTU_MATERIAL_ANISOTROPIC_DAMAGE_HH_ */
