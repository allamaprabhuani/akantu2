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
#include "material_cohesive_damage.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_COHESIVE_DAMAGE_EXTRINSIC_HH_
#define AKANTU_MATERIAL_COHESIVE_DAMAGE_EXTRINSIC_HH_

namespace akantu {

/**
 * Cohesive material linear damage for extrinsinc case
 *
 * parameters in the material files :
 *   - sigma_c   : critical stress sigma_c  (default: 0)
 *   - G_c       : fracture energy for mode I (default: 0)
 *   - k         : cohesive stiffness
 */
template <Int dim> class MaterialCohesiveDamageExtrinsic : public MaterialCohesiveDamage {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialCohesiveDamageExtrinsic(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  /// check stress for cohesive elements' insertion
  void checkInsertion(bool check_only = false) override;

  void computeLambda(GhostType ghost_type = _not_ghost) override;

protected:
  void computeLambdaOnQuad(ElementType type, GhostType ghost_type);

  /// constitutive law
  void computeTraction(ElementType el_type,
                       GhostType ghost_type = _not_ghost) override;

//  /// compute the traction for a given quadrature point
//  template <typename Args> inline void computeTractionOnQuad(Args && args);

  [[nodiscard]] bool needLambda() const override { return false; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// augmented lagrange multiplier
  FacetInternalField<Real> & lambda;

  /// cohesive damage
  FacetInternalField<Real> & czm_damage;

  /// stress at insertion
  CohesiveInternalField<Real> & insertion_stress;

  /// stress at insertion
  CohesiveInternalField<Real> & is_new_crack;

};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

} // namespace akantu

#include "material_cohesive_damage_extrinsic_inline_impl.hh"

#endif /* AKANTU_MATERIAL_COHESIVE_DAMAGE_EXTRINSIC_HH_ */
