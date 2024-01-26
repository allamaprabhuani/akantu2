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
#include "material_cohesive.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_COHESIVE_DAMAGE_HH_
#define AKANTU_MATERIAL_COHESIVE_DAMAGE_HH_

namespace akantu {

/**
 * Cohesive material linear damage
 *
 * parameters in the material files :
 *   - sigma_c   : critical stress sigma_c  (default: 0)
 *   - G_c       : fracture energy for mode I (default: 0)
 *   - k         : cohesive stiffness
 */
class MaterialCohesiveDamage : public MaterialCohesive {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialCohesiveDamage(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material parameters
  void initMaterial() override;

  /// these functions should be put in a separate class
  Real stiffness(Real d);
  Real augmented_stiffness(Real d);
  Real augmented_compliance(Real d);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// cohesive stiffness
  Real k;

  /// mode I fracture energy
  Real G_c;

  /// parameter for softening function
  Real a;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* AKANTU_MATERIAL_COHESIVE_DAMAGE_HH_ */
