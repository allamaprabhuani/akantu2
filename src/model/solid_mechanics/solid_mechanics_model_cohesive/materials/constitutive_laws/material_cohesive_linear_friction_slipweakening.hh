/**
 * @file   material_cohesive_linear_slipweakening.hh
 *
 * @author Mathias Lebihain <mathias.lebihain@epfl.ch>
 *
 * @date creation: Fri Mar 06 2020
 * @date last modification: Fri Mar 06 2020
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading with
 * Mohr-Coulomb insertion criterion and slip-weakening friction for
 * extrinsic type
 *
 * @section LICENSE
 *
 * Copyright (©)  2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */

#include "material_cohesive_linear_friction.hh"

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_COHESIVE_LINEAR_FRICTION_SLIPWEAKENING_HH__
#define __AKANTU_MATERIAL_COHESIVE_LINEAR_FRICTION_SLIPWEAKENING_HH__
/* -------------------------------------------------------------------------- */

namespace akantu {

/**
 * Cohesive material linear with slip-weakening friction force
 *
 * parameters in the material files :
 *   - mu_d   : dynamic friction coefficient
 */
template <UInt spatial_dimension>
class MaterialCohesiveLinearFrictionSlipWeakening
    : public MaterialCohesiveLinearFriction<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  using MaterialParent = MaterialCohesiveLinearFriction<spatial_dimension>;

public:
  MaterialCohesiveLinearFrictionSlipWeakening(SolidMechanicsModel & model,
                                              const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material parameters
  void initMaterial() override;

protected:
  /// constitutive law
  void computeTraction(const Array<Real> & normal, ElementType el_type,
                       GhostType ghost_type = _not_ghost) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// dynamic friction coefficient
  RandomInternalField<Real, CohesiveInternalField> mu_dynamic;
};

} // namespace akantu

#endif /* __AKANTU_MATERIAL_COHESIVE_LINEAR_FRICTION_SLIPWEAKENING_HH__ */
