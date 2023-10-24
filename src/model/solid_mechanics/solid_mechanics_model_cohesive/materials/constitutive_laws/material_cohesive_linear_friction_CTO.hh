/**
 * @file   material_cohesive_linear_friction_CTO.hh
 *
 * @author Emil Gallyamov <emil.gallyamov@gmail.com>
 *
 * @date creation: Thu Sept 21 2023
 * @date last modification: Fri Dec 11 2020
 *
 * @brief  Addition to the existing linear friction cohesive material
 * by adding the consistent tangent operator when in sliping mode
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include "material_cohesive_linear_friction.hh"

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_MATERIAL_COHESIVE_LINEAR_FRICTION_CTO_HH_
#define AKANTU_MATERIAL_COHESIVE_LINEAR_FRICTION_CTO_HH_

/* -------------------------------------------------------------------------- */

namespace akantu {

/**
 * Cohesive material linear with friction force
 *
 * parameters in the material files :
 *   - mu   : friction coefficient
 *   - penalty_for_friction : Penalty parameter for the friction behavior
 */
template <UInt spatial_dimension>
class MaterialCohesiveLinearFrictionCTO
    : public MaterialCohesiveLinearFriction<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  using MaterialParent = MaterialCohesiveLinearFriction<spatial_dimension>;

public:
  MaterialCohesiveLinearFrictionCTO(SolidMechanicsModel & model,
                                    const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material parameters
  void initMaterial() override;

protected:
  /// compute tangent stiffness matrix
  void computeTangentTraction(ElementType el_type, Array<Real> & tangent_matrix,
                              const Array<Real> & normal,
                              GhostType ghost_type) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// basis alligned with cohesive elements
  CohesiveInternalField<Real> basis;
};

} // namespace akantu

#endif /* AKANTU_MATERIAL_COHESIVE_LINEAR_FRICTION_CTO_HH_ */
