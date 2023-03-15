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

#include "material_cohesive_linear.hh"

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_MATERIAL_COHESIVE_LINEAR_FRICTION_HH_
#define AKANTU_MATERIAL_COHESIVE_LINEAR_FRICTION_HH_

/* -------------------------------------------------------------------------- */

namespace akantu {

/**
 * Cohesive material linear with friction force
 *
 * parameters in the material files :
 *   - mu   : friction coefficient
 *   - penalty_for_friction : Penalty parameter for the friction behavior
 */
template <Int dim>
class MaterialCohesiveLinearFriction : public MaterialCohesiveLinear<dim> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  using MaterialParent = MaterialCohesiveLinear<dim>;

public:
  MaterialCohesiveLinearFriction(SolidMechanicsModel & model,
                                 const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material parameters
  void initMaterial() override;

protected:
  /// constitutive law
  void computeTraction(ElementType el_type,
                       GhostType ghost_type = _not_ghost) override;

  /// compute tangent stiffness matrix
  void computeTangentTraction(ElementType el_type, Array<Real> & tangent_matrix,
                              GhostType ghost_type) override;

  decltype(auto) getArguments(ElementType el_type,
                              GhostType ghost_type = _not_ghost) {
    using namespace tuple;
    return zip_append(MaterialParent::getArguments(el_type, ghost_type),
                      "residual_sliding"_n =
                          residual_sliding(el_type, ghost_type),
                      "friction_force"_n =
                          make_view<dim>(friction_force(el_type, ghost_type)),
                      "previous_residual_sliding"_n =
                          residual_sliding.previous(el_type, ghost_type));
  }

  template <typename Args> void computeTractionOnQuad(Args && args);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// maximum value of the friction coefficient
  Real mu_max;

  /// penalty parameter for the friction law
  Real friction_penalty;

  /// history parameter for the friction law
  CohesiveInternalField<Real> residual_sliding;

  /// friction force
  CohesiveInternalField<Real> friction_force;
};

} // namespace akantu

#endif /* AKANTU_MATERIAL_COHESIVE_LINEAR_FRICTION_HH_ */
