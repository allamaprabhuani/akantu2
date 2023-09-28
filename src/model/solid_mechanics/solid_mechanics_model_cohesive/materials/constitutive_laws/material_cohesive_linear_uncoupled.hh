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
#ifndef AKANTU_MATERIAL_COHESIVE_LINEAR_UNCOUPLED_HH_
#define AKANTU_MATERIAL_COHESIVE_LINEAR_UNCOUPLED_HH_

/* -------------------------------------------------------------------------- */

namespace akantu {

/**
 * Cohesive material linear with two different laws for mode I and
 * mode II, for extrinsic case
 *
 * parameters in the material files :
 *  - roughness : define the interaction between mode I and mode II (default: 0)
 */
template <Int dim>
class MaterialCohesiveLinearUncoupled : public MaterialCohesiveLinear<dim> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  //  typedef MaterialCohesiveLinear<spatial_dimension> MaterialParent;
public:
  MaterialCohesiveLinearUncoupled(SolidMechanicsModel & model,
                                  const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// constitutive law
  void computeTraction(ElementType el_type,
                       GhostType ghost_type = _not_ghost) override;

  /// compute tangent stiffness matrix
  void computeTangentTraction(ElementType el_type, Array<Real> & tangent_matrix,
                              GhostType ghost_type) override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// parameter to tune the interaction between mode II and mode I
  Real R;

  /// maximum normal displacement
  CohesiveInternalField<Real> & delta_n_max;

  /// maximum tangential displacement
  CohesiveInternalField<Real> & delta_t_max;

  /// damage associated to normal tractions
  CohesiveInternalField<Real> & damage_n;

  /// damage associated to shear tractions
  CohesiveInternalField<Real> & damage_t;
};

} // namespace akantu

#endif /* AKANTU_MATERIAL_COHESIVE_LINEAR_UNCOUPLED_HH_ */
