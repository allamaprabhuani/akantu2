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
 * Cohesive material linear damage for extrinsic case
 *
 * parameters in the material files :
 *   - sigma_c   : critical stress sigma_c  (default: 0)
 *   - beta      : weighting parameter for sliding and normal opening (default:
 * 0)
 *   - G_cI      : fracture energy for mode I (default: 0)
 *   - G_cII     : fracture energy for mode II (default: 0)
 *   - penalty   : stiffness in compression to prevent penetration
 */
template <Int dim> class MaterialCohesiveDamage : public MaterialCohesive {
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

  /// assemble stiffness
  void assembleStiffnessMatrix(GhostType ghost_type) override;

  /// assemble residual
  void assembleInternalForces(GhostType ghost_type = _not_ghost) override;

protected:
  void computeLambdaOnQuad(ElementType type, GhostType ghost_type);

  inline decltype(auto) getArguments(ElementType element_type,
                                     GhostType ghost_type) {
    using namespace tuple;
    return zip_append(
        MaterialCohesive::getArguments<dim>(element_type, ghost_type),
        "lambda"_n = make_view<dim>(this->lambda(element_type, ghost_type)),
        "err_opening"_n =
            make_view<dim>(this->err_openings(element_type, ghost_type)));
  }

  /// constitutive law
  void computeTraction(ElementType el_type,
                       GhostType ghost_type = _not_ghost) override;

  /// compute tangent stiffness matrix
  /// WARNING : override was removed, not sure it should have
  void computeTangentTraction(ElementType el_type,
                              Array<Real> & tangent_matrix_uu,
                              Array<Real> & tangent_matrix_ll,
                              GhostType ghost_type);

  /// compute the traction for a given quadrature point
  template <typename Args> inline void computeTractionOnQuad(Args && args);

  template <class Derived, class Args>
  inline void
  computeTangentTractionOnQuad(Eigen::MatrixBase<Derived> & tangent_uu,
                               Eigen::MatrixBase<Derived> & tangent_ll,
                               Args && args);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// cohesive stiffness
  Real k;

  /// mode I fracture energy
  Real G_c;

  /// augmented lagrange multiplier
  CohesiveInternalField<Real> & lambda;

  /// target opening
  CohesiveInternalField<Real> & err_openings;

  /// cohesive damage
  CohesiveInternalField<Real> & czm_damage;

  Vector<Real, dim> normal_opening;
  Real normal_opening_norm{0.};
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

} // namespace akantu

#include "material_cohesive_damage_inline_impl.hh"

#endif /* AKANTU_MATERIAL_COHESIVE_DAMAGE_HH_ */
