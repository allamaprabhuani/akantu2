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

  void updateInternalParameters() override;

protected:
  /// constitutive law
  void computeTraction(ElementType el_type,
                       GhostType ghost_type = _not_ghost) override;

  /// compute tangent stiffness matrix
  void computeTangentTraction(ElementType el_type, Array<Real> & tangent_matrix,
                              GhostType ghost_type) override;

  /**
   * Scale insertion traction sigma_c according to the volume of the
   * two elements surrounding a facet
   *
   * see the article: F. Zhou and J. F. Molinari "Dynamic crack
   * propagation with cohesive elements: a methodology to address mesh
   * dependency" International Journal for Numerical Methods in
   * Engineering (2004)
   */
  void scaleInsertionTraction();

  inline decltype(auto) getArguments(ElementType element_type,
                                     GhostType ghost_type) {
    using namespace tuple;
    return zip_append(
        MaterialCohesive::getArguments<dim>(element_type, ghost_type),
        "sigma_c"_n = this->sigma_c_eff(element_type, ghost_type),
        "delta_c"_n = this->delta_c_eff(element_type, ghost_type),
        "insertion_stress"_n =
            make_view<dim>(this->insertion_stress(element_type, ghost_type)));
  }

  /// compute the traction for a given quadrature point
  template <typename Args> inline void computeTractionOnQuad(Args && args);

  template <class Derived, class Args>
  inline void computeTangentTractionOnQuad(Eigen::MatrixBase<Derived> & tangent,
                                           Args && args);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get sigma_c_eff
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(InsertionTraction, sigma_c_eff, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// cohesive stiffness
  Real k;

  /// mode I fracture energy
  Real G_c;


  /// variable defining if we are recomputing the last loading step
  /// after load_reduction
  bool recompute;

  /// critical effective stress
  RandomInternalField<Real, CohesiveInternalField> sigma_c_eff;

  /// effective critical displacement (each element can have a
  /// different value)
  CohesiveInternalField<Real> delta_c_eff;

  /// stress at insertion
  CohesiveInternalField<Real> insertion_stress;

  /// variable saying if there should be penalty contact also after
  /// breaking the cohesive elements
  bool contact_after_breaking;

  /// insertion of cohesive element when stress is high enough just on
  /// one quadrature point
  bool max_quad_stress_insertion;

  Vector<Real, dim> normal_opening;
  Vector<Real, dim> tangential_opening;
  Real normal_opening_norm{0.};
  Real tangential_opening_norm{0.};
  bool penetration{false};
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

} // namespace akantu

#include "material_cohesive_damage_inline_impl.hh"

#endif /* AKANTU_MATERIAL_COHESIVE_DAMAGE_HH_ */
