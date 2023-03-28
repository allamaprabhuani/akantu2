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

#ifndef AKANTU_MATERIAL_COHESIVE_LINEAR_HH_
#define AKANTU_MATERIAL_COHESIVE_LINEAR_HH_

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
template <Int dim> class MaterialCohesiveLinear : public MaterialCohesive {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialCohesiveLinear(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material parameters
  void initMaterial() override;

  void updateInternalParameters() override;

  /// check stress for cohesive elements' insertion
  void checkInsertion(bool check_only = false) override;

  /// compute effective stress norm for insertion check
  template <class D1, class D2, class D3, class D4>
  Real computeEffectiveNorm(const Eigen::MatrixBase<D1> & stress,
                            const Eigen::MatrixBase<D2> & normal,
                            const Eigen::MatrixBase<D3> & tangent,
                            const Eigen::MatrixBase<D4> & normal_stress) const;

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
  /// beta parameter
  Real beta;

  /// beta square inverse to compute effective norm
  Real beta2_inv;

  /// mode I fracture energy
  Real G_c;

  /// kappa parameter
  Real kappa;

  /// constitutive law scalar to compute delta
  Real beta2_kappa2;

  /// constitutive law scalar to compute traction
  Real beta2_kappa;

  /// penalty coefficient
  Real penalty;

  /// reference volume used to scale sigma_c
  Real volume_s;

  /// weibull exponent used to scale sigma_c
  Real m_s;

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

#include "material_cohesive_linear_inline_impl.hh"

#endif /* AKANTU_MATERIAL_COHESIVE_LINEAR_HH_ */
