/**
 * @file   material_cohesive_linear.hh
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue May 08 2012
 * @date last modification: Tue Jul 29 2014
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading with
 * random stress definition for extrinsic type
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */

#include "material_cohesive.hh"

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_COHESIVE_LINEAR_HH__
#define __AKANTU_MATERIAL_COHESIVE_LINEAR_HH__

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/**
 * Cohesive material linear damage for extrinsic case
 *
 * parameters in the material files :
 *   - sigma_c   : critical stress sigma_c  (default: 0)
 *   - beta      : weighting parameter for sliding and normal opening (default: 0)
 *   - G_cI      : fracture energy for mode I (default: 0)
 *   - G_cII     : fracture energy for mode II (default: 0)
 *   - penalty   : stiffness in compression to prevent penetration
 */
template<UInt spatial_dimension>
class MaterialCohesiveLinear : public MaterialCohesive {
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
  virtual void initMaterial();

  /// check stress for cohesive elements' insertion
  virtual void checkInsertion(bool check_only = false);

  /// compute effective stress norm for insertion check
  Real computeEffectiveNorm(const Matrix<Real> & stress,
			    const Vector<Real> & normal,
			    const Vector<Real> & tangent,
			    Vector<Real> & normal_stress) const;

protected:

  /// constitutive law
  virtual void computeTraction(const Array<Real> & normal,
			       ElementType el_type,
			       GhostType ghost_type = _not_ghost);

  /// check delta_max for cohesive elements in case of no convergence
  /// in the solveStep (only for extrinsic-implicit)
  virtual void checkDeltaMax(GhostType ghost_type = _not_ghost);

  /// reset variables when convergence is reached (only for
  /// extrinsic-implicit)
  void resetVariables(GhostType ghost_type = _not_ghost);

  /// compute tangent stiffness matrix
  virtual void computeTangentTraction(const ElementType & el_type,
				      Array<Real> & tangent_matrix,
				      const Array<Real> & normal,
				      GhostType ghost_type);

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

  /// compute the traction for a given quadrature point
  inline void computeTractionOnQuad(Vector<Real> & traction, Real & delta_max,
    const Real & delta_max_prev, const Real & delta_c,
    const Vector<Real> & insertion_stress, const Real & sigma_c,
    Vector<Real> & opening, Vector<Real> & opening_prec, const Vector<Real> & normal,
    Vector<Real> & normal_opening, Vector<Real> & tangential_opening,
    Real & normal_opening_norm, Real & tangential_opening_norm, Real & damage,
    bool & penetration, bool & reduction_penalty, Real & current_penalty,
    Vector<Real> & contact_traction, Vector<Real> & contact_opening);

  inline void computeTangentTractionOnQuad(Matrix<Real> & tangent,
    Real & delta_max, const Real & delta_c,
    const Real & sigma_c,
    Vector<Real> & opening, const Vector<Real> & normal,
    Vector<Real> & normal_opening, Vector<Real> & tangential_opening,
    Real & normal_opening_norm, Real & tangential_opening_norm, Real & damage,
    bool & penetration, bool & reduction_penalty, Real & current_penalty,
    Vector<Real> & contact_opening);

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

  /// delta of the previous step
  CohesiveInternalField<Real> opening_prec;

  /**
   * variable that defines if the penalty parameter for compression
   * has to be decreased for problems of convergence in the solution
   * loops
   */
  CohesiveInternalField<bool> reduction_penalty;

  /// variable saying if there should be penalty contact also after
  /// breaking the cohesive elements
  bool contact_after_breaking;

  /// insertion of cohesive element when stress is high enough just on
  /// one quadrature point
  bool max_quad_stress_insertion;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

__END_AKANTU__

#include "material_cohesive_linear_inline_impl.cc"

#endif /* __AKANTU_MATERIAL_COHESIVE_LINEAR_HH__ */
