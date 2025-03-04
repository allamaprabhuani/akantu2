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
#include "aka_common.hh"
#include "material_cohesive.hh"

/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_COHESIVE_EXPONENTIAL_HH_
#define AKANTU_MATERIAL_COHESIVE_EXPONENTIAL_HH_

/* -------------------------------------------------------------------------- */

namespace akantu {

/**
 * Cohesive material Exponential damage
 *
 * parameters in the material files :
 *   - sigma_c   : critical stress sigma_c  (default: 0)
 *   - beta      : weighting parameter for sliding and normal opening (default:
 * 0)
 *   - delta_c   : critical opening (default: 0)
 */
template <Int spatial_dimension>
class MaterialCohesiveExponential : public MaterialCohesive {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialCohesiveExponential(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// Initialization
  void initMaterial() override;

  /// constitutive law
  void computeTraction(ElementType el_type,
                       GhostType ghost_type = _not_ghost) override;

  /// compute the tangent stiffness matrix for an element type
  void computeTangentTraction(ElementType el_type, Array<Real> & tangent_matrix,
                              GhostType ghost_type = _not_ghost) override;

private:
  template <class D1, class D2, class D3>
  void computeCoupledTraction(Eigen::MatrixBase<D1> & tract,
                              const Eigen::MatrixBase<D2> & normal, Real delta,
                              const Eigen::MatrixBase<D3> & opening,
                              Real & delta_max_new, Real delta_max);

  template <class D1, class D2, class D3>
  void computeCompressiveTraction(Eigen::MatrixBase<D1> & tract,
                                  const Eigen::MatrixBase<D2> & normal,
                                  Real delta_n,
                                  const Eigen::MatrixBase<D3> & opening);

  template <class D1, class D2, class D3>
  void computeCoupledTangent(Eigen::MatrixBase<D1> & tangent,
                             const Eigen::MatrixBase<D2> & normal, Real delta,
                             const Eigen::MatrixBase<D3> & opening,
                             Real delta_max_new);

  template <class D1, class D2>
  void computeCompressivePenalty(Eigen::MatrixBase<D1> & tangent,
                                 const Eigen::MatrixBase<D2> & normal,
                                 Real delta_n);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// beta parameter
  Real beta;

  /// contact penalty = initial slope ?
  bool exp_penalty;

  /// Ratio of contact tangent over the initial exponential tangent
  Real contact_tangent;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

// #include "material_cohesive_exponential_inline_impl.hh"

} // namespace akantu

#endif /* AKANTU_MATERIAL_COHESIVE_EXPONENTIAL_HH_ */
