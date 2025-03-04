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
#include "resolution.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_RESOLUTION_PENALTY_HH_
#define AKANTU_RESOLUTION_PENALTY_HH_

namespace akantu {
class ResolutionPenalty : public Resolution {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ResolutionPenalty(ContactMechanicsModel & model, const ID & id = "");

  ~ResolutionPenalty() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// initialize the resolution
  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Methods for stiffness computation                                        */
  /* ------------------------------------------------------------------------ */
protected:
  /// local computaion of stiffness matrix due to stick state
  void computeStickModuli(const ContactElement & element,
                          Matrix<Real> & stiffness);

  /// local computation of stiffness matrix due to slip state
  void computeSlipModuli(const ContactElement & element,
                         Matrix<Real> & stiffness);

public:
  /// local computation of tangent moduli due to normal traction
  void computeNormalModuli(const ContactElement & element,
                           Matrix<Real> & stiffness) override;

  /// local computation of tangent moduli due to tangential traction
  void computeTangentialModuli(const ContactElement & element,
                               Matrix<Real> & stiffness) override;

  /* ------------------------------------------------------------------------ */
  /* Methods for force computation                                            */
  /* ------------------------------------------------------------------------ */
public:
  /// local computation of normal force due to normal contact
  void computeNormalForce(const ContactElement & element,
                          Vector<Real> & force) override;

  /// local computation of tangential force due to frictional traction
  void computeTangentialForce(const ContactElement & element,
                              Vector<Real> & force) override;

protected:
  /// local computation of normal traction due to penetration
  virtual Real computeNormalTraction(const Real & gap) const;

  /// local computation of trial tangential traction due to friction
  template <typename D>
  void computeTrialTangentialTraction(const ContactElement & element,
                                      const Matrix<Real> & covariant_basis,
                                      Eigen::MatrixBase<D> & traction);

  /// local computation of tangential traction due to stick
  template <typename D1, typename D2>
  void
  computeStickTangentialTraction(const ContactElement & unused,
                                 Eigen::MatrixBase<D1> & traction_trial,
                                 Eigen::MatrixBase<D2> & traction_tangential);

  /// local computation of tangential traction due to slip
  template <typename D1, typename D2>
  void
  computeSlipTangentialTraction(const ContactElement & element,
                                const Matrix<Real> & covariant_basis,
                                Eigen::MatrixBase<D1> & traction_trial,
                                Eigen::MatrixBase<D2> & traction_tangential);

  /// local computation of tangential traction due to friction
  template <typename D>
  void computeTangentialTraction(const ContactElement & element,
                                 const Matrix<Real> & covariant_basis,
                                 Eigen::MatrixBase<D> & traction_tangential);

public:
  void beforeSolveStep() override;

  void afterSolveStep(bool converged = true) override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// penalty parameter for normal traction
  Real epsilon_n;

  /// penalty parameter for tangential traction
  Real epsilon_t;
};

} // namespace akantu

#endif /* AKANTU_RESOLUTION_PENALTY_HH_  */
