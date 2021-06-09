/**
 * @file   contact_resolution_penalty.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Mon Jan 29 2018
 *
 * @brief  Linear Penalty Resolution for Contact Mechanics Model
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "resolution.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_RESOLUTION_PENALTY_HH__
#define __AKANTU_RESOLUTION_PENALTY_HH__

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
  void computeStickModuli(const ContactElement &, Matrix<Real> &);

  /// local computation of stiffness matrix due to slip state
  void computeSlipModuli(const ContactElement &, Matrix<Real> &);

public:
  /// local computation of tangent moduli due to normal traction
  void computeNormalModuli(const ContactElement &, Matrix<Real> &) override;

  /// local computation of tangent moduli due to tangential traction
  void computeTangentialModuli(const ContactElement &, Matrix<Real> &) override;

  /* ------------------------------------------------------------------------ */
  /* Methods for force computation                                            */
  /* ------------------------------------------------------------------------ */
public:
  /// local computation of normal force due to normal contact
  void computeNormalForce(const ContactElement &, Vector<Real> &) override;

  /// local computation of tangential force due to frictional traction
  void computeTangentialForce(const ContactElement &, Vector<Real> &) override;

protected:
  /// local computation of normal traction due to penetration
  Real computeNormalTraction(Real &);

  /// local computation of trial tangential traction due to friction
  void computeTrialTangentialTraction(const ContactElement &,
                                      const Matrix<Real> &, Vector<Real> &);

  /// local computation of tangential traction due to stick
  void computeStickTangentialTraction(const ContactElement &, Vector<Real> &,
                                      Vector<Real> &);

  /// local computation of tangential traction due to slip
  void computeSlipTangentialTraction(const ContactElement &,
                                     const Matrix<Real> &, Vector<Real> &,
                                     Vector<Real> &);

  /// local computation of tangential traction due to friction
  void computeTangentialTraction(const ContactElement &, const Matrix<Real> &,
                                 Vector<Real> &);

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

#endif /* __AKANTU_RESOLUTION_PENALTY_HH__  */
