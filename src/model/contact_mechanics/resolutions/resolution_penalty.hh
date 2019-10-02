/**
 * @file   contact_resolution_penalty.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Mon Jan 29 2018
 *
 * @brief  Material isotropic thermo-elastic
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
  
  /// local computation of stifnness matrix due to normal stress
  void computeNormalStiffness(Matrix<Real> & , Vector<Real> & ,
			      Array<Real>  & , Array<Real>  & ,
			      Matrix<Real> & , Real &);
  
  /// local computation of stiffness matrix due to frictional stress 
  void computeFrictionalStiffness(Vector<Real> &, Array<Real>  &,
				  Array<Real> &, Real & );

  /// local computation of direct stiffness matrix due to friction,
  /// this matrix is common for both stick and slip part
  Array<Real> computeCommonModuli(Array<Real> &, Array<Real> &,
				   Array<Real> &, Vector<Real> &, ContactElement &);

  /// local computaion of stiffness matrix due to stick state
  Matrix<Real> computeStickModuli(Array<Real> &, Array<Real> &, Matrix<Real> &);

  /// local computation of stiffness matrix due to slip state 
  Matrix<Real> computeSlipModuli(Array<Real> &, Array<Real> &,
				 Matrix<Real> &, ContactElement &);

  /* ------------------------------------------------------------------------ */
  /* Methods for stiffness computation                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  /// local computation of tangent moduli due to normal traction
  void computeNormalModuli(Matrix<Real> & , Array<Real>  & , Array<Real>  & ,
			   Vector<Real> & , ContactElement & ) override;

  /// local computation of tangent moduli due to frictional traction
  void computeFrictionalModuli(Matrix<Real> & , Array<Real>  & , Array<Real>  & ,
			       Array<Real>  & , Array<Real>  & , Matrix<Real> &,
			       Vector<Real> &, ContactElement & ) override;

  /* ------------------------------------------------------------------------ */
  /* Methods for force computation                                            */
  /* ------------------------------------------------------------------------ */
public:
  /// local computation of normal force due to normal contact
  void computeNormalForce(ContactElement &, Vector<Real> &) override;
  
  /// local computation of tangential force due to frictional traction 
  void computeTangentialForce(ContactElement &, Vector<Real> &) override;

protected:
  /// local computation of trial tangential traction due to friction
  void computeTrialTangentialTraction(ContactElement &, Vector<Real> &) override;

  /// local computation of tangential traction due to stick 
  void computeStickTangentialTraction(ContactElement &, Vector<Real> &, Vector<Real> &) override;

  /// local computation of tangential traction due to slip
  void computeSlipTangentialTraction(ContactElement &, Vector<Real> &, Vector<Real> &) override;
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// penalty parameter for normal stress
  Real epsilon_n;

  /// penalty parameter for tangential stress
  Real epsilon_t;
};
  

} // akantu


#endif /* __AKANTU_RESOLUTION_PENALTY_HH__  */
