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

protected:
  /// initialize the resolution
  void initialize();
  
  /// local computation of stifnness matrix due to normal stress
  void computeNormalStiffness(Vector<Real>  & n,
			      Array<Real>  & n_alpha,
			      Array<Real>  & d_alpha,
			      Matrix<Real>  & surface_matrix,
			      Real & gap);
  
  /// local computation of stiffness matrix due to frictional stress 
  void computeFrictionalStiffness(Vector<Real> & n,
				  Array<Real>  & n_alpha,
				  Array<Real> & d_alpha,
				  Real & gap);

  /// computation of frictional tractions
  void computeFrictionalTraction(Vector<Real> & tractions);
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// local computation of tangent moduli
  void computeTangentModuli(Vector<Real> & /* N */,
			    Array<Real> & /* N_alpha */,
			    Array<Real> & /* T_alpha */,
			    Array<Real> & /* D_alpha */,
			    Matrix<Real> & /* A matrix */,
			    Real & /* gap */
			    ) override;
  
  /// local computation of normal force
  void computeNormalForce(Vector<Real> &  /* force vector  */,
			  Vector<Real> &   /* n   */,
			  Real & /* gap */) override;
  
  /// local computation of friction force
  void computeFrictionForce(Vector<Real> & /* force vector  */,
			    Array<Real>  & /* D_alpha  */,
			    Real & /* gap */) override;
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// penalty parameter
  Real epsilon;
  
};
  

} // akantu


#endif /* __AKANTU_RESOLUTION_PENALTY_HH__  */
