/**
 * @file   contact_2d_explicit.hh
 *
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 *
 * @date   Fri Nov 19 14:23:18 2010
 *
 * @brief  Interface for 2d explicit contact class
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_CONTACT_2D_EXPLICIT_HH__
#define __AKANTU_CONTACT_2D_EXPLICIT_HH__

/* -------------------------------------------------------------------------- */

#include "aka_common.hh"
#include "contact.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */

class Contact2dExplicit : public Contact{
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Contact2dExplicit(const SolidMechanicsModel & model,
		    const ContactType & type,
		    const ContactID & id = "contact",
		    const MemoryID & memory_id = 0);

  virtual ~Contact2dExplicit();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void solveContact();

private:

  /// Project back nodes on penetrated segments
  void projectNodesOnSegments(PenetrationList & pen_list, Array<UInt> & nodes_index);

  /// Decompose velocities prior impact to compute normal velocities
  void computeNormalVelocities(PenetrationList & pen_list, Array<UInt> & nodes_index, Array<Real> & vel_norm);
  
  /// Decompose velocities prior impact to compute friction velocities
  void computeFrictionVelocities(PenetrationList & pen_list, Array<UInt> & nodes_index, Array<Real> & vel_norm, Array<Real> & vel_fric);

  /// Update velocities adding normal and friction components
  void updatePostImpactVelocities(PenetrationList & pen_list,Array<UInt> & nodes_index, Array<Real> & vel_norm, Array<Real> & vel_fric);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(FrictionCoefficient, friction_coefficient, const Real &);

  /// Set friction coefficient  (default value = 0)
  AKANTU_SET_MACRO(FrictionCoefficient, friction_coefficient, Real);

  /// Set coefficient of restitution (default value = 0)
  AKANTU_SET_MACRO(CoefficientOfRestitution, coefficient_of_restitution, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// friction coefficient
  Real friction_coefficient;

private:
  /// Coefficient of restitution used to compute post impact velocities
  Real coefficient_of_restitution;

};

__END_AKANTU__

#endif /* __AKANTU_CONTACT_2D_EXPLICIT_HH__ */
