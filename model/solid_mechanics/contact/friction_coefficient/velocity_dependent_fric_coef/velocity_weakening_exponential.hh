/**
 * @file   velocity_weakening_exponential.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Jun 20 14:46:57 2011
 *
 * @brief  implementation of a exponential velocity weakening friction coeff.
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

#ifndef __AKANTU_VELOCITY_WEAKENING_EXPONENTIAL_HH__
#define __AKANTU_VELOCITY_WEAKENING_EXPONENTIAL_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "velocity_dependent_fric_coef.hh"
#include "historic_velocity_fric_coef.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class VelocityWeakeningExponential : public VelocityDependentFricCoef, HistoricVelocityFricCoef {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  VelocityWeakeningExponential(ContactRigid & contact,
			       const Surface & master_surface,
			       const Real static_friction_coefficient,
			       const Real dynamic_friction_coefficient,
			       const Real power);
  
  VelocityWeakeningExponential(ContactRigid & contact,
			       const Surface & master_surface,
			       const Real static_friction_coefficient,
			       const Real dynamic_friction_coefficient,
			       const Real power,
			       const Real beta);

  virtual ~VelocityWeakeningExponential();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// no initialization of variables for friction coefficient computation
  virtual void initializeComputeFricCoef();

  /// fill table with friction coefficient
  __aka_inline__ Real computeFricCoef(UInt impactor_node_index);

  /// compute the alpha parameter
  __aka_inline__ void computeAlpha();

  /// function to print the contain of the class
  //virtual void printself(std::ostream & stream, int indent = 0) const;
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  // use instant_velocity
  bool instant_velocity;

protected:
  /// static friction coefficient
  Real static_friction_coefficient;

  /// dynamic friction coefficient
  Real dynamic_friction_coefficient;

  /// power of friction dissipation
  Real power;

  /// exponential parameter
  Real alpha;
};


/* -------------------------------------------------------------------------- */
/* __aka_inline__ functions                                                           */
/* -------------------------------------------------------------------------- */
__END_AKANTU__

#include "contact_rigid.hh"

__BEGIN_AKANTU__

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "velocity_weakening_exponential_inline_impl.cc"
#endif
/*
/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const VelocityWeakeningExponential & _this)
{
  _this.printself(stream);
  return stream;
}
*/

__END_AKANTU__

#endif /* __AKANTU_VELOCITY_WEAKENING_EXPONENTIAL_HH__ */
