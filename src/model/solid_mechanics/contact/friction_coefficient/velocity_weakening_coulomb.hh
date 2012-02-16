/**
 * @file   velocity_weakening_coulomb.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Wed May 11 14:28:18 2011
 *
 * @brief  implementation of a velocity weakening constant friction coefficient
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

#ifndef __AKANTU_VELOCITY_WEAKENING_COULOMB_HH__
#define __AKANTU_VELOCITY_WEAKENING_COULOMB_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "friction_coefficient.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class VelocityWeakeningCoulomb : public FrictionCoefficient {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  VelocityWeakeningCoulomb(ContactRigid & contact,
			   const Surface & master_surface);

  VelocityWeakeningCoulomb(ContactRigid & contact,
			   const Surface & master_surface,
			   const Real static_friction_coefficient,
			   const Real dynamic_friction_coefficient);
  
  virtual ~VelocityWeakeningCoulomb();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// set properties
  virtual void setParam(const std::string & key, const std::string & value);
  
  /// no initialization of variables for friction coefficient computation
  virtual void initializeComputeFricCoef() {};

  /// fill table with friction coefficient
  inline virtual Real computeFricCoef(UInt impactor_node_index);
  
  /// function to print the contain of the class
  //virtual void printself(std::ostream & stream, int indent = 0) const;
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// static friction coefficient
  Real static_friction_coefficient;

  /// dynamic friction coefficient
  Real dynamic_friction_coefficient;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
__END_AKANTU__

#include "contact_rigid.hh"

__BEGIN_AKANTU__

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "velocity_weakening_coulomb_inline_impl.cc"
#endif
/*
/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const VelocityWeakeningCoulomb & _this)
{
  _this.printself(stream);
  return stream;
}
*/

__END_AKANTU__

#endif /* __AKANTU_VELOCITY_WEAKENING_COULOMB_HH__ */
