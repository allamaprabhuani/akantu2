/**
 * @file   rice_kuwano.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Fri Nov  4 13:18:56 2011
 *
 * @brief  implementation of a modified rice-kuwano friction law (velocity
 * weakening-strenghening friction coefficient)
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

#ifndef __AKANTU_RICE_KUWANO_HH__
#define __AKANTU_RICE_KUWANO_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "velocity_dependent_fric_coef.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__


class RiceKuwano : public VelocityDependentFricCoef {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  RiceKuwano(ContactRigid & contact,
	     const Surface & master_surface,
	     const Real static_friction_coefficient,
	     const Real dynamic_friction_coefficient,
	     const Real reference_velocity,
	     const Real alpha);
  
  virtual ~RiceKuwano();
  
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
  //  virtual void printself(std::ostream & stream, int indent = 0) const;
  
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

  /// reference velocity Vw
  Real reference_velocity;

  /// linear parameter
  Real alpha;
};

/* -------------------------------------------------------------------------- */
/* __aka_inline__ functions                                                           */
/* -------------------------------------------------------------------------- */
__END_AKANTU__

#include "contact_rigid.hh"

__BEGIN_AKANTU__

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "rice_kuwano_inline_impl.cc"
#endif

/*
/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const RiceKuwano & _this)
{
  _this.printself(stream);
  return stream;
}
*/

__END_AKANTU__


#endif /* __AKANTU_RICE_KUWANO_HH__ */
