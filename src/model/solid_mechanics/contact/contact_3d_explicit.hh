/**
 * @file   contact_3d_explicit.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Oct 26 18:13:05 2010
 *
 * @brief  Structure that solves contact for 3 dimensions within an explicit time scheme
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

#ifndef __AKANTU_CONTACT_3D_EXPLICIT_HH__
#define __AKANTU_CONTACT_3D_EXPLICIT_HH__

/* -------------------------------------------------------------------------- */

#include "aka_common.hh"
#include "contact.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class Contact3dExplicit : public Contact {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  Contact3dExplicit(const SolidMechanicsModel & model,
		    const ContactType & type,
		    const ContactID & id = "contact",
		    const MemoryID & memory_id = 0);
  
  virtual ~Contact3dExplicit();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// solve the contact
  void solveContact();
  
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
  
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "contact_3d_explicit_inline_impl.cc"

/// standard output stream operator
//inline std::ostream & operator <<(std::ostream & stream, const Contact3dExplicit & _this)
//{
//  _this.printself(stream);
//  return stream;
//}


__END_AKANTU__

#endif /*__AKANTU_CONTACT_3D_EXPLICIT_HH__ */
