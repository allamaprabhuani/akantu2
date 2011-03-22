/**
 * @file   friction_coefficient.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Feb 22 10:54:20 2011
 *
 * @brief  interface for friction coefficient computation
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

#ifndef __AKANTU_FRICTION_COEFFICIENT_HH__
#define __AKANTU_FRICTION_COEFFICIENT_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_vector.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class ContactRigid;
/* -------------------------------------------------------------------------- */

class FrictionCoefficient {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  FrictionCoefficient(ContactRigid & contact,
		      const Surface & master_surface);
  virtual ~FrictionCoefficient();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// read properties
  virtual void setParam(const std::string & key, const std::string & value);

  /// initialize variables needed to compute the friction coefficient
  virtual void initializeComputeFricCoef() = 0;

  /// add an impactor surface to this master surface
  virtual void addImpactorSurface(const Surface & impactor_surface) {};

  /// remove an impactor surface of this master surface
  virtual void removeImpactorSurface(const Surface & impactor_surface) {};

  /// fill table with friction coefficient
  void computeFrictionCoefficient(Vector<Real> & fric_coef);
  
  /// compute the friction coefficient for a given node
  virtual Real computeFricCoef(UInt impactor_node_index) = 0;

  /// function to print the contain of the class
  //virtual void printself(std::ostream & stream, int indent = 0) const;
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(MasterSurface, master_surface, const Surface &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// associated contact
  ContactRigid & contact;

  /// associated master surface
  const Surface master_surface;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

/*
#include "friction_coefficient_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const FrictionCoefficient & _this)
{
  _this.printself(stream);
  return stream;
}
*/

__END_AKANTU__

#endif /* __AKANTU_FRICTION_COEFFICIENT_HH__ */
