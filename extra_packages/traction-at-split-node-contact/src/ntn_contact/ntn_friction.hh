/**
 * @file   ntn_friction.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Jan 22 2016
 *
 * @brief  implementation of friction for node to node contact
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

#ifndef __AST_NTN_FRICTION_HH__
#define __AST_NTN_FRICTION_HH__

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_base_friction.hh"
#include "ntn_friclaw_coulomb.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template <template<class> class FrictionLaw = NTNFricLawCoulomb, 
	  class Regularisation = NTNFricRegNoRegularisation>
class NTNFriction : public FrictionLaw<Regularisation> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NTNFriction(NTNBaseContact * contact,
	      const FrictionID & id = "friction",
	      const MemoryID & memory_id = 0);
  virtual ~NTNFriction() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// apply the friction force
  virtual void applyFrictionTraction();

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  // virtual void addDumpFieldToDumper(const std::string & dumper_name,
  // 				    const std::string & field_id);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
template <template<class> class FrictionLaw, class Regularisation>
inline std::ostream & operator <<(std::ostream & stream, 
				  const NTNFriction<FrictionLaw,
				                    Regularisation> & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#include "ntn_friction_tmpl.hh"

#endif /* __AST_NTN_FRICTION_HH__ */
