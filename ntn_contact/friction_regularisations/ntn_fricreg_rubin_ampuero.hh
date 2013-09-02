/**
 * @file   ntn_fricreg_rubin_ampuero.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Sep  2 16:25:32 2013
 *
 * @brief  regularisation that regularizes the contact pressure
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
#ifndef __AST_NTN_FRICREG_RUBIN_AMPUERO_HH__
#define __AST_NTN_FRICREG_RUBIN_AMPUERO_HH__

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_fricreg_no_regularisation.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
using namespace akantu;

class NTNFricRegRubinAmpuero : public NTNFricRegNoRegularisation {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  NTNFricRegRubinAmpuero(NTNBaseContact * contact,
			 const FrictionID & id = "rubin_ampuero",
			 const MemoryID & memory_id = 0);
  virtual ~NTNFricRegRubinAmpuero() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void registerSynchronizedArray(SynchronizedArrayBase & array);
  virtual void dumpRestart(const std::string & file_name) const;
  virtual void readRestart(const std::string & file_name);

  virtual void setToSteadyState();

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  virtual void addDumpFieldToDumper(const std::string & dumper_name,
				    const std::string & field_id);
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  virtual void setParam(const std::string & param, UInt node, Real value);
  virtual void setParam(const std::string & param, Real value);

protected:
  /// get the contact pressure (the norm: scalar value)
  virtual const SynchronizedArray<Real> & internalGetContactPressure();

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  SynchronizedArray<Real> t_star;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "ntn_fricreg_rubin_ampuero_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, 
				  const NTNFricRegRubinAmpuero & _this)
{
  _this.printself(stream);
  return stream;
}

__END_SIMTOOLS__

#endif /* __AST_NTN_FRICREG_RUBIN_AMPUERO_HH__ */
