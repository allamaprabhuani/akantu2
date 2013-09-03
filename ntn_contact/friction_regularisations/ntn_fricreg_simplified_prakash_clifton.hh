/**
 * @file   ntn_fricreg_simplified_prakash_clifton.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Sep  2 16:25:32 2013
 *
 * @brief  regularisation that regularizes the frictional strength with one parameter
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
#ifndef __AST_NTN_FRICREG_SIMPLIFIED_PRAKASH_CLIFTON_HH__
#define __AST_NTN_FRICREG_SIMPLIFIED_PRAKASH_CLIFTON_HH__

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_fricreg_no_regularisation.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
using namespace akantu;

class NTNFricRegSimplifiedPrakashClifton : public NTNFricRegNoRegularisation {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  NTNFricRegSimplifiedPrakashClifton(NTNBaseContact * contact,
			 const FrictionID & id = "simplified_prakash_clifton",
			 const MemoryID & memory_id = 0);
  virtual ~NTNFricRegSimplifiedPrakashClifton() {};
  
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

protected:
  /// compute frictional strength according to friction law
  virtual void computeFrictionalStrength();

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

  /// get the frictional strength array
  virtual SynchronizedArray<Real> & internalGetFrictionalStrength() {
    return this->spc_internal;
  };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  SynchronizedArray<Real> t_star;

  // to get the incremental frictional strength
  SynchronizedArray<Real> spc_internal;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "ntn_fricreg_simplified_prakash_clifton_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, 
				  const NTNFricRegSimplifiedPrakashClifton & _this)
{
  _this.printself(stream);
  return stream;
}

__END_SIMTOOLS__

#endif /* __AST_NTN_FRICREG_SIMPLIFIED_PRAKASH_CLIFTON_HH__ */
