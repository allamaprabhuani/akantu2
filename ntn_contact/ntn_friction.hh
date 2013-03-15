/**
 * @file   ntn_friction.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Wed Feb 22 11:40:08 2012
 *
 * @brief  implementation of friction for node to node contact
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

#ifndef __AST_NTN_FRICTION_HH__
#define __AST_NTN_FRICTION_HH__

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_contact.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
using namespace akantu;

/* -------------------------------------------------------------------------- */
class NTNFriction : protected Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  NTNFriction(NTNContact & contact,
	      const FrictionID & id = "friction",
	      const MemoryID & memory_id = 0);
  virtual ~NTNFriction() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// compute friction traction
  virtual void computeFrictionTraction();

  /// compute stick traction (friction traction needed to stick the nodes)
  void computeStickTraction();

  /// apply the friction force
  void applyFrictionTraction();

  /// register Syncronizedarrays for sync
  virtual void registerSyncronizedArray(SyncronizedArrayBase & array);
  
  /// dump restart file
  virtual void dumpRestart(const std::string & file_name) const;

  /// read restart file
  virtual void readRestart(const std::string & file_name);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  /// compute frictional strength according to friction law
  virtual void computeFrictionalStrength() = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Contact, contact, const NTNContact &)

  AKANTU_GET_MACRO(IsSticking,                 is_sticking, const SyncronizedArray<bool> &)
  AKANTU_GET_MACRO(FrictionalStrength, frictional_strength, const SyncronizedArray<Real> &)
  AKANTU_GET_MACRO(FrictionTraction,     friction_traction, const SyncronizedArray<Real> &)

  // set friction coefficient to all nodes
  //void setMu(Real mu);

  // set friction coefficient only to node (global index)
  //void setMu(UInt node, Real mu);

protected:
  void setInternalArray(SyncronizedArray<Real> & array, Real value);
  void setInternalArray(SyncronizedArray<Real> & array, UInt node, Real value);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  FrictionID id;

  NTNContact & contact;

  // if node is sticking
  SyncronizedArray<bool> is_sticking;
  // frictional strength
  SyncronizedArray<Real> frictional_strength;
  // friction force
  SyncronizedArray<Real> friction_traction;

  // friction coefficient
  //SyncronizedArray<Real> mu;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "ntn_friction_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const NTNFriction & _this)
{
  _this.printself(stream);
  return stream;
}

__END_SIMTOOLS__

#endif /* __AST_NTN_FRICTION_HH__ */
