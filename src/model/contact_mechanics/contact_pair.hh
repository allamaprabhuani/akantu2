/**
 * @file   contact_pair.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Tue Sep 10 2018
 * @date last modification: Mon Sep 10 2018
 *
 * @brief  Model of Contact Mechanics
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_CONTACT_PAIR_HH__
#define __AKANTU_CONTACT_PAIR_HH__

/* -------------------------------------------------------------------------- */
#include "contact_facet.hh"


namespace akantu {

class ContactPair {

  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructors                                                  */
  /* ------------------------------------------------------------------------ */
public:
  ContactPair();

  ~ContactPair();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the node value of slave
  AKANTU_GET_MACRO(Slave, slave, UInt);

  /// get the master facet
  AKANTU_GET_MACRO(Master, master, ContactFacet);
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// slave node
  UInt slave;
  
  /// master node or element
  ContactFacet master;

  /// gap
  Real gap;

  /// area
  Real area;

  /// penalty;
  Real penalty;

  /// multiplier
  Real lambda;
  
};

} // akantu




#endif /* __AKANTU_CONTACT_PAIR_HH__ */
