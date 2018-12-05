/**
 * @file contact_element.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Wed Sep 12 2018
 * @date last modification: Tue Oct 23 2018
 *
 * @brief  Mother class for all detection algorithms
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
#include "aka_common.hh"


/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONTACT_ELEMENT_HH__
#define __AKANTU_CONTACT_ELEMENT_HH__

/* -------------------------------------------------------------------------- */

namespace akantu {

template <typename S, typename M>
class ContactElement {

  /* ------------------------------------------------------------------------ */
  /* Constructor/ Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  ContactElement(S slave, M master)
    : slave(slave), master(master) {

  }
  
  ~ContactElement() = default;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // sets the value of normal vector
  AKANTU_SET_MACRO(Normal, Vector<Real>, normal);

  // sets the value of tangent vector
  AKANTU_SET_MACRO(Tangent, Vector<Real>, tangent);

  // sets the value of gap
  AKANTU_SET_MACRO(Gap, Real, gap);

  // sets the value of normal vector
  AKANTU_SET_MACRO(Patch, Array<M>, patch);
    
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */ 
public:  

  /// slave node
  S node;
  
  /// master element/node
  M master;

  /// normalized normal direction
  Vector<Real> normal;

  /// normalized tangent direction
  Vector<Real> tangent;

  /// penetration gap between slave and master 
  Real gap;

  /// an array of master nodes/elements around slave node
  Array<M> patch;
};



} // akantu

#endif /* __AKANTU_CONTACT_ELEMENT_HH__ */
