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

using SlaveType = UInt;
using MasterType = Element;
  
class ContactElement {

  /* ------------------------------------------------------------------------ */
  /* Constructor/ Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  ContactElement() = default;

  ContactElement(const SlaveType & slave, const MasterType & master)
    : slave(slave), master(master) {}
  
  ~ContactElement() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  inline UInt getNbNodes() const {
    auto nb_master_nodes = Mesh::getNbNodesPerElement(master.type);
    return nb_master_nodes + 1;
  }
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */ 
public:
  /// slave node
  SlaveType slave;
  
  /// master element/node
  MasterType master;
};

} // akantu

#endif /* __AKANTU_CONTACT_ELEMENT_HH__ */
