/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "element.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_CONTACT_ELEMENT_HH_
#define AKANTU_CONTACT_ELEMENT_HH_

/* -------------------------------------------------------------------------- */

namespace akantu {

using SlaveType = Idx;
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

  bool operator==(const ContactElement & other) const {
    return slave == other.slave and master == other.master;
  }

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

} // namespace akantu

#endif /* AKANTU_CONTACT_ELEMENT_HH_ */
