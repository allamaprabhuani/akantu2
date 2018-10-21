/**
 * @file   contact_element.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Tue Sep 10 2018
 * @date last modification: Mon Sep 20 2018
 *
 * @brief  mother class for contact facets
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
#ifndef __AKANTU_CONTACT_FACET_HH__
#define __AKANTU_CONTACT_FACET_HH__

/* -------------------------------------------------------------------------- */
#include "mesh.hh"
#include "model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

class ContactFacet :
    public Element {

  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructors                                                  */
  /* ------------------------------------------------------------------------ */
public:
  ContactFacet();

  //ContactFacet(SolidMechanicsModel &, ElementType, UInt, GhostType);

  ~ContactFacet() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected :
  

public:
  ///
  inline UInt numNodes() const {
    return Mesh::getNbNodesPerElement(type);
  }
  ///
  void normal();
  ///


  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
private:

  UInt                * connectivity;
  //SolidMechanicsModel * model;
  
  
};


} // akantu

#endif
