/**
 * @file   contact_rigid.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Jan 17 14:10:05 2011
 *
 * @brief Structure that  solves contact for for a  rigid master surface without
 * friction within an explicit time scheme
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

#ifndef __AKANTU_CONTACT_RIGID_HH__
#define __AKANTU_CONTACT_RIGID_HH__

/* -------------------------------------------------------------------------- */

#include "aka_common.hh"
#include "contact.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class ContactRigid : public Contact {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  ContactRigid(const SolidMechanicsModel & model,
	       const ContactType & type,
	       const ContactID & id = "contact",
	       const MemoryID & memory_id = 0);
  
  virtual ~ContactRigid();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// solve the contact
  void solveContact();
  
  /// avoid adhesion by delocking contact nodes that have tensile contact force
  void avoidAdhesion();

  /// add friction forces
  void addFriction();

  /// find nodes that stick
  void addSticking();

  /// function to print the contain of the class
  //virtual void printself(std::ostream & stream, int indent = 0) const;

private:
  /// solve penetration of penetrating nodes
  //void solvePenetration(const PenetrationList & penet_list);

  /// solve penetration to the closest facet
  void solvePenetrationClosestProjection(const PenetrationList & penet_list);

  /// projects the impactor to the projected position
  void projectImpactor(const PenetrationList & penet_list, 
		       const UInt impactor_index, 
		       const ElementType facet_type, 
		       const UInt facet_offset);

  void lockImpactorNode(const PenetrationList & penet_list, 
			const UInt impactor_index, 
			const ElementType facet_type, 
			const UInt facet_offset);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the vector containing the active impactor nodes
  AKANTU_GET_MACRO(ActiveImpactorNodes, active_impactor_nodes, const Vector<UInt> *);

  /// get the vector that indicates if an impactor node sticks
  AKANTU_GET_MACRO(NodeIsSticking, node_is_sticking, const Vector<bool> *);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// spatial dimension of mesh
  UInt spatial_dimension;
  
  /// the mesh
  const Mesh & mesh;

  /// the normal to the master surface
  Vector<Int> * master_normals;

  /// list of active impactor nodes
  Vector<UInt> * active_impactor_nodes;  
  
  /// show if node is sticking
  Vector<bool> * node_is_sticking;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "contact_rigid_inline_impl.cc"

/// standard output stream operator
//inline std::ostream & operator <<(std::ostream & stream, const ContactRigid & _this)
//{
//  _this.printself(stream);
//  return stream;
//}


__END_AKANTU__

#endif /*__AKANTU_CONTACT_RIGID_HH__ */
