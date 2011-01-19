/**
 * @file   contact_rigid_no_friction.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Jan 17 14:10:05 2011
 *
 * @brief  Structure that solves contact for for a rigid master surface without friction within an explicit time scheme
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONTACT_RIGID_NO_FRICTION_HH__
#define __AKANTU_CONTACT_RIGID_NO_FRICTION_HH__

/* -------------------------------------------------------------------------- */

#include "aka_common.hh"
#include "contact.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class ContactRigidNoFriction : public Contact {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  ContactRigidNoFriction(const SolidMechanicsModel & model,
			 const ContactType & type,
			 const ContactID & id = "contact",
			 const MemoryID & memory_id = 0);
  
  virtual ~ContactRigidNoFriction();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// solve the contact
  void solveContact();
  
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

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// spatial dimension of mesh
  UInt spatial_dimension;
  
  /// the mesh
  const Mesh & mesh;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "contact_rigid_no_friction_inline_impl.cc"

/// standard output stream operator
//inline std::ostream & operator <<(std::ostream & stream, const ContactRigidNoFriction & _this)
//{
//  _this.printself(stream);
//  return stream;
//}


__END_AKANTU__

#endif /*__AKANTU_CONTACT_RIGID_NO_FRICTION_HH__ */
