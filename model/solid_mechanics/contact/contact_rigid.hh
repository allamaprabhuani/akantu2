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
  /* Active Contact Nodes Class (per Master Surface)                          */
  /* ------------------------------------------------------------------------ */
public:

  class ImpactorInformationPerMaster {
  public:
    ImpactorInformationPerMaster(const Surface master_id, const UInt spatial_dimension);
    virtual ~ImpactorInformationPerMaster();
  public:
    /// master surface id
    Surface master_id;

    /// impactor surfaces 
    std::vector<Surface> * impactor_surfaces;

    /// list of active impactor nodes
    Vector<UInt> * active_impactor_nodes;  
  
    /// the normal to the master surface element for each active impactor node
    Vector<Int> * master_normals;
    
    /// show if node is sticking
    Vector<bool> * node_is_sticking;
  };
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// solve the contact
  virtual void solveContact();
  
  /// avoid adhesion by delocking contact nodes that have tensile contact force
  virtual void avoidAdhesion();

  /// add friction forces
  virtual void addFriction();

  /// find nodes that stick
  virtual void addSticking();

  /// add a new master surface
  virtual void addMasterSurface(const Surface & master_surface);

  /// add an impactor surface to a master surface
  virtual void addImpactorSurfaceToMasterSurface(const Surface & impactor_surface,
						 const Surface & master_surface);

  /// remove a master surface
  virtual void removeMasterSurface(const Surface & master_surface);

  /// remove an impactor surface from master surface
  virtual void removeImpactorSurfaceFromMasterSurface(const Surface & impactor_surface,
						      const Surface & master_surface);

  /// function to print the contain of the class
  //virtual void printself(std::ostream & stream, int indent = 0) const;

private:
  /// solve penetration of penetrating nodes
  //void solvePenetration(const PenetrationList & penet_list);

  /// solve penetration to the closest facet
  void solvePenetrationClosestProjection(const UInt master_index, 
					 const PenetrationList & penet_list);

  /// projects the impactor to the projected position
  void projectImpactor(const PenetrationList & penet_list, 
		       const UInt impactor_index, 
		       const ElementType facet_type, 
		       const UInt facet_offset);

  void lockImpactorNode(const UInt master_index,
			const PenetrationList & penet_list, 
			const UInt impactor_index, 
			const ElementType facet_type, 
			const UInt facet_offset);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the vector full of impactor information objects
  AKANTU_GET_MACRO(ImpactorsInformation, impactors_information, const std::vector<ImpactorInformationPerMaster *> &);

/// get the vector containing the active impactor nodes
  //AKANTU_GET_MACRO(ActiveImpactorNodes, active_impactor_nodes, const Vector<UInt> *);
  
  /// get the vector that indicates if an impactor node sticks
  //AKANTU_GET_MACRO(NodeIsSticking, node_is_sticking, const Vector<bool> *);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// spatial dimension of mesh
  UInt spatial_dimension;
  
  /// the mesh
  const Mesh & mesh;

  /// list of impactor nodes info for each master surface
  std::vector<ImpactorInformationPerMaster *> impactors_information;

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
