/**
 * @file   contact_rigid.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date   Wed Jan 19 14:45:41 2011
 *
 * @brief  Structure that  solves contact for for a  rigid master surface without
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
#include "friction_coefficient.hh"

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
    Array<UInt> * active_impactor_nodes;  

    // the offset of the associated master surface elemet
    //Array<UInt> * master_element_offset;

    /// the element type of the associated master surface element
    std::vector<ElementType> * master_element_type;
  
    /// the normal to the master surface element for each active impactor node
    Array<Real> * master_normals;
    
    /// show if node is sticking
    Array<bool> * node_is_sticking;

    /// friction force for each
    Array<Real> * friction_forces;
    
    /// stick position for regularized friction
    Array<Real> * stick_positions;

    /// residual forces without friction force
    Array<Real> * residual_forces;
    
    /// velocities before predictor computation
    Array<Real> * previous_velocities;

    /// friction resistance 
    Array<Real> * friction_resistances;
  };

  typedef std::map<Surface, ImpactorInformationPerMaster *> SurfaceToImpactInfoMap;
  typedef std::map<Surface, FrictionCoefficient *> SurfaceToFricCoefMap;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// solve the contact
  virtual void solveContact();
  
  /// avoid adhesion by delocking contact nodes that have tensile contact force
  virtual void avoidAdhesion();

  /// alternative way for addFriction and addSticking
  virtual void addRegularizedFriction(const Real & regularizer);

  /// friction predictor 
  virtual void frictionPredictor();

  /// friction corrector
  virtual void frictionCorrector();

  /// reset stick positions to current positions
  virtual void setStickPositionsToCurrentPositions(const Surface master);

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

  /// set a friction coefficient for a given master surface
  virtual void setFrictionCoefficient(FrictionCoefficient * fric_coefficient);

  /// remove a friction coefficient for a given master surface
  virtual void removeFrictionCoefficient(const Surface master);

  /// put contact information into a map which can be used for restart
  virtual void getRestartInformation(std::map<std::string, ArrayBase* > & map_to_fill, 
				     Surface master);

  /// put contact information into a map which can be used for restart
  virtual void setRestartInformation(std::map<std::string, ArrayBase* > & restart_map, 
				     Surface master);

  /// compute friction based on the simplified prakash clifton regularization (cochard & rice 2000)
  void setSimplifiedPrakashCliftonFriction(Real v_star, Real length);
  void setSimplifiedPrakashCliftonFriction(Real t_star);
  void unsetSimplifiedPrakashCliftonFriction();
  void setPrakashCliftonToSteadyState(const Surface master);

  /// function to print the contain of the class
  //virtual void printself(std::ostream & stream, int indent = 0) const;

private:
  /// solve penetration of penetrating nodes
  //void solvePenetration(const PenetrationList & penet_list);

  /// solve penetration to the closest facet
  void solvePenetrationClosestProjection(const Surface master, 
					 const PenetrationList & penet_list);

  /// find if the penetrated node is not already in list for this master with this normal
  bool isAlreadyActiveImpactor(const Surface master, 
			       const PenetrationList & penet_list, 
			       const UInt impactor_index, 
			       const ElementType facet_type, 
			       const UInt facet_offset);

  /// projects the impactor to the projected position
  void projectImpactor(const PenetrationList & penet_list, 
		       const UInt impactor_index, 
		       const ElementType facet_type, 
		       const UInt facet_offset);

  void lockImpactorNode(const Surface master,
			const PenetrationList & penet_list, 
			const UInt impactor_index, 
			const ElementType facet_type, 
			const UInt facet_offset);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the vector full of impactor information objects
  AKANTU_GET_MACRO(ImpactorsInformation, impactors_information, const SurfaceToImpactInfoMap &);

  /// get spatial dimension
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// spatial dimension of mesh
  UInt spatial_dimension;
  
  /// the mesh
  const Mesh & mesh;

  // simplified prakash clifton regularization
  bool prakash;
  bool dedontney;

  // parameters needed for fomulation of cochard and rice 2000
  Real ref_velocity;
  Real characterstic_length;
  Real t_star;

  /// list of impactor nodes info for each master surface
  SurfaceToImpactInfoMap impactors_information;

  /// list of friction coefficient objects for each master surface
  SurfaceToFricCoefMap friction_coefficient;
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
