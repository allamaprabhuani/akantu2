/**
 * @file   contact_search_3d_explicit.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Oct 26 18:43:27 2010
 *
 * @brief Structure that finds contact  for 3 dimensions within an explicit time
 * scheme
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONTACT_SEARCH_3D_EXPLICIT_HH__
#define __AKANTU_CONTACT_SEARCH_3D_EXPLICIT_HH__

/* -------------------------------------------------------------------------- */

#include "contact_search.hh"
#include "contact_neighbor_structure.hh"
#include "regular_grid_neighbor_structure.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class ContactSearch3dExplicit : public ContactSearch {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  ContactSearch3dExplicit(Contact & contact,
			  const ContactNeighborStructureType & neighbors_structure_type,
			  const ContactSearchType & type,
			  const ContactSearchID & id = "search_contact");

  //virtual ~ContactSearch3dExplicit();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// build the penetration list
  void findPenetration(const Surface & master_surface, PenetrationList & penetration_list);

  /// function to print the contain of the class
  //virtual void printself(std::ostream & stream, int indent = 0) const;

private:
  /// find the closest master node
  void findClosestMasterNodes(const Surface & master_surface, 
			      Vector<UInt> * closest_master_nodes, 
			      Vector<bool> * has_closest_master_node);

  /// compute the square of the distance between two nodes
  inline Real computeSquareDistanceBetweenNodes(const UInt node_1, const UInt node_2);

  /// test if impactor node is inside and in the projection area
  void checkPenetrationSituation(const UInt impactor_node, 
				 const UInt surface_element, 
				 const ElementType type,
				 bool & is_inside,
				 bool & is_in_projection_area);

  /// test if impactor node is inside and in the projection area for segment_2
  void checkPenetrationSituationSegment2(const UInt impactor_node, 
					 const UInt surface_element, 
					 bool & is_inside, 
					 bool & is_in_projection_area);

  /// test if impactor node is inside and in the projection area for triangle_3
  void checkPenetrationSituationTriangle3(const UInt impactor_node, 
					  const UInt surface_element, 
					  bool & is_inside, 
					  bool & is_in_projection_area);

  /// compute the normal, the gap and the projected position for impactor
  void computeComponentsOfProjection(const UInt impactor_node,
				     const UInt surface_element,
				     const ElementType type,
				     Real * normal,
				     Real & gap,
				     Real * projected_position);

  /// compute the normal, the gap and the projected position for impactor for segment_2
  void computeComponentsOfProjectionSegment2(const UInt impactor_node,
					     const UInt surface_element,
					     Real * normal,
					     Real & gap,
					     Real * projected_position);

  /// compute the normal, the gap and the projected position for impactor for triangle_3
  void computeComponentsOfProjectionTriangle3(const UInt impactor_node,
					      const UInt surface_element,
					      Real * normal,
					      Real & gap,
					      Real * projected_position);
  
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

#include "contact_search_3d_explicit_inline_impl.cc"

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const ContactSearch3dExplicit & _this)
// {
//   _this.printself(stream);
//   return stream;
// }

__END_AKANTU__

#endif /* __AKANTU_CONTACT_SEARCH_3D_EXPLICIT_HH__ */

