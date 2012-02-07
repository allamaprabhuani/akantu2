/**
 * @file   contact_search_2d_explicit.hh
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @date   Wed Nov  3 15:09:36 2010
 *
 * @brief  contact search class for contact in 2d (only line1 elements)
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

#ifndef __AKANTU_CONTACT_SEARCH_2D_EXPLICIT_HH__
#define __AKANTU_CONTACT_SEARCH_2D_EXPLICIT_HH__
#define PEN_TOL 1.E-10
#define PROJ_TOL 1.E-06

/* -------------------------------------------------------------------------- */

#include "aka_common.hh"
#include "contact_search.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/// type of possible intersections
enum InterType {
  _not_def = 0,
  _no = 1,
  _yes = 2,
  _node_1 = 3,
  _node_2 = 4,
  _seg_1 = 5,
  _seg_2 = 6
};
/* -------------------------------------------------------------------------- */


class ContactSearch2dExplicit : public ContactSearch {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  ContactSearch2dExplicit(Contact & contact,
			  const ContactNeighborStructureType & neighbors_structure_type,
			  const ContactSearchType & type,
			  const ContactSearchID & id = "search_contact_2d");

  virtual ~ContactSearch2dExplicit();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void findPenetration(const Surface & master_surface, PenetrationList & penetration_list);

private:
  template <typename T> __aka_inline__ Int getSign(T v);

  // __aka_inline__ Real F_LINE(UInt node1, UInt node2, UInt node3);

  InterType Detect_Intersection(UInt node1, UInt node2, UInt node3, Real *vec_surf, Real *vec_dist, Real gap);

  bool checkProjectionAdjacentFacet(PenetrationList & pen_list, UInt facet, UInt c_facet, UInt i_node, Real old_proj, ElementType el_type);


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

};


/* -------------------------------------------------------------------------- */
/* __aka_inline__ functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "contact_search_2d_explicit_inline_impl.cc"
#endif

/// standard output stream operator
// __aka_inline__ std::ostream & operator <<(std::ostream & stream, const ContactSearch2dExplicit & _this)
// {
//   _this.printself(stream);
//   return stream;
// }


__END_AKANTU__



#endif /* __AKANTU_CONTACT_SEARCH_2D_EXPLICIT_HH__ */
