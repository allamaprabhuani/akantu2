/**
 * @file   contact_search_2d_explicit.hh
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @date   Wed Nov  3 15:09:36 2010
 *
 * @brief  contact search class for contact in 2d (only line1 elements)
 *
 * @section LICENSE
 *
 * <insert license here>
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
enum InterType {
  _not_def = 0,
  _no = 1,
  _yes = 2,
  _node_1 = 3,
  _node_2 = 4,
  _seg_1 = 5,
  _seg_2 = 6,
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
  PenetrationList * findPenetration(const Surface & master_surface);
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  template <typename T> inline Int getSign(T v);

  // inline Real F_LINE(UInt node1, UInt node2, UInt node3);

  InterType Detect_Intersection(UInt node1, UInt node2, UInt node3, Real *vec_surf, Real *vec_dist, Real gap);

  bool checkProjectionAdjacentFacet(PenetrationList * const pen_list, UInt facet, UInt c_facet, UInt i_node, Real * x1, Real * x2, Real * x3, Real proj, ElementType el_type);
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "contact_search_2d_explicit_inline_impl.cc"

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const ContactSearch2dExplicit & _this)
// {
//   _this.printself(stream);
//   return stream;
// }


__END_AKANTU__



#endif /* __AKANTU_CONTACT_SEARCH_2D_EXPLICIT_HH__ */
