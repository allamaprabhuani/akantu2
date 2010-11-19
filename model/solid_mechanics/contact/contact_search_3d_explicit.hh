/**
 * @file   contact_search_3d_explicit.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Oct 26 18:43:27 2010
 *
 * @brief  Structure that finds contact for 3 dimensions within an explicit time scheme
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONTACT_SEARCH_3D_EXPLICIT_HH__
#define __AKANTU_CONTACT_SEARCH_3D_EXPLICIT_HH__

/* -------------------------------------------------------------------------- */

#include "contact_search.hh"


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
  PenetrationList * findPenetration(const Surface & master_surface);
    
  /// function to print the contain of the class
  //virtual void printself(std::ostream & stream, int indent = 0) const;
  
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
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "contact_search_3d_explicit_inline_impl.cc"

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const ContactSearch3dExplicit & _this)
// {
//   _this.printself(stream);
//   return stream;
// }

__END_AKANTU__

#endif /* __AKANTU_CONTACT_SEARCH_3D_EXPLICIT_HH__ */

