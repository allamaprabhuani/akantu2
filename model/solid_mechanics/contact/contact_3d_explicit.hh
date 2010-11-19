/**
 * @file   contact_3d_explicit.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Oct 26 18:13:05 2010
 *
 * @brief  Structure that solves contact for 3 dimensions within an explicit time scheme
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONTACT_3D_EXPLICIT_HH__
#define __AKANTU_CONTACT_3D_EXPLICIT_HH__

/* -------------------------------------------------------------------------- */

#include "aka_common.hh"
#include "contact.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class Contact3dExplicit : public Contact {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  Contact3dExplicit(const SolidMechanicsModel & model,
		    const ContactType & type,
		    const ContactID & id = "contact",
		    const MemoryID & memory_id = 0);
  
  virtual ~Contact3dExplicit();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// solve the contact
  void solveContact();
  
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

//#include "contact_3d_explicit_inline_impl.cc"

/// standard output stream operator
//inline std::ostream & operator <<(std::ostream & stream, const Contact3dExplicit & _this)
//{
//  _this.printself(stream);
//  return stream;
//}


__END_AKANTU__

#endif /*__AKANTU_CONTACT_3D_EXPLICIT_HH__ */
