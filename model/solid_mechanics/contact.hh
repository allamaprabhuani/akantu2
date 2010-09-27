/**
 * @file   contact.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Sep 27 09:47:27 2010
 *
 * @brief  Interface for contact classes
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONTACT_HH__
#define __AKANTU_CONTACT_HH__

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class Contact : public Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Contact(const SolidMechanicsModel & model, const ContactID & id = "contact", const MemoryID & memory_id = 0) :
    Memory(memory_id), id(id), model(model) {
    AKANTU_DEBUG_IN();

    AKANTU_DEBUG_OUT();
  };

  virtual ~Contact();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// update the internal structures
  virtual void updateContact() = 0;

  /// solve the contact
  virtual void solveContact() = 0;

  /// function to print the contain of the class
  //  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(ID, id, const ContactID & id);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// id of the contact class
  ContactID id;

  /// Associated model
  SolidMechanicsModel & model;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "contact_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const Contact & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_CONTACT_HH__ */
