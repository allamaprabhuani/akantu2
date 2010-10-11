/**
 * @file   contact.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author David Kammer <david.kammer@epfl.ch>
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
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
#include "contact_search.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

// namespace akantu {
//   class ContactSearch;
// }

class Contact : public Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  Contact(const SolidMechanicsModel & model,
	  ContactSearch & contact_search,
	  const ContactID & id = "contact",
	  const MemoryID & memory_id = 0);

public:
  virtual ~Contact();

  static Contact * newContact(const SolidMechanicsModel & model,
			      const ContactType & contact_type,
			      const ContactSearchType & contact_search_type,
			      const ContactNeighborStructureType & contact_neighbor_structure_type,
			      const ContactID & id = "contact");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// update the internal structures
  virtual void initContact();

  /// check if the neighbor structure need an update
  virtual void checkAndUpdate();

  /// update the internal structures
  virtual void updateContact();

  /// solve the contact
  virtual void solveContact() = 0;

  /// add a new master surface
  void addMasterSurface(const Surface & master_surface);

  /// remove a master surface
  void removeMasterSurface(const Surface & master_surface);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(ID, id, const ContactID &);

  AKANTU_GET_MACRO(Model, model, const SolidMechanicsModel &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// id of the contact class
  ContactID id;

  /// associated model
  const SolidMechanicsModel & model;

  /// contact search object
  ContactSearch * contact_search;

  /// list of master surfaces
  std::vector<Surface> master_surfaces;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "contact_inline_impl.cc"

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const Contact & _this)
// {
//   _this.printself(stream);
//   return stream;
//}


__END_AKANTU__

#endif /* __AKANTU_CONTACT_HH__ */
