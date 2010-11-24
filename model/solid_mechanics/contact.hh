
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
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONTACT_HH__
#define __AKANTU_CONTACT_HH__

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
//#include "contact_search.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {
  class ContactSearch;
  class PenetrationList;
}


__BEGIN_AKANTU__

class Contact : public Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  Contact(const SolidMechanicsModel & model,
	  const ContactType & type,
	  const ContactID & id = "contact",
	  const MemoryID & memory_id = 0);

public:
  virtual ~Contact();

  static Contact * newContact(const SolidMechanicsModel & model,
			      const ContactType & contact_type,
			      const ContactSearchType & contact_search_type,
			      const ContactNeighborStructureType & contact_neighbor_structure_type,
			      const ContactID & id = "contact",
			      const MemoryID & memory_id = 0);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// update the internal structures
  virtual void initContact(bool add_surfaces_flag = true);

  /// initiate the contact search structure
  virtual void initSearch();

  /// initialize all neighbor structures
  virtual void initNeighborStructure();

  /// initialize one neighbor structure
  virtual void initNeighborStructure(const Surface & master_surface);

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

  AKANTU_GET_MACRO(FrictionCoefficient, friction_coefficient, const Real &);

  AKANTU_GET_MACRO(ContactSearch, * contact_search, const ContactSearch &);

  AKANTU_GET_MACRO(SurfaceToNodesOffset, surface_to_nodes_offset, const Vector<UInt> &);

  AKANTU_GET_MACRO(SurfaceToNodes, surface_to_nodes, const Vector<UInt> &);

  AKANTU_GET_MACRO(MasterSurfaces, master_surfaces, const std::vector<Surface> &);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(NodeToElementsOffset, node_to_elements_offset, const Vector<UInt> &);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(NodeToElements, node_to_elements, const Vector<UInt> &);

  AKANTU_GET_MACRO(Type, type, const ContactType &);

  void setContactSearch(ContactSearch & contact_search) {
    this->contact_search = &contact_search;
  }

  /// Set friction coefficient  (default value = 0)
  AKANTU_SET_MACRO(FrictionCoefficient, friction_coefficient, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id of the contact class
  ContactID id;

  /// associated model
  const SolidMechanicsModel & model;

  /// friction coefficient
  Real friction_coefficient;

  /// contact search object
  ContactSearch * contact_search;

  /// list of master surfaces
  std::vector<Surface> master_surfaces;

  /// type of contact object
  ContactType type;

private:
  /// offset of nodes per surface
  Vector<UInt> surface_to_nodes_offset;

  /// list of surface nodes @warning Node can occur multiple time
  Vector<UInt> surface_to_nodes;

  /// offset of surface elements per surface node
  ByElementTypeUInt node_to_elements_offset;

  /// list of surface elements id (elements can occur multiple times)
  ByElementTypeUInt node_to_elements;

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
