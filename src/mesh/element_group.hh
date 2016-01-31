/**
 * @file   element_group.hh
 *
 * @author Dana Christen <dana.christen@gmail.com>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri May 03 2013
 * @date last modification: Fri Sep 05 2014
 *
 * @brief  Stores information relevent to the notion of domain boundary and surfaces.
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_ELEMENT_GROUP_HH__
#define __AKANTU_ELEMENT_GROUP_HH__

#include <set>
#include "aka_common.hh"
#include "aka_memory.hh"
#include "element_type_map.hh"
#include "node_group.hh"
#include "dumpable.hh"
__BEGIN_AKANTU__

class Mesh;
class Element;

/* -------------------------------------------------------------------------- */
class ElementGroup : private Memory, public Dumpable {

   /* ------------------------------------------------------------------------ */
   /* Constructors/Destructors                                                 */
   /* ------------------------------------------------------------------------ */
 public:
  ElementGroup(const std::string & name,
	       const Mesh & mesh,
	       NodeGroup & node_group,
               UInt dimension = _all_dimensions,
	       const std::string & id = "element_group",
	       const MemoryID & memory_id = 0);

  /* ------------------------------------------------------------------------ */
  /* Type definitions                                                         */
  /* ------------------------------------------------------------------------ */
public:
  typedef ElementTypeMapArray<UInt> ElementList;
  typedef Array<UInt> NodeList;

  /* ------------------------------------------------------------------------ */
  /* Node iterator                                                            */
  /* ------------------------------------------------------------------------ */
  typedef NodeGroup::const_node_iterator const_node_iterator;
  inline const_node_iterator node_begin() const;
  inline const_node_iterator node_end() const;

  /* ------------------------------------------------------------------------ */
  /* Element iterator                                                         */
  /* ------------------------------------------------------------------------ */
  typedef ElementList::type_iterator type_iterator;
  inline type_iterator firstType(UInt dim = _all_dimensions,
				 const GhostType & ghost_type = _not_ghost,
				 const ElementKind & kind = _ek_regular) const;

  inline type_iterator lastType(UInt dim = _all_dimensions,
				const GhostType & ghost_type = _not_ghost,
				const ElementKind & kind = _ek_regular) const;

  typedef Array<UInt>::const_iterator<UInt> const_element_iterator;
  inline const_element_iterator element_begin(const ElementType & type,
					      const GhostType & ghost_type = _not_ghost) const;
  inline const_element_iterator element_end(const ElementType & type,
					    const GhostType & ghost_type = _not_ghost) const;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// empty the element group
  void empty();

  /// append another group to this group
  /// BE CAREFUL: it doesn't conserve the element order
  void append(const ElementGroup & other_group);

  /// add an element to the group. By default the it does not add the nodes to the group
  inline void add(const Element & el, bool add_nodes = false, bool check_for_duplicate = true);

  /// \todo fix the default for add_nodes : make it coherent with the other method
  inline void add(const ElementType & type, UInt element,
		  const GhostType & ghost_type = _not_ghost,
		  bool add_nodes = true, bool check_for_duplicate = true);

  inline void addNode(UInt node_id, bool check_for_duplicate = true);

  inline void removeNode(UInt node_id);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;


  /// fill the elements based on the underlying node group.
  virtual void fillFromNodeGroup();

  // sort and remove duplicated values
  void optimize();
private:
  inline void addElement(const ElementType & elem_type,
			 UInt elem_id,
			 const GhostType & ghost_type);


  friend class GroupManager;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Elements, elements, UInt);
  AKANTU_GET_MACRO(Elements, elements, const ElementTypeMapArray<UInt> &);
  AKANTU_GET_MACRO(Nodes, node_group.getNodes(), const Array<UInt> &);
  AKANTU_GET_MACRO(NodeGroup, node_group, const NodeGroup &);
  AKANTU_GET_MACRO_NOT_CONST(NodeGroup, node_group, NodeGroup &);
  AKANTU_GET_MACRO(Dimension, dimension, UInt);
  AKANTU_GET_MACRO(Name, name, std::string);
  inline UInt getNbNodes() const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// Mesh to which this group belongs
  const Mesh & mesh;

  /// name of the group
  std::string name;

  /// list of elements composing the group
  ElementList elements;

  /// sub list of nodes which are composing the elements
  NodeGroup & node_group;

  /// group dimension
  UInt dimension;

  /// empty arry for the iterator to work when an element type not present
  Array<UInt> empty_elements;
};

/// standard output stream operator
inline std::ostream & operator << (std::ostream & stream, const ElementGroup &_this) {
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#include "element.hh"
#include "element_group_inline_impl.cc"

#endif /* __AKANTU_ELEMENT_GROUP_HH__ */
