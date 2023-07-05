/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "dumpable.hh"
#include "element_type_map.hh"
#include "node_group.hh"
/* -------------------------------------------------------------------------- */
#include <set>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_ELEMENT_GROUP_HH_
#define AKANTU_ELEMENT_GROUP_HH_

namespace akantu {
class Mesh;
class Element;
class NewNodesEvent;
} // namespace akantu

namespace akantu {

/* -------------------------------------------------------------------------- */
class AKANTU_EXPORT ElementGroup : public Dumpable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ElementGroup(const std::string & name, const Mesh & mesh,
               NodeGroup & node_group, Int dimension = _all_dimensions,
               const std::string & id = "element_group");

  ElementGroup(const ElementGroup & /*unused*/);

  /* ------------------------------------------------------------------------ */
  /* Type definitions                                                         */
  /* ------------------------------------------------------------------------ */
public:
  using ElementList = ElementTypeMapArray<Idx>;

  /* ------------------------------------------------------------------------ */
  /* Element iterator                                                         */
  /* ------------------------------------------------------------------------ */
  using type_iterator = ElementList::type_iterator;

  template <typename... pack>
  [[nodiscard]] inline decltype(auto) elementTypes(pack &&... _pack) const {
    return elements.elementTypes(_pack...);
  }

  [[nodiscard]] inline auto begin(ElementType type,
                                  GhostType ghost_type = _not_ghost) const;
  [[nodiscard]] inline auto end(ElementType type,
                                GhostType ghost_type = _not_ghost) const;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// empty the element group
  void clear();
  void clear(ElementType type, GhostType ghost_type = _not_ghost);

  [[nodiscard]] bool empty() const;

  /// append another group to this group
  /// BE CAREFUL: it doesn't conserve the element order
  void append(const ElementGroup & other_group);

  /// add an element to the group. By default the it does not add the nodes to
  /// the group
  inline void add(const Element & el, bool add_nodes = false,
                  bool check_for_duplicate = true);

  /// \todo fix the default for add_nodes : make it coherent with the other
  /// method
  inline void add(ElementType type, Idx element,
                  GhostType ghost_type = _not_ghost, bool add_nodes = true,
                  bool check_for_duplicate = true);

  inline void addNode(Idx node_id, bool check_for_duplicate = true);

  inline void removeNode(Idx node_id);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /// fill the elements based on the underlying node group.
  virtual void fillFromNodeGroup();

  // sort and remove duplicated values
  void optimize();

  /// change the dimension if needed
  void addDimension(Int dimension);

  void onNodesAdded(const Array<Idx> & new_nodes, const NewNodesEvent & event);

private:
  inline void addElement(ElementType elem_type, Idx elem_id,
                         GhostType ghost_type);

  friend class GroupManager;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  const Array<Idx> & getElements(ElementType type,
                                 GhostType ghost_type = _not_ghost) const;
  AKANTU_GET_MACRO_AUTO(Elements, elements);
  AKANTU_GET_MACRO_AUTO_NOT_CONST(Elements, elements);

  template <class... Args> auto size(Args &&... pack) const {
    return elements.size(std::forward<Args>(pack)...);
  }

  decltype(auto) getElementsIterable(ElementType type,
                                     GhostType ghost_type = _not_ghost) const;

  //  AKANTU_GET_MACRO(Nodes, node_group.getNodes(), const Array<UInt> &);

  AKANTU_GET_MACRO_AUTO(NodeGroup, node_group);
  AKANTU_GET_MACRO_AUTO_NOT_CONST(NodeGroup, node_group);

  AKANTU_GET_MACRO_AUTO(Dimension, dimension);
  AKANTU_GET_MACRO_AUTO(Name, name);

  inline Int getNbNodes() const;

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
  Int dimension{_all_dimensions};

  /// empty arry for the iterator to work when an element type not present
  Array<Idx> empty_elements;
};

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const ElementGroup & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "element.hh"
#include "element_group_inline_impl.hh"

#endif /* AKANTU_ELEMENT_GROUP_HH_ */
