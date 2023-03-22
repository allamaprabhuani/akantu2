/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_array.hh"
#include "aka_common.hh"
#include "dumpable.hh"
#include "mesh_filter.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NODE_GROUP_HH_
#define AKANTU_NODE_GROUP_HH_

namespace akantu {

class NodeGroup : public Dumpable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NodeGroup(const std::string & name, const Mesh & mesh,
            const std::string & id = "node_group");
  ~NodeGroup() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  using const_node_iterator = Array<Idx>::const_scalar_iterator;

  /// empty the node group
  void clear();

  /// returns treu if the group is empty \warning this changed beahavior if you
  /// want to empty the group use clear
  [[nodiscard]] bool empty() const;

  /// iterator to the beginning of the node group
  [[nodiscard]] inline auto begin() const;
  /// iterator to the end of the node group
  [[nodiscard]] inline auto end() const;

  /// add a node and give the local position through an iterator
  inline auto add(Idx node, bool check_for_duplicate = true);

  /// remove a node
  inline void remove(Idx node);

  [[nodiscard]] inline decltype(auto) find(Idx node) const {
    return node_group.find(node);
  }

  /// remove duplicated nodes
  void optimize();

  /// append a group to current one
  void append(const NodeGroup & other_group);

  /// apply a filter on current node group
  template <typename T> void applyNodeFilter(T & filter);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO_AUTO_NOT_CONST(Nodes, node_group);
  AKANTU_GET_MACRO_AUTO(Nodes, node_group);
  AKANTU_GET_MACRO_AUTO(Name, name);

  /// give the number of nodes in the current group
  [[nodiscard]] inline Idx size() const;

  // UInt * storage() { return node_group.data(); };

  friend class GroupManager;
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// name of the group
  std::string name;

  /// list of nodes in the group
  Array<Idx> node_group;

  /// reference to the mesh in question
  // const Mesh & mesh;
};

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const NodeGroup & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "node_group_inline_impl.hh"

#endif /* AKANTU_NODE_GROUP_HH_ */
