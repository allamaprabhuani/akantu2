/**
 * @file   boundary_inline_impl.cc
 *
 * @author Dana Christen <dana.christen@gmail.com>
 *
 * @date   Wed Mar 06 09:30:00 2013
 *
 * @brief  Stores information relevent to the notion of domain boundary and surfaces.
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
inline void SubBoundary::addNode(UInt node_id) {
  nodes.push_back(node_id);
}

/* -------------------------------------------------------------------------- */
inline void SubBoundary::addElement(ElementType elem_type, UInt elem_id) {

  if(!(elements.exists(elem_type))) {
    elements.alloc(0,1, elem_type);
  }
  elements(elem_type).push_back(elem_id);
}

/* -------------------------------------------------------------------------- */
inline UInt SubBoundary::getNbNodes() const {
  return nodes.getSize();
}

/* -------------------------------------------------------------------------- */
inline SubBoundary::nodes_const_iterator SubBoundary::nodes_begin() const {
  return nodes_const_iterator(nodes.begin());
}

/* -------------------------------------------------------------------------- */
inline SubBoundary::nodes_const_iterator SubBoundary::nodes_end() const {
  return nodes_const_iterator(nodes.end());
}

/* -------------------------------------------------------------------------- */
inline SubBoundary::nodes_const_iterator::nodes_const_iterator(const SubBoundary::NodeList::const_iterator<UInt> & ext_iter)
: iter(ext_iter)
{}

/* -------------------------------------------------------------------------- */
inline SubBoundary::nodes_const_iterator::nodes_const_iterator(const nodes_const_iterator & other)
: iter(other.iter)
{}

/* -------------------------------------------------------------------------- */
inline SubBoundary::nodes_const_iterator SubBoundary::nodes_const_iterator::operator++(int) {
  nodes_const_iterator tmp(*this);
  operator++();
  return tmp;
}

/* -------------------------------------------------------------------------- */
inline SubBoundary::nodes_const_iterator & SubBoundary::nodes_const_iterator::operator++() {
  ++iter;
  return *this;
}

/* -------------------------------------------------------------------------- */
inline const UInt & SubBoundary::nodes_const_iterator::operator*() const {
  return *iter;
}

/* -------------------------------------------------------------------------- */
inline bool SubBoundary::nodes_const_iterator::operator==(const nodes_const_iterator & other) const {
  return (this->iter == other.iter);
}

/* -------------------------------------------------------------------------- */
inline bool SubBoundary::nodes_const_iterator::operator!=(const nodes_const_iterator & other) const {
  return (this->iter != other.iter);
}

/* -------------------------------------------------------------------------- */
inline SubBoundary::nodes_const_iterator & SubBoundary::nodes_const_iterator::operator=(const nodes_const_iterator & other) {
  if(&other != this) {
    this->iter = other.iter;
  }
  return *this;
}

/* -------------------------------------------------------------------------- */
inline void SubBoundary::registerField(const std::string field_name,
                                       DumperIOHelper::Field * field) {

  addDumpFieldToDumper(field_name, field);
}
