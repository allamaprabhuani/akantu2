/**
 * @file   mesh_utils_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Aug 20 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Mesh utils inline functions
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
inline bool MeshUtils::hasElement(const Array<UInt> & connectivity,
                                  const Element & el,
                                  const Vector<UInt> & nodes) {

  UInt nb_nodes_per_element = connectivity.getNbComponent();

  const Vector<UInt> el_nodes(connectivity.storage() +
                                  el.element * nb_nodes_per_element,
                              nb_nodes_per_element);
  UInt * el_nodes_end = el_nodes.storage() + nb_nodes_per_element;

  UInt n = 0;

  while (n < nodes.size() &&
         std::find(el_nodes.storage(), el_nodes_end, nodes[n]) != el_nodes_end)
    ++n;

  return (n == nodes.size());
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
inline bool
MeshUtils::removeElementsInVector(const std::vector<Element> & elem_to_remove,
                                  std::vector<Element> & elem_list) {
  if (elem_list.size() <= elem_to_remove.size())
    return true;

  auto el_it = elem_to_remove.begin();
  auto el_last = elem_to_remove.end();
  std::vector<Element>::iterator el_del;

  UInt deletions = 0;

  for (; el_it != el_last; ++el_it) {
    el_del = std::find(elem_list.begin(), elem_list.end(), *el_it);

    if (el_del != elem_list.end()) {
      elem_list.erase(el_del);
      ++deletions;
    }
  }

  AKANTU_DEBUG_ASSERT(deletions == 0 || deletions == elem_to_remove.size(),
                      "Not all elements have been erased");

  return deletions == 0;
}
