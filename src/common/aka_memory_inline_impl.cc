/**
 * @file   aka_memory_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 15 00:22:48 2010
 *
 * @brief  Implementation of the inline functions of the class Memory
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
template<class T> inline Vector<T> & Memory::alloc(const ID & name,
						    UInt size,
						    UInt nb_component) {
  handeld_vectors_id.push_back(name);
  return static_memory->smalloc<T>(memory_id, name,
				   size, nb_component);
}

/* -------------------------------------------------------------------------- */
template<class T> inline Vector<T> & Memory::alloc(const ID & name,
						   UInt size,
						   UInt nb_component,
						   const T & init_value) {
  handeld_vectors_id.push_back(name);
  return static_memory->smalloc<T>(memory_id, name,
				   size, nb_component, init_value);
}

/* -------------------------------------------------------------------------- */
inline void Memory::dealloc(const ID & name) {
  AKANTU_DEBUG(dblAccessory, "Deleting the vector " << name);
  static_memory->sfree(memory_id, name);
  handeld_vectors_id.remove(name);
}

/* -------------------------------------------------------------------------- */
template<class T> inline Vector<T> & Memory::getVector(const ID & name) {
  return static_cast< Vector<T> & >(const_cast<VectorBase &>(static_memory->getVector(memory_id, name)));
}
