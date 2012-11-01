/**
 * @file   aka_static_memory_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Nov 25 17:00:17 2011
 *
 * @brief  template part of the StaticMemory
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
template<typename T>
Vector<T> & StaticMemory::smalloc(const MemoryID & memory_id,
                                  const ID & name,
                                  UInt size,
                                  UInt nb_component) {
  AKANTU_DEBUG_IN();

  MemoryMap::iterator memory_it;
  memory_it = memories.find(memory_id);

  if(memory_it == memories.end()){
    memories[memory_id] = VectorMap();
    memory_it = memories.find(memory_id);
  }

  if((memory_it->second).find(name) != (memory_it->second).end()) {
    AKANTU_DEBUG_ERROR("The vector \"" << name << "\" is already registred in the memory " << memory_id);
  }

  Vector<T> * tmp_vect = new Vector<T>(size, nb_component, name);
  (memory_it->second)[name] = tmp_vect;

  AKANTU_DEBUG_OUT();
  return *tmp_vect;
}

/* -------------------------------------------------------------------------- */
template<typename T>
Vector<T> & StaticMemory::smalloc(const MemoryID & memory_id,
                                  const ID & name,
                                  UInt size,
                                  UInt nb_component,
                                  const T & init_value) {
  AKANTU_DEBUG_IN();

  MemoryMap::iterator memory_it;
  memory_it = memories.find(memory_id);

  if(memory_it == memories.end()){
    memories[memory_id] = VectorMap();
    memory_it = memories.find(memory_id);
  }

  if((memory_it->second).find(name) != (memory_it->second).end()) {
    AKANTU_DEBUG_ERROR("The vector \"" << name << "\" is already registred in the memory " << memory_id);
  }

  Vector<T> * tmp_vect = new Vector<T>(size, nb_component, init_value, name);
  (memory_it->second)[name] = tmp_vect;

  AKANTU_DEBUG_OUT();
  return *tmp_vect;
}

/* -------------------------------------------------------------------------- */
