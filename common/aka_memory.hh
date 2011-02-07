/**
 * @file   aka_memory.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jun 15 09:30:23 2010
 *
 * @brief  wrapper for the static_memory, all object which wants
 * to access the ststic_memory as to inherit from the class memory
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
#ifndef __AKANTU_MEMORY_HH__
#define __AKANTU_MEMORY_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_static_memory.hh"
#include "aka_vector.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__


class Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Memory(MemoryID memory_id = 0);

  virtual ~Memory();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// malloc
  template<class T>
  inline Vector<T> & alloc(const VectorID & name,
			   UInt size,
			   UInt nb_component);

  /// malloc
  template<class T>
  inline Vector<T> & alloc(const VectorID & name,
			   UInt size,
			   UInt nb_component,
			   const T & init_value);

  /* ------------------------------------------------------------------------ */
  /// free an array
  inline void dealloc(const VectorID & name);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(MemoryID, memory_id, const MemoryID &);

  template<class T>
  inline Vector<T> & getVector(const VectorID & name);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// the static memory instance
  StaticMemory * static_memory;

protected:
  /// the id registred in the static memory
  MemoryID memory_id;

  /// list of allocated vectors id
  std::list<VectorID> handeld_vectors_id;
};


/* -------------------------------------------------------------------------- */
/* Inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "aka_memory_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MEMORY_HH__ */
