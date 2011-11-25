/**
 * @file aka_static_memory.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jun 10 14:19:25 2010
 *
 * @brief Memory management
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
 * @section DESCRIPTION
 *
 * The class handling the memory, allocation/reallocation/desallocation
 * The objects can register their array and ask for allocation or realocation
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_STATIC_MEMORY_HH__
#define __AKANTU_STATIC_MEMORY_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_vector.hh"

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

typedef std::map<ID, VectorBase *> VectorMap;
typedef std::map<MemoryID, VectorMap> MemoryMap;

/**
 * @class StaticMemory
 * @brief Class for memory management common to all objects (this class as to
 * be accessed by an interface class memory)
 */
class StaticMemory {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
private:

  /// Default constructor
  StaticMemory() {};

public:

  virtual ~StaticMemory();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// Get the global instance of the StaticMemory
  static StaticMemory * getStaticMemory();

  static bool isInstantiated() { return is_instantiated; };

  /// remove a reference on the static memory
  void destroy();

  /// access to an Vector
  __aka_inline__ const VectorBase & getVector(const MemoryID & memory_id,
				      const ID & name) const;

  /// get all vectors of a memory
  __aka_inline__ const VectorMap & getMemory(const MemoryID & memory_id) const;

  /* ------------------------------------------------------------------------ */
  /* Class Methods                                                            */
  /* ------------------------------------------------------------------------ */
public:
  /**
   * Allocation of an array of type
   *
   * @param memory_id the id of the memory accessing to the static memory
   * @param name name of the array (for example connectivity)
   * @param size number of size (for example number of nodes)
   * @param nb_component number of component (for example spatial dimension)
   * @param type the type code of the array to be allocated
   *
   * @return pointer to an array of size nb_tupes * nb_component * sizeof(T)
   */
  template<typename T>
  Vector<T> & smalloc(const MemoryID & memory_id, const ID & name,
		      UInt size, UInt nb_component);

  template<typename T>
  Vector<T> & smalloc(const MemoryID & memory_id,
		      const ID & name,
		      UInt size,
		      UInt nb_component,
		      const T & init_value);
  /**
   * free the memory associated to the array name
   *
   * @param memory_id the id of the memory accessing to the static memory
   * @param name the name of an existing array
   */
  void sfree(const MemoryID & memory_id, const ID & name);


  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// is the static memory instantiated
  static bool is_instantiated;

  /// unique instance of the StaticMemory
  static StaticMemory * single_static_memory;

  /// map of all allocated arrays, indexed by their names
  MemoryMap memories;

  /// number of references on the static memory
  static UInt nb_reference;
};

#include "aka_static_memory_tmpl.hh"
#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "aka_static_memory_inline_impl.cc"
#endif

/* -------------------------------------------------------------------------- */
/* __aka_inline__ functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream, const StaticMemory & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AKANTU_STATIC_MEMORY_HH__ */
