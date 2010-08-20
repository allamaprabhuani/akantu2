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
 * <insert license here>
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

  virtual ~Memory() {};

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

};


/* -------------------------------------------------------------------------- */
/* Inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "aka_memory_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MEMORY_HH__ */
