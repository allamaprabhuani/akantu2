/**
 * @file static_memory.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jun 10 14:19:25 2010
 *
 * @brief Memory management
 *
 * @section LICENSE
 *
 * <insert lisence here>
 *
 * @section DESCRIPTION
 *
 * The class handling the memory, allocation/reallocation/desallocation
 * The objects can register their array and ask for allocation or realocation
 */

/* -------------------------------------------------------------------------- */

#ifndef __MYFEM_STATIC_MEMORY_HH__
#define __MYFEM_STATIC_MEMORY_HH__

/* -------------------------------------------------------------------------- */
#include <map>

/* -------------------------------------------------------------------------- */
#include "common.hh"
#include "vector.hh"

/* -------------------------------------------------------------------------- */
__BEGIN_MYFEM__

typedef std::map<VectorID, VectorBase *> VectorMap;
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

  ~StaticMemory();


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// Get the global instance of the StaticMemory
  static StaticMemory * getStaticMemory();

public:

  /// access to an Vector
  inline const VectorBase & getVector(const MemoryID & memory_id,
				      const VectorID & name) const;

  /// get all vectors of a memory
  inline const VectorMap & getMemory(const MemoryID & memory_id) const;

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
   * @param nb_component number of component (for example spatial dimention)
   * @param type the type code of the array to be allocated
   *
   * @return pointer to an array of size nb_tupes * nb_component * sizeof(T)
   */
  template<typename T>
  Vector<T> & smalloc(const MemoryID & memory_id, const VectorID & name,
		      unsigned int size, unsigned int nb_component);

  /**
   * free the memory associated to the array name
   *
   * @param memory_id the id of the memory accessing to the static memory
   * @param name the name of an existing array
   */
  void sfree(const MemoryID & memory_id, const VectorID & name);


  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// unique instance of the StaticMemory
  static StaticMemory * single_static_memory;

  /// map of all allocated arrays, indexed by their names
  MemoryMap memories;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline const VectorMap & StaticMemory::getMemory(const MemoryID & memory_id) const {
  MYFEM_DEBUG_IN();
  MemoryMap::const_iterator memory_it;
  memory_it = memories.find(memory_id);

  if(memory_it == memories.end()) {
    MYFEM_DEBUG_ERROR("StaticMemory as no memory with ID " << memory_id);
  }
  MYFEM_DEBUG_OUT();
  return memory_it->second;
}

/* -------------------------------------------------------------------------- */
inline const VectorBase & StaticMemory::getVector(const MemoryID & memory_id,
						  const VectorID & name) const {
  MYFEM_DEBUG_IN();

  const VectorMap & vectors = getMemory(memory_id);

  VectorMap::const_iterator vectors_it;
  vectors_it = vectors.find(name);
  if(vectors_it == vectors.end()) {
    MYFEM_DEBUG_ERROR("StaticMemory as no array named " << name
		      << " for the Memory " << memory_id);
  }

  MYFEM_DEBUG_OUT();
  return *(vectors_it->second);
}

/* -------------------------------------------------------------------------- */
template<typename T> Vector<T> & StaticMemory::smalloc(const MemoryID & memory_id,
						       const VectorID & name,
						       unsigned int size,
						       unsigned int nb_component) {
  MYFEM_DEBUG_IN();

  MemoryMap::iterator memory_it;
  memory_it = memories.find(memory_id);

  if(memory_it == memories.end()){
    memories[memory_id] = VectorMap();
    memory_it = memories.find(memory_id);
  }

  if((memory_it->second).find(name) != (memory_it->second).end()) {
    MYFEM_DEBUG_ERROR("The vector " << name << " is already registred in the memory " << memory_id);
  }

  (memory_it->second)[name] = new Vector<T>(size, nb_component, name);

  MYFEM_DEBUG_OUT();
  return static_cast<Vector<T> &>(*(memory_it->second)[name]);
}

/* -------------------------------------------------------------------------- */

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream, const StaticMemory & _this)
{
  _this.printself(stream);
  return stream;
}


__END_MYFEM__

#endif /* __MYFEM_STATIC_MEMORY_HH__ */
