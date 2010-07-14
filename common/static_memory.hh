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

#ifndef __AKANTU_STATIC_MEMORY_HH__
#define __AKANTU_STATIC_MEMORY_HH__

/* -------------------------------------------------------------------------- */
#include <map>

/* -------------------------------------------------------------------------- */
#include "common.hh"
#include "vector.hh"

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

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

#include "static_memory_inline_impl.cc"

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream, const StaticMemory & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AKANTU_STATIC_MEMORY_HH__ */
