/**
 * @file static_memory.hpp
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

#ifndef __MYFEM_STATIC_MEMORY__
#define __MYFEM_STATIC_MEMORY__

/* -------------------------------------------------------------------------- */
#include <map>

/* -------------------------------------------------------------------------- */
#include "common.hpp"
#include "array_info.hpp"

/* -------------------------------------------------------------------------- */
__BEGIN_MYFEM__

typedef unsigned int MemoryID;

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
  /* Accesors                                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /// Get the global instance of the StaticMemory
  static StaticMemory * getStaticMemory();

protected:

  /// access to an ArrayInfo
  inline ArrayInfo & getArrayInfo(MemoryID memory_id, std::string name);

  /* ------------------------------------------------------------------------ */
  /* Class Methods                                                            */
  /* ------------------------------------------------------------------------ */
public:
  /**
   * Allocation of an array of type
   *
   * @param memory_id the id of the memory accessing to the static memory
   * @param name name of the array (for example connectivity)
   * @param nb_tupes number of tuples (for example number of nodes)
   * @param nb_component number of component (for example spatial dimention)
   * @param type the type code of the array to be allocated
   *
   * @return pointer to an array of size nb_tupes * nb_component * sizeof(T)
   */
  void* smalloc(MemoryID memory_id, std::string name,
		unsigned int nb_tupes, unsigned int nb_component,
		TypeCode type);

  /// smalloc only for integers
  inline int* smalloc_int(MemoryID memory_id, std::string name,
			  unsigned int nb_tupes, unsigned int nb_component
			  ) {
    return static_cast<int *>(smalloc(memory_id, name, nb_tupes, nb_component, _integer));
  }

  /// smalloc only for single precision real
  inline float* smalloc_float(MemoryID memory_id, std::string name,
			      unsigned int nb_tupes, unsigned int nb_component
			      ) {
    return static_cast<float *>(smalloc(memory_id, name, nb_tupes, nb_component, _float));
  }

  /// smalloc only for double precision real
  inline double* smalloc_double(MemoryID memory_id, std::string name,
				unsigned int nb_tupes, unsigned int nb_component
				) {
    return static_cast<double *>(smalloc(memory_id, name, nb_tupes, nb_component, _double));
  }


  /**
   * Reallocation of a pre-allocated array
   *
   * @param memory_id the id of the memory accessing to the static memory
   * @param name the name of an existing array
   * @param nb_tuples the new number of tuples
   *
   * @return pointer to the reallocated array
   */
  void * srealloc(MemoryID memory_id, std::string name,
			 unsigned int nb_tuples);

  /// realloc for integers array
  inline int* srealloc_int(MemoryID memory_id, std::string name,
			   unsigned int nb_tuples) {
    MYFEM_DEBUG_ASSERT(getArrayInfo(memory_id, name).getType() == _integer,
		       "Call to srealloc_int on a non integer array ("
		       << memory_id << ":" << name
		       << ")");
    return static_cast<int *>(srealloc(memory_id, name, nb_tuples));
  }

  /// realloc for integers array
  inline int* srealloc_float(MemoryID memory_id, std::string name,
			     unsigned int nb_tuples) {
    MYFEM_DEBUG_ASSERT(getArrayInfo(memory_id, name).getType() == _float,
		       "Call to srealloc_int on a non integer array ("
		       << memory_id << ":" << name
		       << ")");
    return static_cast<int *>(srealloc(memory_id, name, nb_tuples));
  }

  /// realloc for integers array
  inline int* srealloc_double(MemoryID memory_id, std::string name,
			      unsigned int nb_tuples) {
    MYFEM_DEBUG_ASSERT(getArrayInfo(memory_id, name).getType() == _double,
		       "Call to srealloc_int on a non integer array ("
		       << memory_id << ":" << name
		       << ")");
    return static_cast<int *>(srealloc(memory_id, name, nb_tuples));
  }

  /**
   * free the memory associated to the array name
   *
   * @param memory_id the id of the memory accessing to the static memory
   * @param name the name of an existing array
   */
  void sfree(MemoryID memory_id, std::string name);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// unique instance of the StaticMemory
  static StaticMemory * single_static_memory;

  /// map of all allocated arrays, indexed by their names
  std::map< MemoryID, std::map<std::string, ArrayInfo> > array_info_maps;

};

/* -------------------------------------------------------------------------- */

__END_MYFEM__

#endif // __MYFEM_STATIC_MEMORY__
