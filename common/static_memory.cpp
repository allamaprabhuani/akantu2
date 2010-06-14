/**
 * @file   static_memory.cpp
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jun 10 14:42:37 2010
 *
 * @brief Memory management
 *
 * @section LICENSE
 *
 * <insert lisence here>
 *
 */

/* -------------------------------------------------------------------------- */
#include <stdexcept>
#include <sstream>

/* -------------------------------------------------------------------------- */
#include "static_memory.hpp"

/* -------------------------------------------------------------------------- */
__BEGIN_MYFEM__

StaticMemory* StaticMemory::single_static_memory = NULL;

/* -------------------------------------------------------------------------- */
StaticMemory::~StaticMemory() {
  delete single_static_memory;
}

/* -------------------------------------------------------------------------- */
StaticMemory * StaticMemory::getStaticMemory() {
  if(!single_static_memory) {
    single_static_memory = new StaticMemory();
  }

  return single_static_memory;
}


/* -------------------------------------------------------------------------- */
inline ArrayInfo & StaticMemory::getArrayInfo(MemoryID memory_id,
					      std::string name) {
  std::map< MemoryID, std::map<std::string, ArrayInfo> >::iterator maps_it;
  maps_it = array_info_maps.find(memory_id);

  if(maps_it == array_info_maps.end()) {
    MYFEM_DEBUG_ERROR("StaticMemory as no memory with ID " << memory_id);
  }

  std::map<std::string, ArrayInfo>::iterator arrays_it;
  arrays_it = maps_it->second.find(name);
  if(arrays_it == maps_it->second.end()) {
    MYFEM_DEBUG_ERROR("StaticMemory as no array named " << name
		      << " for the Memory " << memory_id);
  }

  return arrays_it->second;
}

/* -------------------------------------------------------------------------- */
void * StaticMemory::smalloc(MemoryID memory_id, std::string name,
			     unsigned int nb_tupes, unsigned int nb_component,
			     TypeCode type) {

  void *tmp_array = malloc(nb_component * nb_tupes * MYFEM_SIZEOF(type));

  MYFEM_DEBUG_ASSERT(tmp_array != NULL, "Cannot allocate "
		     << nb_component * nb_tupes * MYFEM_SIZEOF(type)
		     << " bytes for the array " << memory_id << ":" << name);
  if (!tmp_array) return NULL;

  MYFEM_DEBUG_INFO("Allocated "
		   << nb_tuples * nb_component * MYFEM_SIZEOF[type]
		   << " bytes for array " << memory_id << ":" << name)

  std::map< MemoryID, std::map<std::string, ArrayInfo> >::iterator maps_it;
  maps_it = array_info_maps.find(memory_id);

  if(maps_it == array_info_maps.end()){
    array_info_maps[memory_id] = std::map<std::string, ArrayInfo>();
    maps_it = array_info_maps.find(memory_id);
  }

  (maps_it->second)[name] = ArrayInfo(name, nb_tupes, nb_component,
				      type, tmp_array);

  return tmp_array;
}

/* -------------------------------------------------------------------------- */
void * StaticMemory::srealloc(MemoryID memory_id, std::string name,
			      unsigned int nb_tuples) {
  ArrayInfo & array_info = getArrayInfo(memory_id, name);

  unsigned int nb_reserved = array_info.getNbReservedTuples();
  unsigned int nb_component = array_info.getNbComponent();
  unsigned int size_of_type =  MYFEM_SIZEOF(array_info.getType());

  void * cur_ptr = array_info.getAdresse();

  if(nb_tuples <= nb_reserved) {
    if(nb_reserved - nb_tuples > MYFEM_MIN_ALLOCATION) { /// Normally there are no
                                                         /// allocation problem when reducing an array
      MYFEM_DEBUG_INFO("Freeing "
		       << (nb_reserved - nb_tuples)*nb_component*size_of_type
		       << " bytes from array " << memory_id << ":" << name)
      cur_ptr = realloc(cur_ptr, nb_tuples * nb_component * size_of_type);
      array_info.setNbReservedTuples(nb_tuples);
      array_info.setAdresse(cur_ptr);
    }
    array_info.setNbTuples(nb_tuples);
    return cur_ptr;
  }

  unsigned int tuples_to_alloc = (nb_tuples - nb_reserved < MYFEM_MIN_ALLOCATION) ?
    nb_reserved + MYFEM_MIN_ALLOCATION : nb_tuples;

  void * tmp_ptr = realloc(cur_ptr, tuples_to_alloc * nb_component * size_of_type);

  MYFEM_DEBUG_ASSERT(tmp_ptr != NULL, "Cannot allocate "
		     << tuples_to_alloc * nb_component * MYFEM_SIZEOF(type)
		     << " bytes for the array " << memory_id << ":" << name);

  MYFEM_DEBUG_INFO("Allocating "
		   << (tuples_to_alloc - nb_reserved) * nb_component * MYFEM_SIZEOF[type]
		   << " more bytes for array " << memory_id << ":" << name)

  if (!tmp_ptr) return NULL;

  array_info.setNbReservedTuples(nb_tuples);
  array_info.setNbTuples(nb_tuples);
  array_info.setAdresse(tmp_ptr);
  return tmp_ptr;
}

/* -------------------------------------------------------------------------- */

void StaticMemory::sfree(MemoryID memory_id,
			 std::string name) {
  std::map< MemoryID, std::map<std::string, ArrayInfo> >::iterator maps_it;
  maps_it = array_info_maps.find(memory_id);

  if(maps_it == array_info_maps.end()) {
    MYFEM_DEBUG_ERROR("StaticMemory as no memory with ID " << memory_id);
  }

  std::map<std::string, ArrayInfo>::iterator arrays_it;
  arrays_it = maps_it->second.find(name);
  if(arrays_it == maps_it->second.end()) {
    MYFEM_DEBUG_ERROR("StaticMemory as no array named " << name
		      << " for the Memory " << memory_id);
  }

  free(arrays_it->second.getAdresse());
  MYFEM_DEBUG_INFO("Freeing "
		   << arrays_it->second.getMemorySize()
		   << " bytes from array " << memory_id << ":" << name);

  maps_it->second.erase(arrays_it);
}
/* -------------------------------------------------------------------------- */


__END_MYFEM__
