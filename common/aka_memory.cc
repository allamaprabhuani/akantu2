/**
 * @file   aka_memory.cpp
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jun 15 11:15:53 2010
 *
 * @brief  static memory wrapper
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_memory.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
Memory::Memory(MemoryID memory_id) {
  static_memory = StaticMemory::getStaticMemory();
  this->memory_id = memory_id;
}

/* -------------------------------------------------------------------------- */
Memory::~Memory() {
  if(StaticMemory::isInstantiated()) {
    std::list<VectorID>::iterator it;
    for (it = handeld_vectors_id.begin(); it != handeld_vectors_id.end(); ++it) {
      AKANTU_DEBUG(dblAccessory, "Deleting the vector " << *it);
      static_memory->sfree(memory_id, *it);
    }
  }

  handeld_vectors_id.clear();
  static_memory->destroy();
  static_memory = NULL;
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__


