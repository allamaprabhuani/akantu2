/**
 * @file   memory.cpp
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jun 15 11:15:53 2010
 *
 * @brief  static memory wrapper
 *
 * @section LICENSE
 *
 * <insert lisence here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "memory.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_MYFEM__

/* -------------------------------------------------------------------------- */
Memory::Memory(MemoryID memory_id) {
  static_memory = StaticMemory::getStaticMemory();
  this->memory_id = memory_id;
}

/* -------------------------------------------------------------------------- */

__END_MYFEM__


