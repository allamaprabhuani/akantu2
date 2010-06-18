/**
 * @file   static_memory.cc
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
#include "static_memory.hh"

/* -------------------------------------------------------------------------- */
__BEGIN_MYFEM__

StaticMemory* StaticMemory::single_static_memory = NULL;

/* -------------------------------------------------------------------------- */
StaticMemory::~StaticMemory() {
  delete single_static_memory;
}

/* -------------------------------------------------------------------------- */
void StaticMemory::printself(std::ostream & stream, int indent) const{
  std::string space = "";
  for(int i = 0; i < indent; i++, space += MYFEM_INDENT);
  stream << space << "StaticMemory [" << std::endl;
  stream << space << " + nb memories : " << memories.size() << std::endl;

  int tot_size = 0;
  MemoryMap::const_iterator memory_it;
  for(memory_it = memories.begin(); memory_it != memories.end(); ++memory_it) {
    int mem_size = 0;

    stream << space << MYFEM_INDENT << "Memory [" << std::endl;
    stream << space << MYFEM_INDENT << " + memory id   : " << memory_it->first << std::endl;
    stream << space << MYFEM_INDENT << " + nb vectors  : " << (memory_it->second).size() << std::endl;
    VectorMap::const_iterator vector_it;
    for(vector_it = (memory_it->second).begin();
	vector_it != (memory_it->second).end();
	++vector_it) {
      (vector_it->second)->printself(stream, indent + 2);
      mem_size += (vector_it->second)->getMemorySize();
    }
    stream << space << MYFEM_INDENT << " + total size  : " << mem_size<< std::endl;
    stream << space << MYFEM_INDENT << "]" << std::endl;
    tot_size += mem_size;
  }
  stream << space << " + total size  : " << tot_size << std::endl;
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
StaticMemory * StaticMemory::getStaticMemory() {
  if(!single_static_memory) {
    single_static_memory = new StaticMemory();
  }

  return single_static_memory;
}


/* -------------------------------------------------------------------------- */
void StaticMemory::sfree(const MemoryID & memory_id,
			 const VectorID & name) {
  MYFEM_DEBUG_IN();

  VectorMap & vectors = const_cast<VectorMap &>(getMemory(memory_id));

  VectorMap::iterator vector_it;
  vector_it = vectors.find(name);
  if(vector_it == vectors.end()) {
    MYFEM_DEBUG_ERROR("StaticMemory as no array named " << name
		      << " for the Memory " << memory_id);
  }

  vectors.erase(vector_it);

  MYFEM_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

__END_MYFEM__
