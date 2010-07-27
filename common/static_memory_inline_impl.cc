/**
 * @file   static_memory_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 15 00:32:05 2010
 *
 * @brief  Implementation of inline functions of the class StaticMemory
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
inline const VectorMap & StaticMemory::getMemory(const MemoryID & memory_id) const {
  AKANTU_DEBUG_IN();
  MemoryMap::const_iterator memory_it;
  memory_it = memories.find(memory_id);

  if(memory_it == memories.end()) {
    AKANTU_DEBUG_ERROR("StaticMemory as no memory with ID " << memory_id);
  }
  AKANTU_DEBUG_OUT();
  return memory_it->second;
}

/* -------------------------------------------------------------------------- */
inline const VectorBase & StaticMemory::getVector(const MemoryID & memory_id,
						  const VectorID & name) const {
  AKANTU_DEBUG_IN();

  const VectorMap & vectors = getMemory(memory_id);

  VectorMap::const_iterator vectors_it;
  vectors_it = vectors.find(name);
  if(vectors_it == vectors.end()) {
    AKANTU_DEBUG_ERROR("StaticMemory as no array named " << name
		      << " for the Memory " << memory_id);
  }

  AKANTU_DEBUG_OUT();
  return *(vectors_it->second);
}

/* -------------------------------------------------------------------------- */
template<typename T> Vector<T> & StaticMemory::smalloc(const MemoryID & memory_id,
						       const VectorID & name,
						       UInt size,
						       UInt nb_component) {
  AKANTU_DEBUG_IN();

  MemoryMap::iterator memory_it;
  memory_it = memories.find(memory_id);

  if(memory_it == memories.end()){
    memories[memory_id] = VectorMap();
    memory_it = memories.find(memory_id);
  }

  if((memory_it->second).find(name) != (memory_it->second).end()) {
    AKANTU_DEBUG_ERROR("The vector \"" << name << "\" is already registred in the memory " << memory_id);
  }

  (memory_it->second)[name] = new Vector<T>(size, nb_component, name);

  AKANTU_DEBUG_OUT();
  return static_cast<Vector<T> &>(*(memory_it->second)[name]);
}

/* -------------------------------------------------------------------------- */
template<typename T> Vector<T> & StaticMemory::smalloc(const MemoryID & memory_id,
						       const VectorID & name,
						       UInt size,
						       UInt nb_component,
						       const T & init_value) {
  AKANTU_DEBUG_IN();

  MemoryMap::iterator memory_it;
  memory_it = memories.find(memory_id);

  if(memory_it == memories.end()){
    memories[memory_id] = VectorMap();
    memory_it = memories.find(memory_id);
  }

  if((memory_it->second).find(name) != (memory_it->second).end()) {
    AKANTU_DEBUG_ERROR("The vector \"" << name << "\" is already registred in the memory " << memory_id);
  }

  (memory_it->second)[name] = new Vector<T>(size, nb_component, init_value, name);

  AKANTU_DEBUG_OUT();
  return static_cast<Vector<T> &>(*(memory_it->second)[name]);
}

/* -------------------------------------------------------------------------- */
