/**
 * @file   aka_memory_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 15 00:22:48 2010
 *
 * @brief  Implementation of the inline functions of the class Memory
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
template<class T> inline Vector<T> & Memory::alloc(const VectorID & name,
						    UInt size,
						    UInt nb_component) {
  handeld_vectors_id.push_back(name);
  return static_memory->smalloc<T>(memory_id, name,
				   size, nb_component);
}

/* -------------------------------------------------------------------------- */
template<class T> inline Vector<T> & Memory::alloc(const VectorID & name,
						   UInt size,
						   UInt nb_component,
						   const T & init_value) {
  handeld_vectors_id.push_back(name);
  return static_memory->smalloc<T>(memory_id, name,
				   size, nb_component, init_value);
}

/* -------------------------------------------------------------------------- */
inline void Memory::dealloc(const VectorID & name) {
  AKANTU_DEBUG(dblAccessory, "Deleting the vector " << name);
  static_memory->sfree(memory_id, name);
  handeld_vectors_id.remove(name);
}

/* -------------------------------------------------------------------------- */
template<class T> inline Vector<T> & Memory::getVector(const VectorID & name) {
  return static_cast< Vector<T> & >(const_cast<VectorBase &>(static_memory->getVector(memory_id, name)));
}
