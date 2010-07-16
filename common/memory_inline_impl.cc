/**
 * @file   memory_inline_impl.cc
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
template<class T> inline Vector<T> & Memory::malloc(const VectorID & name,
						    UInt size,
						    UInt nb_component) {
  return static_memory->smalloc<T>(memory_id, name,
				   size, nb_component);
}

/* -------------------------------------------------------------------------- */
inline void Memory::free(const VectorID & name) {
  static_memory->sfree(memory_id, name);
}

