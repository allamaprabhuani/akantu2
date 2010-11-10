/**
 * @file   aka_vector_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 15 00:09:33 2010
 *
 * @brief  Inline functions of the classes Vector<T> and VectorBase
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
/* Inline Functions Vector<T>                                                 */
/* -------------------------------------------------------------------------- */
template <class T> inline T & Vector<T>::at(UInt i, UInt j) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(size > 0,
		      "The vector is empty");
  AKANTU_DEBUG_ASSERT((i < size) && (j < nb_component),
		      "The value at position [" << i << "," << j
		      << "] is out of range");

  AKANTU_DEBUG_OUT();
  return values[i*nb_component + j];
}

template <class T> inline const T & Vector<T>::get(UInt i, UInt j) const{
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(size > 0,
		      "The vector is empty");
  AKANTU_DEBUG_ASSERT((i < size) && (j < nb_component),
		      "The value at position [" << i << "," << j
		      << "] is out of range");

  AKANTU_DEBUG_OUT();
  return values[i*nb_component + j];
}

/* -------------------------------------------------------------------------- */
template <class T> inline void Vector<T>::push_back(const T & value) {
  AKANTU_DEBUG_IN();
  UInt pos = size;

  resize(size+1);
  /// @todo see if with std::uninitialized_fill it allow to build vector of objects
  for (UInt i = 0; i < nb_component; ++i) {
    values[pos*nb_component + i] = value;
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T> inline void Vector<T>::push_back(const T new_elem[]) {
  AKANTU_DEBUG_IN();
  UInt pos = size;

  resize(size+1);
  for (UInt i = 0; i < nb_component; ++i) {
    values[pos*nb_component + i] = new_elem[i];
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T> inline void Vector<T>::erase(UInt i){
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT((size > 0),
		      "The vector is empty");
  AKANTU_DEBUG_ASSERT((i < size),
		      "The element at position [" << i << "] is out of range");


  if(i != (size - 1)) {
    for (UInt j = 0; j < nb_component; ++j) {
      values[i*nb_component + j] = values[(size-1)*nb_component + j];
    }
  }

  resize(size - 1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/* Inline Functions VectorBase                                                */
/* -------------------------------------------------------------------------- */

inline UInt VectorBase::getMemorySize() const {
 return allocated_size * nb_component * size_of_type;
}
