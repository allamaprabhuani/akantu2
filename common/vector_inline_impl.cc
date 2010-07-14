/**
 * @file   vector_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 15 00:09:33 2010
 *
 * @brief  Inline functions of the classes Vector<T> and VectorBase
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
/* Inline Functions Vector<T>                                                 */
/* -------------------------------------------------------------------------- */
template <class T> inline T & Vector<T>::at(unsigned int i, unsigned int j) {
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
template <class T> inline void Vector<T>::push_back(const T new_elem[]) {
  AKANTU_DEBUG_IN();
  unsigned int pos = size;

  resize(size+1);
  for (unsigned int i = 0; i < nb_component; ++i) {
    values[pos*nb_component + i] = new_elem[i];
  }
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template <class T> inline void Vector<T>::erase(unsigned int i){
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT((size > 0),
		      "The vector is empty");
  AKANTU_DEBUG_ASSERT((i < size),
		      "The element at position [" << i << "] is out of range");


  if(i != (size - 1)) {
    for (unsigned int j = 0; j < nb_component; ++j) {
      values[i*nb_component + j] = values[(size-1)*nb_component + j];
    }
  }

  resize(size - 1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <class T> void Vector<T>::resize (unsigned int new_size) {
  AKANTU_DEBUG_IN();
  if(new_size <= allocated_size) {
    if(allocated_size - new_size > AKANTU_MIN_ALLOCATION) {
      AKANTU_DEBUG_INFO("Freeing " << (allocated_size - size)*nb_component*sizeof(T) << " bytes");

      /// Normally there are no allocation problem when reducing an array
      values = static_cast<T*>(realloc(values, new_size * nb_component * sizeof(T)));
      allocated_size = size;
    }

    size = new_size;
    AKANTU_DEBUG_OUT();
    return;
  }

  unsigned int size_to_alloc = (new_size - allocated_size < AKANTU_MIN_ALLOCATION) ?
    allocated_size + AKANTU_MIN_ALLOCATION : new_size;

  T *tmp_ptr = static_cast<T*>(realloc(values, size_to_alloc * nb_component * sizeof(T)));
  AKANTU_DEBUG_ASSERT(tmp_ptr != NULL,
		     "Cannot allocate " << size_to_alloc * nb_component * sizeof(T) << " bytes");
  if (!tmp_ptr) return;

  AKANTU_DEBUG_INFO("Allocating " << (size_to_alloc - allocated_size)*nb_component*sizeof(T) << " bytes");

  allocated_size = size_to_alloc;
  size = new_size;
  values = tmp_ptr;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T> void Vector<T>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(int i = 0; i < indent; i++, space += AKANTU_INDENT);

  int real_size = allocated_size * nb_component * size_of_type;

  stream << space << "Vector<" << typeid(T).name() << "> [" << std::endl;
  stream << space << " + id             : " << this->id << std::endl;
  stream << space << " + size           : " << this->size << std::endl;
  stream << space << " + nb_component   : " << this->nb_component << std::endl;
  stream << space << " + allocated size : " << this->allocated_size << std::endl;
  stream << space << " + memory size    : " << real_size << "B" << std::endl;
  stream << space << " + adresse        : " << std::hex << this->values << std::dec << std::endl;
  if(AKANTU_DEBUG_TEST(dblDump)) {
    stream << space << " + values         : {";
    for (unsigned int i = 0; i < this->size; ++i) {
      stream << "{";
      for (unsigned int j = 0; j < this->nb_component; ++j) {
	stream << this->values[i*nb_component + j];
	if(j != this->nb_component - 1) stream << ", ";
      }
      stream << "}";
      if(i != this->size - 1) stream << ", ";
    }
    stream << "}" << std::endl;
  }
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
/* Inline Functions VectorBase                                                */
/* -------------------------------------------------------------------------- */

inline unsigned int VectorBase::getMemorySize() const {
 return allocated_size * nb_component * size_of_type;
}
