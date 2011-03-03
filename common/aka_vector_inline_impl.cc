/**
 * @file   aka_vector_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 15 00:09:33 2010
 *
 * @brief  Inline functions of the classes Vector<T> and VectorBase
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
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
template<typename T>
template<typename Ret>
inline Vector<T>::iterator<Ret> Vector<T>::begin() {
  return iterator<Ret>(values, nb_component);
}

/* -------------------------------------------------------------------------- */
template<typename T>
template<typename Ret>
inline Vector<T>::iterator<Ret> Vector<T>::end() {
  return iterator<Ret>(values + nb_component * size, nb_component);
}


/* -------------------------------------------------------------------------- */
/* Inline Functions VectorBase                                                */
/* -------------------------------------------------------------------------- */

inline UInt VectorBase::getMemorySize() const {
 return allocated_size * nb_component * size_of_type;
}

inline void VectorBase::empty() {
  size = 0;
}


/* -------------------------------------------------------------------------- */
/* Iterators                                                                  */
/* -------------------------------------------------------------------------- */
template<typename T>
template<typename Ret>
Vector<T>::iterator<Ret>::iterator() : offset(0), ret(NULL) {
}

/* -------------------------------------------------------------------------- */
template<typename T>
template<typename Ret>
Vector<T>::iterator<Ret>::iterator(pointer_type data, UInt offset) : offset(offset), ret(new returned_type(data)) {
  AKANTU_DEBUG_ASSERT(offset == ret->size(),
		      "The iterator is not compatible with the type "
		      << typeid(returned_type).name());
}

/* -------------------------------------------------------------------------- */
template<typename T>
template<typename Ret>
Vector<T>::iterator<Ret>::iterator(const Vector<T>::iterator<Ret> & it) {
  if(this != &it) {
    this->offset = it.offset;
    this->ret = new returned_type(it.ret->values);
  }
}

/* -------------------------------------------------------------------------- */
template<typename T>
template<typename Ret>
Vector<T>::iterator<Ret>::~iterator() {
  delete ret;
}

/* -------------------------------------------------------------------------- */
template<typename T>
template<typename Ret>
inline Vector<T>::iterator<Ret> & Vector<T>::iterator<Ret>::operator++() {
  ret->values += offset;
  return *this;
}

/* -------------------------------------------------------------------------- */
template<typename T>
template<typename Ret>
inline Vector<T>::iterator<Ret> & Vector<T>::iterator<Ret>::operator=(const Vector<T>::iterator<Ret> & it) {
  if(this != &it) {
    this->offset = it.offset;
    this->ret = new returned_type(it.ret->values);
  }
  return *this;
}

/* -------------------------------------------------------------------------- */
template<typename T>
template<typename Ret>
inline Vector<T>::iterator<Ret> & Vector<T>::iterator<Ret>::operator+=(const UInt n) {
  ret->values += offset * n;
  return *this;
}

/* -------------------------------------------------------------------------- */
template<typename T>
template<typename Ret>
inline bool Vector<T>::iterator<Ret>::operator==(const iterator & other) {
  return (*this).ret->values == other.ret->values;
}

/* -------------------------------------------------------------------------- */
template<typename T>
template<typename Ret>
inline bool Vector<T>::iterator<Ret>::operator!=(const iterator & other) {
  return (*this).ret->values != other.ret->values;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template<>
template<>
class Vector<Real>::iterator<Matrix> {
public:
  iterator() : ret(NULL), offset(0) {
  };

  iterator(Real * data, UInt offset, UInt m, UInt n) : ret(new Matrix(data, m, n)), offset(offset) {
    AKANTU_DEBUG_ASSERT(offset == n*m,
			"The iterator is not compatible with the type Matrix(" << m << "," << n<< ")");
  };

  iterator(const iterator & it) {
    if(this != &it) {
      offset = it.offset;
      ret = new Matrix(it.ret->values, it.ret->m, it.ret->n);
    }
  };

  ~iterator() { delete ret; };

  inline iterator & operator=(const iterator & it) {
    if(this != &it) {
      offset = it.offset;
      ret = new Matrix(it.ret->values, it.ret->m, it.ret->n);
    }
    return *this;
  };

  inline Matrix & operator*() { return *ret; };
  inline iterator & operator++() { ret->values += offset; return *this; };

  inline iterator & operator+=(const UInt n) { ret->values += n*offset; return *this; };

  inline bool operator==(const iterator & other) { return ret->values == other.ret->values; };
  inline bool operator!=(const iterator & other) { return ret->values != other.ret->values; };

private:
    Matrix * ret;
    UInt offset;
};


/* -------------------------------------------------------------------------- */
template<>
inline Vector<Real>::iterator<Matrix> Vector<Real>::begin(UInt m, UInt n) {
  return iterator<Matrix>(values, nb_component, m, n);
}

/* -------------------------------------------------------------------------- */
template<>
inline Vector<Real>::iterator<Matrix> Vector<Real>::end(UInt m, UInt n) {
  return iterator<Matrix>(values + nb_component * size, nb_component, m, n);
}
