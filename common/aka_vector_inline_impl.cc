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

/* -------------------------------------------------------------------------- */
template <typename T> inline T & Vector<T>::operator()(UInt i, UInt j) {
  AKANTU_DEBUG_ASSERT(size > 0,
		      "The vector is empty");
  AKANTU_DEBUG_ASSERT((i < size) && (j < nb_component),
		      "The value at position [" << i << "," << j
		      << "] is out of range");
  return values[i*nb_component + j];
}

/* -------------------------------------------------------------------------- */
template <typename T> inline const T & Vector<T>::operator()(UInt i, UInt j) const {
  AKANTU_DEBUG_ASSERT(size > 0,
		      "The vector is empty");
  AKANTU_DEBUG_ASSERT((i < size) && (j < nb_component),
		      "The value at position [" << i << "," << j
		      << "] is out of range");
  return values[i*nb_component + j];
}


/* -------------------------------------------------------------------------- */
template <typename T> inline T & Vector<T>::at(UInt i, UInt j) {
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
template <typename T> inline const T & Vector<T>::get(UInt i, UInt j) const{
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
template <typename T> inline void Vector<T>::push_back(const T & value) {
  //  AKANTU_DEBUG_IN();
  UInt pos = size;

  resize(size+1);
  /// @todo see if with std::uninitialized_fill it allow to build vector of objects
  for (UInt i = 0; i < nb_component; ++i) {
    values[pos*nb_component + i] = value;
  }
  //  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename T> inline void Vector<T>::push_back(const T new_elem[]) {
  //  AKANTU_DEBUG_IN();
  UInt pos = size;

  resize(size+1);
  for (UInt i = 0; i < nb_component; ++i) {
    values[pos*nb_component + i] = new_elem[i];
  }
  //  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename T> inline void Vector<T>::erase(UInt i){
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
template <typename T>
Vector<T> & Vector<T>::operator-=(const Vector<T> & vect) {
  AKANTU_DEBUG_ASSERT((size == vect.size) && (nb_component == vect.nb_component),
		      "The too vector don't have the same sizes");

  T * a = values;
  T * b = vect.values;
  for (UInt i = 0; i < size*nb_component; ++i) {
    *a -= *b;
    ++a;++b;
  }

  return *this;
}

/* -------------------------------------------------------------------------- */
template <typename T>
Vector<T> & Vector<T>::operator+=(const Vector<T> & vect) {
  AKANTU_DEBUG_ASSERT((size == vect.size) && (nb_component == vect.nb_component),
		      "The too vector don't have the same sizes");

  T * a = values;
  T * b = vect.values;
  for (UInt i = 0; i < size*nb_component; ++i) {
    *a++ += *b++;
  }

  return *this;
}

/* -------------------------------------------------------------------------- */
template <typename T>
Vector<T> & Vector<T>::operator*=(const T & alpha) {
  T * a = values;
  for (UInt i = 0; i < size*nb_component; ++i) {
    *a++ *= alpha;
  }

  return *this;
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
template<typename Ret, int fps>
Vector<T>::iterator<Ret, fps>::iterator() : offset(0), initial(NULL), ret(NULL) {
}

/* -------------------------------------------------------------------------- */
template<typename T>
template<typename Ret, int fps>
Vector<T>::iterator<Ret, fps>::iterator(pointer_type data, UInt offset) : offset(offset),
								     initial(data),
								     ret(new returned_type(data)) {
  AKANTU_DEBUG_ASSERT(offset == ret->size(),
		      "The iterator is not compatible with the type "
		      << typeid(returned_type).name());
}

/* -------------------------------------------------------------------------- */
template<typename T>
template<typename Ret, int fps>
Vector<T>::iterator<Ret, fps>::iterator(const Vector<T>::iterator<Ret, fps> & it) {
  if(this != &it) {
    this->offset = it.offset;
    this->ret = new returned_type(it.ret->values);
  }
}

/* -------------------------------------------------------------------------- */
template<typename T>
template<typename Ret, int fps>
Vector<T>::iterator<Ret, fps>::~iterator() {
  delete ret;
}

/* -------------------------------------------------------------------------- */
template<typename T>
template<typename Ret, int fps>
inline Vector<T>::iterator<Ret, fps> & Vector<T>::iterator<Ret, fps>::operator++() {
  ret->values += offset;
  return *this;
}

/* -------------------------------------------------------------------------- */
template<typename T>
template<typename Ret, int fps>
inline Vector<T>::iterator<Ret, fps> & Vector<T>::iterator<Ret, fps>::operator=(const Vector<T>::iterator<Ret, fps> & it) {
  if(this != &it) {
    this->offset = it.offset;
    this->ret = new returned_type(it.ret->values);
  }
  return *this;
}

/* -------------------------------------------------------------------------- */
template<typename T>
template<typename Ret, int fps>
inline Vector<T>::iterator<Ret, fps> & Vector<T>::iterator<Ret, fps>::operator+=(const UInt n) {
  ret->values += offset * n;
  return *this;
}

/* -------------------------------------------------------------------------- */
template<typename T>
template<typename Ret, int fps>
inline Ret & Vector<T>::iterator<Ret, fps>::operator[](const UInt n) {
 ret->values = initial + n*offset;
 return *ret;
}

/* -------------------------------------------------------------------------- */
template<typename T>
template<typename Ret, int fps>
inline bool Vector<T>::iterator<Ret, fps>::operator==(const iterator & other) {
  return (*this).ret->values == other.ret->values;
}

/* -------------------------------------------------------------------------- */
template<typename T>
template<typename Ret, int fps>
inline bool Vector<T>::iterator<Ret, fps>::operator!=(const iterator & other) {
  return (*this).ret->values != other.ret->values;
}

/* -------------------------------------------------------------------------- */
/* Specialization : iterator of Matrix                                        */
/* -------------------------------------------------------------------------- */
template<>
template<>
class Vector<Real>::iterator<types::Matrix, 0> {
public:
  iterator() : ret(NULL), initial(NULL), offset(0) {
  };

  iterator(Real * data, UInt offset, UInt m, UInt n) : ret(new types::Matrix(data, m, n)), initial(data), offset(offset) {
    AKANTU_DEBUG_ASSERT(offset == n*m,
			"The iterator is not compatible with the type Matrix(" << m << "," << n<< ")");
  };

  iterator(const iterator & it) {
    if(this != &it) {
      offset = it.offset;
      ret = new types::Matrix(it.ret->values, it.ret->m, it.ret->n);
    }
  };

  ~iterator() { delete ret; };

  inline iterator & operator=(const iterator & it) {
    if(this != &it) {
      offset = it.offset;
      ret = new types::Matrix(it.ret->values, it.ret->m, it.ret->n);
    }
    return *this;
  };

  inline types::Matrix & operator*() { return *ret; };
  inline types::Matrix * operator->() { return ret; };
  inline iterator & operator++() { ret->values += offset; return *this; };
  inline iterator & operator+=(const UInt n) { ret->values += n*offset; return *this; };
  inline types::Matrix & operator[](const UInt n) { ret->values = initial + n*offset; return *ret; };
  inline bool operator==(const iterator & other) { return ret->values == other.ret->values; };
  inline bool operator!=(const iterator & other) { return ret->values != other.ret->values; };

private:
  types::Matrix * ret;
  Real * initial;
  UInt offset;
};


/* -------------------------------------------------------------------------- */
template<>
inline Vector<Real>::iterator<types::Matrix> Vector<Real>::begin(UInt m, UInt n) {
  return iterator<types::Matrix>(values, nb_component, m, n);
}

/* -------------------------------------------------------------------------- */
template<>
inline Vector<Real>::iterator<types::Matrix> Vector<Real>::end(UInt m, UInt n) {
  return iterator<types::Matrix>(values + nb_component * size, nb_component, m, n);
}

/* -------------------------------------------------------------------------- */
/* Specialization : iterator of Vector                                        */
/* -------------------------------------------------------------------------- */
template<typename T>
template<int fps>
class Vector<T>::iterator<types::Vector<T>, fps> {
public:
  typedef types::Vector<T> return_type;
  typedef types::Vector<T> & return_type_ref;
  typedef types::Vector<T> * return_type_ptr;

  iterator() : ret(NULL), initial(NULL), offset(0) {
  };

  iterator(T * data, UInt offset, UInt n) : ret(new types::Vector<T>(data, n)), initial(data), offset(offset) {
    AKANTU_DEBUG_ASSERT(offset == n,
			"The iterator is not compatible with the type Vector(" << n<< ")");
  };

  iterator(const iterator & it) {
    if(this != &it) {
      offset = it.offset;
      ret = new return_type(it.ret->values, it.ret->n);
    }
  };

  ~iterator() { delete ret; };

  inline iterator & operator=(const iterator & it) {
    if(this != &it) {
      offset = it.offset;
      ret = new return_type(it.ret->values, it.ret->n);
    }
    return *this;
  };

  inline return_type_ref operator*() { return *ret; };
  inline return_type_ptr operator->() { return ret; };
  inline return_type_ref operator[](const UInt n) { ret->values = initial + n*offset; return *ret; };

  inline iterator & operator++() { ret->values += offset; return *this; };
  inline iterator & operator+=(const UInt n) { ret->values += n*offset; return *this; };

  inline bool operator==(const iterator & other) { return ret->values == other.ret->values; };
  inline bool operator!=(const iterator & other) { return ret->values != other.ret->values; };

private:
  return_type_ptr ret;
  T * initial;
  UInt offset;
};


/* -------------------------------------------------------------------------- */
template<typename T>
inline Vector<T>::iterator< types::Vector<T> > Vector<T>::begin(UInt n) {
  return iterator< types::Vector<T> >(values, nb_component, n);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline Vector<T>::iterator< types::Vector<T> > Vector<T>::end(UInt n) {
  return iterator< types::Vector<T> >(values + nb_component * size, nb_component, n);
}

/* -------------------------------------------------------------------------- */
/* Specialization : iterator of Scalars                                       */
/* -------------------------------------------------------------------------- */
template<>
template<>
class Vector<Real>::iterator<Real> {
public:
  iterator(Real * current) : current(current), initial(current) {
  };

  iterator(const iterator & it) {
    if(this != &it) {
      current = it.current;
    }
  };

  ~iterator() { };

  inline iterator & operator=(const iterator & it) {
    if(this != &it) {
      current = it.current;
    }
    return *this;
  };

  inline Real & operator*() { return *current; };
  inline iterator & operator++() { current++; return *this; };
  inline iterator & operator+=(const UInt n) { current += n; return *this; };
  inline Real & operator[](const UInt n) { return *(initial + n); };
  inline bool operator==(const iterator & other) { return current == other.current; };
  inline bool operator!=(const iterator & other) { return current != other.current; };

private:
  Real * current;
  Real * initial;
};

/* -------------------------------------------------------------------------- */
template<>
template<>
inline Vector<Real>::iterator<Real> Vector<Real>::begin<Real>() {
  return iterator<Real>(values);
}

/* -------------------------------------------------------------------------- */
template<>
template<> 
inline Vector<Real>::iterator<Real> Vector<Real>::end<Real>() {
  return iterator<Real>(values + nb_component * size);
}
