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
		      "The vector \"" << id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < size) && (j < nb_component),
		      "The value at position [" << i << "," << j
		      << "] is out of range in vector \"" << id << "\"");
  return values[i*nb_component + j];
}

/* -------------------------------------------------------------------------- */
template <typename T> inline const T & Vector<T>::operator()(UInt i, UInt j) const {
  AKANTU_DEBUG_ASSERT(size > 0,
		      "The vector \"" << id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < size) && (j < nb_component),
		      "The value at position [" << i << "," << j
		      << "] is out of range in vector \"" << id << "\"");
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
template<class T>
template<class R, class IR, class Pointer, class Reference>
class Vector<T>::iterator_internal {
public:
  typedef R         value_type;
  typedef Pointer   pointer;
  typedef Reference reference;
  typedef IR        internal_value_type;
  typedef IR*       internal_pointer;
protected:
  iterator_internal(UInt offset, pointer_type data, pointer ret) : offset(offset),
								   initial(data),
								   ret(ret) {
  }

public:
  iterator_internal() : offset(0), initial(NULL), ret(NULL) {};

  iterator_internal(pointer_type data, UInt offset)  :
    offset(offset),
    initial(data),
    ret(new value_type(data)) {
    AKANTU_DEBUG_ASSERT(offset == ret->size(),
			"The iterator_internal is not compatible with the type "
			<< typeid(value_type).name());
  };

  iterator_internal(pointer warped)  : offset(warped->size()),
				       initial(warped->storage()),
				       ret(const_cast<internal_pointer>(warped)) {
  };

  iterator_internal(const iterator_internal & it) {
    if(this != &it) {
      this->offset = it.offset;
      this->initial = it.initial;
      this->ret = new internal_value_type(*it.ret);
    }
  }

  virtual ~iterator_internal() { delete ret; };

  inline iterator_internal & operator=(const iterator_internal & it) {
    if(this != &it) {
      this->offset = it.offset;
      this->initial = it.initial;
      this->ret = new internal_value_type(*it.ret);
    }
    return *this;
  }

  inline reference operator*() { return *ret; };
  inline pointer operator->() { return ret; };
  inline iterator_internal & operator++() { ret->values += offset; return *this; };

  inline iterator_internal & operator+=(const UInt n) {
    ret->values += offset * n;
    return *this;
  }

  inline reference operator[](const UInt n) {
    ret->values = initial + n*offset;
    return *ret;
  }

  inline bool operator==(const iterator_internal & other) {
    return (*this).ret->values == other.ret->values;
  }
  inline bool operator!=(const iterator_internal & other) {
    return (*this).ret->values != other.ret->values;
  }

protected:
  UInt offset;
  pointer_type initial;
  internal_pointer ret;
};


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


// /* -------------------------------------------------------------------------- */
// /* Specialization : iterator of Matrix                                        */
// /* -------------------------------------------------------------------------- */
// template<>
// template<>
// class Vector<Real>::iterator<types::Matrix, 0> {
// public:
//   iterator() : ret(NULL), initial(NULL), offset(0) {
//   };

//   iterator(Real * data, UInt offset, UInt m, UInt n) : ret(new types::Matrix(data, m, n)), initial(data), offset(offset) {
//     AKANTU_DEBUG_ASSERT(offset == n*m,
// 			"The iterator is not compatible with the type Matrix(" << m << "," << n<< ")");
//   };

//   iterator(const iterator & it) {
//     if(this != &it) {
//       offset = it.offset;
//       ret = new types::Matrix(it.ret->values, it.ret->m, it.ret->n);
//     }
//   };

//   ~iterator() { delete ret; };

//   inline iterator & operator=(const iterator & it) {
//     if(this != &it) {
//       offset = it.offset;
//       ret = new types::Matrix(it.ret->values, it.ret->m, it.ret->n);
//     }
//     return *this;
//   };

//   inline types::Matrix & operator*() { return *ret; };
//   inline types::Matrix * operator->() { return ret; };
//   inline iterator & operator++() { ret->values += offset; return *this; };
//   inline iterator & operator+=(const UInt n) { ret->values += n*offset; return *this; };
//   inline types::Matrix & operator[](const UInt n) { ret->values = initial + n*offset; return *ret; };
//   inline bool operator==(const iterator & other) { return ret->values == other.ret->values; };
//   inline bool operator!=(const iterator & other) { return ret->values != other.ret->values; };

// private:
//   types::Matrix * ret;
//   Real * initial;
//   UInt offset;
// };

/* -------------------------------------------------------------------------- */
/* Specialization : iterator of Vector                                        */
/* -------------------------------------------------------------------------- */
// template<typename T>
// template<int fps>
// class Vector<T>::iterator<types::Vector<T>, fps> {
// public:
//   typedef types::Vector<T> return_type;
//   typedef types::Vector<T> & return_type_ref;
//   typedef types::Vector<T> * return_type_ptr;

//   iterator() : ret(NULL), initial(NULL), offset(0) {
//   };

//   iterator(T * data, UInt offset, UInt n) : ret(new types::Vector<T>(data, n)), initial(data), offset(offset) {
//     AKANTU_DEBUG_ASSERT(offset == n,
// 			"The iterator is not compatible with the type Vector(" << n<< ")");
//   };

//   iterator(const iterator & it) {
//     if(this != &it) {
//       offset = it.offset;
//       ret = new return_type(it.ret->values, it.ret->n);
//     }
//   };

//   ~iterator() { delete ret; };

//   inline iterator & operator=(const iterator & it) {
//     if(this != &it) {
//       offset = it.offset;
//       ret = new return_type(it.ret->values, it.ret->n);
//     }
//     return *this;
//   };

//   inline return_type_ref operator*() { return *ret; };
//   inline return_type_ptr operator->() { return ret; };
//   inline return_type_ref operator[](const UInt n) { ret->values = initial + n*offset; return *ret; };

//   inline iterator & operator++() { ret->values += offset; return *this; };
//   inline iterator & operator+=(const UInt n) { ret->values += n*offset; return *this; };

//   inline bool operator==(const iterator & other) { return ret->values == other.ret->values; };
//   inline bool operator!=(const iterator & other) { return ret->values != other.ret->values; };

// private:
//   return_type_ptr ret;
//   T * initial;
//   UInt offset;
// };


/* -------------------------------------------------------------------------- */
template<typename T>
inline Vector<T>::iterator< types::Vector<T> > Vector<T>::begin(UInt n) {
  AKANTU_DEBUG_ASSERT(nb_component == n,
		      "The iterator is not compatible with the type Vector("
		      << n<< ")");
  return iterator< types::Vector<T> >(new types::Vector<T>(values, n));
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline Vector<T>::iterator< types::Vector<T> > Vector<T>::end(UInt n) {
  AKANTU_DEBUG_ASSERT(nb_component == n,
		      "The iterator is not compatible with the type Vector("
		      << n<< ")");
  return iterator< types::Vector<T> >(new types::Vector<T>(values + nb_component * size,
							   n));
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline Vector<T>::const_iterator< types::Vector<T> > Vector<T>::begin(UInt n) const {
  AKANTU_DEBUG_ASSERT(nb_component == n,
		      "The iterator is not compatible with the type Vector("
		      << n<< ")");
  return const_iterator< types::Vector<T> >(new types::Vector<T>(values, n));
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline Vector<T>::const_iterator< types::Vector<T> > Vector<T>::end(UInt n) const {
  AKANTU_DEBUG_ASSERT(nb_component == n,
		      "The iterator is not compatible with the type Vector("
		      << n<< ")");
  return const_iterator< types::Vector<T> >(new types::Vector<T>(values + nb_component * size,
							   n));
}

/* -------------------------------------------------------------------------- */
template<>
inline Vector<Real>::iterator<types::Matrix> Vector<Real>::begin(UInt m, UInt n) {
  AKANTU_DEBUG_ASSERT(nb_component == n*m,
		      "The iterator is not compatible with the type Matrix("
		      << m << "," << n<< ")");
  return iterator< types::Matrix >(new types::Matrix(values, m, n));
}

/* -------------------------------------------------------------------------- */
template<>
inline Vector<Real>::iterator<types::Matrix> Vector<Real>::end(UInt m, UInt n) {
  AKANTU_DEBUG_ASSERT(nb_component == n*m,
		      "The iterator is not compatible with the type Matrix("
		      << m << "," << n<< ")");
  return iterator<types::Matrix>(new types::Matrix(values + nb_component * size, m, n));
}

/* -------------------------------------------------------------------------- */
template<>
inline Vector<Real>::const_iterator<types::Matrix> Vector<Real>::begin(UInt m, UInt n) const {
  AKANTU_DEBUG_ASSERT(nb_component == n*m,
		      "The iterator is not compatible with the type Matrix("
		      << m << "," << n<< ")");
  return const_iterator< types::Matrix >(new types::Matrix(values, m, n));
}

/* -------------------------------------------------------------------------- */
template<>
inline Vector<Real>::const_iterator<types::Matrix> Vector<Real>::end(UInt m, UInt n) const {
  AKANTU_DEBUG_ASSERT(nb_component == n*m,
		      "The iterator is not compatible with the type Matrix("
		      << m << "," << n<< ")");
  return const_iterator<types::Matrix>(new types::Matrix(values + nb_component * size, m, n));
}


// /* -------------------------------------------------------------------------- */
// template<>
// template<>
// inline Vector<Real>::iterator<Real> Vector<Real>::begin<Real>() {
//   return iterator<Real>(values);
// }

// /* -------------------------------------------------------------------------- */
// template<>
// template<>
// inline Vector<Real>::iterator<Real> Vector<Real>::end<Real>() {
//   return iterator<Real>(values + nb_component * size);
// }
