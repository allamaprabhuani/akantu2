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

__END_AKANTU__

#include <memory>

__BEGIN_AKANTU__

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

  resizeUnitialized(size+1);

  /// @todo see if with std::uninitialized_fill it allow to build vector of objects

  std::uninitialized_fill_n(values + pos * nb_component, nb_component, value);

  // for (UInt i = 0; i < nb_component; ++i) {
  //   values[pos*nb_component + i] = value;
  // }

  //  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename T> inline void Vector<T>::push_back(const T new_elem[]) {
  //  AKANTU_DEBUG_IN();
  UInt pos = size;

  resizeUnitialized(size+1);

  T * tmp = values + nb_component * pos;
  std::uninitialized_copy(new_elem, new_elem + nb_component, tmp);

  // for (UInt i = 0; i < nb_component; ++i) {
  //   values[pos*nb_component + i] = new_elem[i];
  // }
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
/* Functions Vector<T>                                                        */
/* -------------------------------------------------------------------------- */
template <class T> Vector<T>::Vector (UInt size,
				      UInt nb_component,
				      const ID & id) :
  VectorBase(id), values(NULL) {
  AKANTU_DEBUG_IN();
  allocate(size, nb_component);

  T val = T();
  std::uninitialized_fill(values, values + size*nb_component, val);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T> Vector<T>::Vector (UInt size,
				      UInt nb_component,
				      const T def_values[],
				      const ID & id) :
  VectorBase(id), values(NULL) {
  AKANTU_DEBUG_IN();
  allocate(size, nb_component);

  T * tmp = values;

  for (UInt i = 0; i < size; ++i) {
    tmp = values + nb_component * i;
    std::uninitialized_copy(def_values, def_values + nb_component, tmp);
    // for (UInt j = 0; j < nb_component; ++j) {
    //   values[i*nb_component + j] = def_values[j];
    // }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T> Vector<T>::Vector (UInt size,
				      UInt nb_component,
				      const T & value,
				      const ID & id) :
  VectorBase(id), values(NULL) {
  AKANTU_DEBUG_IN();
  allocate(size, nb_component);

  std::uninitialized_fill_n(values, size*nb_component, value);

  // for (UInt i = 0; i < nb_component*size; ++i) {
  //   values[i] = value;
  // }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template <class T> Vector<T>::Vector(const Vector<T>& vect, bool deep, const ID & id) {
  AKANTU_DEBUG_IN();
  this->id = (id == "") ? vect.id : id;

  if (deep) {
    allocate(vect.size, vect.nb_component);
    T * tmp = values;
    std::uninitialized_copy(vect.values, vect.values + size * nb_component, tmp);
    // for (UInt i = 0; i < size; ++i) {
    //   for (UInt j = 0; j < nb_component; ++j) {
    // 	values[i*nb_component + j] = vect.values[i*nb_component + j];
    //   }
    // }
  } else {
    this->values = vect.values;
    this->size = vect.size;
    this->nb_component = vect.nb_component;
    this->allocated_size = vect.allocated_size;
    this->size_of_type = vect.size_of_type;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T> Vector<T>::Vector(const std::vector<T>& vect) {
  AKANTU_DEBUG_IN();
  this->id = "";

  allocate(vect.size(), 1);
  T * tmp = values;
  std::uninitialized_copy(&(vect[0]), &(vect[size-1]), tmp);

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template <class T> Vector<T>::~Vector () {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG(dblAccessory, "Freeing "
	       << allocated_size*nb_component*sizeof(T) / 1024.
	       << "kB (" << id <<")");

  if(values)
    free(values);
  size = allocated_size = 0;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T> void Vector<T>::allocate(UInt size,
					    UInt nb_component) {
  AKANTU_DEBUG_IN();
  if (size == 0){
    values = NULL;
  } else {
    values = static_cast<T*>(malloc(nb_component * size * sizeof(T)));
    AKANTU_DEBUG_ASSERT(values != NULL,
			"Cannot allocate "
			<< nb_component * size * sizeof(T) / 1024.
			<< "kB (" << id <<")");
  }

  if (values == NULL) {
    this->size = this->allocated_size = 0;
  } else {
    AKANTU_DEBUG(dblAccessory, "Allocated "
		 << size * nb_component * sizeof(T) / 1024.
		 << "kB (" << id <<")");
    this->size = this->allocated_size = size;
  }

  this->size_of_type = sizeof(T);
  this->nb_component = nb_component;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T>
void Vector<T>::resize(UInt new_size) {
  UInt old_size = size;
  resizeUnitialized(new_size);

  T val = T();
  if(size > old_size)
    std::uninitialized_fill(values + old_size*nb_component, values + size*nb_component, val);
}

/* -------------------------------------------------------------------------- */
template <class T>
void Vector<T>::resizeUnitialized(UInt new_size) {
  //  AKANTU_DEBUG_IN();
  // free some memory
  if(new_size <= allocated_size) {
    if(allocated_size - new_size > AKANTU_MIN_ALLOCATION) {
      AKANTU_DEBUG(dblAccessory, "Freeing "
		   << (allocated_size - size)*nb_component*sizeof(T) / 1024.
		   << "kB (" << id <<")");

      // Normally there are no allocation problem when reducing an array
      T * tmp_ptr = static_cast<T*>(realloc(values, new_size * nb_component * sizeof(T)));
      if(new_size != 0 && tmp_ptr == NULL) {
	AKANTU_DEBUG_ERROR("Cannot free data (" << id << ")"
			   << " [current allocated size : " << allocated_size << " | "
			   << "requested size : " << new_size << "]");
      }
      values = tmp_ptr;
      allocated_size = new_size;
    }

    size = new_size;

    //    AKANTU_DEBUG_OUT();
    return;
  }

  // allocate more memory
  UInt size_to_alloc = (new_size - allocated_size < AKANTU_MIN_ALLOCATION) ?
    allocated_size + AKANTU_MIN_ALLOCATION : new_size;

  T *tmp_ptr = static_cast<T*>(realloc(values, size_to_alloc * nb_component * sizeof(T)));
  AKANTU_DEBUG_ASSERT(tmp_ptr != NULL,
		     "Cannot allocate "
		      << size_to_alloc * nb_component * sizeof(T) / 1024.
		      << "kB");
  if (tmp_ptr == NULL) {
    AKANTU_DEBUG_ERROR("Cannot allocate more data (" << id << ")"
		       << " [current allocated size : " << allocated_size << " | "
		       << "requested size : " << new_size << "]");
  }

  AKANTU_DEBUG(dblAccessory, "Allocating "
	       << (size_to_alloc - allocated_size)*nb_component*sizeof(T) / 1024.
	       << "kB");

  allocated_size = size_to_alloc;
  size = new_size;
  values = tmp_ptr;

  //  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T> void Vector<T>::extendComponentsInterlaced(UInt multiplicator,
							      UInt block_size) {
  AKANTU_DEBUG_IN();

  if (multiplicator == 1) return;

  AKANTU_DEBUG_ASSERT(multiplicator > 1,
  		      "invalid multiplicator");
  AKANTU_DEBUG_ASSERT(nb_component%block_size == 0,
		      "stride must divide actual number of components");

  values = static_cast<T*>(realloc(values, nb_component*multiplicator*size* sizeof(T)));

  UInt new_component = nb_component/block_size * multiplicator;

  for (UInt i = 0,k=size-1; i < size; ++i,--k) {
    for (UInt j = 0; j < new_component; ++j) {
      UInt m = new_component - j -1;
      UInt n = m/multiplicator;
      for (UInt l = 0,p=block_size-1;  l < block_size; ++l,--p) {
	values[k*nb_component*multiplicator+m*block_size+p] =
	  values[k*nb_component+n*block_size+p];
      }
    }
  }

  nb_component = nb_component * multiplicator;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T> Int Vector<T>::find(const T & elem) const {
  AKANTU_DEBUG_IN();
  UInt i = 0;
  for (; (i < size) && (values[i] != elem); ++i);

  AKANTU_DEBUG_OUT();
  return (i == size) ? -1 : (Int) i;
}

/* -------------------------------------------------------------------------- */
template <class T> Int Vector<T>::find(T elem[]) const {
  AKANTU_DEBUG_IN();
  T * it = values;
  UInt i = 0;
  UInt c = 0;
  for (;i < size && (c != nb_component); ++i) {
    c = 0;
    T * cit = it;
    T * celem = elem;
    for(; (c < nb_component) && (*cit == *celem); ++c, ++cit, ++celem);
    it += nb_component;
  }
  AKANTU_DEBUG_OUT();
  return (i == size) ? -1 : (Int) i;
}


/* -------------------------------------------------------------------------- */
template <class T> void Vector<T>::copy(const Vector<T>& vect) {
  AKANTU_DEBUG_IN();

  if(AKANTU_DEBUG_TEST(dblWarning))
    if(vect.nb_component != nb_component) {
      AKANTU_DEBUG(dblWarning, "The two vectors does not have the same number of components");
    }
  //  this->id = vect.id;
  resize((vect.size * vect.nb_component) / nb_component);

  T * tmp = values;
  std::uninitialized_copy(vect.values, vect.values + size * nb_component, tmp);
  //  memcpy(this->values, vect.values, vect.size * vect.nb_component * sizeof(T));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T>
void Vector<T>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  Real real_size = allocated_size * nb_component * size_of_type / 1024.0;

  std::streamsize prec        = stream.precision();
  std::ios_base::fmtflags ff  = stream.flags();

  stream.setf (std::ios_base::showbase);
  stream.precision(2);

  stream << space << "Vector<" << debug::demangle(typeid(T).name()) << "> [" << std::endl;
  stream << space << " + id             : " << this->id << std::endl;
  stream << space << " + size           : " << this->size << std::endl;
  stream << space << " + nb_component   : " << this->nb_component << std::endl;
  stream << space << " + allocated size : " << this->allocated_size << std::endl;
  stream << space << " + memory size    : "
	 << real_size << "kB" << std::endl;
  if(!AKANTU_DEBUG_LEVEL_IS_TEST())
    stream << space << " + address        : " << std::hex << this->values
	   << std::dec << std::endl;

  stream.precision(prec);
  stream.flags(ff);

  if(AKANTU_DEBUG_TEST(dblDump)) {
    stream << space << " + values         : {";
    for (UInt i = 0; i < this->size; ++i) {
      stream << "{";
      for (UInt j = 0; j < this->nb_component; ++j) {
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

inline UInt VectorBase::getMemorySize() const {
 return allocated_size * nb_component * size_of_type;
}

inline void VectorBase::empty() {
  size = 0;
}

/** \todo change definition of iterators
template<class T>
template<class R, class IR, bool value >

template<class T>
template<class R, class IR, bool = TypeTraits<R>::isScalar >
class Vector<T>::iterator_internal : public std::iterator<std::random_access_iterator_tag, R> {
  */
  

/* -------------------------------------------------------------------------- */
/* Iterators                                                                  */
/* -------------------------------------------------------------------------- */
template<class T>
template<class R, class IR, int fps>
class Vector<T>::iterator_internal : public std::iterator<std::random_access_iterator_tag, R> {
  //class Vector<T>::iterator_internal {
public:
  typedef R   value_type;
  typedef R*  pointer;
  typedef R&  reference;
  typedef IR  internal_value_type;
  typedef IR* internal_pointer;
protected:
  iterator_internal(UInt offset, pointer_type data, pointer ret) : offset(offset),
								   initial(data),
								   ret(ret) {
  }

  ~iterator_internal() { delete ret; };

public:
  iterator_internal() : offset(0), initial(NULL), ret(NULL) {};

  iterator_internal(pointer_type data, UInt offset)  :
    offset(offset),
    initial(data),
    ret(new value_type(data)) {
    AKANTU_DEBUG_ASSERT(offset == ret->size(),
			"The iterator_internal is not compatible with the type "
			<< typeid(value_type).name());
  }

  iterator_internal(pointer warped)  : offset(warped->size()),
				       initial(warped->storage()),
				       ret(const_cast<internal_pointer>(warped)) {
  }

  iterator_internal(const iterator_internal & it) {
    if(this != &it) {
      this->offset = it.offset;
      this->initial = it.initial;
      this->ret = new internal_value_type(*it.ret);
    }
  }

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

  inline pointer_type getCurrentStorage() const {
    return ret->storage();
  }

protected:
  UInt offset;
  pointer_type initial;
  internal_pointer ret;
};


// /* -------------------------------------------------------------------------- */
// template<typename T>
// template<typename Ret>
// inline Vector<T>::iterator<Ret> Vector<T>::begin() {
//   return iterator<Ret>(values, nb_component);
// }

// /* -------------------------------------------------------------------------- */
// template<typename T>
// template<typename Ret>
// inline Vector<T>::iterator<Ret> Vector<T>::end() {
//   return iterator<Ret>(values + nb_component * size, nb_component);
// }

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

/* -------------------------------------------------------------------------- */
template<>
inline Vector<Real>::iterator< types::Matrix >
Vector<Real>::begin_reinterpret(UInt m, UInt n,
			       UInt size,
			       UInt nb_component) {
  AKANTU_DEBUG_ASSERT(nb_component * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << nb_component <<
		      ") are not compatible with the one of this vector");

  AKANTU_DEBUG_ASSERT(nb_component == n*m,
		      "The iterator is not compatible with the type Matrix("
		      << m << "," << n<< ")");

  return iterator< types::Matrix >(new types::Matrix(values, m, n));
}

/* -------------------------------------------------------------------------- */
template<>
inline Vector<Real>::iterator< types::Matrix >
Vector<Real>::end_reinterpret(UInt m, UInt n,
			      UInt size,
			      UInt nb_component) {
  AKANTU_DEBUG_ASSERT(nb_component * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << nb_component <<
		      ") are not compatible with the one of this vector");

  AKANTU_DEBUG_ASSERT(nb_component == n*m,
		      "The iterator is not compatible with the type Matrix("
		      << m << "," << n<< ")");

  return iterator<types::Matrix>(new types::Matrix(values + nb_component * size, m, n));
}


/* -------------------------------------------------------------------------- */
template<>
inline Vector<Real>::const_iterator< types::Matrix >
Vector<Real>::begin_reinterpret(UInt m, UInt n,
			       UInt size,
			       UInt nb_component) const {
  AKANTU_DEBUG_ASSERT(nb_component * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << nb_component <<
		      ") are not compatible with the one of this vector");

  AKANTU_DEBUG_ASSERT(nb_component == n*m,
		      "The iterator is not compatible with the type Matrix("
		      << m << "," << n<< ")");

  return const_iterator< types::Matrix >(new types::Matrix(values, m, n));
}

/* -------------------------------------------------------------------------- */
template<>
inline Vector<Real>::const_iterator< types::Matrix >
Vector<Real>::end_reinterpret(UInt m, UInt n,
			      UInt size,
			      UInt nb_component) const {
  AKANTU_DEBUG_ASSERT(nb_component * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << nb_component <<
		      ") are not compatible with the one of this vector");

  AKANTU_DEBUG_ASSERT(nb_component == n*m,
		      "The iterator is not compatible with the type Matrix("
		      << m << "," << n<< ")");

  return const_iterator<types::Matrix>(new types::Matrix(values + nb_component * size, m, n));
}

// template<typename T>
// template<int fps>
// inline typename Vector<T>::template iterator_internal<T> & Vector<T>::iterator_internal<T>::operator++() {
//   ret += offset;
//   return *this;
// }

// inline iterator_internal & operator+=(const UInt n) {
//   ret->values += offset * n;
//   return *this;
// }

/* -------------------------------------------------------------------------- */
/**
 * Specialization for scalar types
 */

/* If gcc <  4.4 some problem with detection of specialization  of a template by
 * an other template so I have copied the code, that is really ugly just to help
 * the compiler. This part of code should  be removed if version of gcc before *
 * 4.4 are note supported anymore. Problems  for know appears on Mac OS with gcc
 * version is 4.2.1
 */
#if (__GNUC__  < 4)  ||				\
  ((__GNUC__ ==  4) && (__GNUC_MINOR__ < 4))

#define specialize_internal_iterator_for_scalar(T)			\
  template <>								\
  template <>								\
  class Vector<T>::iterator_internal<T,T,0> {				\
  public:								\
  typedef T   value_type;						\
  typedef T*  pointer;							\
  typedef T&  reference;						\
  typedef T   internal_value_type;					\
  typedef T*  internal_pointer;						\
  protected:								\
  iterator_internal(UInt offset, pointer data, pointer ret) :		\
    initial(data),							\
    ret(ret) { }							\
  									\
  public:								\
  iterator_internal() : initial(NULL), ret(NULL) {};			\
  									\
  iterator_internal(pointer data, UInt offset)  : initial(data),	\
						  ret(data) {		\
    AKANTU_DEBUG_ASSERT(offset == 1,					\
			"The iterator is not compatible with the type "	\
			<< typeid(value_type).name());			\
  };									\
  									\
  iterator_internal(const iterator_internal & it) {			\
    if(this != &it) {							\
      this->initial = it.initial;					\
      this->ret = new internal_value_type(*it.ret);			\
    }									\
  }									\
  									\
  virtual ~iterator_internal() { };					\
  									\
  inline iterator_internal & operator=(const iterator_internal & it) {	\
    if(this != &it) {							\
      this->initial = it.initial;					\
      this->ret = it.ret;						\
    }									\
    return *this;							\
  }									\
  									\
  inline reference operator*() { return *ret; };			\
  inline pointer operator->() { return ret; };				\
  inline iterator_internal & operator++() { ret += offset; return *this; }; \
  									\
  inline iterator_internal & operator+=(const UInt n) {			\
    ret += n;								\
    return *this;							\
  }									\
  									\
  inline reference operator[](const UInt n) {				\
    ret = initial + n;							\
    return *ret;							\
  }									\
  									\
  inline bool operator==(const iterator_internal & other) {		\
    return (*this).ret == other.ret;					\
  }									\
  inline bool operator!=(const iterator_internal & other) {		\
    return (*this).ret != other.ret;					\
  }									\
  									\
  inline iterator_internal & operator-(size_t n) {			\
    ret -= n;                                                           \
    return *this;							\
  }									\
									\
  inline size_t operator-(const iterator_internal & b) {		\
    return ret - b.getCurrentStorage();					\
  }									\
                                                                        \
  inline pointer getCurrentStorage() const {				\
    return ret;								\
  }									\
  									\
  protected:								\
  pointer initial;							\
  pointer ret;								\
  }

specialize_internal_iterator_for_scalar(Real);
specialize_internal_iterator_for_scalar(UInt);

#else
template <typename T>
template <int fps>
class Vector<T>::iterator_internal<T,T,fps> : public std::iterator<std::random_access_iterator_tag, T>{
public:
  typedef T   value_type;
  typedef T*  pointer;
  typedef T&  reference;
  typedef T   internal_value_type;
  typedef T*  internal_pointer;
protected:
  iterator_internal(UInt offset, pointer data, pointer ret) :
    initial(data),
    ret(ret) { }

public:
  iterator_internal() : initial(NULL), ret(NULL) {};

  iterator_internal(pointer data, UInt offset)  :
    initial(data),
    ret(data) {
    AKANTU_DEBUG_ASSERT(offset == 1,
			"The iterator_internal is not compatible with the type "
			<< typeid(value_type).name());
  };

  iterator_internal(const iterator_internal & it) {
    if(this != &it) {
      this->initial = it.initial;
      this->ret = new internal_value_type(*it.ret);
    }
  }

  virtual ~iterator_internal() { };

  inline iterator_internal & operator=(const iterator_internal & it) {
    if(this != &it) {
      this->initial = it.initial;
      this->ret = it.ret;
    }
    return *this;
  }

  inline reference operator*() { return *ret; };
  inline pointer operator->() { return ret; };
  inline iterator_internal & operator++() { ++ret; return *this; };
  inline iterator_internal & operator--() { --ret; return *this; };

  inline iterator_internal & operator+=(const UInt n) {
    ret += n;
    return *this;
  }

  inline reference operator[](const UInt n) {
    ret = initial + n;
    return *ret;
  }

  inline bool operator==(const iterator_internal & other) {
    return (*this).ret == other.ret;
  }
  inline bool operator!=(const iterator_internal & other) {
    return (*this).ret != other.ret;
  }

  inline bool operator<(const iterator_internal & other) {
    return ret < other.ret;
  }

  inline bool operator<=(const iterator_internal & other) {
    return ret <= other.ret;
  }

  inline bool operator>(const iterator_internal & other) {
    return ret > other.ret;
  }

  inline bool operator>=(const iterator_internal & other) {
    return ret >= other.ret;
  }


  inline iterator_internal & operator-(size_t n) {
    ret -= n;
    return *this;
  }

  inline size_t operator-(const iterator_internal & b) {
    return ret - b.getCurrentStorage();
  }

  inline iterator_internal & operator+(size_t n) {
    ret += n;
    return *this;
  }

  inline pointer getCurrentStorage() const {
    return ret;
  }

protected:
  pointer initial;
  pointer ret;
};
#endif


/* -------------------------------------------------------------------------- */
template<typename T>
template<typename R>
inline void Vector<T>::erase(const iterator<R> & it) {
  T * curr = it.getCurrentStorage();
  UInt pos = (curr - values) / nb_component;
  erase(pos);
}

// inline reference operator[](const UInt n) {
//   ret->values = initial + n*offset;
//   return *ret;
// }


/* -------------------------------------------------------------------------- */
template<typename T>
inline Vector<T>::iterator<T> Vector<T>::begin() {
  return iterator<T>(values, 1);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline Vector<T>::iterator<T> Vector<T>::end() {
  return iterator<T>(values + nb_component * size, 1);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline Vector<T>::const_iterator<T> Vector<T>::begin() const {
  return const_iterator<T>(values, 1);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline Vector<T>::const_iterator<T> Vector<T>::end() const {
  return const_iterator<T>(values + nb_component * size, 1);
}
