/**
 * @file   aka_vector_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Thu Jul 15 00:41:12 2010
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
template <class T, bool is_scal>
inline T & Vector<T, is_scal>::operator()(UInt i, UInt j) {
  AKANTU_DEBUG_ASSERT(size > 0,
		      "The vector \"" << id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < size) && (j < nb_component),
		      "The value at position [" << i << "," << j
		      << "] is out of range in vector \"" << id << "\"");
  return values[i*nb_component + j];
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline const T & Vector<T, is_scal>::operator()(UInt i, UInt j) const {
  AKANTU_DEBUG_ASSERT(size > 0,
		      "The vector \"" << id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < size) && (j < nb_component),
		      "The value at position [" << i << "," << j
		      << "] is out of range in vector \"" << id << "\"");
  return values[i*nb_component + j];
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline void Vector<T, is_scal>::push_back(const T & value) {
  UInt pos = size;

  resizeUnitialized(size+1);

  std::uninitialized_fill_n(values + pos * nb_component, nb_component, value);
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline void Vector<T, is_scal>::push_back(const T new_elem[]) {
  UInt pos = size;

  resizeUnitialized(size+1);

  T * tmp = values + nb_component * pos;
  std::uninitialized_copy(new_elem, new_elem + nb_component, tmp);
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
template<class Ret>
inline void Vector<T, is_scal>::push_back(const Vector<T, is_scal>::iterator<Ret> & it) {
  UInt pos = size;

  resizeUnitialized(size+1);

  T * tmp = values + nb_component * pos;
  T * new_elem = it.data();
  std::uninitialized_copy(new_elem, new_elem + nb_component, tmp);
}

/* -------------------------------------------------------------------------- */




/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline void Vector<T, is_scal>::erase(UInt i){
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT((size > 0),
		      "The vector is empty");
  AKANTU_DEBUG_ASSERT((i < size),
		      "The element at position [" << i << "] is out of range (" << i << ">=" << size << ")");


  if(i != (size - 1)) {
    for (UInt j = 0; j < nb_component; ++j) {
      values[i*nb_component + j] = values[(size-1)*nb_component + j];
    }
  }

  resize(size - 1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Vector<T, is_scal> & Vector<T, is_scal>::operator-=(const Vector<T, is_scal> & vect) {
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
template <class T, bool is_scal>
Vector<T, is_scal> & Vector<T, is_scal>::operator+=(const Vector<T, is_scal> & vect) {
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
template <class T, bool is_scal>
Vector<T, is_scal> & Vector<T, is_scal>::operator*=(const T & alpha) {
  T * a = values;
  for (UInt i = 0; i < size*nb_component; ++i) {
    *a++ *= alpha;
  }

  return *this;
}

/* -------------------------------------------------------------------------- */
/* Functions Vector<T, is_scal>                                               */
/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Vector<T, is_scal>::Vector (UInt size,
			    UInt nb_component,
			    const ID & id) :
  VectorBase(id), values(NULL) {
  AKANTU_DEBUG_IN();
  allocate(size, nb_component);

  if(!is_scal) {
    T val = T();
    std::uninitialized_fill(values, values + size*nb_component, val);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Vector<T, is_scal>::Vector (UInt size,
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
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Vector<T, is_scal>::Vector (UInt size,
			    UInt nb_component,
			    const T & value,
			    const ID & id) :
  VectorBase(id), values(NULL) {
  AKANTU_DEBUG_IN();
  allocate(size, nb_component);

  std::uninitialized_fill_n(values, size*nb_component, value);

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Vector<T, is_scal>::Vector(const Vector<T, is_scal> & vect,
			   bool deep,
			   const ID & id) {
  AKANTU_DEBUG_IN();
  this->id = (id == "") ? vect.id : id;

  if (deep) {
    allocate(vect.size, vect.nb_component);
    T * tmp = values;
    std::uninitialized_copy(vect.values, vect.values + size * nb_component, tmp);
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
template <class T, bool is_scal>
Vector<T, is_scal>::Vector(const std::vector<T>& vect) {
  AKANTU_DEBUG_IN();
  this->id = "";

  allocate(vect.size(), 1);
  T * tmp = values;
  std::uninitialized_copy(&(vect[0]), &(vect[size-1]), tmp);

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Vector<T, is_scal>::~Vector () {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG(dblAccessory, "Freeing "
	       << allocated_size*nb_component*sizeof(T) / 1024.
	       << "kB (" << id <<")");

  if(values){
    if(!is_scal)
      for (UInt i = 0; i < size * nb_component; ++i) {
	T * obj = values+i;
	obj->~T();
      }
    free(values);
  }
  size = allocated_size = 0;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
void Vector<T, is_scal>::allocate(UInt size,
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
template <class T, bool is_scal>
void Vector<T, is_scal>::resize(UInt new_size) {
  UInt old_size = size;

  T * old_values = values;
  if(new_size < size) {
    for (UInt i = new_size * nb_component; i < size * nb_component; ++i) {
      T * obj = old_values+i;
      obj->~T();
    }
  }

  resizeUnitialized(new_size);


  T val = T();
  if(size > old_size)
    std::uninitialized_fill(values + old_size*nb_component, values + size*nb_component, val);
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
void Vector<T, is_scal>::resizeUnitialized(UInt new_size) {
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
template <class T, bool is_scal>
void Vector<T, is_scal>::extendComponentsInterlaced(UInt multiplicator,
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
template <class T, bool is_scal>
Int Vector<T, is_scal>::find(const T & elem) const {
  AKANTU_DEBUG_IN();
  UInt i = 0;
  for (; (i < size) && (values[i] != elem); ++i);

  AKANTU_DEBUG_OUT();
  return (i == size) ? -1 : (Int) i;
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
Int Vector<T, is_scal>::find(T elem[]) const {
  AKANTU_DEBUG_IN();
  T * it = values;
  UInt i = 0;
  for (;i < size; ++i) {
    if(*it == elem[0]) {
      T * cit = it;
      UInt c = 0;
      for(; (c < nb_component) && (*cit == elem[c]); ++c, ++cit);
      if(c == nb_component) {
	AKANTU_DEBUG_OUT();
	return i;
      }
    }
    it += nb_component;
  }
  return -1;
}


/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
void Vector<T, is_scal>::copy(const Vector<T, is_scal>& vect) {
  AKANTU_DEBUG_IN();

  if(AKANTU_DEBUG_TEST(dblWarning))
    if(vect.nb_component != nb_component) {
      AKANTU_DEBUG(dblWarning, "The two vectors do not have the same number of components");
    }
  //  this->id = vect.id;
  resize((vect.size * vect.nb_component) / nb_component);

  T * tmp = values;
  std::uninitialized_copy(vect.values, vect.values + size * nb_component, tmp);
  //  memcpy(this->values, vect.values, vect.size * vect.nb_component * sizeof(T));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
void Vector<T, is_scal>::printself(std::ostream & stream, int indent) const {
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

/* -------------------------------------------------------------------------- */
/* Iterators                                                                  */
/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
template<class R, class IR, bool is_r_scal>
class Vector<T, is_scal>::iterator_internal {
public:
  typedef R                               value_type;
  typedef R*                              pointer;
  typedef R&                              reference;
  typedef IR                              internal_value_type;
  typedef IR*                             internal_pointer;
  typedef std::ptrdiff_t                  difference_type;
  typedef std::random_access_iterator_tag iterator_category;

public:
  iterator_internal() : _offset(0), initial(NULL), ret(NULL) {};

  iterator_internal(pointer_type data, UInt _offset)  :
    _offset(_offset),
    initial(data),
    ret(new value_type(data)) {
    AKANTU_DEBUG_ASSERT(_offset == ret->size(),
			"The iterator_internal is not compatible with the type "
			<< typeid(value_type).name());
  }

  iterator_internal(pointer warped)  : _offset(warped->size()),
				       initial(warped->storage()),
				       ret(const_cast<internal_pointer>(warped)) {
  }

  iterator_internal(const iterator_internal & it) {
    if(this != &it) {
      this->_offset = it._offset;
      this->initial = it.initial;
      this->ret = new internal_value_type(*it.ret);
    }
  }

  virtual ~iterator_internal() { delete ret; };

  inline iterator_internal & operator=(const iterator_internal & it) {
    if(this != &it) {
      this->_offset = it._offset;
      this->initial = it.initial;
      if(this->ret) this->ret->shallowCopy(*it.ret);
      else this->ret = new internal_value_type(*it.ret);
    }
    return *this;
  }

  inline reference operator*() { return *ret; };
  inline pointer operator->() { return ret; };
  inline iterator_internal & operator++() { ret->values += _offset; return *this; };
  inline iterator_internal & operator--() { ret->values -= _offset; return *this; };

  inline iterator_internal & operator+=(const UInt n) { ret->values += _offset * n; return *this; }
  inline iterator_internal & operator-=(const UInt n) { ret->values -= _offset * n; return *this; }

  inline reference operator[](const UInt n) { ret->values = initial + n*_offset; return *ret; }

  inline bool operator==(const iterator_internal & other) const { return (*this).ret->storage() == other.ret->storage(); }
  inline bool operator!=(const iterator_internal & other) const { return (*this).ret->storage() != other.ret->storage(); }
  inline bool operator <(const iterator_internal & other) const { return (*this).ret->storage()  < other.ret->storage(); }
  inline bool operator<=(const iterator_internal & other) const { return (*this).ret->storage() <= other.ret->storage(); }
  inline bool operator> (const iterator_internal & other) const { return (*this).ret->storage() >  other.ret->storage(); }
  inline bool operator>=(const iterator_internal & other) const { return (*this).ret->storage() >= other.ret->storage(); }

  inline iterator_internal operator+(difference_type n) { iterator_internal tmp(*this); tmp += n; return tmp; }
  inline iterator_internal operator-(difference_type n) { iterator_internal tmp(*this); tmp -= n; return tmp; }


  inline difference_type operator-(const iterator_internal & b) { return ret->values - b.ret->values; }


  inline pointer_type data() const { return ret->storage(); }
  inline difference_type offset() const { return _offset; }

protected:
  UInt _offset;
  pointer_type initial;
  internal_pointer ret;
};

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::iterator< types::Vector<T> > Vector<T, is_scal>::begin(UInt n) {
  AKANTU_DEBUG_ASSERT(nb_component == n,
		      "The iterator is not compatible with the type Vector("
		      << n<< ")");
  return iterator< types::Vector<T> >(new types::Vector<T>(values, n));
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::iterator< types::Vector<T> > Vector<T, is_scal>::end(UInt n) {
  AKANTU_DEBUG_ASSERT(nb_component == n,
		      "The iterator is not compatible with the type Vector("
		      << n<< ")");
  return iterator< types::Vector<T> >(new types::Vector<T>(values + nb_component * size,
							   n));
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::const_iterator< types::Vector<T> > Vector<T, is_scal>::begin(UInt n) const {
  AKANTU_DEBUG_ASSERT(nb_component == n,
		      "The iterator is not compatible with the type Vector("
		      << n<< ")");
  return const_iterator< types::Vector<T> >(new types::Vector<T>(values, n));
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::const_iterator< types::Vector<T> > Vector<T, is_scal>::end(UInt n) const {
  AKANTU_DEBUG_ASSERT(nb_component == n,
		      "The iterator is not compatible with the type Vector("
		      << n<< ")");
  return const_iterator< types::Vector<T> >(new types::Vector<T>(values + nb_component * size,
							   n));
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::iterator< types::Vector<T> >
Vector<T, is_scal>::begin_reinterpret(UInt n, __attribute__((unused)) UInt size) {
  AKANTU_DEBUG_ASSERT(n * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << n <<
		      ") are not compatible with the one of this vector");
  return iterator< types::Vector<T> >(new types::Vector<T>(values, n));
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::iterator< types::Vector<T> >
Vector<T, is_scal>::end_reinterpret(UInt n, UInt size) {
  AKANTU_DEBUG_ASSERT(n * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << n <<
		      ") are not compatible with the one of this vector");
  return iterator< types::Vector<T> >(new types::Vector<T>(values + n * size, n));
}


/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::const_iterator< types::Vector<T> >
Vector<T, is_scal>::begin_reinterpret(UInt n, __attribute__((unused)) UInt size) const {
  AKANTU_DEBUG_ASSERT(n * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << n <<
		      ") are not compatible with the one of this vector");
  return const_iterator< types::Vector<T> >(new types::Vector<T>(values, n));
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::const_iterator< types::Vector<T> >
Vector<T, is_scal>::end_reinterpret(UInt n, UInt size) const {
  AKANTU_DEBUG_ASSERT(n * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << n <<
		      ") are not compatible with the one of this vector");
  return const_iterator< types::Vector<T> >(new types::Vector<T>(values + n * size, n));
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::iterator< types::Matrix<T> > Vector<T, is_scal>::begin(UInt m, UInt n) {
  AKANTU_DEBUG_ASSERT(nb_component == n*m,
		      "The iterator is not compatible with the type Matrix("
		      << m << "," << n<< ")");
  return iterator< types::Matrix<T> >(new types::Matrix<T>(values, m, n));
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::iterator< types::Matrix<T> > Vector<T, is_scal>::end(UInt m, UInt n) {
  AKANTU_DEBUG_ASSERT(nb_component == n*m,
		      "The iterator is not compatible with the type Matrix("
		      << m << "," << n<< ")");
  return iterator< types::Matrix<T> >(new types::Matrix<T>(values + nb_component * size, m, n));
}


/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::const_iterator< types::Matrix<T> > Vector<T, is_scal>::begin(UInt m, UInt n) const {
  AKANTU_DEBUG_ASSERT(nb_component == n*m,
		      "The iterator is not compatible with the type Matrix("
		      << m << "," << n<< ")");
  return const_iterator< types::Matrix<T> >(new types::Matrix<T>(values, m, n));
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::const_iterator< types::Matrix<T> > Vector<T, is_scal>::end(UInt m, UInt n) const {
  AKANTU_DEBUG_ASSERT(nb_component == n*m,
		      "The iterator is not compatible with the type Matrix("
		      << m << "," << n<< ")");
  return const_iterator< types::Matrix<T> >(new types::Matrix<T>(values + nb_component * size, m, n));
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::iterator< types::Matrix<T> >
Vector<T, is_scal>::begin_reinterpret(UInt m, UInt n, __attribute__((unused)) UInt size) {
  AKANTU_DEBUG_ASSERT(n * m * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << n * m <<
		      ") are not compatible with the one of this vector");
  return iterator< types::Matrix<T> >(new types::Matrix<T>(values, m, n));
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::iterator< types::Matrix<T> >
Vector<T, is_scal>::end_reinterpret(UInt m, UInt n, UInt size) {
  AKANTU_DEBUG_ASSERT(n * m * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << n * m <<
		      ") are not compatible with the one of this vector");
  return iterator< types::Matrix<T> >(new types::Matrix<T>(values + n * m * size, m, n));
}


/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::const_iterator< types::Matrix<T> >
Vector<T, is_scal>::begin_reinterpret(UInt m, UInt n, __attribute__((unused)) UInt size) const {
  AKANTU_DEBUG_ASSERT(n * m * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << n * m <<
		      ") are not compatible with the one of this vector");
  return const_iterator< types::Matrix<T> >(new types::Matrix<T>(values, m, n));
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::const_iterator< types::Matrix<T> >
Vector<T, is_scal>::end_reinterpret(UInt m, UInt n, UInt size) const {
  AKANTU_DEBUG_ASSERT(n * m * size == this->nb_component * this->size,
		      "The new values for size (" << size
		      << ") and nb_component (" << n * m <<
		      ") are not compatible with the one of this vector");
  return const_iterator< types::Matrix<T> >(new types::Matrix<T>(values + n * m * size, m, n));
}


/* -------------------------------------------------------------------------- */
/**
 * Specialization for scalar types
 */
template <class T, bool is_scal>
template <class R, class IR>
class Vector<T, is_scal>::iterator_internal<R, IR, true> {
public:
  typedef R                               value_type;
  typedef R*                              pointer;
  typedef R&                              reference;
  typedef IR                              internal_value_type;
  typedef IR*                             internal_pointer;
  typedef std::ptrdiff_t                  difference_type;
  typedef std::random_access_iterator_tag iterator_category;

public:
  iterator_internal(pointer data = NULL, __attribute__ ((unused)) UInt _offset = 1) : _offset(_offset), ret(data), initial(data) { };
  iterator_internal(const iterator_internal & it) {
    if(this != &it) { this->ret = it.ret; this->initial = it.initial; }
  }

  virtual ~iterator_internal() { };

  inline iterator_internal & operator=(const iterator_internal & it)
  { if(this != &it) { this->ret = it.ret; this->initial = it.initial; } return *this; }

  inline reference operator*() { return *ret; };
  inline pointer operator->() { return ret; };
  inline iterator_internal & operator++() { ++ret; return *this; };
  inline iterator_internal & operator--() { --ret; return *this; };

  inline iterator_internal & operator+=(const UInt n) { ret += n; return *this; }
  inline iterator_internal & operator-=(const UInt n) { ret -= n; return *this; }

  inline reference operator[](const UInt n) { ret = initial + n; return *ret; }

  inline bool operator==(const iterator_internal & other) const { return ret == other.ret; }
  inline bool operator!=(const iterator_internal & other) const { return ret != other.ret; }
  inline bool operator< (const iterator_internal & other) const { return ret <  other.ret; }
  inline bool operator<=(const iterator_internal & other) const { return ret <= other.ret; }
  inline bool operator> (const iterator_internal & other) const { return ret >  other.ret; }
  inline bool operator>=(const iterator_internal & other) const { return ret >= other.ret; }

  inline iterator_internal operator-(difference_type n) { return iterator_internal(ret - n); }
  inline iterator_internal operator+(difference_type n) { return iterator_internal(ret + n); }

  inline difference_type operator-(const iterator_internal & b) { return ret - b.ret; }

  inline pointer data() const { return ret; }
  inline difference_type offset() const { return _offset; }
protected:
  difference_type _offset;
  pointer ret;
  pointer initial;
};


/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
template<typename R>
inline void Vector<T, is_scal>::erase(const iterator<R> & it) {
  T * curr = it.getCurrentStorage();
  UInt pos = (curr - values) / nb_component;
  erase(pos);
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::iterator<T> Vector<T, is_scal>::begin() {
  AKANTU_DEBUG_ASSERT(nb_component == 1, "this iterator cannot be used on a vector which has nb_component != 1");
  return iterator<T>(values);
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::iterator<T> Vector<T, is_scal>::end() {
  AKANTU_DEBUG_ASSERT(nb_component == 1, "this iterator cannot be used on a vector which has nb_component != 1");
  return iterator<T>(values + size);
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::const_iterator<T> Vector<T, is_scal>::begin() const {
  AKANTU_DEBUG_ASSERT(nb_component == 1, "this iterator cannot be used on a vector which has nb_component != 1");
  return const_iterator<T>(values);
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline Vector<T, is_scal>::const_iterator<T> Vector<T, is_scal>::end() const {
  AKANTU_DEBUG_ASSERT(nb_component == 1, "this iterator cannot be used on a vector which has nb_component != 1");
  return const_iterator<T>(values + size);
}
