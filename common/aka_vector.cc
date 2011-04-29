/**
 * @file   aka_vector.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jun 17 15:14:24 2010
 *
 * @brief  class vector
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
#include "aka_common.hh"
#include "aka_vector.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* Functions VectorBase                                                       */
/* -------------------------------------------------------------------------- */
VectorBase::VectorBase(const VectorID & id) :
  id(id), allocated_size(0), size(0), nb_component(1), size_of_type(0) {
}

/* -------------------------------------------------------------------------- */
VectorBase::~VectorBase() {
}

/* -------------------------------------------------------------------------- */
void VectorBase::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  stream << space << "VectorBase [" << std::endl;
  stream << space << " + size             : " << size << std::endl;
  stream << space << " + nb component     : " << nb_component << std::endl;
  stream << space << " + allocated size   : " << allocated_size << std::endl;
  Real mem_size = (allocated_size * nb_component * size_of_type) / 1024.;
  stream << space << " + size of type     : " << size_of_type << "B" << std::endl;
  stream << space << " + memory allocated : " << mem_size << "kB" << std::endl;
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
/* Functions Vector<T>                                                        */
/* -------------------------------------------------------------------------- */
template <class T> Vector<T>::Vector (UInt size,
				      UInt nb_component,
				      const VectorID & id) :
  VectorBase(id), values(NULL) {
  AKANTU_DEBUG_IN();
  allocate(size, nb_component);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T> Vector<T>::Vector (UInt size,
				      UInt nb_component,
				      const T def_values[],
				      const VectorID & id) :
  VectorBase(id), values(NULL) {
  AKANTU_DEBUG_IN();
  allocate(size, nb_component);

  for (UInt i = 0; i < size; ++i) {
    for (UInt j = 0; j < nb_component; ++j) {
      values[i*nb_component + j] = def_values[j];
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T> Vector<T>::Vector (UInt size,
				      UInt nb_component,
				      const T & value,
				      const VectorID & id) :
  VectorBase(id), values(NULL) {
  AKANTU_DEBUG_IN();
  allocate(size, nb_component);

  for (UInt i = 0; i < nb_component*size; ++i) {
    values[i] = value;
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template <class T> Vector<T>::Vector(const Vector<T>& vect, bool deep) {
  AKANTU_DEBUG_IN();
  this->id = vect.id;

  if (deep) {
    allocate(vect.size, vect.nb_component);
    for (UInt i = 0; i < size; ++i) {
      for (UInt j = 0; j < nb_component; ++j) {
	values[i*nb_component + j] = vect.values[i*nb_component + j];
      }
    }
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
template <class T> void Vector<T>::resize(UInt new_size) {
  //  AKANTU_DEBUG_IN();
  /// free some memory
  if(new_size <= allocated_size) {
    if(allocated_size - new_size > AKANTU_MIN_ALLOCATION) {
      AKANTU_DEBUG(dblAccessory, "Freeing "
		   << (allocated_size - size)*nb_component*sizeof(T) / 1024.
		   << "kB (" << id <<")");

      /// Normally there are no allocation problem when reducing an array
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

  /// allocate more memory
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
							      UInt stride) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(multiplicator > 1,
		      "you can only extend a vector by changing the number of components");
  AKANTU_DEBUG_ASSERT(nb_component%stride == 0,
		      "stride must divide actual number of components");

  values = static_cast<T*>(realloc(values, nb_component*multiplicator*size* sizeof(T)));

  UInt strided_component = nb_component/stride;
  UInt new_component = strided_component * multiplicator;

  for (UInt i = 0,k=size-1; i < size; ++i,--k) {
    for (UInt j = 0; j < new_component; ++j) {
      UInt m = new_component - j -1;
      UInt n = m*strided_component/new_component;
      for (UInt l = 0,p=stride-1;  l < stride; ++l,--p) {
	values[k*new_component*stride+m*stride+p] = values[k*strided_component*stride+n*stride+p];
      }
    }
  }

  nb_component = new_component*stride;

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
template <> Int Vector<Real>::find(const Real & elem) const {
  AKANTU_DEBUG_IN();
  UInt i = 0;
  Real epsilon = std::numeric_limits<Real>::epsilon();
  for (; (i < size) && (fabs(values[i] - elem) <= epsilon); ++i);

  AKANTU_DEBUG_OUT();
  return (i == size) ? -1 : (Int) i;
}

/* -------------------------------------------------------------------------- */
template <class T> void Vector<T>::copy(const Vector<T>& vect) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(vect.nb_component == nb_component,
		      "The two vectors does not have the same number of components");

  //  this->id = vect.id;
  resize(vect.size);

  memcpy(this->values, vect.values, vect.size * nb_component * sizeof(T));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T> void Vector<T>::printself(std::ostream & stream, int indent) const {
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
  if(!AKANTU_DEBUG_LEVEL_IS_TEST())
  stream << space << " + memory size    : "
	 << real_size << "kB" << std::endl;
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
template class Vector<char>;
template class Vector<Int>;
template class Vector<UInt>;
template class Vector<UInt64>;
template class Vector<Real>;
template class Vector<bool>;

__END_AKANTU__
