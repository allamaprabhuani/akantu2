/**
 * @file   vector.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jun 17 15:14:24 2010
 *
 * @brief  class vector
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "common.hh"
#include "vector.hh"

__BEGIN_AKANTU__


/* -------------------------------------------------------------------------- */
/* Functions VectorBase                                                       */
/* -------------------------------------------------------------------------- */
VectorBase::VectorBase() {
  allocated_size = 0;
  size = 0;
  nb_component = 1;
  size_of_type=0;
  id = "";
};

/* -------------------------------------------------------------------------- */
void VectorBase::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(int i = 0; i < indent; i++, space += AKANTU_INDENT);
  stream << indent << "VectorBase [" << std::endl;
  stream << indent << " + size             : " << size << std::endl;
  stream << indent << " + nb component     : " << nb_component << std::endl;
  stream << indent << " + allocated size   : " << allocated_size << std::endl;
  double mem_size = (allocated_size * nb_component * size_of_type) / 1024.;
  stream << indent << " + size of type     : " << size_of_type << "B" << std::endl;
  stream << indent << " + memory allocated : " << mem_size << "kB" << std::endl;
  stream << indent << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
/* Functions Vector<T>                                                        */
/* -------------------------------------------------------------------------- */
template <class T> Vector<T>::Vector (unsigned int size,
				      unsigned int nb_component,
				      const VectorID & id) {
  AKANTU_DEBUG_IN();
  this->id = id;
  this->values = NULL;
  allocate(size, nb_component);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T> Vector<T>::Vector (unsigned int size,
				      unsigned int nb_component,
				      const T def_values[],
				      const VectorID & id) {
  AKANTU_DEBUG_IN();
  this->id = id;
  this->values = NULL;
  allocate(size, nb_component);

  for (unsigned int i = 0; i < size; ++i) {
    for (unsigned int j = 0; j < nb_component; ++j) {
      values[i*nb_component + j] = def_values[j];
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T> Vector<T>::Vector (unsigned int size,
				      unsigned int nb_component,
				      const T& value,
				      const VectorID & id) {
  AKANTU_DEBUG_IN();
  this->id = id;
  this->values = NULL;
  allocate(size, nb_component);

  for (unsigned int i = 0; i < nb_component*size; ++i) {
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
    for (unsigned int i = 0; i < size; ++i) {
      for (unsigned int j = 0; j < nb_component; ++j) {
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
template <class T> void Vector<T>::allocate(unsigned int size,
					    unsigned int nb_component) {
  AKANTU_DEBUG_IN();
  if (size == 0){
    values = NULL;
  } else {
    values = static_cast<T*>(malloc(nb_component * size * sizeof(T)));
    AKANTU_DEBUG_ASSERT(values != NULL,
			"Cannot allocate " << nb_component * size * sizeof(T) << " bytes");
  }

  if (values == NULL) {
    this->size = this->allocated_size = 0;
  } else {
    AKANTU_DEBUG_INFO("Allocated " << size * nb_component * sizeof(T) << " bytes");
    this->size = this->allocated_size = size;
  }

  this->size_of_type = sizeof(T);
  this->nb_component = nb_component;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T> Vector<T>::~Vector () {
  free(values);
  size = allocated_size = 0;
}

/* -------------------------------------------------------------------------- */

template class Vector<int>;
template class Vector<long int>;
template class Vector<float>;
template class Vector<double>;

__END_AKANTU__
