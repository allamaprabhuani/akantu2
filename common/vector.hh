/**
 * @file   vector.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jun 17 10:04:55 2010
 *
 * @brief  class of vectors
 *
 * @section LICENSE
 *
 * <insert lisence here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_VECTOR_HH__
#define __AKANTU_VECTOR_HH__

/* -------------------------------------------------------------------------- */
#include <typeinfo>

/* -------------------------------------------------------------------------- */
#include "common.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/// class that afford to store vectors in static memory
class VectorBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  //  VectorBase();

  virtual ~VectorBase() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// get the amount of space allocated in bytes
  inline unsigned int getMemorySize() const;

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {};

  /* ------------------------------------------------------------------------ */
  /* Accesors                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO_SCALAR(AllocatedSize, allocated_size, unsigned int);

  AKANTU_GET_MACRO_SCALAR(Size, size, unsigned int);

  AKANTU_GET_MACRO_SCALAR(NbComponent, nb_component, unsigned int);
  

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// the size allocated
  unsigned int allocated_size;

  /// the size used
  unsigned int size;

  /// number of components
  unsigned int nb_component;

  /// size of the stored type
  unsigned int size_of_type;

  /// id of the vector
  VectorID id;
};



/* -------------------------------------------------------------------------- */
template<class T> class Vector : public VectorBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /// Allocation of a new vector
  Vector(unsigned int size, unsigned int nb_component,
	 const VectorID & id = "");

  /// Allocation of a new vector with a default value
  Vector(unsigned int size, unsigned int nb_component,
	 const T def_value[], const VectorID & id = "");

  /// Copy constructor (deep copy if deep=true) \todo to implement 
  Vector(const Vector<T>& vect, bool deep = true);

  virtual ~Vector();



  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// get jth componemt of the ith tuple
  inline T & at(unsigned int i, unsigned int j = 0);

  /// add an element at the and of the vector
  inline void push_back(const T new_elem[]);

  /**
   * remove an element and move the last one in the hole
   * /!\ change the order in the vector
   */
  inline void erase(unsigned int i);

  /// change the size of the vector and allocate more memory if needed
  void resize(unsigned int size);

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;


  Vector<T>& operator = (const Vector<T>& vect);
  
protected:
  /// perform the allocation for the constructors
  inline void allocate(unsigned int size, unsigned int nb_component);

  /* ------------------------------------------------------------------------ */
  /* Accesors                                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
public:
  /// array of values
  T *values; // /!\ very dangerous

};



/* -------------------------------------------------------------------------- */
/* Inline Functions Vector<T>                                                 */
/* -------------------------------------------------------------------------- */

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
				      const T def_value[],
				      const VectorID & id) {
  AKANTU_DEBUG_IN();
  this->id = id;
  this->values = NULL;
  allocate(size, nb_component);

  for (unsigned int i = 0; i < size; ++i) {
    for (unsigned int j = 0; j < nb_component; ++j) {
      values[i*nb_component + j] = def_value[j];
    }
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
template <class T> Vector<T>& Vector<T>::operator= (const Vector<T>& vect) {
  AKANTU_DEBUG_IN();
  //  this->id = std::string(vect.id);
  allocate(vect.size, vect.nb_component);
  for (unsigned int i = 0; i < size; ++i) {
    for (unsigned int j = 0; j < nb_component; ++j) {
      values[i*nb_component + j] = vect.values[i*nb_component + j];
    }
  }
  
  AKANTU_DEBUG_OUT();
  return *this;
}
  
/* -------------------------------------------------------------------------- */
template <class T> inline void Vector<T>::allocate(unsigned int size,
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
template <class T> inline T & Vector<T>::at(unsigned int i, unsigned int j) {
  AKANTU_DEBUG_IN();
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
  stream << space << " + id             : " << id << std::endl;
  stream << space << " + size           : " << size << std::endl;
  stream << space << " + nb_component   : " << nb_component << std::endl;
  stream << space << " + allocated size : " << allocated_size << std::endl;
  stream << space << " + memory size    : " << real_size << "B" << std::endl;
  stream << space << " + adresse        : " << std::hex << values << std::dec << std::endl;
  if(AKANTU_DEBUG_TEST(dblDump)) {
    stream << space << " + values         : {";
    for (unsigned int i = 0; i < size; ++i) {
      stream << "{";
      for (unsigned int j = 0; j < nb_component; ++j) {
	stream << values[i*nb_component + j];
	if(j != nb_component - 1) stream << ", ";
      }
      stream << "}";
      if(i != size - 1) stream << ", ";
    }
    stream << "}" << std::endl;
  }
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
template <class T>
inline std::ostream & operator<<(std::ostream & stream, const Vector<T> & _this)
{
  _this.printself(stream);
  return stream;
}

/* -------------------------------------------------------------------------- */
/* Inline Functions VectorBase                                                */
/* -------------------------------------------------------------------------- */

inline unsigned int VectorBase::getMemorySize() const {
 return allocated_size * nb_component * size_of_type;
}

/* -------------------------------------------------------------------------- */
inline std::ostream & operator<<(std::ostream & stream, const VectorBase & _this)
{
  _this.printself(stream);
  return stream;
}



__END_AKANTU__


#endif /* __AKANTU_VECTOR_HH__ */
