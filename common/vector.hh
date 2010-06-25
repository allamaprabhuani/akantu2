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
  AKANTU_GET_MACRO(AllocatedSize, allocated_size, unsigned int);

  AKANTU_GET_MACRO(Size, size, unsigned int);

  AKANTU_GET_MACRO(NbComponent, nb_component, unsigned int);

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

  Vector(unsigned int size, unsigned int nb_component, const VectorID & id = "");

  virtual ~Vector();

private:
  Vector(const Vector<T>& vect) {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// get jth componemt of the ith tuple
  inline T & at(unsigned int i, unsigned int j = 0);

  void resize(unsigned int size);

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

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

template <class T> Vector<T>::Vector (unsigned int size,
				      unsigned int nb_component,
				      const VectorID & id) {
  this->id = id;
  values = static_cast<T*>(malloc(nb_component * size * sizeof(T)));
  AKANTU_DEBUG_ASSERT(values != NULL,
		     "Cannot allocate " << nb_component * size * sizeof(T) << " bytes");

  if (values == NULL) {
    this->size = this->allocated_size = 0;
  } else {
    AKANTU_DEBUG_INFO("Allocated " << size * nb_component * sizeof(T) << " bytes");
    this->size = this->allocated_size = size;
  }

  this->size_of_type = sizeof(T);
  this->nb_component = nb_component;
}

/* -------------------------------------------------------------------------- */
template <class T> Vector<T>::~Vector () {
  free(values);
  size = allocated_size = 0;
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
}

/* -------------------------------------------------------------------------- */
template <class T> void Vector<T>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(int i = 0; i < indent; i++, space += AKANTU_INDENT);

  int size = allocated_size * nb_component * size_of_type;

  stream << space << "Vector<" << typeid(T).name() << "> [" << std::endl;
  stream << space << " + id             : " << id << std::endl;
  stream << space << " + size           : " << size << std::endl;
  stream << space << " + nb_component   : " << nb_component << std::endl;
  stream << space << " + allocated size : " << allocated_size << std::endl;
  stream << space << " + memory size    : " << size << "B" << std::endl;
  stream << space << " + adresse        : " << std::hex << values << std::dec << std::endl;
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
