/**
 * @file   vector.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jun 17 10:04:55 2010
 *
 * @brief  class of vectors
 *
 * @section LICENSE
 *
 * <insert license here>
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
  VectorBase();

  virtual ~VectorBase() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// get the amount of space allocated in bytes
  inline unsigned int getMemorySize() const;

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                 */
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
	 const T def_values[], const VectorID & id = "");

  /// Allocation of a new vector with a default value
  Vector(unsigned int size, unsigned int nb_component,
	 const T & value, const VectorID & id = "");

  /// Copy constructor (deep copy if deep=true) \todo to implement
  Vector(const Vector<T>& vect, bool deep = true);

  virtual ~Vector();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// get jth componemt of the ith tuple
  inline T & at(unsigned int i, unsigned int j = 0);

  /// add an  element at  the and  of the vector  with the  value value  for all
  /// component
  inline void push_back(const T& value);

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

  //  Vector<T>& operator=(const Vector<T>& vect);

protected:
  /// perform the allocation for the constructors
  void allocate(unsigned int size, unsigned int nb_component);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  unsigned int getSize() const{ return this->size; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
public:
  /// array of values
  T *values; // /!\ very dangerous

};

#include "vector_inline_impl.cc"

/* -------------------------------------------------------------------------- */
/* Inline Functions Vector<T>                                                 */
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
inline std::ostream & operator<<(std::ostream & stream, const VectorBase & _this)
{
  _this.printself(stream);
  return stream;
}



__END_AKANTU__


#endif /* __AKANTU_VECTOR_HH__ */
