/**
 * @file   aka_vector.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jun 17 10:04:55 2010
 *
 * @brief  class of vectors
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_VECTOR_HH__
#define __AKANTU_VECTOR_HH__

/* -------------------------------------------------------------------------- */
#include <typeinfo>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/// class that afford to store vectors in static memory
class VectorBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  VectorBase(const VectorID & id = "");

  virtual ~VectorBase();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// get the amount of space allocated in bytes
  inline UInt getMemorySize() const;

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(AllocatedSize, allocated_size, UInt);

  AKANTU_GET_MACRO(Size, size, UInt);

  AKANTU_GET_MACRO(NbComponent, nb_component, UInt);

  AKANTU_GET_MACRO(ID, id, const VectorID &);
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id of the vector
  VectorID id;

  /// the size allocated
  UInt allocated_size;

  /// the size used
  UInt size;

  /// number of components
  UInt nb_component;

  /// size of the stored type
  UInt size_of_type;
};



/* -------------------------------------------------------------------------- */
template<class T> class Vector : public VectorBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /// Allocation of a new vector
  Vector(UInt size = 0, UInt nb_component = 1,
	 const VectorID & id = "");

  /// Allocation of a new vector with a default value
  Vector(UInt size, UInt nb_component,
  	 const T def_values[], const VectorID & id = "");

  /// Allocation of a new vector with a default value
  Vector(UInt size, UInt nb_component,
	 const T & value, const VectorID & id = "");

  /// Copy constructor (deep copy if deep=true) \todo to implement
  Vector(const Vector<T>& vect, bool deep = true);

  virtual ~Vector();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// get jth componemt of the ith tuple in read-only
  inline const T & get(UInt i, UInt j = 0) const;
  /// get jth componemt of the ith tuple
  inline T & at(UInt i, UInt j = 0);

  /// add an  element at  the and  of the vector  with the  value value  for all
  /// component
  inline void push_back(const T& value);

  /// add an element at the and of the vector
  inline void push_back(const T new_elem[]);

  /**
   * remove an element and move the last one in the hole
   * /!\ change the order in the vector
   */
  inline void erase(UInt i);

  /// change the size of the vector and allocate more memory if needed
  void resize(UInt size);

  /// change the number of components by interlacing data
  void extendComponentsInterlaced(UInt multiplicator,UInt stride);

  
  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  //  Vector<T>& operator=(const Vector<T>& vect);

  /// search elem in the vector, return  the position of the first occurrence or
  /// -1 if not found
  Int find(const T & elem) const;


protected:
  /// perform the allocation for the constructors
  void allocate(UInt size, UInt nb_component = 1);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  UInt getSize() const{ return this->size; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
public:
  /// array of values
  T *values; // /!\ very dangerous

};

#include "aka_vector_inline_impl.cc"

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
