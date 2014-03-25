/**
 * @file   aka_array.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Jun 18 11:48:28 2010
 *
 * @brief  class of vectors
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

#ifndef __AKANTU_VECTOR_HH__
#define __AKANTU_VECTOR_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <typeinfo>
//#include <cstring>

#include <vector>

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__


/// class that afford to store vectors in static memory
class ArrayBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ArrayBase(const ID & id = "");

  virtual ~ArrayBase();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// get the amount of space allocated in bytes
  inline UInt getMemorySize() const;

  /// set the size to zero without freeing the allocated space
  inline void empty();

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(AllocatedSize, allocated_size, UInt);

  AKANTU_GET_MACRO(Size, size, UInt);

  AKANTU_GET_MACRO(NbComponent, nb_component, UInt);

  AKANTU_GET_MACRO(ID, id, const ID &);

  AKANTU_GET_MACRO(Tag, tag, const std::string &);
  AKANTU_SET_MACRO(Tag, tag, const std::string &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id of the vector
  ID id;

  /// the size allocated
  UInt allocated_size;

  /// the size used
  UInt size;

  /// number of components
  UInt nb_component;

  /// size of the stored type
  UInt size_of_type;

  /// User defined tag
  std::string tag;
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */
template<typename T, bool is_scal>
class Array : public ArrayBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  typedef T                  value_type;
  typedef value_type &       reference;
  typedef value_type *       pointer_type;
  typedef const value_type & const_reference;

  /// Allocation of a new vector
  inline Array(UInt size = 0, UInt nb_component = 1,
		const ID & id = "");

  /// Allocation of a new vector with a default value
  Array(UInt size, UInt nb_component,
  	 const value_type def_values[], const ID & id = "");

  /// Allocation of a new vector with a default value
  Array(UInt size, UInt nb_component,
	 const_reference value, const ID & id = "");

  /// Copy constructor (deep copy if deep=true)
  Array(const Array<value_type, is_scal>& vect, bool deep = true, const ID & id = "");

  /// Copy constructor (deep copy)
  Array(const std::vector<value_type> & vect);


  virtual inline ~Array();

  /* ------------------------------------------------------------------------ */
  /* Iterator                                                                 */
  /* ------------------------------------------------------------------------ */
  /// \todo protected: does not compile with intel  check why
public:
  template <class R, class IR = R, bool issame = is_same<IR, T>::value >
  class iterator_internal;
public:
  /* ------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------ */
// #if defined(AKANTU_CORE_CXX11)
//   template<class R> using const_iterator = iterator_internal<const R, R>;
// #else
  template<typename R = T>  class const_iterator;
// #endif

// #if defined(AKANTU_CORE_CXX11)
//   template<class R> using iterator = iterator_internal<R>;
// #else
  template<typename R = T>  class iterator;
// #endif
  typedef iterator< T > scalar_iterator;
  typedef const_iterator< T > const_scalar_iterator;

  typedef iterator< Vector<T> > vector_iterator;
  typedef const_iterator< Vector<T> > const_vector_iterator;

  typedef iterator< Matrix<T> > matrix_iterator;
  typedef const_iterator< Matrix<T> > const_matrix_iterator;

  /*! Get an iterator that behaves like a pointer T * to the
   *  first entry in the member array values
   *  @return a scalar_iterator
   */
  inline iterator<T> begin();
  /*! Get an iterator that behaves like a pointer T * that points *past* the
   *  last entry in the member array values
   *  @return a scalar_iterator
   */
  inline iterator<T> end();
  /*! Get a const iterator that behaves like a pointer T * to the
   *  first entry in the member array values
   *  @return a const_scalar_iterator
   */
  inline const_iterator<T> begin() const;
  /*! Get a const iterator that behaves like a pointer T * that points *past* the
   *  last entry in the member array values
   *  @return a const_scalar_iterator
   */
  inline const_iterator<T> end() const;

  /*! Get an iterator that behaves like a pointer akantu::Vector<T> * to the
   *  first tuple of the array.
   *  @param n vector size. Has to be equal to nb_component. This unfortunate
   *  redundancy is necessary to distinguish it from ::begin() which it
   *  overloads. If compiled in debug mode, an incorrect value of n will result
   *  in an exception being thrown. Optimized code will fail in an unpredicted
   *  manner.
   *  @return a vector_iterator
   */
  inline vector_iterator begin(UInt n);
  /*! Get an iterator that behaves like a pointer akantu::Vector<T> * pointing
   *  *past* the last tuple of the array.
   *  @param n vector size. see Array::begin(UInt n) for more
   *  @return a vector_iterator
   */
  inline vector_iterator end(UInt n);
  /*! Get a const iterator that behaves like a pointer akantu::Vector<T> * to the
   *  first tuple of the array.
   *  @param n vector size. see Array::begin(UInt n) for more
   *  @return a vector_iterator
   */
  inline const_vector_iterator begin(UInt n) const;
  /*! Get a const iterator that behaves like a pointer akantu::Vector<T> * pointing
   *  *past* the last tuple of the array.
   *  @param n vector size. see Array::begin(UInt n) for more
   *  @return a const_vector_iterator
   */
  inline const_vector_iterator end(UInt n) const;

  /*! Get an iterator that behaves like a pointer akantu::Matrix<T> * to the
   *  first tuple of the array.
   *  @param m number of rows
   *  @param n number of columns. m times n has to equal nb_component.
   *  If compiled in debug mode, an incorrect combination of m and n will result
   *  in an exception being thrown. Optimized code will fail in an unpredicted
   *  manner.
   *  @return a matrix_iterator
   */
  inline matrix_iterator begin(UInt m, UInt n);
  /*! Get an iterator that behaves like a pointer akantu::Matrix<T> * pointing
   *  *past* the last tuple of the array.
   *  @param m number of rows
   *  @param n number of columns. See Array::begin(UInt m, UInt n)
   *  @return a matrix_iterator
   */
  inline matrix_iterator end(UInt m, UInt n);
  /*! Get a const iterator that behaves like a pointer akantu::Matrix<T> * to the
   *  first tuple of the array.
   *  @param m number of rows
   *  @param n number of columns. See Array::begin(UInt m, UInt n)
   *  @return a matrix_iterator
   */
  inline const_matrix_iterator begin(UInt m, UInt n) const;
  /*! Get a const iterator that behaves like a pointer akantu::Matrix<T> * pointing
   *  *past* the last tuple of the array.
   *  @param m number of rows
   *  @param n number of columns. See Array::begin(UInt m, UInt n)
   *  @return a const_matrix_iterator
   */
  inline const_matrix_iterator end(UInt m, UInt n) const;

  /*! Get an iterator that behaves like a pointer akantu::Vector<T> * to the
   *  first tuple of the array.
   *
   *  The reinterpret iterators allow to iterate over an array in any way that
   *  preserves the number of entries of the array. This can for instance be use
   *  full if the shape of the data in an array is not initially known.
   *  @param n vector size.
   *  @param size number of tuples in array. n times size must match the number
   *  of entries of the array. If compiled in debug mode, an incorrect
   *  combination of n and size will result
   *  in an exception being thrown. Optimized code will fail in an unpredicted
   *  manner.
   *  @return a vector_iterator
   */
  inline vector_iterator begin_reinterpret(UInt n, UInt size);
  /*! Get an iterator that behaves like a pointer akantu::Vector<T> * pointing
   *  *past* the last tuple of the array.
   *  @param n vector size.
   *  @param size number of tuples in array. See Array::begin_reinterpret(UInt n, UInt size)
   *  @return a vector_iterator
   */
  inline vector_iterator end_reinterpret(UInt n, UInt size);
  /*! Get a const iterator that behaves like a pointer akantu::Vector<T> * to the
   *  first tuple of the array.
   *  @param n vector size.
   *  @param size number of tuples in array. See Array::begin_reinterpret(UInt n, UInt size)
   *  @return a const_vector_iterator
   */
  inline const_vector_iterator begin_reinterpret(UInt n, UInt size) const;
  /*! Get a const iterator that behaves like a pointer akantu::Vector<T> * pointing
   *  *past* the last tuple of the array.
   *  @param n vector size.
   *  @param size number of tuples in array. See Array::begin_reinterpret(UInt n, UInt size)
   *  @return a const_vector_iterator
   */
  inline const_vector_iterator end_reinterpret(UInt n, UInt size) const;


  /*! Get an iterator that behaves like a pointer akantu::Matrix<T> * to the
   *  first tuple of the array.
   *
   *  The reinterpret iterators allow to iterate over an array in any way that
   *  preserves the number of entries of the array. This can for instance be use
   *  full if the shape of the data in an array is not initially known.
   *  @param m number of rows
   *  @param n number of columns
   *  @param size number of tuples in array. m times n times size must match the number
   *  of entries of the array. If compiled in debug mode, an incorrect
   *  combination of m, n and size will result
   *  in an exception being thrown. Optimized code will fail in an unpredicted
   *  manner.
   *  @return a matrix_iterator
   */
  inline matrix_iterator begin_reinterpret(UInt m, UInt n, UInt size);
  /*! Get an iterator that behaves like a pointer akantu::Matrix<T> * pointing
   *  *past* the last tuple of the array.
   *  @param m number of rows
   *  @param n number of columns
   *  @param size number of tuples in array. See Array::begin_reinterpret(UInt m, UInt n, UInt size)
   *  @return a matrix_iterator
   */
  inline matrix_iterator end_reinterpret(UInt m, UInt n, UInt size);
  /*! Get a const iterator that behaves like a pointer akantu::Matrix<T> * to the
   *  first tuple of the array.
   *  @param m number of rows
   *  @param n number of columns
   *  @param size number of tuples in array. See Array::begin_reinterpret(UInt m, UInt n, UInt size)
   *  @return a const_matrix_iterator
   */
  inline const_matrix_iterator begin_reinterpret(UInt m, UInt n, UInt size) const;
  /*! Get a const iterator that behaves like a pointer akantu::Matrix<T> * pointing
   *  *past* the last tuple of the array.
   *  @param m number of rows
   *  @param n number of columns
   *  @param size number of tuples in array. See Array::begin_reinterpret(UInt m, UInt n, UInt size)
   *  @return a const_matrix_iterator
   */
  inline const_matrix_iterator end_reinterpret(UInt m, UInt n, UInt size) const;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /*! append a tuple to the array with the value value for all
   * components
   * @param value the new last tuple or the array will contain nb_component copies of value
   */
  inline void push_back(const_reference value);

  /*! append a tuple to the array
   *  @param new_elem a C-array containing the values to be copied to the end of the array */
  inline void push_back(const value_type new_elem[]);

  /*! append a tuple to the array
   * @param it an iterator to the tuple to be copied to the end of the array */
  template<typename Ret>
  inline void push_back(const iterator<Ret> & it);

  /*! append a matrix or a vector to the array
   *  @param new_elem a reference to a Matrix<T> or Vector<T> */
  template<template<typename> class C>
  inline void push_back(const C<T> & new_elem);

  /**
   * erase an element. If the erased element is not the last of the array, the
   * last element is moved into the hole in order to maintain contiguity. This
   * may invalidate existing iterators (For instance an iterator obtained by
   * Array::end() is no longer correct) and will change the order of the
   * elements.
   * @param i index of element to erase
   */
  inline void erase(UInt i);


  /// ask Nico, clarify
  template<typename R>
  inline iterator<R> erase(const iterator<R> & it);


  /*! change the size of the array and allocate or free memory if needed. If the
   *  size increases, the new tuples are filled with zeros
   *  @param size new number of tuples contained in the array */
  void resize(UInt size);

  /// change the number of components by interlacing data
  /// ask Nico, clarify
  void extendComponentsInterlaced(UInt multiplicator, UInt stride);

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  //  Array<T, is_scal>& operator=(const Array<T, is_scal>& vect);

  /// search elem in the vector, return  the position of the first occurrence or
  /// -1 if not found
  /// @param elem the element to look for
  /// @return index of the first occurrence of elem or -1 if elem is not present
  Int find(const_reference elem) const;\
  /// @see Array::find(const_reference elem) const
  Int find(T elem[]) const;

  /// set all entries of the array to 0
  inline void clear() { std::fill_n(values, size*nb_component, T()); }

  /// set all entries of the array to the value t
  /// @param t value to fill the array with
  inline void set(T t) { std::fill_n(values, size*nb_component, t); }

  /// set all tuples of the array to a given vector or matrix
  /// @param vm Matrix or Vector to fill the array with
  template<template<typename> class C>
  inline void set(const C<T> & vm);

  /*! copy the content of another array. This overwrites the current content.
   *  @param other Array to copy into this array. It has to have the same
   *  nb_component as this. If compiled in debug mode, an incorrect other will
   *  result in an exception being thrown. Optimised code may result in
   *  unpredicted behaviour. */
  void copy(const Array<T, is_scal> & other);

  /// give the address of the memory allocated for this vector
  T * storage() const { return values; };

protected:
  /// perform the allocation for the constructors
  void allocate(UInt size, UInt nb_component = 1);

  /// resize without initializing the memory
  void resizeUnitialized(UInt new_size);

  /* ------------------------------------------------------------------------ */
  /* Operators                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /*! Subtract another array entry by entry from this array in place. Both arrays must
   *  have the same size and nb_component. If the arrays have different shapes,
   *  code compiled in debug mode will throw an expeption and optimised code
   *  will behave in an unpredicted manner
   *  @param other array to subtract from this
   *  @return reference to modified this */
  Array<T, is_scal> & operator-=(const Array<T, is_scal> & other);
  /*! Add another array entry by entry to this array in place. Both arrays must
   *  have the same size and nb_component. If the arrays have different shapes,
   *  code compiled in debug mode will throw an expeption and optimised code
   *  will behave in an unpredicted manner
   *  @param other array to add to this
   *  @return reference to modified this */
  Array<T, is_scal> & operator+=(const Array<T, is_scal> & other);
  /*! Multiply all entries of this array by a scalar in place
   *  @param alpha scalar multiplicant
   *  @return reference to modified this */
  array<T, is_scal> & operator*=(const T & alpha);

  /*! Compare this array element by element to another.
   *  @param other array to compare to
   *  @return true it all element are equal and arrays have the same shape, else false */
  bool operator==(const Array<T, is_scal> & other) const;
  /// @see Array::operator==(const Array<T, is_scal> & other) const
  bool operator!=(const Array<T, is_scal> & other) const;

  /// return a reference to the j-th entry of the i-th tuple
  inline reference operator()(UInt i, UInt j = 0);
  /// return a const reference to the j-th entry of the i-th tuple
  inline const_reference operator()(UInt i, UInt j = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the number of tuples contained in the array
  UInt getSize() const{ return this->size; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
public:
  /// array of values
  T * values; // /!\ very dangerous

};

__END_AKANTU__

#include "aka_types.hh"

__BEGIN_AKANTU__

#include "aka_array_tmpl.hh"

/* -------------------------------------------------------------------------- */
/* Inline Functions Array<T, is_scal>                                         */
/* -------------------------------------------------------------------------- */
template <typename T, bool is_scal>
inline std::ostream & operator<<(std::ostream & stream, const Array<T, is_scal> & _this)
{
  _this.printself(stream);
  return stream;
}

/* -------------------------------------------------------------------------- */
/* Inline Functions ArrayBase                                                 */
/* -------------------------------------------------------------------------- */
inline std::ostream & operator<<(std::ostream & stream, const ArrayBase & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__


#endif /* __AKANTU_VECTOR_HH__ */
