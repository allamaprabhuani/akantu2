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

  inline iterator<T> begin();
  inline iterator<T> end();
  inline const_iterator<T> begin() const;
  inline const_iterator<T> end() const;

  inline vector_iterator begin(UInt n);
  inline vector_iterator end(UInt n);
  inline const_vector_iterator begin(UInt n) const;
  inline const_vector_iterator end(UInt n) const;

  inline matrix_iterator begin(UInt m, UInt n);
  inline matrix_iterator end(UInt m, UInt n);
  inline const_matrix_iterator begin(UInt m, UInt n) const;
  inline const_matrix_iterator end(UInt m, UInt n) const;

  /// /!\ to use with caution
  inline vector_iterator begin_reinterpret(UInt n, UInt size);
  inline vector_iterator end_reinterpret(UInt n, UInt size);
  inline const_vector_iterator begin_reinterpret(UInt n, UInt size) const;
  inline const_vector_iterator end_reinterpret(UInt n, UInt size) const;


  inline matrix_iterator begin_reinterpret(UInt m, UInt n, UInt size);
  inline matrix_iterator end_reinterpret(UInt m, UInt n, UInt size);
  inline const_matrix_iterator begin_reinterpret(UInt m, UInt n, UInt size) const;
  inline const_matrix_iterator end_reinterpret(UInt m, UInt n, UInt size) const;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// add an  element at  the end  of the vector  with the  value value  for all
  /// component
  inline void push_back(const_reference value);

  /// add an element at the end of the vector
  inline void push_back(const value_type new_elem[]);

  template<typename Ret>
  inline void push_back(const iterator<Ret> & it);

  /// push_back for Vector<T> and Matrix<T>
  template<template<typename> class C>
  inline void push_back(const C<T> & new_elem);

  /**
   * remove an element and move the last one in the hole
   * /!\ change the order in the vector
   */
  inline void erase(UInt i);

  template<typename R>
  inline iterator<R> erase(const iterator<R> & it);


  /// change the size of the vector and allocate more memory if needed
  void resize(UInt size);

  /// change the number of components by interlacing data
  void extendComponentsInterlaced(UInt multiplicator, UInt stride);

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  //  Array<T, is_scal>& operator=(const Array<T, is_scal>& vect);

  /// search elem in the vector, return  the position of the first occurrence or
  /// -1 if not found
  Int find(const_reference elem) const;
  Int find(T elem[]) const;

  /// set a vvector to 0
  inline void clear() { std::fill_n(values, size*nb_component, T()); }

  /// set a vector to the value t
  inline void set(T t) { std::fill_n(values, size*nb_component, t); }

  template<template<typename> class C>
  inline void set(const C<T> & vm);

  /// copy the content of an other vector
  void copy(const Array<T, is_scal> & vect);

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

  Array<T, is_scal> & operator-=(const Array<T, is_scal> & vect);
  Array<T, is_scal> & operator+=(const Array<T, is_scal> & vect);
  Array<T, is_scal> & operator*=(const T & alpha);

  bool operator==(const Array<T, is_scal> & vect) const;
  bool operator!=(const Array<T, is_scal> & vect) const;

  inline reference operator()(UInt i, UInt j = 0);
  inline const_reference operator()(UInt i, UInt j = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
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
