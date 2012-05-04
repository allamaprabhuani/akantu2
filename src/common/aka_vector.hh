/**
 * @file   aka_vector.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jun 17 10:04:55 2010
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
#include <vector>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class Matrix;

/// class that afford to store vectors in static memory
class VectorBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  VectorBase(const ID & id = "");

  virtual ~VectorBase();

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


namespace types {
  class Matrix;
  template<typename T> class Vector;
}
/* -------------------------------------------------------------------------- */
template<typename T> class Vector : public VectorBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  typedef T                  value_type;
  typedef value_type &       reference;
  typedef value_type *       pointer_type;
  typedef const value_type & const_reference;

  /// Allocation of a new vector
  Vector(UInt size = 0, UInt nb_component = 1,
	 const ID & id = "");

  /// Allocation of a new vector with a default value
  Vector(UInt size, UInt nb_component,
  	 const value_type def_values[], const ID & id = "");

  /// Allocation of a new vector with a default value
  Vector(UInt size, UInt nb_component,
	 const_reference value, const ID & id = "");

  /// Copy constructor (deep copy if deep=true)
  Vector(const Vector<value_type>& vect, bool deep = true, const ID & id = "");

  /// Copy constructor (deep copy)
  Vector(const std::vector<value_type>& vect);


  virtual ~Vector();

  /* ------------------------------------------------------------------------ */
  /* Iterator                                                                 */
  /* ------------------------------------------------------------------------ */
  /// \todo protected: does not compile with intel  check why 
public:
  template <class R, class IR = R, int fps = 0>
  class iterator_internal;
public:
  /* ------------------------------------------------------------------------ */
  // template<typename R, int fps = 0> class iterator : public iterator_internal<R> {};
  // template<typename R, int fps = 0> class const_iterator : public iterator_internal<const R> {};
  /* ------------------------------------------------------------------------ */
  template<typename R, int fps = 0>
  class iterator : public iterator_internal<R> {
  public:
    typedef iterator_internal<R> parent;
    typedef typename parent::pointer pointer;
  public:
    iterator() : parent() {};
    iterator(pointer_type data, UInt offset) : parent(data, offset) {};
    iterator(pointer warped) : parent(warped) {};
    iterator(const iterator & it) : parent(it) {};
  };

  /* ------------------------------------------------------------------------ */
  template<typename R, int fps = 0>
  class const_iterator : public iterator_internal<const R, R> {
  public:
    typedef iterator_internal<const R, R> parent;
    typedef typename parent::pointer pointer;
  public:
    const_iterator() : parent() {};
    const_iterator(pointer_type data, UInt offset) : parent(data, offset) {};
    const_iterator(pointer warped) : parent(warped) {};
    const_iterator(const const_iterator & it) : parent(it) {};
  };

  /// specialization of the previous iterators for the scalar types
  /* ------------------------------------------------------------------------ */
  template<int fps>
  class iterator<T, fps> : public iterator_internal<T> {
  public:
    typedef iterator_internal<T> parent;
    typedef typename parent::pointer pointer;
  public:
    iterator() : parent() {};
    iterator(pointer_type data, UInt offset) : parent(data, offset) {};
    iterator(pointer warped) : parent(warped) {};
    iterator(const iterator & it) : parent(it) {};
  };

  /* -------------------------------------------------------------------------- */
  template<int fps>
  class const_iterator<T, fps> : public iterator_internal<const T, T> {
  public:
    typedef iterator_internal<const T, T> parent;
    typedef typename parent::pointer pointer;
  public:
    const_iterator() : parent() {};
    const_iterator(pointer_type data, UInt offset) : parent(data, offset) {};
    const_iterator(pointer warped) : parent(warped) {};
    const_iterator(const const_iterator & it) : parent(it) {};
  };


  /* ------------------------------------------------------------------------ */
  // template<typename Ret> inline iterator<Ret> begin();
  // template<typename Ret> inline iterator<Ret> end();
  // template<typename Ret> inline const_iterator<Ret> begin() const;
  // template<typename Ret> inline const_iterator<Ret> end() const;

  inline iterator<T> begin();
  inline iterator<T> end();
  inline const_iterator<T> begin() const;
  inline const_iterator<T> end() const;

  inline iterator< types::Vector<T> > begin(UInt n);
  inline iterator< types::Vector<T> > end(UInt n);
  inline const_iterator< types::Vector<T> > begin(UInt n) const;
  inline const_iterator< types::Vector<T> > end(UInt n) const;

  inline iterator< types::Matrix > begin(UInt m, UInt n);
  inline iterator< types::Matrix > end(UInt m, UInt n);
  inline const_iterator< types::Matrix > begin(UInt m, UInt n) const;
  inline const_iterator< types::Matrix > end(UInt m, UInt n) const;

  inline iterator< types::Matrix > begin_reinterpret(UInt m, UInt n,
						     UInt size, UInt nb_component);
  inline iterator< types::Matrix > end_reinterpret(UInt m, UInt n, UInt size, UInt nb_component);
  inline const_iterator< types::Matrix > begin_reinterpret(UInt m, UInt n,
							   UInt size, UInt nb_component) const;
  inline const_iterator< types::Matrix > end_reinterpret(UInt m, UInt n,
							 UInt size, UInt nb_component) const;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// get jth componemt of the ith tuple in read-only
  inline const_reference get(UInt i, UInt j = 0) const;
  /// get jth component of the ith tuple
  inline reference at(UInt i, UInt j = 0);

  /// add an  element at  the end  of the vector  with the  value value  for all
  /// component
  inline void push_back(const_reference value);

  /// add an element at the end of the vector
  inline void push_back(const value_type new_elem[]);

  /**
   * remove an element and move the last one in the hole
   * /!\ change the order in the vector
   */
  inline void erase(UInt i);

  template<typename R>
  inline void erase(const iterator<R> & it);


  /// change the size of the vector and allocate more memory if needed
  void resize(UInt size);

  /// change the number of components by interlacing data
  void extendComponentsInterlaced(UInt multiplicator, UInt stride);

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  //  Vector<T>& operator=(const Vector<T>& vect);

  /// search elem in the vector, return  the position of the first occurrence or
  /// -1 if not found
  Int find(const_reference elem) const;
  Int find(T elem[]) const;

  /// set a vvector to 0
  inline void clear() { std::fill_n(values, size*nb_component, T()); };

  /// copy the content of an other vector
  void copy(const Vector<T> & vect);

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

  Vector<T> & operator-=(const Vector<T> & vect);
  Vector<T> & operator+=(const Vector<T> & vect);
  Vector<T> & operator*=(const T & alpha);

  inline reference operator()(UInt i, UInt j = 0);
  inline const_reference operator()(UInt i, UInt j = 0) const;

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
  T * values; // /!\ very dangerous

};

__END_AKANTU__

#include "aka_types.hh"

__BEGIN_AKANTU__

#include "aka_vector_tmpl.hh"


/* -------------------------------------------------------------------------- */
/* Inline Functions Vector<T>                                                 */
/* -------------------------------------------------------------------------- */
template <typename T>
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
