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
#include <typeinfo>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class Matrix;


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
	 const VectorID & id = "");

  /// Allocation of a new vector with a default value
  Vector(UInt size, UInt nb_component,
  	 const value_type def_values[], const VectorID & id = "");

  /// Allocation of a new vector with a default value
  Vector(UInt size, UInt nb_component,
	 const_reference value, const VectorID & id = "");

  /// Copy constructor (deep copy if deep=true) \todo to implement
  Vector(const Vector<value_type>& vect, bool deep = true, const VectorID & id = "");

  virtual ~Vector();

  /* ------------------------------------------------------------------------ */
  /* Iterator                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  template<typename R, int for_part_spec = 0>
  class iterator {
  public:
    typedef R returned_type;
    typedef returned_type & returned_type_ref;
    typedef returned_type * returned_type_ptr;

    iterator();
    iterator(pointer_type data, UInt offset);
    iterator(const iterator & it);

    ~iterator();

    inline iterator & operator=(const iterator & it);

    inline returned_type_ref operator*() { return *ret; };
    inline returned_type_ptr operator->() { return ret; };
    inline iterator & operator++();

    inline iterator & operator+=(const UInt n);
    inline returned_type_ref operator[](const UInt n);

    inline bool operator==(const iterator & other);
    inline bool operator!=(const iterator & other);

  private:
    UInt offset;
    pointer_type initial;
    returned_type_ptr ret;
  };

  template<typename Ret> inline iterator<Ret> begin();
  template<typename Ret> inline iterator<Ret> end();

  inline iterator<types::Vector<T> > begin(UInt n);
  inline iterator<types::Vector<T> > end(UInt n);

  inline iterator<types::Matrix> begin(UInt m, UInt n);
  inline iterator<types::Matrix> end(UInt m, UInt n);

  inline reference operator()(UInt i, UInt j = 0);
  inline const_reference operator()(UInt i, UInt j = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// get jth componemt of the ith tuple in read-only
  inline const_reference get(UInt i, UInt j = 0) const;
  /// get jth componemt of the ith tuple
  inline reference at(UInt i, UInt j = 0);

  /// add an  element at  the and  of the vector  with the  value value  for all
  /// component
  inline void push_back(const_reference value);

  /// add an element at the and of the vector
  inline void push_back(const value_type new_elem[]);

  /**
   * remove an element and move the last one in the hole
   * /!\ change the order in the vector
   */
  inline void erase(UInt i);

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

  /// set a vvector to 0
  inline void clear() { memset(values, 0, size*nb_component*sizeof(T)); };

  /// copy the content of an other vector
  void copy(const Vector<T> & vect);

protected:
  /// perform the allocation for the constructors
  void allocate(UInt size, UInt nb_component = 1);

  /* ------------------------------------------------------------------------ */
  /* Operators                                                                */
  /* ------------------------------------------------------------------------ */
public:

  Vector<T> & operator-=(const Vector<T> & vect);
  Vector<T> & operator+=(const Vector<T> & vect);
  Vector<T> & operator*=(const T & alpha);

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




__END_AKANTU__

#include "aka_types.hh"

__BEGIN_AKANTU__

#include "aka_vector_inline_impl.cc"

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
