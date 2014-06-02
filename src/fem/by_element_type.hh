/**
 * @file   by_element_type.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Aug 31 11:09:48 2011
 *
 * @brief  storage class by element type
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either dataversion 3 of the License, or (at your option) any
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

#ifndef __AKANTU_BY_ELEMENT_TYPE_HH__
#define __AKANTU_BY_ELEMENT_TYPE_HH__

#include "aka_common.hh"
#include "aka_array.hh"
#include "aka_memory.hh"

__BEGIN_AKANTU__


// Forward declaration of the ByElementTypeVector class
template<class Stored, typename SupportType = ElementType>
class ByElementType;

/* -------------------------------------------------------------------------- */
/* ByElementTypeBase                                                          */
/* -------------------------------------------------------------------------- */
/// Common non templated base class for the ByElementType class
class ByElementTypeBase {
public:
  virtual ~ByElementTypeBase() {};
};

/* -------------------------------------------------------------------------- */
/* ByElementType                                                              */
/* -------------------------------------------------------------------------- */

template<class Stored, typename SupportType>
class ByElementType : public ByElementTypeBase {
private:
  void operator=(const ByElementType &) {};
public:
  ByElementType();
  ~ByElementType();

  inline static std::string printType(const SupportType & type, const GhostType & ghost_type);

  /*! Tests whether a type is present in the object
   *  @param type the type to check for
   *  @param ghost_type optional: by default, the data map for non-ghost
   *         elements is searched
   *  @return true if the type is present. */
  inline bool exists(const SupportType & type, const GhostType & ghost_type = _not_ghost) const;

  /*! get the stored data corresponding to a type
   *  @param type the type to check for
   *  @param ghost_type optional: by default, the data map for non-ghost
   *         elements is searched
   *  @return stored data corresponding to type. */
  inline const Stored & operator()(const SupportType & type,
				   const GhostType & ghost_type = _not_ghost) const;

  /*! get the stored data corresponding to a type
   *  @param type the type to check for
   *  @param ghost_type optional: by default, the data map for non-ghost
   *         elements is searched
   *  @return stored data corresponding to type. */
  inline Stored & operator()(const SupportType & type,
			     const GhostType & ghost_type = _not_ghost);

  /*! insert data of a new type (not yet present) into the map. THIS METHOD IS
   *  NOT ARRAY SAFE, when using ByElementTypeArray, use setArray instead
   *  @param data to insert
   *  @param type type of data (if this type is already present in the map,
   *         an exception is thrown).
   *  @param ghost_type optional: by default, the data map for non-ghost
   *         elements is searched
   *  @return stored data corresponding to type. */
  inline Stored & operator()(const Stored & insert,
			     const SupportType & type,
			     const GhostType & ghost_type = _not_ghost);

  /// print helper
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Element type Iterator                                                    */
  /* ------------------------------------------------------------------------ */
  /*! iterator allows to iterate over type-data pairs of the map. The interface
   *  expects the SupportType to be ElementType. */
  typedef std::map<SupportType, Stored> DataMap;
  class type_iterator : private std::iterator<std::forward_iterator_tag, const SupportType> {
  public:
    typedef const SupportType   value_type;
    typedef const SupportType*  pointer;
    typedef const SupportType&  reference;
  protected:
    typedef typename ByElementType<Stored>::DataMap::const_iterator DataMapIterator;
  public:
    type_iterator(DataMapIterator & list_begin,
		  DataMapIterator & list_end,
		  UInt dim,
		  ElementKind ek);

    type_iterator(const type_iterator & it);
    type_iterator() {}

    inline reference operator*();
    inline reference operator*() const;
    inline type_iterator & operator++();
    type_iterator operator++(int);
    inline bool operator==(const type_iterator & other) const;
    inline bool operator!=(const type_iterator & other) const;
    type_iterator & operator=(const type_iterator & other);

  private:
    DataMapIterator list_begin;
    DataMapIterator list_end;
    UInt dim;
    ElementKind kind;
  };

  /*! Get an iterator to the beginning of a subset datamap. This method expects
   *  the SupportType to be ElementType.
   *  @param dim optional: iterate over data of dimension dim (e.g. when
   *         iterating over (surface) facets of a 3D mesh, dim would be 2).
   *         by default, all dimensions are considered.
   *  @param ghost_type optional: by default, the data map for non-ghost
   *         elements is iterated over.
   *  @param kind optional: the kind of element to search for (see
   *         aka_common.hh), by default all kinds are considered
   *  @return an iterator to the first stored data matching the filters
   *          or an iterator to the end of the map if none match*/
  inline type_iterator firstType(UInt dim = _all_dimensions,
				 GhostType ghost_type = _not_ghost,
				 ElementKind kind = _ek_not_defined) const;
  /*! Get an iterator to the end of a subset datamap. This method expects
   *  the SupportType to be ElementType.
   *  @param dim optional: iterate over data of dimension dim (e.g. when
   *         iterating over (surface) facets of a 3D mesh, dim would be 2).
   *         by default, all dimensions are considered.
   *  @param ghost_type optional: by default, the data map for non-ghost
   *         elements is iterated over.
   *  @param kind optional: the kind of element to search for (see
   *         aka_common.hh), by default all kinds are considered
   *  @return an iterator to the last stored data matching the filters
   *          or an iterator to the end of the map if none match */
  inline type_iterator lastType(UInt dim = _all_dimensions,
				GhostType ghost_type = _not_ghost,
				ElementKind kind = _ek_not_defined) const;


protected:
  /*! Direct access to the underlying data map. for internal use by daughter
   *  classes only
   *  @param ghost_type whether to return the data map or the ghost_data map
   *  @return the raw map */
  inline DataMap & getData(GhostType ghost_type);
  /*! Direct access to the underlying data map. for internal use by daughter
   *  classes only
   *  @param ghost_type whether to return the data map or the ghost_data map
   *  @return the raw map */
  inline const DataMap & getData(GhostType ghost_type) const;

  /* ------------------------------------------------------------------------ */
protected:
  DataMap data;
  DataMap ghost_data;
};

/* -------------------------------------------------------------------------- */
/* Some typedefs                                                              */
/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType = ElementType>
class ByElementTypeArray : public ByElementType<Array<T> *, SupportType>, public Memory {
protected:
  typedef ByElementType<Array<T> *, SupportType> parent;
  typedef typename parent::DataMap DataMap;
private:
private:
  /// standard assigment (copy) operator
  void operator=(const ByElementType<T, SupportType> &) {};

public:
  typedef typename parent::type_iterator type_iterator;

  /*! Constructor
   *  @param id optional: identifier (string)
   *  @param parent_id optional: parent identifier. for organizational purposes
   *         only
   *  @param memory_id optional: choose a specific memory, defaults to memory 0
   */
  ByElementTypeArray(const ID & id = "by_element_type_array", const ID & parent_id = "no_parent",
                     const MemoryID & memory_id = 0) :
    parent(), Memory(parent_id + ":" + id, memory_id) {};

  /*! allocate memory for a new array
   *  @param size number of tuples of the new array
   *  @param nb_component tuple size
   *  @param type the type under which the array is indexed in the map
   *  @param ghost_type whether to add the field to the data map or the
   *         ghost_data map
   *  @return a reference to the allocated array */
  inline Array<T> & alloc(UInt size,
			  UInt nb_component,
			  const SupportType & type,
			  const GhostType & ghost_type);

  /*! allocate memory for a new array in both the data and the ghost_data map
   *  @param size number of tuples of the new array
   *  @param nb_component tuple size
   *  @param type the type under which the array is indexed in the map*/
  inline void alloc(UInt size,
		    UInt nb_component,
		    const SupportType & type);

  /* get a reference to the array of certain type
   * @param type data filed under type is returned
   * @param ghost_type optional: by default the non-ghost map is searched
   * @return a reference to the array */
  inline const Array<T> & operator()(const SupportType & type,
				      const GhostType & ghost_type = _not_ghost) const;

  /* get a reference to the array of certain type
   * @param type data filed under type is returned
   * @param ghost_type optional: by default the non-ghost map is searched
   * @return a const reference to the array */
  inline Array<T> & operator()(const SupportType & type,
                               const GhostType & ghost_type = _not_ghost);

  /*! insert data of a new type (not yet present) into the map.
   *  @param type type of data (if this type is already present in the map,
   *         an exception is thrown).
   *  @param ghost_type optional: by default, the data map for non-ghost
   *         elements is searched
   *  @param vect the vector to include into the map
   *  @return stored data corresponding to type. */
  inline void setArray(const SupportType & type,
			const GhostType & ghost_type,
			const Array<T> & vect);
  /*! frees all memory related to the data*/
  inline void free();

  /*! deletes and reorders entries in the stored arrays
   *  @param new_numbering a ByElementTypeArray of new indices. UInt(-1) indicates
   *         deleted entries. */
  inline void onElementsRemoved(const ByElementTypeArray<UInt> & new_numbering);

  /// text output helper
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /*! set the id
   *  @param id the new name
   */
  inline void setID(const ID & id) { this->id = id; }

private:
  ByElementTypeArray operator=(__attribute__((unused)) const ByElementTypeArray & other) {};
};

/// to store data Array<Real> by element type
typedef ByElementTypeArray<Real> ByElementTypeReal;
/// to store data Array<Int> by element type
typedef ByElementTypeArray<Int>  ByElementTypeInt;
/// to store data Array<UInt> by element type
typedef ByElementTypeArray<UInt, ElementType> ByElementTypeUInt;

/// Map of data of type UInt stored in a mesh
typedef std::map<std::string, Array<UInt> *> UIntDataMap;
typedef ByElementType<UIntDataMap, ElementType> ByElementTypeUIntDataMap;

__END_AKANTU__

#endif /* __AKANTU_BY_ELEMENT_TYPE_HH__ */
