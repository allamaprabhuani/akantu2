/**
 * @file   element_type_map.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 31 2011
 * @date last modification: Fri Oct 02 2015
 *
 * @brief  storage class by element type
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "aka_array.hh"
#include "aka_memory.hh"
#include "aka_named_argument.hh"
#include "element.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ELEMENT_TYPE_MAP_HH__
#define __AKANTU_ELEMENT_TYPE_MAP_HH__

namespace akantu {
class FEEngine;
} // namespace akantu

namespace akantu {

namespace {
  DECLARE_NAMED_ARGUMENT(all_ghost_types);
  DECLARE_NAMED_ARGUMENT(default_value);
  DECLARE_NAMED_ARGUMENT(element_kind);
  DECLARE_NAMED_ARGUMENT(ghost_type);
  DECLARE_NAMED_ARGUMENT(nb_component);
  DECLARE_NAMED_ARGUMENT(with_nb_element);
  DECLARE_NAMED_ARGUMENT(with_nb_nodes_per_element);
  DECLARE_NAMED_ARGUMENT(spatial_dimension);
} // namespace

template <class Stored, typename SupportType = ElementType>
class ElementTypeMap;

/* -------------------------------------------------------------------------- */
/* ElementTypeMapBase */
/* -------------------------------------------------------------------------- */
/// Common non templated base class for the ElementTypeMap class
class ElementTypeMapBase {
public:
  virtual ~ElementTypeMapBase() = default;
};

/* -------------------------------------------------------------------------- */
/* ElementTypeMap */
/* -------------------------------------------------------------------------- */

template <class Stored, typename SupportType>
class ElementTypeMap : public ElementTypeMapBase {

public:
  ElementTypeMap();
  ~ElementTypeMap();

  inline static std::string printType(const SupportType & type,
                                      const GhostType & ghost_type);

  /*! Tests whether a type is present in the object
   *  @param type the type to check for
   *  @param ghost_type optional: by default, the data map for non-ghost
   *         elements is searched
   *  @return true if the type is present. */
  inline bool exists(const SupportType & type,
                     const GhostType & ghost_type = _not_ghost) const;

  /*! get the stored data corresponding to a type
   *  @param type the type to check for
   *  @param ghost_type optional: by default, the data map for non-ghost
   *         elements is searched
   *  @return stored data corresponding to type. */
  inline const Stored &
  operator()(const SupportType & type,
             const GhostType & ghost_type = _not_ghost) const;

  /*! get the stored data corresponding to a type
   *  @param type the type to check for
   *  @param ghost_type optional: by default, the data map for non-ghost
   *         elements is searched
   *  @return stored data corresponding to type. */
  inline Stored & operator()(const SupportType & type,
                             const GhostType & ghost_type = _not_ghost);

  /*! insert data of a new type (not yet present) into the map. THIS METHOD IS
   *  NOT ARRAY SAFE, when using ElementTypeMapArray, use setArray instead
   *  @param data to insert
   *  @param type type of data (if this type is already present in the map,
   *         an exception is thrown).
   *  @param ghost_type optional: by default, the data map for non-ghost
   *         elements is searched
   *  @return stored data corresponding to type. */
  inline Stored & operator()(const Stored & insert, const SupportType & type,
                             const GhostType & ghost_type = _not_ghost);

  /// print helper
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Element type Iterator                                                    */
  /* ------------------------------------------------------------------------ */
  /*! iterator allows to iterate over type-data pairs of the map. The interface
   *  expects the SupportType to be ElementType. */
  typedef std::map<SupportType, Stored> DataMap;
  class type_iterator
      : private std::iterator<std::forward_iterator_tag, const SupportType> {
  public:
    using value_type = const SupportType;
    using pointer = const SupportType *;
    using reference = const SupportType &;

  protected:
    using DataMapIterator =
        typename ElementTypeMap<Stored>::DataMap::const_iterator;

  public:
    type_iterator(DataMapIterator & list_begin, DataMapIterator & list_end,
                  UInt dim, ElementKind ek);

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

  /// helper class to use in range for constructions
  class ElementTypesIteratorHelper {
  public:
    using Container = ElementTypeMap<Stored, SupportType>;
    using iterator = typename Container::type_iterator;

    ElementTypesIteratorHelper(const Container & container, UInt dim,
                               GhostType ghost_type, ElementKind kind)
        : container(std::cref(container)), dim(dim), ghost_type(ghost_type),
          kind(kind) {}

    template <typename... pack>
    ElementTypesIteratorHelper(const Container & container, use_named_args_t,
                               pack &&... _pack)
        : ElementTypesIteratorHelper(
              container, OPTIONAL_NAMED_ARG(spatial_dimension, _all_dimensions),
              OPTIONAL_NAMED_ARG(ghost_type, _not_ghost),
              OPTIONAL_NAMED_ARG(element_kind, _ek_regular)) {}

    ElementTypesIteratorHelper(const ElementTypesIteratorHelper &) = default;
    ElementTypesIteratorHelper &
    operator=(const ElementTypesIteratorHelper &) = default;
    ElementTypesIteratorHelper &
    operator=(ElementTypesIteratorHelper &&) = default;

    iterator begin() {
      return container.get().firstType(dim, ghost_type, kind);
    }
    iterator end() { return container.get().lastType(dim, ghost_type, kind); }

  private:
    std::reference_wrapper<const Container> container;
    UInt dim;
    GhostType ghost_type;
    ElementKind kind;
  };

private:
  ElementTypesIteratorHelper
  elementTypesImpl(UInt dim = _all_dimensions,
                   GhostType ghost_type = _not_ghost,
                   ElementKind kind = _ek_regular) const;

  template <typename... pack>
  ElementTypesIteratorHelper
  elementTypesImpl(const use_named_args_t & /*unused*/, pack &&... _pack) const;

  template <typename... pack>
  using named_argument_test =
      std::conjunction<std::bool_constant<(sizeof...(pack) > 0)>,
                       is_named_argument<pack>...>;

public:
  template <typename... pack>
  std::enable_if_t<named_argument_test<pack...>{}, ElementTypesIteratorHelper>
  elementTypes(pack &&... _pack) const {
    return elementTypesImpl(use_named_args,
                            std::forward<decltype(_pack)>(_pack)...);
  }

  template <typename... pack>
  std::enable_if_t<not named_argument_test<pack...>{},
                   ElementTypesIteratorHelper>
  elementTypes(pack &&... _pack) const {
    return elementTypesImpl(std::forward<decltype(_pack)>(_pack)...);
  }

public:
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
template <typename T, typename SupportType>
class ElementTypeMapArray : public ElementTypeMap<Array<T> *, SupportType>,
                            public Memory {
public:
  using type = T;
  using array_type = Array<T>;

protected:
  using parent = ElementTypeMap<Array<T> *, SupportType>;
  using DataMap = typename parent::DataMap;

public:
  using type_iterator = typename parent::type_iterator;

  /// standard assigment (copy) operator
  void operator=(const ElementTypeMapArray &) = delete;
  ElementTypeMapArray(const ElementTypeMapArray &) = delete;

  /*! Constructor
   *  @param id optional: identifier (string)
   *  @param parent_id optional: parent identifier. for organizational purposes
   *         only
   *  @param memory_id optional: choose a specific memory, defaults to memory 0
   */
  ElementTypeMapArray(const ID & id = "by_element_type_array",
                      const ID & parent_id = "no_parent",
                      const MemoryID & memory_id = 0)
      : parent(), Memory(parent_id + ":" + id, memory_id), name(id){};

  /*! allocate memory for a new array
   *  @param size number of tuples of the new array
   *  @param nb_component tuple size
   *  @param type the type under which the array is indexed in the map
   *  @param ghost_type whether to add the field to the data map or the
   *         ghost_data map
   *  @return a reference to the allocated array */
  inline Array<T> & alloc(UInt size, UInt nb_component,
                          const SupportType & type,
                          const GhostType & ghost_type,
                          const T & default_value = T());

  /*! allocate memory for a new array in both the data and the ghost_data map
   *  @param size number of tuples of the new array
   *  @param nb_component tuple size
   *  @param type the type under which the array is indexed in the map*/
  inline void alloc(UInt size, UInt nb_component, const SupportType & type,
                    const T & default_value = T());

  /* get a reference to the array of certain type
   * @param type data filed under type is returned
   * @param ghost_type optional: by default the non-ghost map is searched
   * @return a reference to the array */
  inline const Array<T> &
  operator()(const SupportType & type,
             const GhostType & ghost_type = _not_ghost) const;

  /// access the data of an element, this combine the map and array accessor
  inline const T & operator()(const Element & element) const;

  /// access the data of an element, this combine the map and array accessor
  inline T & operator()(const Element & element);

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
  inline void setArray(const SupportType & type, const GhostType & ghost_type,
                       const Array<T> & vect);
  /*! frees all memory related to the data*/
  inline void free();

  /*! set all values in the ElementTypeMap to zero*/
  inline void clear();

  /*! deletes and reorders entries in the stored arrays
   *  @param new_numbering a ElementTypeMapArray of new indices. UInt(-1)
   * indicates
   *         deleted entries. */
  inline void
  onElementsRemoved(const ElementTypeMapArray<UInt> & new_numbering);

  /// text output helper
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /*! set the id
   *  @param id the new name
   */
  inline void setID(const ID & id) { this->id = id; }

  ElementTypeMap<UInt>
  getNbComponents(UInt dim = _all_dimensions, GhostType ghost_type = _not_ghost,
                  ElementKind kind = _ek_not_defined) const {

    ElementTypeMap<UInt> nb_components;

    for (auto & type : this->elementTypes(dim, ghost_type, kind)) {
      UInt nb_comp = (*this)(type, ghost_type).getNbComponent();
      nb_components(type, ghost_type) = nb_comp;
    }

    return nb_components;
  }

  /* ------------------------------------------------------------------------ */
  /* more evolved allocators                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the arrays in accordance to a functor
  template <class Func>
  void initialize(const Func & f, const T & default_value = T());

  /// initialize with sizes and number of components in accordance of a mesh
  /// content
  template <typename... pack>
  void initialize(const Mesh & mesh, pack &&... _pack);

  /// initialize with sizes and number of components in accordance of a fe
  /// engine content (aka integration points)
  template <typename... pack>
  void initialize(const FEEngine & fe_engine, pack &&... _pack);

  /* ------------------------------------------------------------------------ */
  /* Accesssors                                                               */
  /* ------------------------------------------------------------------------ */
public:
  /// get the name of the internal field
  AKANTU_GET_MACRO(Name, name, ID);

  /// name of the elment type map: e.g. connectivity, grad_u
  ID name;
};

/// to store data Array<Real> by element type
using ElementTypeMapReal = ElementTypeMapArray<Real>;
/// to store data Array<Int> by element type
using ElementTypeMapInt = ElementTypeMapArray<Int>;
/// to store data Array<UInt> by element type
using ElementTypeMapUInt = ElementTypeMapArray<UInt, ElementType>;

/// Map of data of type UInt stored in a mesh
using UIntDataMap = std::map<std::string, Array<UInt> *>;
using ElementTypeMapUIntDataMap = ElementTypeMap<UIntDataMap, ElementType>;

} // namespace akantu

#endif /* __AKANTU_ELEMENT_TYPE_MAP_HH__ */
