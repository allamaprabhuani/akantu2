/**
 * Copyright (©) 2011-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
#include "aka_named_argument.hh"
#include "element.hh"
/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_ELEMENT_TYPE_MAP_HH_
#define AKANTU_ELEMENT_TYPE_MAP_HH_

namespace akantu {
class FEEngine;
class IntegrationPoint;
} // namespace akantu

namespace akantu {

namespace {
  DECLARE_NAMED_ARGUMENT(all_ghost_types);
  DECLARE_NAMED_ARGUMENT(default_value);
  DECLARE_NAMED_ARGUMENT(element_kind);
  DECLARE_NAMED_ARGUMENT(ghost_type);
  DECLARE_NAMED_ARGUMENT(nb_component);
  DECLARE_NAMED_ARGUMENT(nb_component_functor);
  DECLARE_NAMED_ARGUMENT(with_nb_element);
  DECLARE_NAMED_ARGUMENT(with_nb_nodes_per_element);
  DECLARE_NAMED_ARGUMENT(spatial_dimension);
  DECLARE_NAMED_ARGUMENT(do_not_default);
  DECLARE_NAMED_ARGUMENT(element_filter);
} // namespace

template <class Stored, typename SupportType = ElementType>
class ElementTypeMap;

/* -------------------------------------------------------------------------- */
/* ElementTypeMapBase */
/* -------------------------------------------------------------------------- */
/// Common non templated base class for the ElementTypeMap class
class ElementTypeMapBase {
public:
  ElementTypeMapBase() = default;
  ElementTypeMapBase(const ElementTypeMapBase &) = default;
  ElementTypeMapBase(ElementTypeMapBase &&) = default;
  ElementTypeMapBase & operator=(const ElementTypeMapBase &) = default;
  ElementTypeMapBase & operator=(ElementTypeMapBase &&) = default;
  virtual ~ElementTypeMapBase() = default;
};

/* -------------------------------------------------------------------------- */
/* ElementTypeMap */
/* -------------------------------------------------------------------------- */

template <class Stored, typename SupportType>
class ElementTypeMap : public ElementTypeMapBase {
public:
  using value_type = Stored;

  ElementTypeMap() = default;

  inline static auto printType(SupportType type, GhostType ghost_type)
      -> std::string;

  /*! Tests whether a type is present in the object
   *  @param type the type to check for
   *  @param ghost_type optional: by default, the data map for non-ghost
   *         elements is searched
   *  @return true if the type is present. */
  [[nodiscard]] inline auto exists(SupportType type,
                                   GhostType ghost_type = _not_ghost) const
      -> bool;

  /*! get the stored data corresponding to a type
   *  @param type the type to check for
   *  @param ghost_type optional: by default, the data map for non-ghost
   *         elements is searched
   *  @return stored data corresponding to type. */
  inline auto operator()(SupportType type,
                         GhostType ghost_type = _not_ghost) const
      -> const Stored &;

  /*! get the stored data corresponding to a type
   *  @param type the type to check for
   *  @param ghost_type optional: by default, the data map for non-ghost
   *         elements is searched
   *  @return stored data corresponding to type. */
  inline auto operator()(SupportType type, GhostType ghost_type = _not_ghost)
      -> Stored &;

  /*! insert data of a new type (not yet present) into the map. THIS METHOD IS
   *  NOT ARRAY SAFE, when using ElementTypeMapArray, use setArray instead
   *  @param data to insert
   *  @param type type of data (if this type is already present in the map,
   *         an exception is thrown).
   *  @param ghost_type optional: by default, the data map for non-ghost
   *         elements is searched
   *  @return stored data corresponding to type. */
  template <typename U>
  inline auto operator()(U && insertee, SupportType type,
                         GhostType ghost_type = _not_ghost) -> Stored &;

public:
  /// print helper
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Element type Iterator                                                    */
  /* ------------------------------------------------------------------------ */
  /*! iterator allows to iterate over type-data pairs of the map. The interface
   *  expects the SupportType to be ElementType. */
  using DataMap = std::map<SupportType, Stored>;

  /// helper class to use in range for constructions
  class type_iterator {
  public:
    using value_type = const SupportType;
    using pointer = const SupportType *;
    using reference = SupportType;
    using iterator_category = std::input_iterator_tag;
    using difference_type = Int;

  protected:
    using DataMapIterator =
        typename ElementTypeMap<Stored>::DataMap::const_iterator;

  public:
    type_iterator(DataMapIterator & list_begin, DataMapIterator & list_end,
                  Int dim, ElementKind ek);
    type_iterator() = default;

    inline auto operator*() -> reference;
    inline auto operator*() const -> reference;
    inline auto operator++() -> type_iterator &;
    auto operator++(int) -> type_iterator;
    inline auto operator==(const type_iterator & other) const -> bool;
    inline auto operator!=(const type_iterator & other) const -> bool;

  private:
    DataMapIterator list_begin;
    DataMapIterator list_end;
    Int dim{_all_dimensions};
    ElementKind kind{_ek_not_defined};
  };

  /// helper class to use in range for constructions
  class ElementTypesIteratorHelper {
  public:
    using Container = ElementTypeMap<Stored, SupportType>;
    using iterator = typename Container::type_iterator;

    ElementTypesIteratorHelper(const Container & container, Int dim,
                               GhostType ghost_type, ElementKind kind)
        : container(container), dim(dim), ghost_type(ghost_type), kind(kind) {}

    auto begin() -> iterator;
    auto end() -> iterator;

  private:
    std::reference_wrapper<const Container> container;
    Int dim{_all_dimensions};
    GhostType ghost_type{_not_ghost};
    ElementKind kind{_ek_not_defined};
  };

protected:
  [[nodiscard]] virtual auto
  elementTypesImpl(Int dim = _all_dimensions, GhostType ghost_type = _not_ghost,
                   ElementKind kind = _ek_not_defined) const
      -> ElementTypesIteratorHelper;

public:
  /*!
   * \param _pack
   * \parblock
   *  represent optional parameters:
   * \li \c _spatial_dimension filter for elements of given spatial
   * dimension
   * \li \c _ghost_type filter for a certain ghost_type
   * \li \c _element_kind filter for elements of given kind
   * \endparblock
   */
  template <typename... pack>
  [[nodiscard]] auto elementTypes(pack &&... _pack) const
      -> std::enable_if_t<are_named_argument<pack...>::value,
                          ElementTypesIteratorHelper> {
    return elementTypesImpl(
        OPTIONAL_NAMED_ARG(spatial_dimension, _all_dimensions),
        OPTIONAL_NAMED_ARG(ghost_type, _not_ghost),
        OPTIONAL_NAMED_ARG(element_kind, _ek_not_defined));
  }

  [[nodiscard]] auto elementTypes(Int dim, GhostType ghost_type = _not_ghost,
                                  ElementKind kind = _ek_not_defined) const {
    return elementTypesImpl(dim, ghost_type, kind);
  }

  [[nodiscard]] auto elementTypes(GhostType ghost_type) const {
    return elementTypesImpl(_all_dimensions, ghost_type, _ek_not_defined);
  }

  /*! Direct access to the underlying data map. for internal use by daughter
   *  classes only
   *  @param ghost_type whether to return the data map or the ghost_data map
   *  @return the raw map */
  [[nodiscard]] inline auto getData(GhostType ghost_type) -> DataMap &;
  /*! Direct access to the underlying data map. for internal use by daughter
   *  classes only
   *  @param ghost_type whether to return the data map or the ghost_data map
   *  @return the raw map */
  [[nodiscard]] inline auto getData(GhostType ghost_type) const
      -> const DataMap &;

  /* ------------------------------------------------------------------------ */
protected:
  DataMap data;
  DataMap ghost_data;
};

/* -------------------------------------------------------------------------- */
/* Some typedefs                                                              */
/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
class ElementTypeMapArray
    : public ElementTypeMap<std::unique_ptr<Array<T>>, SupportType> {
public:
  using value_type = T;
  using array_type = Array<T>;

protected:
  using parent = ElementTypeMap<std::unique_ptr<Array<T>>, SupportType>;
  using DataMap = typename parent::DataMap;

public:
  using type_iterator = typename parent::type_iterator;

  /// standard assigment (copy) operator
  auto operator=(const ElementTypeMapArray & other) -> ElementTypeMapArray &;
  ElementTypeMapArray(const ElementTypeMapArray & other);

  /// explicit copy
  void copy(const ElementTypeMapArray & other);

  /*! Constructor
   *  @param id optional: identifier (string)
   *  @param parent_id optional: parent identifier. for organizational purposes
   *         only
   */
  ElementTypeMapArray(const ID & id = "by_element_type_array",
                      const ID & parent_id = "no_parent")
      : parent(), id(parent_id + ":" + id), name(id){};

  /*! allocate memory for a new array
   *  @param size number of tuples of the new array
   *  @param nb_component tuple size
   *  @param type the type under which the array is indexed in the map
   *  @param ghost_type whether to add the field to the data map or the
   *         ghost_data map
   *  @param default_value the default value to use to fill the array
   *  @return a reference to the allocated array */
  inline auto alloc(Int size, Int nb_component, SupportType type,
                    GhostType ghost_type, const T & default_value = T())
      -> Array<T> &;

  /*! allocate memory for a new array in both the data and the ghost_data map
   *  @param size number of tuples of the new array
   *  @param nb_component tuple size
   *  @param type the type under which the array is indexed in the map*/
  inline void alloc(Int size, Int nb_component, SupportType type,
                    const T & default_value = T());

  /* get a reference to the array of certain type
   * @param type data filed under type is returned
   * @param ghost_type optional: by default the non-ghost map is searched
   * @return a reference to the array */
  inline auto operator()(SupportType type,
                         GhostType ghost_type = _not_ghost) const
      -> const Array<T> &;

  /// access the data of an element, this combine the map and array accessor
  inline auto operator()(const Element & element, Int component = 0) const
      -> const T &;

  /// access the data of an element, this combine the map and array accessor
  inline auto operator()(const Element & element, Int component = 0) -> T &;

  /// access the data of an element, this combine the map and array accessor
  inline auto operator()(const IntegrationPoint & point,
                         Int component = 0) const -> const T &;

  /// access the data of an element, this combine the map and array accessor
  inline auto operator()(const IntegrationPoint & point, Int component = 0)
      -> T &;

  /// access the data of an element, this combine the map and array accessor
  inline decltype(auto) get(const Element & element);
  inline decltype(auto) get(const Element & element) const;

  inline decltype(auto) get(const IntegrationPoint & point);
  inline decltype(auto) get(const IntegrationPoint & point) const;

  template <typename... Ns,
            std::enable_if_t<
                std::conjunction_v<std::is_integral<std::decay_t<Ns>>...> and
                sizeof...(Ns) >= 1> * = nullptr>
  inline decltype(auto) get(const Element & element, Ns &&... ns);
  template <typename... Ns,
            std::enable_if_t<
                std::conjunction_v<std::is_integral<std::decay_t<Ns>>...> and
                sizeof...(Ns) >= 1> * = nullptr>
  inline decltype(auto) get(const Element & element, Ns &&... ns) const;

  /* get a reference to the array of certain type
   * @param type data filed under type is returned
   * @param ghost_type optional: by default the non-ghost map is searched
   * @return a const reference to the array */
  inline auto operator()(SupportType type, GhostType ghost_type = _not_ghost)
      -> Array<T> &;

  /*! insert data of a new type (not yet present) into the map.
   *  @param type type of data (if this type is already present in the map,
   *         an exception is thrown).
   *  @param ghost_type optional: by default, the data map for non-ghost
   *         elements is searched
   *  @param vect the vector to include into the map
   *  @return stored data corresponding to type. */
  inline void setArray(SupportType type, GhostType ghost_type, Array<T> & vect);
  /*! frees all memory related to the data*/
  inline void free();

  inline void clear();

  [[nodiscard]] inline bool empty() const;

  /*! set all values in the ElementTypeMap to zero*/
  inline void zero() { this->set(T()); }

  /*! set all values in the ElementTypeMap to value */
  template <typename ST> inline void set(const ST & value);

  /*! deletes and reorders entries in the stored arrays
   *  @param new_numbering a ElementTypeMapArray of new indices. UInt(-1)
   * indicates
   *         deleted entries. */
  inline void onElementsRemoved(const ElementTypeMapArray<Int> & new_numbering);

  /// text output helper
  void printself(std::ostream & stream, int indent = 0) const override;

  /*! set the id
   *  @param id the new name
   */
  inline void setID(const ID & id) { this->id = id; }

  /// return the id
  [[nodiscard]] [[nodiscard]] inline auto getID() const -> ID {
    return this->id;
  }

  [[nodiscard]] auto getNbComponents(
      Int dim = _all_dimensions, GhostType requested_ghost_type = _not_ghost,
      ElementKind kind = _ek_not_defined) const -> ElementTypeMap<Int> {
    ElementTypeMap<Int> nb_components;
    auto all_ghost_types = requested_ghost_type == _casper;
    for (auto ghost_type : ghost_types) {
      if ((not(ghost_type == requested_ghost_type)) and (not all_ghost_types)) {
        continue;
      }

      for (const auto & type : this->elementTypes(dim, ghost_type, kind)) {
        auto nb_comp = (*this)(type, ghost_type).getNbComponent();
        nb_components(type, ghost_type) = nb_comp;
      }
    }
    return nb_components;
  }

  /* ------------------------------------------------------------------------ */
  /* more evolved allocators                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the arrays in accordance to a functor
  template <class Func>
  void initialize(const Func & f, const T & default_value, bool do_not_default);

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

  /**
   * get the size of the ElementTypeMapArray<T>
   * @param[in] _pack
   * \parblock
   * optional arguments can be any of:
   * \li \c _spatial_dimension the dimension to consider (default:
   * _all_dimensions)
   * \li \c _ghost_type  (default: _not_ghost)
   * \li \c _element_kind (default: _ek_not_defined)
   * \li \c _all_ghost_types (default: false)
   * \endparblock
   **/
  template <typename... pack> auto size(pack &&... _pack) const -> Int;

  [[nodiscard]] auto isNodal() const -> bool { return is_nodal; }
  void isNodal(bool is_nodal) { this->is_nodal = is_nodal; }

private:
  [[nodiscard]] auto sizeImpl(Int spatial_dimension, GhostType ghost_type,
                              ElementKind kind) const -> Int;

private:
  ID id;

protected:
  /// name of the element type map: e.g. connectivity, grad_u
  ID name;

  /// Is the data stored by node of the element
  bool is_nodal{false};
};

/// to store data Array<Real> by element type
using ElementTypeMapReal = ElementTypeMapArray<Real>;
/// to store data Array<Int> by element type
using ElementTypeMapInt = ElementTypeMapArray<Int>;
/// to store data Array<UInt> by element type
using ElementTypeMapUInt = ElementTypeMapArray<UInt>;
/// to store data Array<Idx> by element type
using ElementTypeMapIdx = ElementTypeMapArray<Idx>;

} // namespace akantu

// namespace std {
// template <class Stored, typename SupportType>
// struct iterator_traits<
//     ::akantu::template ElementTypeMap<Stored, SupportType>::type_iterator> {
// private:
//   using iterator_type =
//       typename ::akantu::ElementTypeMap<Stored, SupportType>::type_iterator;

// public:
//   using iterator_category = typename iterator_type::iterator_category;
//   using value_type = typename iterator_type::value_type;
//   using difference_type = typename iterator_type::difference_type;
//   using pointer = typename iterator_type::pointer;
//   using reference = typename iterator_type::reference;
// };
// } // namespace std

#endif /* AKANTU_ELEMENT_TYPE_MAP_HH_ */
