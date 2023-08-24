/**
 * @file   dumper_field.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Tue Jan 19 2016
 *
 * @brief  Common interface for fields
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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
#include "aka_array.hh"
#include "element_type_map.hh"
#include "internal_field.hh"
#include "support.hh"
/* -------------------------------------------------------------------------- */
#include <sstream>
#include <typeindex>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_DUMPER_FIELD_HH_
#define AKANTU_DUMPER_FIELD_HH_

namespace akantu {

namespace dumper {

  namespace details {
    // primary template handles types that do not support pre-increment:
    template <class, class = void> constexpr bool has_getNbComponent_member{};

    // specialization recognizes types that do support pre-increment:
    template <class T>
    constexpr bool has_getNbComponent_member<
        T, std::void_t<decltype(std::declval<T &>().getNbComponent(1))>> = true;

    template <class T>
    constexpr bool has_getNbComponent_member<
        T, std::void_t<decltype(std::declval<T &>().getNbComponent(
               1, _not_defined))>> = true;

    template <class T, class Function> struct function_with_type_scalar {
      using type = typename decltype(std::declval<Function &>().operator()(
          VectorProxy<T>(nullptr, 1), _not_defined))::Scalar;
    };

    template <class T, class Function> struct function_scalar {
      using type = typename decltype(std::declval<Function &>().operator()(
          VectorProxy<T>(nullptr, 1)))::Scalar;
    };

    template <class T, class Function>
    using function_with_type_scalar_t =
        typename function_with_type_scalar<T, Function>::type;

    template <class T, class Function>
    using function_scalar_t = typename function_scalar<T, Function>::type;

    template <class T> struct is_array : public std::false_type {};
    template <class T> struct is_array<Array<T>> : public std::true_type {};
    template <class T>
    struct is_array<const Array<T>> : public std::true_type {};

    template <class T> constexpr auto is_array_v = is_array<T>::value;

    template <class T>
    struct is_element_type_map_array : public std::false_type {};
    template <class T>
    struct is_element_type_map_array<ElementTypeMapArray<T>>
        : public std::true_type {};
    template <class T>
    struct is_element_type_map_array<const ElementTypeMapArray<T>>
        : public std::true_type {};

    template <class T>
    constexpr auto is_element_type_map_array_v =
        is_element_type_map_array<T>::value;

  } // namespace details

  /// Field interface
  class FieldBase : public PropertiesManager,
                    public std::enable_shared_from_this<FieldBase> {
  public:
    FieldBase(const SupportBase & support, std::type_index type,
              FieldType field_type)
        : support(support), type_(type), field_type(field_type) {}

    virtual ~FieldBase() = default;

    [[nodiscard]] std::type_index type() const { return type_; }
    AKANTU_GET_MACRO(FieldType, field_type, FieldType);
    AKANTU_GET_MACRO(Support, support, const SupportBase &);

  protected:
    void setFieldType(FieldType type) { field_type = type; }

  private:
    const SupportBase & support;
    std::type_index type_;
    FieldType field_type{FieldType::_not_defined};
  };

  inline bool is_nodal_field(const FieldBase & field) {
    auto type = field.getFieldType();
    return (type == FieldType::_node_array_function or
            type == FieldType::_node_array);
  }

  inline bool is_quadrature_points_field(const FieldBase & field) {
    auto type = field.getFieldType();
    return (type == FieldType::_internal_field_function or
            type == FieldType::_internal_field);
  }

  inline bool is_elemental_field(const FieldBase & field) {
    auto type = field.getFieldType();
    return (type == FieldType::_element_map_array_function or
            type == FieldType::_element_map_array) or
           is_quadrature_points_field(field);
  }

  /* ------------------------------------------------------------------------ */
  class FieldNodeArrayBase : public FieldBase {
  public:
    FieldNodeArrayBase(const SupportBase & support, std::type_index type,
                       FieldType field_type = FieldType::_node_array)
        : FieldBase(support, type, field_type) {}

    [[nodiscard]] virtual const void * data() const = 0;
    [[nodiscard]] virtual Int size() const = 0;
    [[nodiscard]] virtual Int getNbComponent() const = 0;
  };

  /* ------------------------------------------------------------------------ */
  template <class T>
  class FieldNodeArrayTemplateBase : public FieldNodeArrayBase {
  public:
    using FieldNodeArrayBase::FieldNodeArrayBase;
    [[nodiscard]] virtual const Array<T> & getArray() const = 0;
    [[nodiscard]] virtual const Array<T> & getArray() = 0;

    auto getSharedPointer() {
      return std::dynamic_pointer_cast<FieldNodeArrayTemplateBase<T>>(
          this->shared_from_this());
    }
  };

  /* ------------------------------------------------------------------------ */
  template <typename T, class Array_ = const Array<T> &>
  class FieldNodeArray : public FieldNodeArrayTemplateBase<T> {
  public:
    FieldNodeArray(Array_ && array, const SupportBase & support)
        : FieldNodeArrayTemplateBase<T>(support, typeid(T)),
          array(std::forward<Array_>(array)) {}

    [[nodiscard]] const void * data() const override { return array.data(); }

    [[nodiscard]] Int size() const override { return array.size(); }
    [[nodiscard]] Int getNbComponent() const override {
      return array.getNbComponent();
    }

    [[nodiscard]] const Array<T> & getArray() const override { return array; }
    [[nodiscard]] const Array<T> & getArray() override { return array; }

  private:
    Array_ array;
  };

  /* ------------------------------------------------------------------------ */
  template <class T, class Function>
  class FieldFunctionNodeArray : public FieldNodeArrayTemplateBase<
                                     details::function_scalar_t<T, Function>> {
    using OutT = details::function_scalar_t<T, Function>;

  public:
    FieldFunctionNodeArray(
        const std::shared_ptr<FieldNodeArrayTemplateBase<T>> & array_in,
        const SupportBase & support,
        Function && function) // NOLINT
        : FieldNodeArrayTemplateBase<OutT>(support, typeid(T),
                                           FieldType::_node_array_function),
          function(std::forward<Function>(function)), array_in(array_in) {
      update();
    }

    FieldFunctionNodeArray(const Array<T> & array_in,
                           const SupportBase & support,
                           Function && function) // NOLINT
        : FieldFunctionNodeArray(make_field(array_in, support), support,
                                 std::forward<Function>(function)) {}

    [[nodiscard]] const void * data() const override {
      return array_out->data();
    }

    template <class OutT = T>
    [[nodiscard]] std::enable_if<not std::is_const_v<OutT>> * data() {
      update();
      return array_out->data();
    }

    [[nodiscard]] Int size() const override { return array_in->size(); }
    [[nodiscard]] Int getNbComponent() const override {
      if constexpr (details::has_getNbComponent_member<Function>) {
        return function.getNbComponent(array_in->getNbComponent());
      } else {
        return array_in->getNbComponent();
      }
    }

    [[nodiscard]] const Array<OutT> & getArray() const override {
      return *array_out;
    }
    [[nodiscard]] const Array<OutT> & getArray() override {
      update();
      return *array_out;
    }

  private:
    void update() {
      if (not array_out) {
        array_out = std::make_unique<Array<OutT>>(array_in->size(),
                                                  this->getNbComponent());
      } else {
        array_out->resize(array_in->size());
      }

      for (auto && [out, in] :
           zip(make_view(*array_out, array_out->getNbComponent()),
               make_view(array_in->getArray(), array_in->getNbComponent()))) {
        auto && res = function(in);
        out = res;
      }
    }

  private:
    Function function;
    std::shared_ptr<FieldNodeArrayTemplateBase<T>> array_in;
    std::unique_ptr<Array<OutT>> array_out;
  };

  /* ------------------------------------------------------------------------ */
  class FieldElementMapArrayBase : public FieldBase {
  public:
    using ElementTypesIteratorHelper =
        ElementTypeMapArray<Idx, ElementType>::ElementTypesIteratorHelper;

    explicit FieldElementMapArrayBase(
        const SupportBase & support, std::type_index type,
        FieldType field_type = FieldType::_element_map_array)
        : FieldBase(support, type, field_type),
          support_element(dynamic_cast<const SupportElements &>(support)) {}

    [[nodiscard]] virtual const FieldNodeArrayBase &
    array(ElementType type, GhostType ghost_type = _not_ghost) const {
      return *fields.at({type, ghost_type});
    }

    [[nodiscard]] virtual FieldNodeArrayBase &
    array(ElementType type, GhostType ghost_type = _not_ghost) {
      return *fields.at({type, ghost_type});
    }

    [[nodiscard]] virtual ElementTypesIteratorHelper elementTypes() const {
      return this->support_element.elementTypes();
    };

    [[nodiscard]] Int getNbComponent(const ElementType & type) const {
      return this->array(type).getNbComponent();
    }

    [[nodiscard]] auto size() const {
      Int total_element{0};
      Int total_size{0};

      for (auto && type : this->elementTypes()) {
        auto && array = this->array(type);
        total_element += Int(array.size());
        total_size += Int(array.size()) * array.getNbComponent();
      }
      return std::make_pair(total_element, total_size);
    }

  protected:
    const SupportElements & support_element;
    std::map<std::pair<ElementType, GhostType>,
             std::shared_ptr<FieldNodeArrayBase>>
        fields;
  };

  /* ------------------------------------------------------------------------ */
  template <typename T>
  class FieldElementMapArrayTemplateBase : public FieldElementMapArrayBase {
  public:
    using FieldElementMapArrayBase::FieldElementMapArrayBase;

    auto getSharedPointer() {
      return std::dynamic_pointer_cast<FieldElementMapArrayTemplateBase<T>>(
          this->shared_from_this());
    }

    [[nodiscard]] decltype(auto) arrayTyped(ElementType type,
                                            GhostType ghost_type = _not_ghost) {
      return (aka::as_type<FieldNodeArrayTemplateBase<T>>(
          *this->fields.at({type, ghost_type})));
    }

    [[nodiscard]] decltype(auto)
    arrayTyped(ElementType type, GhostType ghost_type = _not_ghost) const {
      return (aka::as_type<FieldNodeArrayTemplateBase<T>>(
          *this->fields.at({type, ghost_type})));
    }
  };

  /* ------------------------------------------------------------------------ */
  template <typename T,
            class ElementTypeMapArray_ = const ElementTypeMapArray<T> &,
            std::enable_if_t<details::is_element_type_map_array_v<
                std::decay_t<ElementTypeMapArray_>>> * = nullptr>
  class FieldElementMapArray : public FieldElementMapArrayTemplateBase<T> {
  public:
    FieldElementMapArray(ElementTypeMapArray_ && map_array,
                         const SupportBase & support,
                         FieldType field_type = FieldType::_element_map_array)
        : FieldElementMapArrayTemplateBase<T>(support, typeid(T), field_type),
          map_array(std::forward<ElementTypeMapArray_>(map_array)) {
      createByTypesNodalFields();
    }

  private:
    void createByTypesNodalFields() {
      for (auto ghost_type : ghost_types) {
        for (auto && type : this->support_element.elementTypes(ghost_type)) {
          auto name = std::to_string(type);
          if (ghost_type != _not_ghost) {
            name += ":" + std::to_string(ghost_type);
          }

          auto field = std::make_shared<FieldNodeArray<T>>(
              map_array(type, ghost_type), this->getSupport());
          field->addProperty("name", name);
          this->fields.emplace(std::pair(type, ghost_type), std::move(field));
        }
      }
    }

  private:
    ElementTypeMapArray_ map_array;
  };

  /* ------------------------------------------------------------------------ */
  template <class T, class Function>
  class FieldFunctionElementMapArray
      : public FieldElementMapArrayTemplateBase<
            details::function_with_type_scalar_t<T, Function>> {
    using OutT = details::function_with_type_scalar_t<T, Function>;

    template <class ElementFunction> struct ElementTypeMapArrayFunctor {
      ElementTypeMapArrayFunctor(ElementType type,
                                 const ElementFunction & function)
          : type(type), function(function) {}

      [[nodiscard]] Int getNbComponent(Int nb_component) const {
        if constexpr (details::has_getNbComponent_member<Function>) {
          return function.getNbComponent(nb_component, type);
        } else {
          return nb_component;
        }
      }

      template <class Derived>
      [[nodiscard]] decltype(auto)
      operator()(const Eigen::MatrixBase<Derived> & in) const {
        return function(in, type);
      }

      ElementType type;
      const ElementFunction & function;
    };

  public:
    FieldFunctionElementMapArray(
        const std::shared_ptr<FieldElementMapArrayTemplateBase<T>> &
            map_array_in,
        const SupportBase & support,
        Function && function, // NOLINT
        FieldType field_type = FieldType::_element_map_array_function)
        : FieldElementMapArrayTemplateBase<OutT>(support, typeid(T),
                                                 field_type),
          function(std::forward<Function>(function)),
          map_array_in(map_array_in) {
      createByTypesNodalFields();
    }

    template <class ElementTypeMapArray_ = const ElementTypeMapArray<T> &,
              std::enable_if_t<details::is_element_type_map_array_v<
                  std::decay_t<ElementTypeMapArray_>>> * = nullptr>
    FieldFunctionElementMapArray(
        ElementTypeMapArray_ && map_array_in, const SupportBase & support,
        Function && function, // NOLINT
        FieldType field_type = FieldType::_element_map_array_function)
        : FieldFunctionElementMapArray(
              make_field(std::forward<ElementTypeMapArray_>(map_array_in),
                         support),
              support, std::forward<Function>(function), field_type) {}

    [[nodiscard]] decltype(auto) arrayTyped(ElementType type,
                                            GhostType ghost_type) const {
      return (aka::as_type<FieldFunctionNodeArray<T, decltype(function)>>(
          *this->fields.at({type, ghost_type})));
    }

  private:
    void createByTypesNodalFields() {
      for (auto ghost_type : ghost_types) {
        for (auto && type : aka::as_type<SupportElements>(this->support_element)
                                .elementTypes(ghost_type)) {
          auto name = std::to_string(type);
          if (ghost_type != _not_ghost) {
            name += ":" + std::to_string(ghost_type);
          }

          auto field = make_field(
              map_array_in->arrayTyped(type, ghost_type).getSharedPointer(),
              this->getSupport(),
              ElementTypeMapArrayFunctor(type, this->function));
          field->addProperty("name", name);
          field->removeProperty("data_location");
          this->fields.emplace(std::pair(type, ghost_type), std::move(field));
        }
      }
    }

  private:
    Function function;
    std::shared_ptr<FieldElementMapArrayTemplateBase<T>> map_array_in;
  };

  struct toVTKConnectivity {
    template <class Derived>
    Vector<Idx> operator()(const Eigen::MatrixBase<Derived> & connectivity,
                           ElementType /*type*/) const {
      return connectivity;
    }
  };

  struct ElementGroupConnectivityFunctor {
    ElementGroupConnectivityFunctor(
        const ElementTypeMapArray<Idx> & connectivities)
        : connectivities(connectivities) {}

    [[nodiscard]] Int getNbComponent(Int /*nb_component*/,
                                     ElementType type) const {
      return connectivities(type).getNbComponent();
    }

    template <class Derived>
    Vector<Idx> operator()(const Eigen::MatrixBase<Derived> & element,
                           ElementType type) const {
      Element el{type, element[0], _not_ghost};
      return rewrite(connectivities.get(el), type);
    }

  private:
    toVTKConnectivity rewrite;
    const ElementTypeMapArray<Idx> & connectivities;
  };

  struct ElementGroupNodesFunctor {
    ElementGroupNodesFunctor(const Array<Real> & nodes) : nodes(nodes) {}

    [[nodiscard]] Int getNbComponent(Int /*nb_component*/
    ) const {
      return nodes.getNbComponent();
    }

    template <class Derived>
    Vector<Real> operator()(const Eigen::MatrixBase<Derived> & node) const {
      return make_view(nodes, nodes.getNbComponent()).begin()[node[0]];
    }

  private:
    toVTKConnectivity rewrite;
    const Array<Real> & nodes;
  };

  using FieldConnectivity =
      FieldFunctionElementMapArray<Idx, toVTKConnectivity>;

  template <class T> class FieldInternalField : public FieldElementMapArray<T> {
  public:
    FieldInternalField(InternalField<T> & map_array_in,
                       const SupportBase & support) // NOLINT
        : FieldElementMapArray<T>(map_array_in, support,
                                  FieldType::_internal_field) {}
  };

  template <class T, class Function>
  class FieldFunctionInternalField
      : public FieldFunctionElementMapArray<T, Function> {
  public:
    FieldFunctionInternalField(InternalField<T> & map_array_in,
                               const SupportBase & support,
                               Function && function) // NOLINT
        : FieldFunctionElementMapArray<T, Function>(
              map_array_in, support, std::forward<Function>(function),
              FieldType::_internal_field_function) {}
  };

  /* ------------------------------------------------------------------------ */
  template <class Array_, std::enable_if_t<dumper::details::is_array_v<
                              std::decay_t<Array_>>> * = nullptr>
  auto make_field(Array_ && array, const SupportBase & support) {
    using T = typename std::decay_t<Array_>::value_type;
    return std::make_shared<FieldNodeArray<T, Array_>>(
        std::forward<Array_>(array), support);
  }

  /* ------------------------------------------------------------------------ */
  template <class Function, class Array_,
            std::enable_if_t<
                dumper::details::is_array_v<std::decay_t<Array_>>> * = nullptr>
  auto make_field(Array_ && array, const SupportBase & support,
                  Function && function) {
    using T = typename std::decay_t<Array_>::value_type;
    return std::make_shared<FieldFunctionNodeArray<T, Function>>(
        std::forward<Array_>(array), support, std::forward<Function>(function));
  }

  /* ------------------------------------------------------------------------ */
  template <typename T, class Function>
  auto make_field(const std::shared_ptr<FieldNodeArrayTemplateBase<T>> & array,
                  const SupportBase & support, Function && function) {
    return std::make_shared<FieldFunctionNodeArray<T, Function>>(
        array, support, std::forward<Function>(function));
  }

  /* ------------------------------------------------------------------------ */
  template <class ElementTypeMapArray_,
            std::enable_if_t<dumper::details::is_element_type_map_array_v<
                std::decay_t<ElementTypeMapArray_>>> * = nullptr>
  auto make_field(ElementTypeMapArray_ && array, const SupportBase & support) {
    using T = typename std::decay_t<ElementTypeMapArray_>::value_type;
    return std::make_shared<FieldElementMapArray<T, ElementTypeMapArray_>>(
        std::forward<ElementTypeMapArray_>(array), support);
  }

  /* ------------------------------------------------------------------------ */
  template <class Function, class ElementTypeMapArray_,
            std::enable_if_t<dumper::details::is_element_type_map_array_v<
                std::decay_t<ElementTypeMapArray_>>> * = nullptr>
  auto make_field(ElementTypeMapArray_ && array, const SupportBase & support,
                  Function && function) {
    using T = typename std::decay_t<ElementTypeMapArray_>::value_type;
    return std::make_shared<FieldFunctionElementMapArray<T, Function>>(
        std::forward<ElementTypeMapArray_>(array), support,
        std::forward<Function>(function));
  }

  /* ------------------------------------------------------------------------ */
  template <typename T, class Function>
  auto
  make_field(const std::shared_ptr<FieldElementMapArrayTemplateBase<T>> & array,
             const SupportBase & support, Function && function) {
    return std::make_shared<FieldFunctionElementMapArray<T, Function>>(
        array, support, std::forward<Function>(function));
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  auto make_field(InternalField<T> & array, const SupportBase & support) {
    return std::make_shared<FieldInternalField<T>>(array, support);
  }

  /* ------------------------------------------------------------------------ */
  template <typename T, class Function>
  auto make_field(InternalField<T> & array, const SupportBase & support,
                  Function && function) {
    return std::make_shared<FieldFunctionInternalField<T, Function>>(
        array, support, std::forward<Function>(function));
  }

} // namespace dumper
} // namespace akantu

#endif /* AKANTU_DUMPER_FIELD_HH_ */
