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
#include "dumper_functors.hh"
#include "dumper_internal_types.hh"
/* -------------------------------------------------------------------------- */
#include <sstream>
#include <typeindex>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_DUMPER_FIELD_HH_
#define AKANTU_DUMPER_FIELD_HH_

namespace akantu {
namespace dumper {

  /// Field interface
  class FieldBase : public PropertiesManager,
                    public std::enable_shared_from_this<FieldBase> {
  public:
    FieldBase(const SupportBase & support, std::type_index type,
              FieldType field_type)
        : support(support), type_(type), field_type(field_type) {}

    [[nodiscard]] std::type_index type() const { return type_; }
    AKANTU_GET_MACRO(FieldType, field_type, FieldType);
    AKANTU_GET_MACRO(Support, support, const SupportBase &);

    [[nodiscard]] inline bool is_visualization() const {
      return this->is_usage_field(FieldUsageType::_visualisation);
    }
    [[nodiscard]] inline bool is_restart() const {
      return this->is_usage_field(FieldUsageType::_restart);
    }
    [[nodiscard]] inline bool is_internal() const {
      return this->is_usage_field(FieldUsageType::_internal);
    }

    [[nodiscard]] inline bool is_nodal_field() const {
      return (field_type == FieldType::_node_array_function or
              field_type == FieldType::_node_array);
    }

    [[nodiscard]] inline bool is_quadrature_points_field() const {
      return (field_type == FieldType::_internal_field_function or
              field_type == FieldType::_internal_field);
    }

    [[nodiscard]] inline bool is_elemental_field() const {
      return (field_type == FieldType::_element_map_array_function or
              field_type == FieldType::_element_map_array) or
             this->is_quadrature_points_field();
    }

    virtual void update() {}

  protected:
    void setFieldType(FieldType type) { field_type = type; }
    inline bool is_usage_field(FieldUsageType usage_type) const {
      if (not this->hasProperty("field_usage")) {
        return false;
      }

      FieldUsageType usage = this->getProperty<FieldUsageType>("field_usage");
      return (usage & usage_type);
    }

  private:
    const SupportBase & support;
    std::type_index type_;
    FieldType field_type{FieldType::_not_defined};
  };

  /* ------------------------------------------------------------------------ */
  class FieldArrayBase : public FieldBase {
  public:
    FieldArrayBase(const SupportBase & support, std::type_index type,
                   FieldType field_type = FieldType::_node_array)
        : FieldBase(support, type, field_type) {}

    [[nodiscard]] virtual const void * data() const = 0;
    [[nodiscard]] virtual Int size() const = 0;
    [[nodiscard]] virtual Int getNbComponent() const = 0;
  };

  /* ------------------------------------------------------------------------ */
  template <class T> class FieldArrayTemplateBase : public FieldArrayBase {
  public:
    using FieldArrayBase::FieldArrayBase;
    [[nodiscard]] virtual const Array<T> & getArray() const = 0;
    [[nodiscard]] virtual const Array<T> & getArray() = 0;

    auto getSharedPointer() {
      return std::dynamic_pointer_cast<FieldArrayTemplateBase<T>>(
          this->shared_from_this());
    }
  };

  /* ------------------------------------------------------------------------ */
  template <typename T, class Array_ = const Array<T> &>
  class FieldArray : public FieldArrayTemplateBase<T> {
  public:
    FieldArray(Array_ && array, const SupportBase & support)
        : FieldArrayTemplateBase<T>(support, typeid(T)),
          array(std::forward<Array_>(array)) {}

    [[nodiscard]] const void * data() const override { return array.data(); }

    [[nodiscard]] Int getNbComponent() const override {
      // auto factor = array.size() / this->size();
      return array.getNbComponent();
    }

    [[nodiscard]] const Array<T> & getArray() const override { return array; }
    [[nodiscard]] const Array<T> & getArray() override { return array; }

  private:
    Array_ array;
  };

  /* ------------------------------------------------------------------------ */
  template <class T, class Function>
  class FieldFunctionArray
      : public FieldArrayTemplateBase<
            details::function_return_scalar_t<T, Function>> {
    using OutT = details::function_return_scalar_t<T, Function>;

  public:
    FieldFunctionArray(
        const std::shared_ptr<FieldArrayTemplateBase<T>> & array_in,
        const SupportBase & support,
        Function && function) // NOLINT
        : FieldArrayTemplateBase<OutT>(support, typeid(T),
                                       FieldType::_node_array_function),
          function(std::forward<Function>(function)), array_in(array_in) {
      // update();
    }

    FieldFunctionArray(const Array<T> & array_in, const SupportBase & support,
                       Function && function) // NOLINT
        : FieldFunctionArray(make_field(array_in, support), support,
                             std::forward<Function>(function)) {}

    [[nodiscard]] const void * data() const override {
      return array_out->data();
    }

    template <class OutT = T>
    [[nodiscard]] std::enable_if<not std::is_const_v<OutT>> * data() {
      update();

      return array_out->data();
    }

    [[nodiscard]] Int getNbComponent() const override {
      // auto old_size = this->array_in->size();
      // auto new_size = this->size();
      // auto factor = old_size / new_size;
      if constexpr (details::has_getNbComponent_member<Function>) {
        return function.getNbComponent(this->array_in->getNbComponent());
      } else {
        return this->array_in->getNbComponent();
      }
    }

    [[nodiscard]] const Array<OutT> & getArray() const override {
      return *array_out;
    }
    [[nodiscard]] const Array<OutT> & getArray() override {
      update();
      return *array_out;
    }

    [[nodiscard]] Release getRelease() const override { return last_release; }

    void update() override {
      if (unchanged() and array_out) {
        return;
      }

      if (not array_out) {
        array_out =
            std::make_unique<Array<OutT>>(this->size(), this->getNbComponent());
      } else {
        array_out->resize(this->size());
      }

      auto nb_component_in = array_in->getNbComponent();
      auto nb_component_out = array_out->getNbComponent();
      for (auto && [out, in] :
           zip(make_view(*array_out, nb_component_out),
               make_view(array_in->getArray(), nb_component_in))) {
        auto && res = function(in);
        out = res;
      }

      last_release = getRelease();
    }

  private:
    bool unchanged() { return getRelease() == last_release; }

  private:
    Function function;
    std::shared_ptr<FieldArrayTemplateBase<T>> array_in;
    std::unique_ptr<Array<OutT>> array_out;
    Release last_release;
  };

  /* ------------------------------------------------------------------------ */
  template <class T, class Array_ = const Array<T> &>
  class FieldNodeArray : public FieldArray<T, Array_> {
  public:
    using FieldArray<T, Array_>::FieldArray;

    [[nodiscard]] Int size() const override {
      return aka::as_type<SupportElements>(this->getSupport()).getNbNodes();
    }
  };

  /* ------------------------------------------------------------------------ */
  template <class T, class Array_ = const Array<T> &>
  class FieldElementsArray : public FieldArray<T, Array_> {
  public:
    FieldElementsArray(Array_ && array, const SupportBase & support,
                       const ElementType & type) // NOLINT
        : FieldArray<T, Array_>(std::forward<Array_>(array), support),
          type(type) {}

    [[nodiscard]] Int size() const override {
      return aka::as_type<SupportElements>(this->getSupport())
          .getNbElements(type);
    }

  protected:
    ElementType type;
  };

  /* ------------------------------------------------------------------------ */
  template <class T, class Function>
  class FieldFunctionNodeArray : public FieldFunctionArray<T, Function> {
  public:
    FieldFunctionNodeArray(
        const std::shared_ptr<FieldArrayTemplateBase<T>> & array_in,
        const SupportBase & support, Function && function) // NOLINT
        : FieldFunctionArray<T, Function>(array_in, support,
                                          std::forward<Function>(function)) {
      this->update();
    }

    FieldFunctionNodeArray(const Array<T> & array_in,
                           const SupportBase & support, Function && function)
        : FieldFunctionArray<T, Function>(array_in, support,
                                          std::forward<Function>(function)) {
      this->update();
    }

    [[nodiscard]] Int size() const override {
      return aka::as_type<SupportElements>(this->getSupport()).getNbNodes();
    }
  };

  /* ------------------------------------------------------------------------ */
  template <class T, class Function>
  class FieldFunctionElementsArray
      : public FieldFunctionArray<T, ElementTypeMapArrayFunctor<Function>> {
  public:
    FieldFunctionElementsArray(
        const std::shared_ptr<FieldArrayTemplateBase<T>> & array_in,
        const SupportBase & support, Function && function, // NOLINT
        const ElementType & type)
        : FieldFunctionArray<T, ElementTypeMapArrayFunctor<Function>>(
              array_in, support,
              ElementTypeMapArrayFunctor<Function>(
                  type, std::forward<Function>(function))),
          type(type) {
      this->update();
    }

    FieldFunctionElementsArray(const Array<T> & array_in,
                               const SupportBase & support,
                               Function && function, // NOLINT
                               const ElementType & type)
        : FieldFunctionArray<T, ElementTypeMapArrayFunctor<Function>>(
              array_in, support,
              ElementTypeMapArrayFunctor<Function>(
                  type, std::forward<Function>(function))),
          type(type) {
      this->update();
    }

    [[nodiscard]] Int size() const override {
      return aka::as_type<SupportElements>(this->getSupport())
          .getNbElements(type);
    }

  protected:
    ElementType type;
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

    [[nodiscard]] virtual const FieldArrayBase &
    array(ElementType type, GhostType ghost_type = _not_ghost) const {
      return *fields.at({type, ghost_type});
    }

    [[nodiscard]] virtual FieldArrayBase &
    array(ElementType type, GhostType ghost_type = _not_ghost) {
      return *fields.at({type, ghost_type});
    }

    [[nodiscard]] virtual ElementTypesIteratorHelper elementTypes() const {
      return this->support_element.elementTypes();
    };

    [[nodiscard]] virtual std::pair<Int, Int> size() const = 0;
    [[nodiscard]] virtual Int
    getNbComponent(const ElementType & type) const = 0;

    void update() override {
      for (auto && [_, field] : fields) {
        field->update();
      }
    }

  protected:
    virtual void createByTypesFields() {}

  protected:
    const SupportElements & support_element;
    std::map<std::pair<ElementType, GhostType>, std::shared_ptr<FieldArrayBase>>
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
      return (aka::as_type<FieldArrayTemplateBase<T>>(
          this->array(type, ghost_type)));
    }

    [[nodiscard]] decltype(auto)
    arrayTyped(ElementType type, GhostType ghost_type = _not_ghost) const {
      return (aka::as_type<FieldArrayTemplateBase<T>>(
          this->array(type, ghost_type)));
    }

    [[nodiscard]] std::pair<Int, Int> size() const override {
      Int total_element{0};
      Int total_size{0};

      for (auto && type : this->elementTypes()) {
        auto && array = this->array(type);
        total_element += array.size();
        total_size += Int(array.size()) * array.getNbComponent();
      }
      return std::make_pair(total_element, total_size);
    }

    [[nodiscard]] Int getNbComponent(const ElementType & type) const override {
      return this->array(type).getNbComponent();
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
      createByTypesFields();
    }

    [[nodiscard]] Release getRelease() const override {
      return map_array.getRelease();
    }

  private:
    void createByTypesFields() override {
      for (auto ghost_type : ghost_types) {
        for (auto && type : this->support_element.elementTypes(ghost_type)) {
          auto name = std::to_string(type);
          if (ghost_type != _not_ghost) {
            name += ":" + std::to_string(ghost_type);
          }

          auto field = std::make_shared<FieldElementsArray<T>>(
              map_array(type, ghost_type), this->getSupport(), type);
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
            details::function_with_type_return_scalar_t<T, Function>> {
    using OutT = details::function_with_type_return_scalar_t<T, Function>;

  public:
    FieldFunctionElementMapArray(
        const std::shared_ptr<FieldElementMapArrayTemplateBase<T>> &
            map_array_in,
        const SupportBase & support,
        Function && function, // NOLINT
        FieldType field_type = FieldType::_element_map_array_function,
        bool create_fields = true)
        : FieldElementMapArrayTemplateBase<OutT>(support, typeid(T),
                                                 field_type),
          function(std::forward<Function>(function)),
          map_array_in(map_array_in) {
      if (create_fields) {
        createByTypesFields();
      }
    }

    template <class ElementTypeMapArray_ = const ElementTypeMapArray<T> &,
              std::enable_if_t<details::is_element_type_map_array_v<
                  std::decay_t<ElementTypeMapArray_>>> * = nullptr>
    FieldFunctionElementMapArray(
        ElementTypeMapArray_ && map_array_in, const SupportBase & support,
        Function && function, // NOLINT
        FieldType field_type = FieldType::_element_map_array_function,
        bool create_fields = true)
        : FieldFunctionElementMapArray(
              make_field(std::forward<ElementTypeMapArray_>(map_array_in),
                         support),
              support, std::forward<Function>(function), field_type,
              create_fields) {}

    [[nodiscard]] decltype(auto) arrayTyped(ElementType type,
                                            GhostType ghost_type) const {
      return (aka::as_type<FieldFunctionNodeArray<T, decltype(function)>>(
          this->array(type, ghost_type)));
    }

    [[nodiscard]] Release getRelease() const override {
      return map_array_in->getRelease();
    }

  protected:
    void createByTypesFields() override {
      for (auto ghost_type : ghost_types) {
        for (auto && type : aka::as_type<SupportElements>(this->support_element)
                                .elementTypes(ghost_type)) {
          auto name = std::to_string(type);
          if (ghost_type != _not_ghost) {
            name += ":" + std::to_string(ghost_type);
          }

          auto field =
              std::make_shared<FieldFunctionElementsArray<T, const Function &>>(
                  map_array_in->arrayTyped(type, ghost_type).getSharedPointer(),
                  this->getSupport(), this->function, type);
          field->addProperty("name", name);
          field->removeProperty("data_location");
          this->fields.emplace(std::pair(type, ghost_type), std::move(field));
        }
      }
    }

  protected:
    Function function;
    std::shared_ptr<FieldElementMapArrayTemplateBase<T>> map_array_in;
  };

  /* ------------------------------------------------------------------------ */
  using FieldConnectivity =
      FieldFunctionElementMapArray<Idx, toVTKConnectivity>;

  /* ------------------------------------------------------------------------ */
  template <class T> class FieldInternalField : public FieldElementMapArray<T> {
  public:
    FieldInternalField(InternalField<T> & map_array_in,
                       const SupportBase & support) // NOLINT
        : FieldElementMapArray<T>(map_array_in, support,
                                  FieldType::_internal_field) {}
  };

  /* ------------------------------------------------------------------------ */

  template <class T, class Function>
  class FieldFunctionInternalField
      : public FieldFunctionElementMapArray<T, Function> {
    using parent = FieldFunctionElementMapArray<T, Function>;

  public:
    FieldFunctionInternalField(InternalField<T> & map_array_in,
                               const SupportBase & support,
                               Function && function) // NOLINT
        : parent(map_array_in, support, std::forward<Function>(function),
                 FieldType::_internal_field_function, false) {
      createByTypesFields();
    }

    [[nodiscard]] virtual FieldArrayBase &
    array(ElementType type, GhostType ghost_type = _not_ghost) {
      update_nb_integration_points();
      return *this->fields.at({type, ghost_type});
    }

    void update() override {
      this->update_nb_integration_points();
      parent::update();
    }

  private:
    bool unchanged() { return this->getRelease() == last_release; }

    void createByTypesFields() override {
      update_nb_integration_points();
      parent::createByTypesFields();
    }

    void update_nb_integration_points() {
      if (unchanged()) {
        return;
      }

      if constexpr (details::has_set_nb_integration_points_member<Function>) {
        ElementTypeMap<Int> nb_integration_points_per_elem;
        auto && support = aka::as_type<SupportElements>(this->getSupport());
        auto && connectivities = support.getConnectivities();
        for (auto type : this->map_array_in->elementTypes()) {
          auto nb_elements = connectivities.array(type).size();
          auto nb_integration_points = this->map_array_in->array(type).size();

          nb_integration_points_per_elem(type) =
              nb_integration_points / nb_elements;
        }
        this->function.setNbIntegtrationPoints(nb_integration_points_per_elem);
      }
    }

    Release last_release;
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
  auto make_field(const std::shared_ptr<FieldArrayTemplateBase<T>> & array,
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
