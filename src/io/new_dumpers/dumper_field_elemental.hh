/* -------------------------------------------------------------------------- */
#include "dumper_field_nodal.hh"
/* -------------------------------------------------------------------------- */

#ifndef DUMPER_FIELD_ELEMENTAL_HH_
#define DUMPER_FIELD_ELEMENTAL_HH_

namespace akantu {
namespace dumper {

  /* ------------------------------------------------------------------------ */
  class FieldElementalArrayBase : public FieldArrayBase {
  public:
    FieldElementalArrayBase(const SupportBase & support, std::type_index type,
                            const ElementType & element_type,
                            const GhostType & ghost_type)
        : FieldArrayBase(support, type, FieldType::_element_array),
          element_type(element_type), ghost_type(ghost_type) {}
    [[nodiscard]] Int size() const override {
      return aka::as_type<SupportElements>(this->getSupport())
          .getNbElements(element_type, ghost_type);
    }

    auto getElementType() const { return element_type; }
    auto getGhostType() const { return ghost_type; }

  protected:
    ElementType element_type;
    GhostType ghost_type;
  };

  /* ------------------------------------------------------------------------ */
  template <class T, class Array_ = const Array<T> &>
  class FieldElementalArray
      : public FieldArray<T, Array_, FieldElementalArrayBase> {
  public:
    FieldElementalArray(Array_ && array, const SupportBase & support,
                        const ElementType & element_type,
                        const GhostType & ghost_type) // NOLINT
        : FieldArray<T, Array_, FieldElementalArrayBase>(
              std::forward<Array_>(array), support, element_type, ghost_type) {}
  };

  /* ------------------------------------------------------------------------ */
  template <class T, class Function>
  class FieldFunctionElementalArray
      : public FieldFunctionArray<T, ElementTypeMapArrayFunctor<Function>,
                                  FieldElementalArrayBase> {
  public:
    FieldFunctionElementalArray(
        const std::shared_ptr<
            FieldArrayTemplateBase<T, FieldElementalArrayBase>> & array_in,
        const SupportBase & support, Function && function, // NOLINT
        const ElementType & type, const GhostType & ghost_type)
        : FieldFunctionArray<T, ElementTypeMapArrayFunctor<Function>,
                             FieldElementalArrayBase>(
              array_in, support,
              ElementTypeMapArrayFunctor<Function>(
                  type, ghost_type, std::forward<Function>(function)),
              type, ghost_type) {
      this->update();
    }

    FieldFunctionElementalArray(const Array<T> & array_in,
                                const SupportBase & support,
                                Function && function, // NOLINT
                                const ElementType & type,
                                const GhostType & ghost_type)
        : FieldFunctionArray<T, ElementTypeMapArrayFunctor<Function>>(
              array_in, support,
              ElementTypeMapArrayFunctor<Function>(
                  type, ghost_type, std::forward<Function>(function)),
              type, ghost_type) {
      this->update();
    }
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

    [[nodiscard]] virtual ElementTypesIteratorHelper
    elementTypes(const GhostType & ghost_type) const {
      return this->support_element.elementTypes(ghost_type);
    };

    [[nodiscard]] virtual std::pair<Int, Int>
    size(const GhostType & ghost_type) const = 0;
    [[nodiscard]] virtual Int
    getNbComponent(const ElementType & type,
                   const GhostType & ghost_type) const = 0;

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
      return (aka::as_type<FieldElementalArrayTemplateBase<T>>(
          this->array(type, ghost_type)));
    }

    [[nodiscard]] decltype(auto)
    arrayTyped(ElementType type, GhostType ghost_type = _not_ghost) const {
      return (aka::as_type<FieldElementalArrayTemplateBase<T>>(
          this->array(type, ghost_type)));
    }

    [[nodiscard]] std::pair<Int, Int>
    size(const GhostType & ghost_type) const override {
      Int total_element{0};
      Int total_size{0};

      for (auto && type : this->elementTypes(ghost_type)) {
        auto && array = this->array(type, ghost_type);
        total_element += array.size();
        total_size += Int(array.size()) * array.getNbComponent();
      }
      return std::make_pair(total_element, total_size);
    }

    [[nodiscard]] Int
    getNbComponent(const ElementType & type,
                   const GhostType & ghost_type) const override {
      return this->array(type, ghost_type).getNbComponent();
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

          auto field = std::make_shared<FieldElementalArray<T>>(
              map_array(type, ghost_type), this->getSupport(), type,
              ghost_type);
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
        FieldType field_type = FieldType::_element_map_array,
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
        FieldType field_type = FieldType::_element_map_array,
        bool create_fields = true)
        : FieldFunctionElementMapArray(
              make_field(std::forward<ElementTypeMapArray_>(map_array_in),
                         support),
              support, std::forward<Function>(function), field_type,
              create_fields) {}

    [[nodiscard]] decltype(auto) arrayTyped(ElementType type,
                                            GhostType ghost_type) const {
      return (aka::as_type<FieldFunctionElementalArray<T, decltype(function)>>(
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

          auto field = std::make_shared<
              FieldFunctionElementalArray<T, const Function &>>(
              map_array_in->arrayTyped(type, ghost_type).getSharedPointer(),
              this->getSupport(), this->function, type, ghost_type);
          field->addProperty("name", name);
          if (field->hasProperty("data_location")) {
            field->removeProperty("data_location");
          }
          this->fields.emplace(std::pair(type, ghost_type), std::move(field));
        }
      }
    }

  protected:
    Function function;
    std::shared_ptr<FieldElementMapArrayTemplateBase<T>> map_array_in;
  };

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
                 FieldType::_internal_field, false) {
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
        for (auto ghost_type : ghost_types) {
          for (auto type : this->map_array_in->elementTypes(ghost_type)) {
            auto nb_elements = connectivities.array(type, ghost_type).size();
            auto nb_integration_points =
                this->map_array_in->array(type, ghost_type).size();

            nb_integration_points_per_elem(type, ghost_type) =
                nb_integration_points / nb_elements;
          }
        }
        this->function.setNbIntegtrationPoints(nb_integration_points_per_elem);
      }
    }

    Release last_release;
  };

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

#endif // DUMPER_FIELD_ELEMENTAL_HH_
