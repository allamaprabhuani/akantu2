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
    FieldBase(SupportBase & support, std::type_index type, FieldType field_type)
        : support(support), type_(type), field_type(field_type) {}

    [[nodiscard]] std::type_index type() const { return type_; }
    AKANTU_GET_MACRO_AUTO(FieldType, field_type);
    AKANTU_GET_MACRO_AUTO(Support, support);
    AKANTU_GET_MACRO_AUTO_NOT_CONST(Support, support);

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
      return (field_type == FieldType::_node_array);
    }

    [[nodiscard]] inline bool is_quadrature_points_field() const {
      return (field_type == FieldType::_internal_field);
    }

    [[nodiscard]] inline bool is_elemental_map_field() const {
      return (field_type == FieldType::_element_map_array) or
             this->is_quadrature_points_field();
    }

    // for compute fields to apply the computation
    virtual void update() {}

    // for when the field needs internal memory
    virtual void allocate() {}

    // for compute fields to apply the reverse computation
    virtual void back_propagate() {}

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
    SupportBase & support;
    std::type_index type_;
    FieldType field_type{FieldType::_not_defined};
  };

  /* ------------------------------------------------------------------------ */
  class FieldArrayBase : public FieldBase {
  public:
    FieldArrayBase(SupportBase & support, std::type_index type,
                   FieldType field_type = FieldType::_array)
        : FieldBase(support, type, field_type) {}

    [[nodiscard]] virtual const void * data() const = 0;
    [[nodiscard]] virtual void * data() = 0;
    [[nodiscard]] virtual Int size() const = 0;
    [[nodiscard]] virtual Int getNbComponent() const = 0;
  };

  /* ------------------------------------------------------------------------ */
  template <class T, class Base> class FieldArrayTemplateBase : public Base {
  public:
    using Base::Base;
    [[nodiscard]] virtual Array<T> & getArray() = 0;
    [[nodiscard]] virtual const Array<T> & getArray() const = 0;

    auto getSharedPointer() {
      return std::dynamic_pointer_cast<FieldArrayTemplateBase<T, Base>>(
          this->shared_from_this());
    }
  };

  /* ------------------------------------------------------------------------ */
  template <typename T, class Array_ = Array<T> &, class Base = FieldArrayBase>
  class FieldArray : public FieldArrayTemplateBase<T, Base> {
  public:
    template <class... Args>
    FieldArray(Array_ && array, SupportBase & support, Args &&... args)
        : FieldArrayTemplateBase<T, Base>(
              support, typeid(T), std::forward<decltype(args)>(args)...),
          array(array) {}

    [[nodiscard]] const void * data() const override {
      if (this->size() == 0) {
        return nullptr;
      }
      return array.data();
    }

    [[nodiscard]] void * data() override {
      if (this->size() == 0) {
        return nullptr;
      }
      if constexpr (std::is_const_v<std::remove_reference_t<Array_>>) {
        AKANTU_EXCEPTION("The dumper field contains a const Array, thus cannot "
                         "give a non const data pointer");
      } else {
        return array.data();
      }
    }

    [[nodiscard]] Int getNbComponent() const override {
      // auto factor = array.size() / this->size();
      return array.getNbComponent();
    }

    [[nodiscard]] Array<T> & getArray() override {
      if constexpr (std::is_const_v<std::remove_reference_t<Array_>>) {
        AKANTU_EXCEPTION("The dumper field contains a const Array, thus cannot "
                         "give a non const reference");
      } else {

        return array;
      }
    }
    [[nodiscard]] const Array<T> & getArray() const override { return array; }

  private:
    Array_ array;
  };

  /* ------------------------------------------------------------------------ */
  template <class InT, class OutT, class Base = FieldArrayBase>
  class FieldComputeArray : public FieldArrayTemplateBase<OutT, Base> {
  public:
    template <class... Args>
    FieldComputeArray(
        std::shared_ptr<FieldArrayTemplateBase<InT, Base>> array_in,
        SupportBase & support, Args &&... args)
        : FieldArrayTemplateBase<OutT, Base>(
              support, typeid(OutT), std::forward<decltype(args)>(args)...),
          array_in(array_in) {
      for (auto prop : this->array_in->getPropertiesList()) {
        this->addProperty(prop, this->array_in->getPropertyVariant(prop));
      }
    }

    template <class... Args>
    FieldComputeArray(Array<InT> & array_in, SupportBase & support,
                      Args &&... args)
        : FieldComputeArray(make_field(array_in, support), support,
                            std::forward<decltype(args)>(args)...) {}

    ~FieldComputeArray() override {
      for (auto prop : this->getPropertiesList()) {
        this->array_in->addProperty(prop, this->getPropertyVariant(prop));
      }
    }

    [[nodiscard]] Int getNbComponent() const override {
      return this->array_in->getNbComponent();
    }

    [[nodiscard]] const void * data() const override {
      if (this->size() == 0) {
        return nullptr;
      }
      return array_out->data();
    }

    [[nodiscard]] void * data() override {
      if (this->size() == 0) {
        return nullptr;
      }
      return array_out->data();
    }

    [[nodiscard]] Array<OutT> & getArray() override { return *array_out; }
    [[nodiscard]] const Array<OutT> & getArray() const override {
      return *array_out;
    }

    void allocate() override {
      if (not array_out) {
        array_out =
            std::make_unique<Array<OutT>>(this->size(), this->getNbComponent());
      } else {
        array_out->resize(this->size());
      }
    }

    void update() override {
      if (unchanged() or this->size() == 0) {
        return;
      }

      this->array_in->update();

      allocate();
      compute();

      this->getRelease() = array_in->getRelease();
    }

    void back_propagate() override {
      back_compute();
      this->array_in->back_propagate();
    }

    virtual void compute() = 0;
    virtual void back_compute() { AKANTU_TO_IMPLEMENT(); };

  public:
    [[nodiscard]] virtual bool unchanged() const {
      return (this->getRelease() == array_in->getRelease()) and array_out;
    }

  protected:
    std::shared_ptr<FieldArrayTemplateBase<InT, Base>> array_in;
    std::unique_ptr<Array<OutT>> array_out;
  };

  /* ------------------------------------------------------------------------ */
  struct NoOpFunction {
    template <class T> decltype(auto) operator()(const T & t) { return t; }
    template <class T>
    decltype(auto) operator()(const T & t, ElementType /*type*/,
                              GhostType /*ghost_type*/) const {
      return t;
    }
  };
  /* ------------------------------------------------------------------------ */
  template <class T, class Function, class ReverseFunction,
            class Base = FieldArrayBase>
  class FieldFunctionArray
      : public FieldComputeArray<
            T, details::function_return_scalar_t<T, Function>, Base> {
    using OutT = details::function_return_scalar_t<T, Function>;
    using parent = FieldComputeArray<T, OutT, Base>;

  public:
    template <class... Args>
    FieldFunctionArray(
        std::shared_ptr<FieldArrayTemplateBase<T, Base>> array_in,
        SupportBase & support, Function && function,
        ReverseFunction && reverse_function, Args &&... args)
        : parent(array_in, support, std::forward<decltype(args)>(args)...),
          function(std::forward<decltype(function)>(function)),
          reverse_function(
              std::forward<decltype(reverse_function)>(reverse_function)) {}

    template <class... Args>
    FieldFunctionArray(Array<T> & array_in, SupportBase & support,
                       Function && function,
                       ReverseFunction && reverse_function, Args &&... args)
        : parent(array_in, support, std::forward<decltype(args)>(args)...),
          function(std::forward<decltype(function)>(function)),
          reverse_function(
              std::forward<decltype(reverse_function)>(reverse_function)) {}

    [[nodiscard]] Int getNbComponent() const override {
      if constexpr (details::has_getNbComponent_member<Function>) {
        return function.getNbComponent(this->array_in->getNbComponent());
      } else {
        return this->array_in->getNbComponent();
      }
    }

    void compute() override {
      auto nb_component_in = this->array_in->getNbComponent();
      auto nb_component_out = this->array_out->getNbComponent();
      for (auto && [out, in] :
           zip(make_view(*this->array_out, nb_component_out),
               make_view(const_cast<const FieldArrayTemplateBase<T, Base> &>(
                             *this->array_in)
                             .getArray(),
                         nb_component_in))) {
        auto && res = function(in);
        out = res;
      }
    }

  private:
    Function function;
    ReverseFunction reverse_function;
  };

} // namespace dumper
} // namespace akantu

/* -------------------------------------------------------------------------- */
#include "dumper_field_elemental.hh"
#include "dumper_field_nodal.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
namespace dumper {
  inline auto & field_cast(FieldBase & field,
                           field_type_not_defined_t /*unused*/) {
    AKANTU_EXCEPTION("There are no good reason to call this function");
    return field;
  }
  inline auto & field_cast(const FieldBase & field,
                           field_type_not_defined_t /*unused*/) {
    AKANTU_EXCEPTION("There are no good reason to call this function");
    return field;
  }

  inline decltype(auto) field_cast(FieldBase & field,
                                   field_type_array_t /*unused*/) {
    return aka::as_type<FieldArrayBase>(field);
  }
  inline decltype(auto) field_cast(FieldBase & field,
                                   field_type_element_array_t /*unused*/) {
    return aka::as_type<FieldElementalArrayBase>(field);
  }
  inline decltype(auto) field_cast(FieldBase & field,
                                   field_type_node_array_t /*unused*/) {
    return aka::as_type<FieldNodalArrayBase>(field);
  }
  inline decltype(auto) field_cast(FieldBase & field,
                                   field_type_element_map_array_t /*unused*/) {
    return aka::as_type<FieldElementMapArrayBase>(field);
  }
  inline decltype(auto) field_cast(FieldBase & field,
                                   field_type_internal_field_t /*unused*/) {
    return aka::as_type<FieldElementMapArrayBase>(field);
  }

  inline decltype(auto) field_cast(const FieldBase & field,
                                   field_type_array_t /*unused*/) {
    return aka::as_type<FieldArrayBase>(field);
  }
  inline decltype(auto) field_cast(const FieldBase & field,
                                   field_type_element_array_t /*unused*/) {
    return aka::as_type<FieldElementalArrayBase>(field);
  }
  inline decltype(auto) field_cast(const FieldBase & field,
                                   field_type_node_array_t /*unused*/) {
    return aka::as_type<FieldNodalArrayBase>(field);
  }
  inline decltype(auto) field_cast(const FieldBase & field,
                                   field_type_element_map_array_t /*unused*/) {
    return aka::as_type<FieldElementMapArrayBase>(field);
  }
  inline decltype(auto) field_cast(const FieldBase & field,
                                   field_type_internal_field_t /*unused*/) {
    return aka::as_type<FieldElementMapArrayBase>(field);
  }

} // namespace dumper
} // namespace akantu

#endif /* AKANTU_DUMPER_FIELD_HH_ */
