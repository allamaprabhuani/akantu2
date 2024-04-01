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
      return (field_type == FieldType::_node_array);
    }

    [[nodiscard]] inline bool is_quadrature_points_field() const {
      return (field_type == FieldType::_internal_field);
    }

    [[nodiscard]] inline bool is_elemental_map_field() const {
      return (field_type == FieldType::_element_map_array) or
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
                   FieldType field_type = FieldType::_array)
        : FieldBase(support, type, field_type) {}

    [[nodiscard]] virtual const void * data() = 0;
    [[nodiscard]] virtual Int size() const = 0;
    [[nodiscard]] virtual Int getNbComponent() const = 0;
  };

  /* ------------------------------------------------------------------------ */
  template <class T, class Base> class FieldArrayTemplateBase : public Base {
  public:
    using Base::Base;
    //[[nodiscard]] virtual const Array<T> & getArray() const = 0;
    [[nodiscard]] virtual const Array<T> & getArray() = 0;

    auto getSharedPointer() {
      return std::dynamic_pointer_cast<FieldArrayTemplateBase<T, Base>>(
          this->shared_from_this());
    }
  };

  /* ------------------------------------------------------------------------ */
  template <typename T, class Array_ = const Array<T> &,
            class Base = FieldArrayBase>
  class FieldArray : public FieldArrayTemplateBase<T, Base> {
  public:
    template <class... Args>
    FieldArray(Array_ && array, const SupportBase & support, Args &&... args)
        : FieldArrayTemplateBase<T, Base>(
              support, typeid(T), std::forward<decltype(args)>(args)...),
          array(array) {}

    [[nodiscard]] const void * data() override { return array.data(); }

    [[nodiscard]] Int getNbComponent() const override {
      // auto factor = array.size() / this->size();
      return array.getNbComponent();
    }

    //[[nodiscard]] const Array<T> & getArray() const override { return array; }
    [[nodiscard]] const Array<T> & getArray() override { return array; }

  private:
    Array_ array;
  };

  /* ------------------------------------------------------------------------ */
  template <class InT, class OutT, class Base = FieldArrayBase>
  class FieldComputeArray : public FieldArrayTemplateBase<OutT, Base> {
  public:
    template <class... Args>
    FieldComputeArray(
        const std::shared_ptr<FieldArrayTemplateBase<InT, Base>> & array_in,
        const SupportBase & support, Args &&... args)
        : FieldArrayTemplateBase<OutT, Base>(
              support, typeid(OutT), std::forward<decltype(args)>(args)...),
          array_in(array_in) {}

    template <class... Args>
    FieldComputeArray(const Array<InT> & array_in, const SupportBase & support,
                      Args &&... args)
        : FieldComputeArray(make_field(array_in, support), support,
                            std::forward<decltype(args)>(args)...) {}

    [[nodiscard]] Int getNbComponent() const override {
      return this->array_in->getNbComponent();
    }

    [[nodiscard]] const void * data() {
      update();
      return array_out->data();
    }

    [[nodiscard]] const Array<OutT> & getArray() override {
      update();
      return *array_out;
    }

    [[nodiscard]] Release getRelease() const override { return last_release; }

    void update() override {
      if (unchanged() or this->size() == 0) {
        return;
      }

      if (not array_out) {
        array_out =
            std::make_unique<Array<OutT>>(this->size(), this->getNbComponent());
      } else {
        array_out->resize(this->size());
      }

      compute();

      last_release = getRelease();
    }

    virtual void compute() = 0;

  public:
    [[nodiscard]] virtual bool unchanged() const {
      return (getRelease() == last_release) and array_out;
    }

  protected:
    std::shared_ptr<FieldArrayTemplateBase<InT, Base>> array_in;
    std::unique_ptr<Array<OutT>> array_out;
    Release last_release;
  };

  /* ------------------------------------------------------------------------ */
  template <class T, class Function, class Base = FieldArrayBase>
  class FieldFunctionArray
      : public FieldComputeArray<
            T, details::function_return_scalar_t<T, Function>, Base> {
    using OutT = details::function_return_scalar_t<T, Function>;

    using parent = FieldComputeArray<T, OutT, Base>;

  public:
    template <class... Args>
    FieldFunctionArray(
        const std::shared_ptr<FieldArrayTemplateBase<T, Base>> & array_in,
        const SupportBase & support, Function && function, Args &&... args)
        : parent(array_in, support, std::forward<decltype(args)>(args)...),
          function(std::forward<decltype(function)>(function)) {}

    template <class... Args>
    FieldFunctionArray(const Array<T> & array_in, const SupportBase & support,
                       Function && function, Args &&... args)
        : parent(array_in, support, std::forward<decltype(args)>(args)...),
          function(std::forward<decltype(function)>(function)) {}

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
               make_view(this->array_in->getArray(), nb_component_in))) {
        auto && res = function(in);
        out = res;
      }
    }

  private:
    Function function;
  };

} // namespace dumper
} // namespace akantu

/* -------------------------------------------------------------------------- */
#include "dumper_field_elemental.hh"
#include "dumper_field_nodal.hh"
/* -------------------------------------------------------------------------- */

#endif /* AKANTU_DUMPER_FIELD_HH_ */
