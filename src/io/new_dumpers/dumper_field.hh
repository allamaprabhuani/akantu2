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
#include "element_type_map.hh"
#include "support.hh"
/* -------------------------------------------------------------------------- */
#include <sstream>
#include <typeindex>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_DUMPER_FIELD_HH__
#define __AKANTU_DUMPER_FIELD_HH__

namespace akantu {

namespace dumper {
  enum class FieldType {
    _not_defined,
    _node_array,
    _node_array_function,
    _element_map_array,
    _element_map_array_function
  };

  /// Field interface
  class FieldBase : public PropertiesManager {
  public:
    FieldBase(const SupportBase & support, std::type_index type,
              FieldType field_type)
        : support(support), type_(type), field_type(field_type) {}

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

  /* ------------------------------------------------------------------------ */
  class FieldNodeArrayBase : public FieldBase {
  public:
    FieldNodeArrayBase(const SupportBase & support, std::type_index type,
                       FieldType field_type = FieldType::_node_array)
        : FieldBase(support, type, field_type) {}

    [[nodiscard]] virtual const void * data() const = 0;
    [[nodiscard]] virtual void * data() = 0;

    [[nodiscard]] virtual Int size() const = 0;
    [[nodiscard]] virtual Int getNbComponent() const = 0;
  };

  /* ------------------------------------------------------------------------ */
  template <typename T> class FieldNodeArray : public FieldNodeArrayBase {
  public:
    FieldNodeArray(Array<T> & array, const SupportBase & support)
        : FieldNodeArrayBase(support, typeid(T)), array(array) {}

    [[nodiscard]] const void * data() const override { return array.data(); }
    [[nodiscard]] void * data() override { return array.data(); }

    [[nodiscard]] Int size() const override { return array.size(); }
    [[nodiscard]] Int getNbComponent() const override {
      return array.getNbComponent();
    }

    [[nodiscard]] const Array<T> & getArray() const { return array; }

  private:
    Array<T> & array;
  };

  /* ------------------------------------------------------------------------ */
  template <class T, class Function>
  class FieldFunctionNodeArray : public FieldNodeArrayBase {
  public:
    FieldFunctionNodeArray(const Array<T> & array_in,
                           const SupportBase & support, Function && function)
        : FieldNodeArrayBase(support, typeid(T),
                             FieldType::_node_array_function),
          function(std::forward<Function>(function)), array_in(array_in) {
      update();
    }

    [[nodiscard]] const void * data() const override {
      return array_out->data();
    }
    [[nodiscard]] void * data() override {
      update();
      return array_out->data();
    }

    [[nodiscard]] Int size() const override { return array_in.size(); }
    [[nodiscard]] Int getNbComponent() const override {
      if constexpr (std::is_class_v<Function>) {
        return function.getNbComponent(array_in.getNbComponent());
      } else {
        return array_in.getNbComponent();
      }
    }

    [[nodiscard]] const Array<T> & getArray() const { return *array_out; }
    [[nodiscard]] const Array<T> & getArray() {
      update();
      return *array_out;
    }

  private:
    void update() {
      if (not array_out) {
        array_out =
            std::make_unique<Array<T>>(array_in.size(), this->getNbComponent());
      } else {
        array_out->resize(array_in.size());
      }

      for (auto && [out, in] :
           zip(make_view(*array_out, array_out->getNbComponent()),
               make_view(array_in, array_in.getNbComponent()))) {
        out = function(in);
      }
    }

  private:
    Function function;
    const Array<T> & array_in;
    std::unique_ptr<Array<T>> array_out;
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

    [[nodiscard]] virtual FieldNodeArrayBase &
    array(ElementType type, GhostType ghost_type = _not_ghost) const {
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
             std::unique_ptr<FieldNodeArrayBase>>
        fields;
  };

  /* ------------------------------------------------------------------------ */
  template <typename T>
  class FieldElementMapArray : public FieldElementMapArrayBase {
  public:
    FieldElementMapArray(ElementTypeMapArray<T> & map_array,
                         const SupportBase & support)
        : FieldElementMapArrayBase(support, typeid(T)), map_array(map_array) {
      for (auto ghost_type : ghost_types) {
        for (auto && type : this->support_element.elementTypes(ghost_type)) {
          auto name = std::to_string(type);
          if (ghost_type != _not_ghost) {
            name + ":" + std::to_string(ghost_type);
          }

          auto field =
              std::make_unique<FieldNodeArray<T>>(map_array(type), support);
          field->addProperty("name", name);

          this->fields.emplace(std::pair(type, ghost_type), std::move(field));
        }
      }
    }

    [[nodiscard]] decltype(auto)
    arrayTyped(ElementType type, GhostType ghost_type = _not_ghost) const {
      return (aka::as_type<FieldNodeArray<T>>(
          *this->fields.at({type, ghost_type})));
    }

  private:
    const ElementTypeMapArray<T> & map_array;
  };

  /* ------------------------------------------------------------------------ */
  template <class T, class Function>
  class FieldFunctionElementMapArray : public FieldElementMapArrayBase {
    template <class ElementFunction> struct ElementTypeMapArrayFunctor {
      ElementTypeMapArrayFunctor(ElementType type,
                                 const ElementFunction & function)
          : type(type), function(function) {}

      [[nodiscard]] Int getNbComponent(UInt nb_component) const {
        return function.getNbComponent(nb_component, type);
      }

      [[nodiscard]] decltype(auto) operator()(const Vector<T> & in) const {
        return function(in, type);
      }

      ElementType type;
      const ElementFunction & function;
    };

  public:
    FieldFunctionElementMapArray(ElementTypeMapArray<T> & map_array_in,
                                 const SupportBase & support,
                                 Function && function)
        : FieldElementMapArrayBase(support, typeid(T),
                                   FieldType::_element_map_array_function),
          function(std::forward<Function>(function)),
          map_array_in(map_array_in) {
      for (auto ghost_type : ghost_types) {
        for (auto && type : aka::as_type<SupportElements>(support_element)
                                .elementTypes(ghost_type)) {
          auto name = std::to_string(type);
          if (ghost_type != _not_ghost) {
            name += ":" + std::to_string(ghost_type);
          }

          // if constexpr (std::is_class_v<Function>) {
          this->fields.emplace(
              std::pair(type, ghost_type),
              make_field(map_array_in(type, ghost_type), support,
                         ElementTypeMapArrayFunctor(type, this->function)));
          // } else {
          // auto field = make_field(map_array_in(type), support, [&](auto &&
          // in)
          // {
          //   return this->function(in, type);
          // });
          // field->addProperty("name", name);
          // this->fields.emplace(type, std::move(field));
          //}
        }
      }
    }

    [[nodiscard]] decltype(auto) arrayTyped(ElementType type,
                                            GhostType ghost_type) const {
      return (aka::as_type<FieldFunctionNodeArray<T, decltype(function)>>(
          *this->fields.at({type, ghost_type})));
    }

  private:
    Function function;
    const ElementTypeMapArray<T> & map_array_in;
  };

  struct ConnectivityFunctor {
    Int getNbComponent(Int nb_component, ElementType type) const {
      if (type == _segment_2 or type == _segment_3) {
        return nb_component + 2;
      }
      return nb_component + 1;
    }

    Vector<Idx> operator()(const Vector<Idx> & connectivity,
                           ElementType type) const {
      Int offset{1};
      if (type == _segment_2 or type == _segment_3) {
        offset = 2;
      }

      Vector<Idx> result(connectivity.size() + offset);
      for (auto && [i, c] : enumerate(connectivity)) {
        result(i + offset) = c;
      }

      result[0] = aka_type_to_dumper_type.at(type);

      if (type == _segment_2) {
        result[0] = 2;
      } else if (type == _segment_3) {
        result[0] = 3;
      }

      return result;
    }

  private:
    const std::unordered_map<ElementType, int> aka_type_to_dumper_type{
        {_point_1, 1},        {_segment_2, 2},      {_segment_3, 2},
        {_triangle_3, 4},     {_triangle_6, 36},    {_quadrangle_4, 5},
        {_quadrangle_8, 37},  {_tetrahedron_4, 6},  {_tetrahedron_10, 38},
        {_hexahedron_8, 9},   {_hexahedron_20, 48}, {_pentahedron_6, 8},
        {_pentahedron_15, 40}};
  };

  using FieldConnectivity =
      FieldFunctionElementMapArray<Idx, ConnectivityFunctor>;

  /* ------------------------------------------------------------------------ */
  template <typename T>
  auto make_field(Array<T> & array, const SupportBase & support) {
    return std::make_unique<FieldNodeArray<T>>(array, support);
  }

  /* ------------------------------------------------------------------------ */
  template <typename T, class Function>
  auto make_field(Array<T> & array, const SupportBase & support,
                  Function && function) {
    return std::make_unique<FieldFunctionNodeArray<T, Function>>(
        array, support, std::forward<Function>(function));
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  auto make_field(ElementTypeMapArray<T> & array, const SupportBase & support) {
    return std::make_unique<FieldElementMapArray<T>>(array, support);
  }

  /* ------------------------------------------------------------------------ */
  template <typename T, class Function>
  auto make_field(ElementTypeMapArray<T> & array, const SupportBase & support,
                  Function && function) {
    return std::make_unique<FieldFunctionElementMapArray<T, Function>>(
        array, support, std::forward<Function>(function));
  }
} // namespace dumper
} // namespace akantu

#endif /* __AKANTU_DUMPER_FIELD_HH__ */
