/**
 * @file   support.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Oct 03 2017
 *
 * @brief Describe a support for data
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
#include "aka_common.hh"
#include "internal_field.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include <map>
#include <typeindex>
#include <variant>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SUPPORT_HH__
#define __AKANTU_SUPPORT_HH__

namespace akantu {
namespace dumper {
  class FieldBase;
  template <class T> class FieldNodeArrayTemplateBase;
  template <class T> class FieldElementMapArrayTemplateBase;

  class VariableBase;

  enum class SupportType {
    _mesh,
    _element_group,
  };

  enum class FieldType {
    _not_defined,
    _node_array,
    _node_array_function,
    _element_map_array,
    _element_map_array_function,
    _internal_field,
    _internal_field_function,
  };

  /* ------------------------------------------------------------------------ */
  class PropertiesManager {
  public:
    using Property = std::variant<Int, Real, bool, std::string, FieldType,
                                  std::tuple<ElementType, GhostType>>;

    PropertiesManager() = default;
    PropertiesManager(const PropertiesManager &) = default;
    PropertiesManager(PropertiesManager &&) = default;
    PropertiesManager & operator=(const PropertiesManager &) = default;
    PropertiesManager & operator=(PropertiesManager &&) = default;

    virtual ~PropertiesManager() = default;

    template <class T,
              std::enable_if_t<std::is_integral_v<T> and
                               not std::is_same_v<T, bool>> * = nullptr>
    void addProperty(const std::string & property, const T & value) {
      properties[property] = Int(value);
    }

    template <class T,
              std::enable_if_t<std::is_floating_point_v<T>> * = nullptr>
    void addProperty(const std::string & property, const T & value) {
      properties[property] = Real(value);
    }

    template <class T,
              std::enable_if_t<not((std::is_integral_v<T> and
                                    not std::is_same_v<T, bool>) or
                                   std::is_floating_point_v<T>)> * = nullptr>
    void addProperty(const std::string & property, const T & value) {
      properties[property] = value;
    }

    void removeProperty(const std::string & property) {
      properties.erase(property);
    }

    template <class T, std::enable_if_t<std::is_integral_v<T>> * = nullptr>
    [[nodiscard]] T getProperty(const std::string & property) const {
      return std::get<Int>(properties.at(property));
    }

    template <class T,
              std::enable_if_t<std::is_floating_point_v<T>> * = nullptr>
    [[nodiscard]] T getProperty(const std::string & property) const {
      return std::get<Real>(properties.at(property));
    }

    template <class T,
              std::enable_if_t<not std::disjunction_v<
                  std::is_integral<T>, std::is_floating_point<T>>> * = nullptr>
    [[nodiscard]] const T & getProperty(const std::string & property) const {
      return std::get<T>(properties.at(property));
    }

    [[nodiscard]] bool hasProperty(const std::string & property) const {
      return (properties.find(property) != properties.end());
    }

    [[nodiscard]] virtual ID getName() const {
      return getProperty<std::string>("name");
    }

  private:
    std::map<std::string, Property> properties;
  };

  /* ------------------------------------------------------------------------ */
  class SupportBase : public PropertiesManager {
  public:
    using Fields = std::map<ID, std::shared_ptr<FieldBase>>;
    using Variables = std::map<ID, std::shared_ptr<VariableBase>>;
    using SubSupports = std::map<ID, std::shared_ptr<SupportBase>>;

    explicit SupportBase(SupportType type) : type(type) {}

    AKANTU_GET_MACRO(Type, type, SupportType);

    template <typename T> inline auto & addSubSupport(T & type);
    [[nodiscard]] inline const auto & getSubSupports() const;
    [[nodiscard]] inline const auto & getParentSupport() const;

    inline void addField(const ID & id, std::shared_ptr<FieldBase> && field);

    template <class T> void addField(const ID & id, InternalField<T> & data);
    template <class T>
    void addField(const ID & id, ElementTypeMapArray<T> & data);
    template <class T> void addField(const ID & id, Array<T> & data);

    [[nodiscard]] decltype(auto) getFields() const { return (fields_); }

    [[nodiscard]] decltype(auto) getField(const ID & field) const {
      return fields_.at(field);
    }

    [[nodiscard]] decltype(auto) getField(const ID & field) {
      return fields_.at(field);
    }

    [[nodiscard]] bool hasField(const ID & field) const {
      return (fields_.find(field) != fields_.end());
    }

  protected:
    Fields fields_;
    SubSupports sub_supports_;
    Variables variables_;
    SupportType type;

  private:
    SupportBase * parent{nullptr};
  };

  /* ------------------------------------------------------------------------ */
  template <class Inner> class Support : public SupportBase {
    //[[nodiscard]] ID getName() const override { return "not implemented"; };
  };

  /* ------------------------------------------------------------------------ */
  class SupportElements {
  public:
    using ElementTypesIteratorHelper =
        ElementTypeMapArray<Idx, ElementType>::ElementTypesIteratorHelper;

    [[nodiscard]] virtual ElementTypesIteratorHelper
    elementTypes(GhostType ghost_type = _not_ghost) const = 0;

    SupportElements() = default;
    SupportElements(const SupportElements &) = default;
    SupportElements(SupportElements &&) = default;
    SupportElements & operator=(const SupportElements &) = default;
    SupportElements & operator=(SupportElements &&) = default;

    virtual ~SupportElements() = default;

    [[nodiscard]] auto & getNodes() const { return *nodes; }
    [[nodiscard]] auto & getConnectivities() const { return *connectivities; }

  protected:
    std::shared_ptr<FieldNodeArrayTemplateBase<Real>> nodes;
    std::shared_ptr<FieldElementMapArrayTemplateBase<Idx>> connectivities;
  };

} // namespace dumper
} // namespace akantu

#include "support_tmpl.hh"

namespace akantu {
template <typename T> auto make_support(T & t) {
  return std::make_unique<dumper::Support<T>>(t);
}
} // namespace akantu

#endif /* __AKANTU_SUPPORT_HH__ */
