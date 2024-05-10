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
#include "dumper_types.hh"
#include "internal_field.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include <hdf5.h>
#include <tuple>
#include <type_traits>
#include <variant>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SUPPORT_HH__
#define __AKANTU_SUPPORT_HH__

namespace akantu {
namespace dumper {
  /* ------------------------------------------------------------------------ */
  class PropertiesManager {
    using integral_type = std::common_type_t<size_t, hsize_t, hid_t>;

  public:
    using slabs_type =
        std::vector<std::tuple<integral_type, integral_type, integral_type>>;

    using Property =
        std::variant<integral_type, Real, bool, std::string, FieldType,
                     FieldUsageType, std::tuple<ElementType, GhostType>,
                     Release, ElementTypeMap<integral_type>, slabs_type>;

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
      properties[property] = integral_type(value);
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
      AKANTU_DEBUG_ASSERT(properties.find(property) != properties.end(),
                          "The property " << property << " is not registered");
      properties.erase(property);
    }

    [[nodiscard]] Property
    getPropertyVariant(const std::string & property) const {
      AKANTU_DEBUG_ASSERT(properties.find(property) != properties.end(),
                          "The property " << property << " is not registered");
      return properties.at(property);
    }

    template <class T, std::enable_if_t<std::is_integral_v<T>> * = nullptr>
    [[nodiscard]] T getProperty(const std::string & property) const {
      AKANTU_DEBUG_ASSERT(properties.find(property) != properties.end(),
                          "The property " << property << " is not registered");
      return std::get<integral_type>(properties.at(property));
    }

    template <class T,
              std::enable_if_t<std::is_floating_point_v<T>> * = nullptr>
    [[nodiscard]] T getProperty(const std::string & property) const {
      AKANTU_DEBUG_ASSERT(properties.find(property) != properties.end(),
                          "The property " << property << " is not registered");
      return std::get<Real>(properties.at(property));
    }

    template <class T,
              std::enable_if_t<not std::disjunction_v<
                  std::is_integral<T>, std::is_floating_point<T>>> * = nullptr>
    [[nodiscard]] const T & getProperty(const std::string & property) const {
      AKANTU_DEBUG_ASSERT(properties.find(property) != properties.end(),
                          "The property " << property << " is not registered");
      return std::get<T>(properties.at(property));
    }

    [[nodiscard]] bool hasProperty(const std::string & property) const {
      return (properties.find(property) != properties.end());
    }

    [[nodiscard]] virtual ID getName() const {
      return getProperty<std::string>("name");
    }

    [[nodiscard]] virtual const Release & getRelease() const { return release; }
    [[nodiscard]] virtual Release & getRelease() { return release; }

    decltype(auto) getPropertiesList() const {
      return make_keys_adaptor(properties);
    }

  private:
    std::map<std::string, Property> properties;
    Release release;
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

    inline void addField(const ID & id,
                         const std::shared_ptr<FieldBase> & field,
                         FieldUsageType usage = FieldUsageType::_visualisation);

    template <class Cont,
              std::enable_if_t<not std::is_convertible_v<
                  std::decay_t<Cont>, std::shared_ptr<FieldBase>>> * = nullptr>
    inline void addField(const ID & id, Cont && data,
                         FieldUsageType usage = FieldUsageType::_visualisation);

    template <class Cont, class Func,
              std::enable_if_t<not std::is_convertible_v<
                  std::decay_t<Cont>, std::shared_ptr<FieldBase>>> * = nullptr>
    inline void addField(const ID & id, Cont && data, Func && func,
                         FieldUsageType usage = FieldUsageType::_visualisation);

    [[nodiscard]] decltype(auto) getFields() const { return (fields_); }

    [[nodiscard]] decltype(auto) getField(const ID & field) const {
      AKANTU_DEBUG_ASSERT(fields_.find(field) != fields_.end(),
                          "The field " << field
                                       << " is not registered in the support");
      return fields_.at(field);
    }

    [[nodiscard]] decltype(auto) getField(const ID & field) {
      AKANTU_DEBUG_ASSERT(fields_.find(field) != fields_.end(),
                          "The field " << field
                                       << " is not registered in the support");
      return fields_.at(field);
    }

    [[nodiscard]] bool hasField(const ID & field) const {
      return (fields_.find(field) != fields_.end());
    }

    [[nodiscard]] virtual bool isDistributed() const { return false; }
    [[nodiscard]] virtual const Communicator & getCommunicator() const {
      return Communicator::getSelfCommunicator();
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

    [[nodiscard]] auto & getNodes() const { return *positions; }
    [[nodiscard]] auto & getConnectivities() const { return *connectivities; }

    [[nodiscard]] virtual Int getNbNodes() const = 0;
    [[nodiscard]] virtual Int getNbGlobalNodes() const = 0;
    [[nodiscard]] virtual Int getNbLocalNodes() const = 0;
    [[nodiscard]] virtual Int
    getNbElements(const ElementType & type,
                  const GhostType & ghost_type = _not_ghost) const = 0;

    [[nodiscard]] virtual Int
    getElementsOffsets(const ElementType & type,
                       const GhostType & ghost_type = _not_ghost) const = 0;

    [[nodiscard]] virtual Int
    getNbGlobalElements(const ElementType & type,
                        const GhostType & ghost_type = _not_ghost) const = 0;

    virtual void updateOffsets() {}

  protected:
    std::shared_ptr<FieldNodalArrayTemplateBase<Real>> positions;
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
