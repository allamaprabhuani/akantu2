/**
 * @file   constitutive_law_selector.hh
 *
 * @author Mohit Pundir <mohit.pundir@ethz.ch>
 *
 * @date creation: Tue May 3 2022
 * @date last modification: Tue May 3 2022
 *
 * @brief  class describing how to choose a constitutive law
 * for a given element
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2014-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "element.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONSTITUTIVE_LAW_SELECTOR_HH__
#define __AKANTU_CONSTITUTIVE_LAW_SELECTOR_HH__

namespace akantu {
class PoissonModel;
} // namespace akantu

/* -------------------------------------------------------------------------- */
namespace akantu {

/**
 * main class to assign same or differentconstitutive laws to
 * different elements
 */
class ConstitutiveLawSelector
    : public std::enable_shared_from_this<ConstitutiveLawSelector> {
public:
  ConstitutiveLawSelector() = default;
  virtual ~ConstitutiveLawSelector() = default;
  virtual inline UInt operator()(const Element & element) {
    if (fallback_selector) {
      return (*fallback_selector)(element);
    }

    return fallback_value;
  }

  inline void setFallback(UInt f) { fallback_value = f; }
  inline void
  setFallback(const std::shared_ptr<ConstitutiveLawSelector> & fallback_selector) {
    this->fallback_selector = fallback_selector;
  }

  inline void setFallback(ConstitutiveLawSelector & fallback_selector) {
    this->fallback_selector = fallback_selector.shared_from_this();
  }

  inline std::shared_ptr<ConstitutiveLawSelector> & getFallbackSelector() {
    return this->fallback_selector;
  }

  inline UInt getFallbackValue() const { return this->fallback_value; }

protected:
  UInt fallback_value{0};
  std::shared_ptr<ConstitutiveLawSelector> fallback_selector;
};

/* -------------------------------------------------------------------------- */
/**
 * class that assigns the first constitutive law to regular elements by default
 */
class DefaultConstitutiveLawSelector : public ConstitutiveLawSelector {
public:
  explicit DefaultConstitutiveLawSelector(
      const ElementTypeMapArray<UInt> & constitutive_law_index)
      : constitutive_law_index(constitutive_law_index) {}

  UInt operator()(const Element & element) override {
    if (not constitutive_law_index.exists(element.type, element.ghost_type)) {
      return ConstitutiveLawSelector::operator()(element);
    }

    const auto & law_indexes =
        constitutive_law_index(element.type, element.ghost_type);
    if (element.element < law_indexes.size()) {
      auto && tmp_phase = law_indexes(element.element);
      if (tmp_phase != UInt(-1)) {
        return tmp_phase;
      }
    }

    return ConstitutiveLawSelector::operator()(element);
  }

private:
  const ElementTypeMapArray<UInt> & constitutive_law_index;
};

/* -------------------------------------------------------------------------- */
/**
 * Use elemental data to assign constitutive laws
 */
template <typename T>
class ElementDataConstitutiveLawSelector : public ConstitutiveLawSelector {
public:
  ElementDataConstitutiveLawSelector(const ElementTypeMapArray<T> & element_data,
                                const PoissonModel & model,
                                UInt first_index = 1)
      : element_data(element_data), model(model), first_index(first_index) {}

  inline T elementData(const Element & element) {
    DebugLevel dbl = debug::getDebugLevel();
    debug::setDebugLevel(dblError);
    T data = element_data(element.type, element.ghost_type)(element.element);
    debug::setDebugLevel(dbl);
    return data;
  }

  inline UInt operator()(const Element & element) override {
    return ConstitutiveLawSelector::operator()(element);
  }

protected:
  /// list of element with the specified data (i.e. tag value)
  const ElementTypeMapArray<T> & element_data;

  /// the model that the constitutive laws belong
  const PoissonModel & model;

  /// first constitutive law index: equal to 1 if none specified
  UInt first_index;
};

/* -------------------------------------------------------------------------- */
/**
 * class to use mesh data information to assign different phasefields
 * where name is the tag value: tag_0, tag_1
 */
template <typename T>
class MeshDataConstitutiveLawSelector : public ElementDataConstitutiveLawSelector<T> {
public:
  MeshDataConstitutiveLawSelector(const std::string & name,
                             const PoissonModel & model,
                             UInt first_index = 1);
};

} // namespace akantu

#endif /* __AKANTU_CONSTITUTIVE_LAW_SELECTOR_HH__ */
