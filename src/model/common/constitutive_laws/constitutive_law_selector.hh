/**
 * @file   constitutive_law_selector.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  ven mar 26 2021
 *
 * @brief class to dispatch cl to elements
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
#include "element_type_map.hh"
/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_CONSTITUTIVE_LAW_SELECTOR_HH
#define AKANTU_CONSTITUTIVE_LAW_SELECTOR_HH

namespace akantu {

/**
 * main class to assign same or different constitutive_laws for different
 * elements
 */
class ConstitutiveLawSelector : public std::enable_shared_from_this<ConstitutiveLawSelector> {
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
 * class that assigns the first constitutive_law to regular elements by default
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

    const auto & mat_indexes = constitutive_law_index(element.type, element.ghost_type);
    if (element.element < mat_indexes.size()) {
      auto && tmp_mat = mat_indexes(element.element);
      if (tmp_mat != UInt(-1)) {
        return tmp_mat;
      }
    }

    return ConstitutiveLawSelector::operator()(element);
  }

private:
  const ElementTypeMapArray<UInt> & constitutive_law_index;
};

} // akantu


#endif /* AKANTU_CONSTITUTIVE_LAW_SELECTOR_HH */
