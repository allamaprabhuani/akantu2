/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
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
 */

/* -------------------------------------------------------------------------- */
#include "element.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_SELECTOR_HH_
#define AKANTU_MATERIAL_SELECTOR_HH_

/* -------------------------------------------------------------------------- */
namespace akantu {

class SolidMechanicsModel;

/**
 * main class to assign same or different materials for different
 * elements
 */
class AKANTU_EXPORT MaterialSelector {
public:
  MaterialSelector() = default;
  virtual ~MaterialSelector() = default;
  virtual inline Int operator()(const Element & element) {
    if (fallback_selector) {
      return (*fallback_selector)(element);
    }

    return fallback_value;
  }

  inline void setFallback(Int f) { fallback_value = f; }
  inline void
  setFallback(const std::shared_ptr<MaterialSelector> & fallback_selector) {
    this->fallback_selector = fallback_selector;
  }

  inline std::shared_ptr<MaterialSelector> & getFallbackSelector() {
    return this->fallback_selector;
  }

  inline Int getFallbackValue() const { return this->fallback_value; }

protected:
  Int fallback_value{0};
  std::shared_ptr<MaterialSelector> fallback_selector;
};

/* -------------------------------------------------------------------------- */
/**
 * class that assigns the first material to regular elements by default
 */
class AKANTU_EXPORT DefaultMaterialSelector : public MaterialSelector {
public:
  explicit DefaultMaterialSelector(
      const ElementTypeMapArray<Idx> & material_index)
      : material_index(material_index) {}

  Int operator()(const Element & element) override {
    if (not material_index.exists(element.type, element.ghost_type)) {
      return MaterialSelector::operator()(element);
    }

    const auto & mat_indexes = material_index(element.type, element.ghost_type);
    if (element.element < mat_indexes.size()) {
      auto && tmp_mat = mat_indexes(element.element);
      if (tmp_mat != Int(-1)) {
        return tmp_mat;
      }
    }

    return MaterialSelector::operator()(element);
  }

private:
  const ElementTypeMapArray<Idx> & material_index;
};

/* -------------------------------------------------------------------------- */
/**
 * Use elemental data to assign materials
 */
template <typename T>
class AKANTU_EXPORT ElementDataMaterialSelector : public MaterialSelector {
public:
  ElementDataMaterialSelector(const ElementTypeMapArray<T> & element_data,
                              const SolidMechanicsModel & model,
                              Int first_index = 1)
      : element_data(element_data), model(model), first_index(first_index) {}

  inline T elementData(const Element & element) {
    DebugLevel dbl = debug::getDebugLevel();
    debug::setDebugLevel(dblError);
    T data = element_data(element);
    debug::setDebugLevel(dbl);
    return data;
  }

  inline Int operator()(const Element & element) override;

protected:
  /// list of element with the specified data (i.e. tag value)
  const ElementTypeMapArray<T> & element_data;

  /// the model that the materials belong
  const SolidMechanicsModel & model;

  /// first material index: equal to 1 if none specified
  Int first_index;
};

/* -------------------------------------------------------------------------- */
/**
 * class to use mesh data information to assign different materials
 * where name is the tag value: tag_0, tag_1
 */
template <typename T>
class AKANTU_EXPORT MeshDataMaterialSelector
    : public ElementDataMaterialSelector<T> {
public:
  MeshDataMaterialSelector(const std::string & name,
                           const SolidMechanicsModel & model,
                           Int first_index = 1);
};

} // namespace akantu

#endif /* AKANTU_MATERIAL_SELECTOR_HH_ */
