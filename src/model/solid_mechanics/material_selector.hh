/**
 * @file   material_selector.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Fri Apr 09 2021
 *
 * @brief  class describing how to choose a material for a given element
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

#ifndef AKANTU_MATERIAL_SELECTOR_HH_
#define AKANTU_MATERIAL_SELECTOR_HH_

/* -------------------------------------------------------------------------- */
namespace akantu {

class SolidMechanicsModel;

using DefaultMaterialSelector = DedfaultConstitutiveLawSelector;
/* -------------------------------------------------------------------------- */
/**
 * Use elemental data to assign materials
 */
template <typename T>
class ElementDataMaterialSelector : public MaterialSelector {
public:
  ElementDataMaterialSelector(const ElementTypeMapArray<T> & element_data,
                              const SolidMechanicsModel & model,
                              UInt first_index = 1)
      : element_data(element_data), model(model), first_index(first_index) {}

  inline T elementData(const Element & element) {
    DebugLevel dbl = debug::getDebugLevel();
    debug::setDebugLevel(dblError);
    T data = element_data(element);
    debug::setDebugLevel(dbl);
    return data;
  }

  inline UInt operator()(const Element & element) override;

protected:
  /// list of element with the specified data (i.e. tag value)
  const ElementTypeMapArray<T> & element_data;

  /// the model that the materials belong
  const SolidMechanicsModel & model;

  /// first material index: equal to 1 if none specified
  UInt first_index;
};

/* -------------------------------------------------------------------------- */
/**
 * class to use mesh data information to assign different materials
 * where name is the tag value: tag_0, tag_1
 */
template <typename T>
class MeshDataMaterialSelector : public ElementDataMaterialSelector<T> {
public:
  MeshDataMaterialSelector(const std::string & name,
                           const SolidMechanicsModel & model,
                           UInt first_index = 1);
};

} // namespace akantu

#endif /* AKANTU_MATERIAL_SELECTOR_HH_ */
