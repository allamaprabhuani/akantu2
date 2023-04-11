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
#include "phasefield_selector.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_PHASEFIELD_SELECTOR_TMPL_HH_
#define AKANTU_PHASEFIELD_SELECTOR_TMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <>
inline Idx ElementDataPhaseFieldSelector<std::string>::operator()(
    const Element & element) {
  try {
    std::string material_name = this->elementData(element);
    return model.getPhaseFieldIndex(material_name);
  } catch (...) {
    return PhaseFieldSelector::operator()(element);
  }
}

/* -------------------------------------------------------------------------- */
template <>
inline Idx
ElementDataPhaseFieldSelector<Idx>::operator()(const Element & element) {
  try {
    return this->elementData(element) - first_index;
  } catch (...) {
    return PhaseFieldSelector::operator()(element);
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
MeshDataPhaseFieldSelector<T>::MeshDataPhaseFieldSelector(
    const std::string & name, const PhaseFieldModel & model, Idx first_index)
    : ElementDataPhaseFieldSelector<T>(model.getMesh().getData<T>(name), model,
                                       first_index) {}

} // namespace akantu

#endif /* AKANTU_PHASEFIELD_SELECTOR_TMPL_HH_ */
