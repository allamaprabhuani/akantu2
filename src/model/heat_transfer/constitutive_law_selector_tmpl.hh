/**
 * @file   constitutive_law_selector_tmpl.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Tue May 3 2022
 * @date last modification: Tue May 3 2022
 *
 * @brief  Implementation of the template ConstitutiveLawSelector
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
#include "constitutive_law_selector.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONSTITUTIVE_LAW_SELECTOR_TMPL_HH__
#define __AKANTU_CONSTITUTIVE_LAW_SELECTOR_TMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <>
inline UInt ElementDataConstitutiveLawSelector<std::string>::operator()(
    const Element & element) {
  try {
    std::string material_name = this->elementData(element);
    return model.getConstitutiveLawIndex(material_name);
  } catch (...) {
    return ConstitutiveLawSelector::operator()(element);
  }
}

/* -------------------------------------------------------------------------- */
template <>
inline UInt
ElementDataConstitutiveLawSelector<UInt>::operator()(const Element & element) {
  try {
    return this->elementData(element) - first_index;
  } catch (...) {
    return ConstitutiveLawSelector::operator()(element);
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
MeshDataConstitutiveLawSelector<T>::MeshDataConstitutiveLawSelector(
    const std::string & name, const PoissonModel & model, UInt first_index)
    : ElementDataConstitutiveLawSelector<T>(model.getMesh().getData<T>(name), model,
                                       first_index) {}

} // namespace akantu

#endif /* __AKANTU_CONSTITUTIVE_LAW_SELECTOR_TMPL_HH__ */
