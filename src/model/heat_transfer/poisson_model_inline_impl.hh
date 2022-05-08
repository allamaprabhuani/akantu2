/**
 * @file   poisson_model_inline_impl.hh
 *
 * @author Mohit Pundir <mohit.pundir@ethz.ch>
 *
 * @date creation: Sat May 07 2022
 * @date last modification: Sat May 07 2022
 *
 * @brief  Poisson model implementation of inline functions
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_named_argument.hh"
#include "constitutive_law_selector.hh"
#include "constitutive_law_selector_tmpl.hh"
#include "poisson_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_POISSON_MODEL_INLINE_IMPL_HH__
#define __AKANTU_POISSON_MODEL_INLINE_IMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline decltype(auto) PoissonModel::getConstitutiveLaws() {
  return make_dereference_adaptor(constitutive_laws);
}

/* -------------------------------------------------------------------------- */
inline decltype(auto) PoissonModel::getConstitutiveLaws() const {
  return make_dereference_adaptor(constitutive_laws);
}

/* -------------------------------------------------------------------------- */
inline ConstitutiveLaw & PoissonModel::getConstitutiveLaw(UInt law_index) {
  AKANTU_DEBUG_ASSERT(law_index < constitutive_laws.size(),
                      "The model " << id << " has no constituitve law no "
                                   << law_index);
  return *constitutive_laws[law_index];
}

/* -------------------------------------------------------------------------- */
inline const ConstitutiveLaw & PoissonModel::getConstitutiveLaw(UInt law_index) const {
  AKANTU_DEBUG_ASSERT(law_index < constitutive_laws.size(),
                      "The model " << id << " has no constituitve law no "
                                   << law_index);
  return *constitutive_laws[law_index];
}

/* -------------------------------------------------------------------------- */
inline ConstitutiveLaw & PoissonModel::getConstitutiveLaw(const std::string & name) {
  std::map<std::string, UInt>::const_iterator it =
      constitutive_laws_names_to_id.find(name);
  AKANTU_DEBUG_ASSERT(it != constitutive_laws_names_to_id.end(),
                      "The model " << id << " has no constituitve law named "
                                   << name);
  return *constitutive_laws[it->second];
}

/* -------------------------------------------------------------------------- */
inline UInt
PoissonModel::getConstitutiveLawIndex(const std::string & name) const {
  auto it = constitutive_laws_names_to_id.find(name);
  AKANTU_DEBUG_ASSERT(it != constitutive_laws_names_to_id.end(),
                      "The model " << id << " has no constituitve law named "
                                   << name);
  return it->second;
}

/* -------------------------------------------------------------------------- */
inline const ConstitutiveLaw &
PoissonModel::getConstitutiveLaw(const std::string & name) const {
  auto it = constitutive_laws_names_to_id.find(name);
  AKANTU_DEBUG_ASSERT(it != constitutive_laws_names_to_id.end(),
                      "The model " << id << " has no constituitve law named "
                                   << name);
  return *constitutive_laws[it->second];
}

/* -------------------------------------------------------------------------- */
} // namespace akantu

#endif /* __AKANTU_POISSON_MODEL_INLINE_IMPL_HH__ */
