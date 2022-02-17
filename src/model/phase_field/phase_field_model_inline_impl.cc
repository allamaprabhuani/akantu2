/**
 * @file   phase_field_model_inline_impl.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Fri Jun 19 2020
 *
 * @brief  Phase field implementation of inline functions
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
#include "phasefield_selector.hh"
#include "phasefield_selector_tmpl.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PHASE_FIELD_MODEL_INLINE_IMPL_CC__
#define __AKANTU_PHASE_FIELD_MODEL_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline decltype(auto) PhaseFieldModel::getPhaseFields() {
  return make_dereference_adaptor(phasefields);
}

/* -------------------------------------------------------------------------- */
inline decltype(auto) PhaseFieldModel::getPhaseFields() const {
  return make_dereference_adaptor(phasefields);
}

/* -------------------------------------------------------------------------- */
inline PhaseField & PhaseFieldModel::getPhaseField(UInt mat_index) {
  AKANTU_DEBUG_ASSERT(mat_index < phasefields.size(),
                      "The model " << id << " has no phasefield no "
                                   << mat_index);
  return *phasefields[mat_index];
}

/* -------------------------------------------------------------------------- */
inline const PhaseField & PhaseFieldModel::getPhaseField(UInt mat_index) const {
  AKANTU_DEBUG_ASSERT(mat_index < phasefields.size(),
                      "The model " << id << " has no phasefield no "
                                   << mat_index);
  return *phasefields[mat_index];
}

/* -------------------------------------------------------------------------- */
inline PhaseField & PhaseFieldModel::getPhaseField(const std::string & name) {
  std::map<std::string, UInt>::const_iterator it =
      phasefields_names_to_id.find(name);
  AKANTU_DEBUG_ASSERT(it != phasefields_names_to_id.end(),
                      "The model " << id << " has no phasefield named "
                                   << name);
  return *phasefields[it->second];
}

/* -------------------------------------------------------------------------- */
inline UInt
PhaseFieldModel::getPhaseFieldIndex(const std::string & name) const {
  auto it = phasefields_names_to_id.find(name);
  AKANTU_DEBUG_ASSERT(it != phasefields_names_to_id.end(),
                      "The model " << id << " has no phasefield named "
                                   << name);
  return it->second;
}

/* -------------------------------------------------------------------------- */
inline const PhaseField &
PhaseFieldModel::getPhaseField(const std::string & name) const {
  auto it = phasefields_names_to_id.find(name);
  AKANTU_DEBUG_ASSERT(it != phasefields_names_to_id.end(),
                      "The model " << id << " has no phasefield named "
                                   << name);
  return *phasefields[it->second];
}

/* -------------------------------------------------------------------------- */
} // namespace akantu

#endif /* __AKANTU_PHASE_FIELD_MODEL_INLINE_IMPL_CC__ */
