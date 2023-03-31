/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef AKANTU_IGFEM_INTERNAL_FIELD_TMPL_HH_
#define AKANTU_IGFEM_INTERNAL_FIELD_TMPL_HH_

namespace akantu {
/* -------------------------------------------------------------------------- */
template <typename T>
IGFEMInternalField<T>::IGFEMInternalField(const ID & id, Material & material)
    : InternalField<T>(
          id, material, material.getModel().getFEEngine("IGFEMFEEngine"),
          dynamic_cast<MaterialIGFEM &>(material).getElementFilter()) {
  this->element_kind = _ek_igfem;
}

/* -------------------------------------------------------------------------- */
template <typename T> IGFEMInternalField<T>::~IGFEMInternalField(){};

} // namespace akantu

#endif /* AKANTU_IGFEM_INTERNAL_FIELD_TMPL_HH_ */
