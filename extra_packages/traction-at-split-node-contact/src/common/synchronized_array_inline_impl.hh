/**
 * Copyright (©) 2016-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "synchronized_array.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {
/* -------------------------------------------------------------------------- */
template <typename T>
template <typename P>
inline void SynchronizedArray<T>::push_back(P && value) {
  AKANTU_DEBUG_ASSERT(deleted_elements.size() == 0,
                      "Cannot push_back element if SynchronizedArray"
                          << " is already modified without synchronization");

  Array<T>::push_back(std::forward<P>(value));
  this->nb_added_elements++;
}

/* -------------------------------------------------------------------------- */
template <typename T> inline void SynchronizedArray<T>::erase(Idx i) {
  AKANTU_DEBUG_ASSERT(nb_added_elements == 0,
                      "Cannot erase element if SynchronizedArray"
                          << " is already modified without synchronization");

  for (Int j = 0; j < this->nb_component; ++j) {
    this->values[i * this->nb_component + j] =
        this->values[(this->size_ - 1) * this->nb_component + j];
  }
  this->size_--;

  this->deleted_elements.push_back(i);
}

} // namespace akantu
