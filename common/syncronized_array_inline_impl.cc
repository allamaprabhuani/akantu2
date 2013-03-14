/**
 * @file   syncronized_array_inline_impl.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Wed Mar 13 07:22:25 2013
 *
 * @brief  inlined methods for the syncronized array
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
template <typename T> inline void SyncronizedArray<T>::push_back(const_reference value) {
  AKANTU_DEBUG_IN();
 
  AKANTU_DEBUG_ASSERT(deleted_elements.size() == 0, 
		      "Cannot push_back element if SyncronizedArray" << 
		      " is already modified without synchronization");

  Array<T>::push_back(value);
  this->nb_added_elements++;
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename T> inline void SyncronizedArray<T>::push_back(const value_type new_elem[]) {
  AKANTU_DEBUG_IN();
 
  AKANTU_DEBUG_ASSERT(deleted_elements.size() == 0, 
		      "Cannot push_back element if SyncronizedArray" << 
		      " is already modified without synchronization");

  Array<T>::push_back(new_elem);
  this->nb_added_elements++;
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename T> inline void SyncronizedArray<T>::erase(UInt i) {
  AKANTU_DEBUG_IN();
  
  AKANTU_DEBUG_ASSERT(nb_added_elements == 0, 
		      "Cannot erase element if SyncronizedArray" << 
		      " is already modified without synchronization");
  
  for (UInt j=0; j<this->nb_component; ++j)
    this->values[i*this->nb_component+j] = this->values[(this->size-1)*this->nb_component+j];
  this->size--;
  
  this->deleted_elements.push_back(i);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/*
template<typename T>
template<typename R>
inline void SyncronizedArray<T>::erase(const iterator<R> & it) {
  T * curr = it.getCurrentStorage();
  UInt pos = (curr - values) / nb_component;
  erase(pos);
}
*/
