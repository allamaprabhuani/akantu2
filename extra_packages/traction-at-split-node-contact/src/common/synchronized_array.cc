/**
 * @file   synchronized_array.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Fri Mar 16 2018
 * @date last modification: Tue Sep 29 2020
 *
 * @brief  implementation of synchronized array function
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
// std
#include <fstream>
#include <iostream>

// simtools
#include "synchronized_array.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class T>
SynchronizedArray<T>::SynchronizedArray(Int size, Int nb_component,
                                        const_reference value, const ID & id,
                                        const_reference default_value,
                                        const std::string & restart_name)
    : SynchronizedArrayBase(), Array<T>(size, nb_component, value, id),
      default_value(default_value), restart_name(restart_name),
      deleted_elements(0), nb_added_elements(size), depending_arrays(0) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T>
void SynchronizedArray<T>::syncElements(SyncChoice sync_choice) {
  AKANTU_DEBUG_IN();

  if (sync_choice == _deleted) {
    for (auto && array : depending_arrays) {
      auto vec_size [[gnu::unused]] =
          array->syncDeletedElements(this->deleted_elements);
      AKANTU_DEBUG_ASSERT(vec_size == this->size_,
                          "Synchronized arrays do not have the same length"
                              << "(may be a double synchronization)");
    }
    this->deleted_elements.clear();
  }

  else if (sync_choice == _added) {
    for (auto && array : depending_arrays) {
      auto vec_size [[gnu::unused]] =
          array->syncAddedElements(this->nb_added_elements);
      AKANTU_DEBUG_ASSERT(vec_size == this->size_,
                          "Synchronized arrays do not have the same length"
                              << "(may be a double synchronization)");
    }
    this->nb_added_elements = 0;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class T>
Int SynchronizedArray<T>::syncDeletedElements(std::vector<Idx> & del_elements) {
  AKANTU_DEBUG_ASSERT(
      nb_added_elements == 0 and deleted_elements.empty(),
      "Cannot sync with a SynchronizedArray if it has already been modified");

  for (auto && el : del_elements) {
    erase(el);
  }

  syncElements(_deleted);

  return this->size_;
}

/* -------------------------------------------------------------------------- */
template <class T>
Int SynchronizedArray<T>::syncAddedElements(Int nb_add_elements) {
  AKANTU_DEBUG_ASSERT(
      nb_added_elements == 0 and deleted_elements.empty(),
      "Cannot sync with a SynchronizedArray if it has already been modified");

  for (UInt i = 0; i < nb_add_elements; ++i) {
    push_back(this->default_value);
  }
  syncElements(_added);

  return this->size_;
}

/* -------------------------------------------------------------------------- */
template <typename T>
void SynchronizedArray<T>::registerDependingArray(
    SynchronizedArrayBase & array) {
  this->depending_arrays.push_back(&array);
  array.syncAddedElements(this->size_);
}

/* -------------------------------------------------------------------------- */
template <typename T>
void SynchronizedArray<T>::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "SynchronizedArray<" << debug::demangle(typeid(T).name())
         << "> [" << std::endl;
  stream << space << " + default_value     : " << this->default_value
         << std::endl;
  stream << space << " + nb_added_elements : " << this->nb_added_elements
         << std::endl;
  stream << space << " + deleted_elements  : ";
  for (auto && deleted_element : deleted_elements) {
    stream << deleted_element << " ";
  }
  stream << std::endl;

  stream << space << " + depending_arrays : ";
  for (auto && depending_array : this->depending_arrays) {
    stream << depending_array->getID() << " ";
  }
  stream << std::endl;

  Array<T>::printself(stream, indent + 1);

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
template <typename T>
void SynchronizedArray<T>::dumpRestartFile(
    const std::string & file_name) const {
  AKANTU_DEBUG_ASSERT(
      nb_added_elements == 0 and deleted_elements.empty(),
      "Restart File for SynchronizedArray "
          << this->id << " should not be dumped as it is not synchronized yet");

  std::stringstream name;
  name << file_name << "-" << this->restart_name << ".rs";

  std::ofstream out_restart;
  out_restart.open(name.str().c_str());

  out_restart << this->size_ << " " << this->nb_component << std::endl;
  Real size_comp = this->size_ * this->nb_component;
  for (Int i = 0; i < size_comp; ++i) {
    out_restart << std::setprecision(12) << this->values[i] << " ";
  }

  out_restart << std::endl;
}

/* -------------------------------------------------------------------------- */
template <typename T>
void SynchronizedArray<T>::readRestartFile(const std::string & file_name) {
  AKANTU_DEBUG_ASSERT(
      nb_added_elements == 0 and deleted_elements.empty(),
      "Restart File for SynchronizedArray "
          << this->id << " should not be read as it is not synchronized yet");

  std::stringstream name;
  name << file_name << "-" << this->restart_name << ".rs";
  std::ifstream infile;
  infile.open(name.str().c_str());

  std::string line;

  // get size and nb_component info
  AKANTU_DEBUG_ASSERT(infile.good(), "Could not read restart file for "
                                         << "SynchronizedArray " << this->id);
  getline(infile, line);
  std::stringstream size_comp(line);
  size_comp >> this->size_;
  size_comp >> this->nb_component;

  // get elements in array
  getline(infile, line);
  std::stringstream data(line);
  for (Int i = 0; i < this->size_ * this->nb_component; ++i) {
    AKANTU_DEBUG_ASSERT(
        !data.eof(),
        "Read SynchronizedArray "
            << this->id
            << " got to the end of the file before having read all data!");
    data >> this->values[i];
  }
}

/* -------------------------------------------------------------------------- */
template class SynchronizedArray<Real>;
template class SynchronizedArray<UInt>;
template class SynchronizedArray<Int>;
template class SynchronizedArray<bool>;

} // namespace akantu
