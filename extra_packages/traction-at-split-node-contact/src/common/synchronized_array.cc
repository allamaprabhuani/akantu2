/**
 * @file   synchronized_array.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Jan 22 2016
 *
 * @brief  implementation of synchronized array function
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
// std
#include <iostream>
#include <fstream>

// simtools
#include "synchronized_array.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<class T>
SynchronizedArray<T>::SynchronizedArray(UInt size,
					UInt nb_component,
					SynchronizedArray<T>::const_reference value,
					const ID & id,
					SynchronizedArray<T>::const_reference default_value,
					const std::string restart_name) :
  SynchronizedArrayBase(),
  Array<T>(size, nb_component, value, id),
  default_value(default_value),
  restart_name(restart_name),
  deleted_elements(0),
  nb_added_elements(size),
  depending_arrays(0) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<class T>
void SynchronizedArray<T>::syncElements(SyncChoice sync_choice) {
  AKANTU_DEBUG_IN();

  if (sync_choice == _deleted) {
    std::vector<SynchronizedArrayBase *>::iterator it;
    for (it=depending_arrays.begin(); it != depending_arrays.end(); ++it) {
      UInt vec_size = (*it)->syncDeletedElements(this->deleted_elements);
      AKANTU_DEBUG_ASSERT(vec_size == this->size,
			  "Synchronized arrays do not have the same length" <<
			  "(may be a double synchronization)");
    }
    this->deleted_elements.clear();
  }

  else if (sync_choice == _added) {
    std::vector<SynchronizedArrayBase *>::iterator it;
    for (it=depending_arrays.begin(); it != depending_arrays.end(); ++it) {
      UInt vec_size = (*it)->syncAddedElements(this->nb_added_elements);
      AKANTU_DEBUG_ASSERT(vec_size == this->size,
			  "Synchronized arrays do not have the same length" <<
			  "(may be a double synchronization)");
    }
    this->nb_added_elements = 0;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<class T>
UInt SynchronizedArray<T>::syncDeletedElements(std::vector<UInt> & del_elements) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(nb_added_elements == 0 && deleted_elements.size() == 0,
		      "Cannot sync with a SynchronizedArray if it has already been modified");

  std::vector<UInt>::const_iterator it;
  for (it = del_elements.begin(); it != del_elements.end(); ++it) {
    erase(*it);
  }
  syncElements(_deleted);

  AKANTU_DEBUG_OUT();
  return this->size;
}

/* -------------------------------------------------------------------------- */
template<class T>
UInt SynchronizedArray<T>::syncAddedElements(UInt nb_add_elements) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(nb_added_elements == 0 && deleted_elements.size() == 0,
		      "Cannot sync with a SynchronizedArray if it has already been modified");

  for (UInt i=0; i<nb_add_elements; ++i) {
    push_back(this->default_value);
  }
  syncElements(_added);

  AKANTU_DEBUG_OUT();
  return this->size;
}

/* -------------------------------------------------------------------------- */
template<typename T>
void SynchronizedArray<T>::registerDependingArray(SynchronizedArrayBase & array) {
  AKANTU_DEBUG_IN();

  this->depending_arrays.push_back(&array);
  array.syncAddedElements(this->size);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<typename T>
void SynchronizedArray<T>::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();

  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "SynchronizedArray<" << debug::demangle(typeid(T).name()) << "> [" << std::endl;
  stream << space << " + default_value     : " << this->default_value << std::endl;
  stream << space << " + nb_added_elements : " << this->nb_added_elements << std::endl;
  stream << space << " + deleted_elements  : ";
  for(std::vector<UInt>::const_iterator it = this->deleted_elements.begin();
      it != this->deleted_elements.end();
      ++it)
    stream << *it << " ";
  stream << std::endl;

  stream << space << " + depending_arrays : ";
  for (std::vector<SynchronizedArrayBase *>::const_iterator it = this->depending_arrays.begin();
       it!=this->depending_arrays.end();
       ++it)
    stream << (*it)->getID() << " ";
  stream << std::endl;

  Array<T>::printself(stream, indent+1);

  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<typename T>
void SynchronizedArray<T>::dumpRestartFile(std::string file_name) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(nb_added_elements == 0 && deleted_elements.size() == 0,
		      "Restart File for SynchronizedArray " << this->id <<
		      " should not be dumped as it is not synchronized yet");

  std::stringstream name;
  name << file_name << "-" << this->restart_name << ".rs";

  std::ofstream out_restart;
  out_restart.open(name.str().c_str());

  out_restart << this->size << " " << this->nb_component << std::endl;
  Real size_comp = this->size * this->nb_component;
  for (UInt i=0; i<size_comp; ++i)
    out_restart << std::setprecision(12) << this->values[i] << " ";
  //    out_restart << std::hex << std::setprecision(12) << this->values[i] << " ";
  out_restart <<  std::endl;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<typename T>
void SynchronizedArray<T>::readRestartFile(std::string file_name) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(nb_added_elements == 0 && deleted_elements.size() == 0,
		      "Restart File for SynchronizedArray " << this->id <<
		      " should not be read as it is not synchronized yet");

  std::stringstream name;
  name << file_name << "-" << this->restart_name << ".rs";
  std::ifstream infile;
  infile.open(name.str().c_str());

  std::string line;

  // get size and nb_component info
  AKANTU_DEBUG_ASSERT(infile.good(), "Could not read restart file for " <<
		      "SynchronizedArray " << this->id);
  getline(infile, line);
  std::stringstream size_comp(line);
  size_comp >> this->size;
  size_comp >> this->nb_component;

  // get elements in array
  getline(infile, line);
  std::stringstream data(line);
  for (UInt i=0; i<this->size * this->nb_component; ++i) {
    AKANTU_DEBUG_ASSERT(!data.eof(), "Read SynchronizedArray " << this->id <<
			" got to the end of the file before having read all data!");
    data >> this->values[i];
    //    data >> std::hex >> this->values[i];
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template class SynchronizedArray<Real>;
template class SynchronizedArray<UInt>;
template class SynchronizedArray<Int>;
template class SynchronizedArray<bool>;

__END_AKANTU__
