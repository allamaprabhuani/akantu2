/**
 * @file   boundary_inline_impl.cc
 *
 * @author Dana Christen <dana.christen@gmail.com>
 *
 * @date   Wed Mar 06 09:30:00 2013
 *
 * @brief  Stores information relevent to the notion of domain boundary and surfaces.
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
inline Boundary::const_iterator Boundary::begin() const {
  return const_iterator(boundaries.begin());
}

/* -------------------------------------------------------------------------- */
inline Boundary::const_iterator Boundary::find(std::string name) const {
  return const_iterator(boundaries.find(name));
}

/* -------------------------------------------------------------------------- */
inline Boundary::const_iterator Boundary::end() const {
  return const_iterator(boundaries.end());
}

/* -------------------------------------------------------------------------- */
inline Boundary::iterator Boundary::begin() {
  return iterator(boundaries.begin());
}

/* -------------------------------------------------------------------------- */
inline Boundary::iterator Boundary::find(std::string name)  {
  return iterator(boundaries.find(name));
}

/* -------------------------------------------------------------------------- */
inline Boundary::iterator Boundary::end() {
  return iterator(boundaries.end());
}


/* -------------------------------------------------------------------------- */
inline const SubBoundary & Boundary::operator()(const std::string & name) const {
  BoundaryList::const_iterator iter = boundaries.find(name);
  if(iter == boundaries.end()) {
    AKANTU_EXCEPTION("No boundary named " << name << "!");
  }
  return *(iter->second);
}

/* -------------------------------------------------------------------------- */
template<typename T, typename container_iterator>
inline Boundary::internal_iterator<T, container_iterator>::internal_iterator(const container_iterator & ext_iter)
{
  iter = ext_iter;
}

/* -------------------------------------------------------------------------- */
template<typename T, typename container_iterator>
inline Boundary::internal_iterator<T, container_iterator>::internal_iterator(const Boundary::internal_iterator<T, container_iterator> & other)
{
  iter = other.iter;
}

/* -------------------------------------------------------------------------- */
template<typename T, typename container_iterator>
inline T & Boundary::internal_iterator<T, container_iterator>::operator*() const {
  return *(iter->second);
}

/* -------------------------------------------------------------------------- */
template<typename T, typename container_iterator>
inline T * Boundary::internal_iterator<T, container_iterator>::operator->() const {
  return iter->second;
}

/* -------------------------------------------------------------------------- */
template<typename T, typename container_iterator>
inline Boundary::internal_iterator<T, container_iterator> Boundary::internal_iterator<T, container_iterator>::operator++(int) {
  const_iterator tmp(*this);
  operator++();
  return tmp;
}

/* -------------------------------------------------------------------------- */
template<typename T, typename container_iterator>
inline Boundary::internal_iterator<T, container_iterator> & Boundary::internal_iterator<T, container_iterator>::operator++() {
  ++iter;
  return *this;
}



/* -------------------------------------------------------------------------- */
template<typename T, typename container_iterator>
inline bool Boundary::internal_iterator<T, container_iterator>::operator==(const internal_iterator<T, container_iterator> & other) const {
  return (this->iter == other.iter);
}

/* -------------------------------------------------------------------------- */
template<typename T, typename container_iterator>
inline bool Boundary::internal_iterator<T, container_iterator>::operator!=(const internal_iterator<T, container_iterator> & other) const {
  return (this->iter != other.iter);
}

/* -------------------------------------------------------------------------- */
template<typename T, typename container_iterator>
inline Boundary::internal_iterator<T, container_iterator> & Boundary::internal_iterator<T, container_iterator>::operator=(const internal_iterator<T, container_iterator> & other) {
  if(&other != this) {
    this->iter = other.iter;
  }
  return *this;
}


