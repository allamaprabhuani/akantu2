/**
 * @file   communication_buffer_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Apr  6 22:20:37 2011
 *
 * @brief
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


/* -------------------------------------------------------------------------- */
template<typename T>
inline CommunicationBuffer & CommunicationBuffer::operator<< (const T & to_pack) {
  T * tmp = reinterpret_cast<T *>(ptr_current);
  *tmp = to_pack;
  ptr_current += sizeof(T);
  return *this;
}


/* -------------------------------------------------------------------------- */
template<typename T>
inline CommunicationBuffer & CommunicationBuffer::operator>> (T & to_unpack) {
  T * tmp = reinterpret_cast<T *>(ptr_current);
  to_unpack = *tmp;
  ptr_current += sizeof(T);
  return *this;
}


/* -------------------------------------------------------------------------- */
/* Specialization                                                             */
/* -------------------------------------------------------------------------- */

/**
 * types::Vector
 */

/* -------------------------------------------------------------------------- */
template<typename T>
inline CommunicationBuffer & CommunicationBuffer::operator<< (const types::Vector<T> & to_pack) {
  UInt size = to_pack.size() * sizeof(T);
  memcpy(ptr_current, to_pack.storage(), size);
  ptr_current += size;
  return *this;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline CommunicationBuffer & CommunicationBuffer::operator>> (types::Vector<T> & to_unpack) {
  UInt size = to_unpack.size() * sizeof(T);
  memcpy(to_unpack.storage(), ptr_current, size);
  ptr_current += size;
  return *this;
}

/**
 * types::Matrix
 */

/* -------------------------------------------------------------------------- */
template<> inline CommunicationBuffer &
CommunicationBuffer::operator<< <types::Matrix> (const types::Matrix & to_pack) {
  UInt size = to_pack.size() * sizeof(Real);
  memcpy(ptr_current, to_pack.storage(), size);
  ptr_current += size;
  return *this;
}

/* -------------------------------------------------------------------------- */
template<> inline CommunicationBuffer &
CommunicationBuffer::operator>> <types::Matrix> (types::Matrix & to_unpack) {
  UInt size = to_unpack.size() * sizeof(Real);
  memcpy(to_unpack.storage(), ptr_current, size);
  ptr_current += size;
  return *this;
}
