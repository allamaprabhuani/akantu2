/**
 * @file   communication_buffer_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Thu Apr 14 18:22:18 2011
 *
 * @brief  CommunicationBuffer inline implementation
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
  T * tmp = reinterpret_cast<T *>(ptr_pack);
  AKANTU_DEBUG_ASSERT(buffer.storage() + buffer.getSize() >= ptr_pack + sizeof(T),
		      "Packing too much data in the CommunicationBuffer");
  *tmp = to_pack;
  ptr_pack += sizeof(T);
  return *this;
}


/* -------------------------------------------------------------------------- */
template<typename T>
inline CommunicationBuffer & CommunicationBuffer::operator>> (T & to_unpack) {
  T * tmp = reinterpret_cast<T *>(ptr_unpack);
  to_unpack = *tmp;
  ptr_unpack += sizeof(T);
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
  AKANTU_DEBUG_ASSERT(buffer.storage() + buffer.getSize() >= ptr_pack + size,
		      "Packing too much data in the CommunicationBuffer");
  memcpy(ptr_pack, to_pack.storage(), size);
  ptr_pack += size;
  return *this;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline CommunicationBuffer & CommunicationBuffer::operator>> (types::Vector<T> & to_unpack) {
  UInt size = to_unpack.size() * sizeof(T);
  memcpy(to_unpack.storage(), ptr_unpack, size);
  ptr_unpack += size;
  return *this;
}

/**
 * types::Matrix
 */

/* -------------------------------------------------------------------------- */
template<typename T>
inline CommunicationBuffer & CommunicationBuffer::operator<< (const types::Matrix<T> & to_pack) {
  UInt size = to_pack.size() * sizeof(Real);
  AKANTU_DEBUG_ASSERT(buffer.storage() + buffer.getSize() >= ptr_pack + size,
		      "Packing too much data in the CommunicationBuffer");
  memcpy(ptr_pack, to_pack.storage(), size);
  ptr_pack += size;
  return *this;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline CommunicationBuffer & CommunicationBuffer::operator>> (types::Matrix<T> & to_unpack) {
  UInt size = to_unpack.size() * sizeof(Real);
  memcpy(to_unpack.storage(), ptr_unpack, size);
  ptr_unpack += size;
  return *this;
}

/* -------------------------------------------------------------------------- */
template<typename T> inline std::string
CommunicationBuffer::extractStream(UInt block_size) {
  std::stringstream str;
  T * ptr = reinterpret_cast<T*>(buffer.values);
  UInt sz = buffer.getSize()/sizeof(T);
  UInt sz_block = block_size/sizeof(T);

  UInt n_block = 0;
  for (UInt i = 0; i < sz; ++i) {
    if (i% sz_block == 0) {
      str << std::endl << n_block << " ";
      ++n_block;
    }
    str << *ptr << " ";
    ++ptr;
  }
  return str.str();
}

/* -------------------------------------------------------------------------- */
inline void CommunicationBuffer::resize(UInt size) {
    buffer.resize(size);
    reset();
#ifndef AKANTU_NDEBUG
    clear();
#endif
}

/* -------------------------------------------------------------------------- */
inline void CommunicationBuffer::clear() {
  buffer.clear();
}

/* -------------------------------------------------------------------------- */
inline void CommunicationBuffer::reset() {
  ptr_pack = buffer.values;
  ptr_unpack = buffer.values;
}
