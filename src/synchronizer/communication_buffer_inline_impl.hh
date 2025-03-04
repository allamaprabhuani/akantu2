/**
 * Copyright (©) 2011-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "communication_buffer.hh"
#include <cstring>
/* -------------------------------------------------------------------------- */
namespace akantu {

// NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic,
//             cppcoreguidelines-pro-type-reinterpret-cast,
//             cppcoreguidelines-narrowing-conversions)
/* -------------------------------------------------------------------------- */
template <bool is_static>
template <typename T, std::enable_if_t<std::is_standard_layout_v<T> and
                                       not aka::is_tensor_v<T>> *>
inline std::size_t
CommunicationBufferTemplated<is_static>::sizeInBuffer(const T & /*unused*/) {
  return sizeof(T);
}

template <bool is_static>
template <typename Tensor, std::enable_if_t<aka::is_tensor_v<Tensor>> *>
inline std::size_t
CommunicationBufferTemplated<is_static>::sizeInBuffer(const Tensor & data) {
  std::size_t size = data.size() * sizeof(typename Tensor::Scalar);
  return size;
}

template <bool is_static>
template <typename T>
inline std::size_t CommunicationBufferTemplated<is_static>::sizeInBuffer(
    const std::vector<T> & data) {
  std::size_t size = data.size() * sizeof(T) + sizeof(size_t);
  return size;
}

template <bool is_static>
inline std::size_t CommunicationBufferTemplated<is_static>::sizeInBuffer(
    const std::string & data) {
  std::size_t size =
      data.size() * sizeof(std::string::value_type) + sizeof(size_t);
  return size;
}

/* -------------------------------------------------------------------------- */
template <bool is_static>
inline void
CommunicationBufferTemplated<is_static>::packResize(std::size_t size) {
  if (not is_static) {
    char * values = buffer.data();
    auto nb_packed = ptr_pack - values;

    if (std::size_t(buffer.size()) > nb_packed + size) {
      return;
    }

    buffer.resize(nb_packed + size);
    ptr_pack = buffer.data() + nb_packed;
    ptr_unpack = buffer.data() + (ptr_unpack - values);
  }
}

/* -------------------------------------------------------------------------- */
template <bool is_static>
template <typename T, std::enable_if_t<std::is_standard_layout_v<T> and
                                       not aka::is_tensor_v<T>> *>
inline CommunicationBufferTemplated<is_static> &
CommunicationBufferTemplated<is_static>::operator<<(const T & to_pack) {
  std::size_t size = sizeInBuffer(to_pack);
  packResize(size);
  AKANTU_DEBUG_ASSERT(
      (buffer.data() + buffer.size()) >= (ptr_pack + size),
      "Packing too much data in the CommunicationBufferTemplated");
  std::memcpy(ptr_pack, reinterpret_cast<const char *>(&to_pack), size);
  ptr_pack += size;
  return *this;
}

/* -------------------------------------------------------------------------- */
template <bool is_static>
template <typename T, std::enable_if_t<std::is_standard_layout_v<T> and
                                       not aka::is_tensor_v<T>> *>
inline CommunicationBufferTemplated<is_static> &
CommunicationBufferTemplated<is_static>::operator>>(T & to_unpack) {
  std::size_t size = sizeInBuffer(to_unpack);

  alignas(alignof(T)) std::array<char, sizeof(T)> aligned_ptr;
  std::memcpy(aligned_ptr.data(), ptr_unpack, size);

  auto * tmp = reinterpret_cast<T *>(aligned_ptr.data());
  AKANTU_DEBUG_ASSERT(
      (buffer.data() + buffer.size()) >= (ptr_unpack + size),
      "Unpacking too much data in the CommunicationBufferTemplated");
  to_unpack = *tmp;
  // memcpy(reinterpret_cast<char *>(&to_unpack), ptr_unpack, size);
  ptr_unpack += size;
  return *this;
}

/* -------------------------------------------------------------------------- */
template <bool is_static>
template <typename Tensor, std::enable_if_t<aka::is_tensor_v<Tensor>> *>
inline CommunicationBufferTemplated<is_static> &
CommunicationBufferTemplated<is_static>::operator<<(const Tensor & to_pack) {
  std::size_t size = sizeInBuffer(to_pack);
  packResize(size);
  AKANTU_DEBUG_ASSERT(
      (buffer.data() + buffer.size()) >= (ptr_pack + size),
      "Packing too much data in the CommunicationBufferTemplated");
  std::memcpy(ptr_pack, to_pack.data(), size);
  ptr_pack += size;
  return *this;
}

/* --------------------------------------------------------------------------
 */
template <bool is_static>
template <typename Tensor, std::enable_if_t<aka::is_tensor_v<Tensor>> *>
inline CommunicationBufferTemplated<is_static> &
CommunicationBufferTemplated<is_static>::operator>>(Tensor & to_unpack) {
  std::size_t size = sizeInBuffer(to_unpack);
  AKANTU_DEBUG_ASSERT(
      (buffer.data() + buffer.size()) >= (ptr_unpack + size),
      "Unpacking too much data in the CommunicationBufferTemplated");
  std::memcpy(to_unpack.data(), ptr_unpack, size);
  ptr_unpack += size;
  return *this;
}

/* --------------------------------------------------------------------------
 */
template <bool is_static>
template <typename T>
inline void CommunicationBufferTemplated<is_static>::packIterable(T & to_pack) {
  operator<<(std::size_t(to_pack.size()));
  auto it = to_pack.begin();
  auto end = to_pack.end();
  for (; it != end; ++it) {
    operator<<(*it);
  }
}

/* --------------------------------------------------------------------------
 */
template <bool is_static>
template <typename T>
inline void
CommunicationBufferTemplated<is_static>::unpackIterable(T & to_unpack) {
  std::size_t size{};
  operator>>(size);

  to_unpack.resize(size);
  auto it = to_unpack.begin();
  auto end = to_unpack.end();
  for (; it != end; ++it) {
    operator>>(*it);
  }
}

/**
 * std::vector<T>
 */
/* --------------------------------------------------------------------------
 */

template <bool is_static>
template <typename T>
inline CommunicationBufferTemplated<is_static> &
CommunicationBufferTemplated<is_static>::operator<<(
    const std::vector<T> & to_pack) {
  packIterable(to_pack);
  return *this;
}

/* --------------------------------------------------------------------------
 */
template <bool is_static>
template <typename T>
inline CommunicationBufferTemplated<is_static> &
CommunicationBufferTemplated<is_static>::operator>>(
    std::vector<T> & to_unpack) {
  unpackIterable(to_unpack);
  return *this;
}

/**
 * std::string
 */
/* --------------------------------------------------------------------------
 */

template <bool is_static>
inline CommunicationBufferTemplated<is_static> &
CommunicationBufferTemplated<is_static>::operator<<(
    const std::string & to_pack) {
  packIterable(to_pack);
  return *this;
}

/* --------------------------------------------------------------------------
 */
template <bool is_static>
inline CommunicationBufferTemplated<is_static> &
CommunicationBufferTemplated<is_static>::operator>>(std::string & to_unpack) {
  unpackIterable(to_unpack);
  return *this;
}

/* --------------------------------------------------------------------------
 */
template <bool is_static>
template <typename T>
inline std::string
CommunicationBufferTemplated<is_static>::extractStream(std::size_t block_size) {
  std::stringstream str;
  auto * ptr = reinterpret_cast<T *>(buffer.data());
  auto sz = buffer.size() / sizeof(T);
  auto sz_block = block_size / sizeof(T);

  std::size_t n_block = 0;
  for (std::size_t i = 0; i < sz; ++i) {
    if (i % sz_block == 0) {
      str << "\n" << n_block << " ";
      ++n_block;
    }
    str << *ptr << " ";
    ++ptr;
  }
  return str.str();
}

/* --------------------------------------------------------------------------
 */
template <bool is_static>
inline void CommunicationBufferTemplated<is_static>::resize(std::size_t size) {
  if (!is_static) {
    buffer.resize(0, 0);
  } else {
    buffer.resize(size, 0);
  }
  reset();
#ifndef AKANTU_NDEBUG
  zero();
#endif
}

/* --------------------------------------------------------------------------
 */
template <bool is_static>
inline void CommunicationBufferTemplated<is_static>::reserve(std::size_t size) {
  char * values = buffer.data();
  std::size_t nb_packed = ptr_pack - values;

  buffer.resize(size);
  ptr_pack = buffer.data() + nb_packed;
  ptr_unpack = buffer.data() + (ptr_unpack - values);
}

/* --------------------------------------------------------------------------
 */
template <bool is_static>
inline void CommunicationBufferTemplated<is_static>::zero() {
  buffer.zero();
}

/* --------------------------------------------------------------------------
 */
template <bool is_static>
inline void CommunicationBufferTemplated<is_static>::reset() {
  ptr_pack = buffer.data();
  ptr_unpack = buffer.data();
}

// NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic,
//           cppcoreguidelines-pro-type-reinterpret-cast,
//           cppcoreguidelines-narrowing-conversions)
} // namespace akantu
