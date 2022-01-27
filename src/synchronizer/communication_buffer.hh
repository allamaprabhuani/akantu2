/**
 * @file   communication_buffer.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Dec 11 2019
 *
 * @brief  Buffer for packing and unpacking data
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
#include "aka_array.hh"
#include "aka_common.hh"
#include "element.hh"
/* -------------------------------------------------------------------------- */
#include <array>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_COMMUNICATION_BUFFER_HH_
#define AKANTU_COMMUNICATION_BUFFER_HH_

namespace akantu {

template <bool is_static = true> class CommunicationBufferTemplated {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  explicit CommunicationBufferTemplated(std::size_t size)
      : buffer(size, 1, char()) {
    ptr_pack = buffer.data();
    ptr_unpack = buffer.data();
  };

  CommunicationBufferTemplated() : CommunicationBufferTemplated(0) {}

  CommunicationBufferTemplated(const CommunicationBufferTemplated & other) =
      delete;
  CommunicationBufferTemplated &
  operator=(const CommunicationBufferTemplated & other) = delete;

  CommunicationBufferTemplated(CommunicationBufferTemplated && other) noexcept =
      default;

  virtual ~CommunicationBufferTemplated() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// reset to "empty"
  inline void reset();

  /// resize the internal buffer do not allocate on dynamic buffers
  inline void resize(std::size_t size);

  /// resize the internal buffer allocate always
  inline void reserve(std::size_t size);

  /// clear buffer context
  inline void zero();

private:
  inline void packResize(std::size_t size);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  [[deprecated("use data instead to be stl compatible")]] inline char *
  storage() {
    return buffer.data();
  };
  [[deprecated("use data instead to be stl compatible")]] inline const char *
  storage() const {
    return buffer.data();
  };

  inline char * data() { return buffer.data(); };
  inline const char * data() const { return buffer.data(); };

  /* ------------------------------------------------------------------------ */
  /* Operators                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// printing tool
  template <typename T>
  inline std::string extractStream(std::size_t block_size);

  /// packing data
  template <typename T,
            std::enable_if_t<std::is_standard_layout<T>::value and
                             not aka::is_tensor<T>::value> * = nullptr>
  inline CommunicationBufferTemplated & operator<<(const T & to_pack);

  template <typename T, std::enable_if_t<aka::is_tensor<T>::value> * = nullptr>
  inline CommunicationBufferTemplated & operator<<(const T & to_pack);

  template <typename T>
  inline CommunicationBufferTemplated &
  operator<<(const std::vector<T> & to_pack);

  inline CommunicationBufferTemplated & operator<<(const std::string & to_pack);

  /// unpacking data
  template <typename T,
            std::enable_if_t<std::is_standard_layout<T>::value and
                             not aka::is_tensor<T>::value> * = nullptr>
  inline CommunicationBufferTemplated & operator>>(T & to_unpack);

  template <typename T, std::enable_if_t<aka::is_tensor<T>::value> * = nullptr>
  inline CommunicationBufferTemplated & operator>>(T & to_pack);

  template <typename T>
  inline CommunicationBufferTemplated & operator>>(std::vector<T> & to_unpack);

  inline CommunicationBufferTemplated & operator>>(std::string & to_unpack);

private:
  template <typename T> inline void packIterable(T & to_pack);
  template <typename T> inline void unpackIterable(T & to_unpack);

  /* ------------------------------------------------------------------------ */
  /* Accessor                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  template <typename T,
            std::enable_if_t<std::is_standard_layout<T>::value and
                             not aka::is_tensor<T>::value> * = nullptr>
  static inline std::size_t sizeInBuffer(const T & data);

  template <typename T, std::enable_if_t<aka::is_tensor<T>::value> * = nullptr>
  static inline std::size_t sizeInBuffer(const T & data);

  template <typename T>
  static inline std::size_t sizeInBuffer(const std::vector<T> & data);

  static inline std::size_t sizeInBuffer(const std::string & data);

  /// return the size in bytes of the stored values
  inline std::size_t getPackedSize() const { return ptr_pack - buffer.data(); };
  /// return the size in bytes of data left to be unpacked
  inline std::size_t getLeftToUnpack() const {
    return buffer.size() - (ptr_unpack - buffer.data());
  };
  /// return the global size allocated
  inline std::size_t size() const { return buffer.size(); };

  /// is the buffer empty
  inline bool empty() const {
    return (getPackedSize() == 0) and (getLeftToUnpack() == 0);
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// current position for packing
  char * ptr_pack;

  /// current position for unpacking
  char * ptr_unpack;

  /// storing buffer
  Array<char> buffer;
};

using CommunicationBuffer = CommunicationBufferTemplated<true>;
using DynamicCommunicationBuffer = CommunicationBufferTemplated<false>;

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "communication_buffer_inline_impl.hh"

#endif /* AKANTU_COMMUNICATION_BUFFER_HH_ */
