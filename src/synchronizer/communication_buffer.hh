/**
 * @file   communication_buffer.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Thu Apr 14 18:22:18 2011
 *
 * @brief  Buffer for packing and unpacking data
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
#include "aka_common.hh"
#include "aka_array.hh"

#ifndef __AKANTU_COMMUNICATION_BUFFER_HH__
#define __AKANTU_COMMUNICATION_BUFFER_HH__

__BEGIN_AKANTU__

template<bool is_static = true>
class CommunicationBufferTemplated {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  CommunicationBufferTemplated(UInt size = 0) : buffer(size, 1) {
    ptr_pack = buffer.storage();
    ptr_unpack = buffer.storage();
  };

  virtual ~CommunicationBufferTemplated() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// reset to "empty"
  inline void reset();

  /// resize the internal buffer
  inline void resize(UInt size);

  /// clear buffer context
  inline void clear();

private:
  inline void packResize(UInt size);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  inline char * storage() { return buffer.values; };

  /* ------------------------------------------------------------------------ */
  /* Operators                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// printing tool
  template <typename T> inline std::string extractStream(UInt packet_size);

  /// packing data
  template<typename T>
  inline CommunicationBufferTemplated & operator<< (const T & to_pack);

  template<typename T>
  inline CommunicationBufferTemplated & operator<< (const Vector<T> & to_pack);

  template<typename T>
  inline CommunicationBufferTemplated & operator<< (const Matrix<T> & to_pack);

  template<typename T>
  inline CommunicationBufferTemplated & operator<< (const std::vector<T> & to_pack);

  /// unpacking data
  template<typename T>
  inline CommunicationBufferTemplated & operator>> (T & to_unpack);

  template<typename T>
  inline CommunicationBufferTemplated & operator>> (Vector<T> & to_unpack);

  template<typename T>
  inline CommunicationBufferTemplated & operator>> (Matrix<T> & to_unpack);

  template<typename T>
  inline CommunicationBufferTemplated & operator>> (std::vector<T> & to_unpack);

  inline CommunicationBufferTemplated & operator<< (const std::string & to_pack);
  inline CommunicationBufferTemplated & operator>> (std::string & to_unpack);

private:

  template<typename T>
  inline void packIterable (T & to_pack);

  template<typename T>
  inline void unpackIterable (T & to_pack);

  /* ------------------------------------------------------------------------ */
  /* Accessor                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  /// return the size in bytes of the stored values
  inline UInt getPackedSize(){ return ptr_pack - buffer.values; };
  /// return the size in bytes of data left to be unpacked
  inline UInt getLeftToUnpack(){ return buffer.getSize() - (ptr_unpack - buffer.values); };
  /// return the global size allocated
  inline UInt getSize(){ return buffer.getSize(); };

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


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "communication_buffer_inline_impl.cc"
#endif

typedef CommunicationBufferTemplated<true> CommunicationBuffer;
typedef CommunicationBufferTemplated<false> DynamicCommunicationBuffer;


__END_AKANTU__

#endif /* __AKANTU_COMMUNICATION_BUFFER_HH__ */
