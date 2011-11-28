/**
 * @file   communication_buffer.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Apr  6 21:23:16 2011
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

#ifndef __AKANTU_COMMUNICATION_BUFFER_HH__
#define __AKANTU_COMMUNICATION_BUFFER_HH__

__BEGIN_AKANTU__

class CommunicationBuffer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  CommunicationBuffer(UInt size = 0) : buffer(size, 1) {
    ptr_pack = buffer.values;
    ptr_unpack = buffer.values;
  };

  virtual ~CommunicationBuffer() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// reset to "empty"
  __aka_inline__ void reset();

  /// resize the internal buffer
  __aka_inline__ void resize(UInt size);

  /// clear buffer context
  __aka_inline__ void clear();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  __aka_inline__ char * storage() { return buffer.values; };

  /* ------------------------------------------------------------------------ */
  /* Operators                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// printing tool
  template <typename T> __aka_inline__ std::string extractStream(UInt packet_size);

  /// packing data
  template<typename T>
  __aka_inline__ CommunicationBuffer & operator<< (const T & to_pack);

  template<typename T>
  __aka_inline__ CommunicationBuffer & operator<< (const types::Vector<T> & to_pack);


  /// unpacking data
  template<typename T>
  __aka_inline__ CommunicationBuffer & operator>> (T & to_unpack);

  template<typename T>
  __aka_inline__ CommunicationBuffer & operator>> (types::Vector<T> & to_unpack);

  /* ------------------------------------------------------------------------ */
  /* Accessor                                                                 */
  /* ------------------------------------------------------------------------ */

  /// return the size in bytes of the stored values
  __aka_inline__ UInt getPackedSize(){return ptr_pack-buffer.values;}; 
  /// return the size in bytes of data left to be unpacked
  __aka_inline__ UInt getLeftToUnpack(){return buffer.getSize()-(ptr_unpack-buffer.values);}; 
  /// return the global size allocated
  __aka_inline__ UInt getSize(){return buffer.getSize();}; 

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// current position for packing
  char * ptr_pack;

  /// current position for unpacking
  char * ptr_unpack;

  /// storing buffer
  Vector<char> buffer;
};


/* -------------------------------------------------------------------------- */
/* __aka_inline__ functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "communication_buffer_inline_impl.cc"
#endif


__END_AKANTU__

#endif /* __AKANTU_COMMUNICATION_BUFFER_HH__ */
