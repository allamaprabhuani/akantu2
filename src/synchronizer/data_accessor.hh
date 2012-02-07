/**
 * @file   data_accessor.hh
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @date   Wed Jun 15 14:43:23 2011
 *
 * @brief  Interface of accessors for pack_unpack system
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


#ifndef __AKANTU_DATA_ACCESSOR_HH__
#define __AKANTU_DATA_ACCESSOR_HH__


/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "communication_buffer.hh"
/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

class DataAccessor {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  DataAccessor();
  virtual ~DataAccessor();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /**
   * @brief get  the number of  data to send  for a given akantu::Element  and a
   * given akantu::SynchronizationTag
   */
  virtual UInt getNbDataToPack(__attribute__((unused)) const Element & element,
                               __attribute__((unused)) SynchronizationTag tag) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /**
   * @brief get  the number of  data to send  for a given
   * akantu::SynchronizationTag
   */
  virtual UInt getNbDataToPack(__attribute__((unused)) SynchronizationTag tag) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /**
   * @brief get the number of data  to receive for a given akantu::Element and a
   * given akantu::SynchronizationTag
   */
  virtual UInt getNbDataToUnpack(__attribute__((unused)) const Element & element,
				 __attribute__((unused)) SynchronizationTag tag) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /**
   * @brief get the number of data  to receive for a given
   * akantu::SynchronizationTag
   */
  virtual UInt getNbDataToUnpack(__attribute__((unused)) SynchronizationTag tag) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /**
   * @brief   pack  the   data  for   a  given   akantu::Element  and   a  given
   * akantu::SynchronizationTag
   */
  virtual void packData(__attribute__((unused)) CommunicationBuffer & buffer,
                        __attribute__((unused)) const Element & element,
                        __attribute__((unused)) SynchronizationTag tag) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /**
   * @brief   pack  the   data  for   a  given  index  and   a  given
   * akantu::SynchronizationTag
   */
  virtual void packData(__attribute__((unused)) CommunicationBuffer & buffer,
			__attribute__((unused)) const UInt index,
                        __attribute__((unused)) SynchronizationTag tag) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /**
   * @brief   unpack  the   data  for   a  given   akantu::Element  and   a  given
   * akantu::SynchronizationTag
   */
  virtual void unpackData(__attribute__((unused)) CommunicationBuffer & buffer,
                          __attribute__((unused)) const Element & element,
                          __attribute__((unused)) SynchronizationTag tag) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /**
   * @brief   unpack  the   data  for   a  given  index  and   a  given
   * akantu::SynchronizationTag
   */
  virtual void unpackData(__attribute__((unused)) CommunicationBuffer & buffer,
                          __attribute__((unused)) const UInt index,
                          __attribute__((unused)) SynchronizationTag tag) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

};


/* -------------------------------------------------------------------------- */
/* __aka_inline__ functions                                                           */
/* -------------------------------------------------------------------------- */

// #include "data_accessor_inline_impl.cc"

// /// standard output stream operator
// __aka_inline__ std::ostream & operator <<(std::ostream & stream, const DataAccessor & _this)
// {
//   _this.printself(stream);
//   return stream;
// }

__END_AKANTU__

#endif /* __AKANTU_DATA_ACCESSOR_HH__ */
