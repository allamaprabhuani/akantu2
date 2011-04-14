/**
 * @file   real_static_communicator.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Apr  7 12:19:34 2011
 *
 * @brief  empty class just for inheritance
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

/* -------------------------------------------------------------------------- */


#ifndef __AKANTU_REAL_STATIC_COMMUNICATOR_HH__
#define __AKANTU_REAL_STATIC_COMMUNICATOR_HH__

__BEGIN_AKANTU__

class CommunicationRequest {
public:
  CommunicationRequest(UInt source, UInt dest);
  virtual ~CommunicationRequest();

  virtual void printself(std::ostream & stream, int indent = 0) const;

  AKANTU_GET_MACRO(Source, source, UInt);
  AKANTU_GET_MACRO(Destination, destination, UInt);
private:
  UInt source;
  UInt destination;
  UInt id;
  static UInt counter;
};

class StaticCommunicator;

class RealStaticCommunicator {
public:
  RealStaticCommunicator(__attribute__ ((unused)) int * argc,
			 __attribute__ ((unused)) char *** argv) {
    prank = -1;
    psize = -1;
  };
  virtual ~RealStaticCommunicator() { };

  friend class StaticCommunicator;
protected:
  Int prank;

  Int psize;
};

__END_AKANTU__

#endif /* __AKANTU_REAL_STATIC_COMMUNICATOR_HH__ */
