/**
 * @file   static_communicator_dummy.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 19 15:34:09 2010
 *
 * @brief  Class handling the parallel communications
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

#ifndef __AKANTU_STATIC_COMMUNICATOR_DUMMY_HH__
#define __AKANTU_STATIC_COMMUNICATOR_DUMMY_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class StaticCommunicatorDummy : public StaticCommunicator {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  virtual ~StaticCommunicatorDummy() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void send(__attribute__ ((unused)) UInt * buffer,
		    __attribute__ ((unused)) Int size,
		    __attribute__ ((unused)) Int receiver,
		    __attribute__ ((unused)) Int tag) {};
  virtual void send(__attribute__ ((unused)) Real * buffer,
		    __attribute__ ((unused)) Int size,
		    __attribute__ ((unused)) Int receiver,
		    __attribute__ ((unused)) Int tag) {};

  virtual void receive(__attribute__ ((unused)) UInt * buffer,
		       __attribute__ ((unused)) Int size,
		       __attribute__ ((unused)) Int sender,
		       __attribute__ ((unused)) Int tag) {};
  virtual void receive(__attribute__ ((unused)) Real * buffer,
		       __attribute__ ((unused)) Int size,
		       __attribute__ ((unused)) Int sender,
		       __attribute__ ((unused)) Int tag) {};

  virtual CommunicationRequest * asyncSend(__attribute__ ((unused)) UInt * buffer,
					   __attribute__ ((unused)) Int size,
					   __attribute__ ((unused)) Int receiver,
					   __attribute__ ((unused)) Int tag) {
    return new CommunicationRequest(0, 0);
  };

  virtual CommunicationRequest * asyncSend(__attribute__ ((unused)) Real * buffer,
					   __attribute__ ((unused)) Int size,
					   __attribute__ ((unused)) Int receiver,
					   __attribute__ ((unused)) Int tag) {
    return new CommunicationRequest(0, 0);
  };

  virtual CommunicationRequest * asyncReceive(__attribute__ ((unused)) UInt * buffer,
					      __attribute__ ((unused)) Int size,
					      __attribute__ ((unused)) Int sender,
					      __attribute__ ((unused)) Int tag) {
    return new CommunicationRequest(0, 0);
  };

  virtual CommunicationRequest * asyncReceive(__attribute__ ((unused)) Real * buffer,
					      __attribute__ ((unused)) Int size,
					      __attribute__ ((unused)) Int sender,
					      __attribute__ ((unused)) Int tag) {
    return new CommunicationRequest(0, 0);
  };

  virtual bool testRequest(__attribute__ ((unused)) CommunicationRequest * request) { return true; };


  virtual void wait(__attribute__ ((unused)) CommunicationRequest * request) {};

  virtual void waitAll(__attribute__ ((unused)) std::vector<CommunicationRequest *> & requests) {};

  virtual void barrier() {};

  virtual void allReduce(__attribute__ ((unused)) Real * values,
			 __attribute__ ((unused)) Int nb_values,
			 __attribute__ ((unused)) const SynchronizerOperation & op) {};
  virtual void allReduce(__attribute__ ((unused)) UInt * values,
			 __attribute__ ((unused)) Int nb_values,
			 __attribute__ ((unused)) const SynchronizerOperation & op) {};

  inline void gather(__attribute__ ((unused)) Real * values,
		     __attribute__ ((unused)) Int nb_values,
		     __attribute__ ((unused)) Int root = 0) {};
  inline void gather(__attribute__ ((unused)) UInt * values,
		     __attribute__ ((unused)) Int nb_values,
		     __attribute__ ((unused)) Int root = 0) {};
  inline void gather(__attribute__ ((unused)) Int * values,
		     __attribute__ ((unused)) Int nb_values,
		     __attribute__ ((unused)) Int root = 0) {};

  inline void gatherv(__attribute__ ((unused)) Real * values,
		      __attribute__ ((unused)) Int * nb_values,
		      __attribute__ ((unused)) Int root = 0) {};
  inline void gatherv(__attribute__ ((unused)) UInt * values,
		      __attribute__ ((unused)) Int * nb_values,
		      __attribute__ ((unused)) Int root = 0) {};
  inline void gatherv(__attribute__ ((unused)) Int * values,
		      __attribute__ ((unused)) Int * nb_values,
		      __attribute__ ((unused)) Int root = 0) {};

  inline void broadcast(__attribute__ ((unused)) Real * values,
			__attribute__ ((unused)) Int nb_values,
			__attribute__ ((unused)) Int root = 0) {};
  inline void broadcast(__attribute__ ((unused)) UInt * values,
			__attribute__ ((unused)) Int nb_values,
			__attribute__ ((unused)) Int root = 0) {};
  inline void broadcast(__attribute__ ((unused)) Int * values,
			__attribute__ ((unused)) Int nb_values,
			__attribute__ ((unused)) Int root = 0) {};


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  virtual Int getNbProc() const { return 1; };
  virtual Int whoAmI() const { return 0; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
};

__END_AKANTU__

#endif /* __AKANTU_STATIC_COMMUNICATOR_DUMMY_HH__ */
