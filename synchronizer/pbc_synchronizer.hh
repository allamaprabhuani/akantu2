/**
 * @file   pbc_synchronizer.hh
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @date   Wed Jun 15 12:56:28 2011
 *
 * @brief  Dofs Synchronizer for periodic boundary condition 
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

#ifndef __AKANTU_PBC_SYNCHRONIZER_HH__
#define __AKANTU_PBC_SYNCHRONIZER_HH__


class PBCSynchronizer : public Synchronizer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  PBCSynchronizer();
  virtual ~PBCSynchronizer();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// get a mesh and create the local dof links due to pbc 
  static PBCSynchronizer * createPBCSynchronizer(Mesh & mesh,
						 SynchronizerID id = "pbc_synch",
						 MemoryID memory_id = 0);
  
  /* ------------------------------------------------------------------------ */
  /* Inherited from Synchronizer                                              */
  /* ------------------------------------------------------------------------ */
  
  /// synchronize ghosts
  void synchronize(GhostSynchronizationTag tag);

  /// asynchronous synchronization of ghosts
  void asynchronousSynchronize(GhostSynchronizationTag tag){}

  /// wait end of asynchronous synchronization of ghosts
  void waitEndSynchronize(GhostSynchronizationTag tag){}

  /// do a all reduce operation
  void allReduce(Real * values, UInt nb_values, const SynchronizerOperation & op){};


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
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "pbc_synchronizer_inline_impl.cc"

// /// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const PBCSynchronizer & _this)
// {
//   _this.printself(stream);
//   return stream;
// }


#endif /* __AKANTU_PBC_SYNCHRONIZER_HH__ */
