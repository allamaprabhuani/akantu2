/**
 * @file   mesh_partition_scotch.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Aug 16 13:17:22 2010
 *
 * @brief  mesh partitioning based on libScotch
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

#ifndef __AKANTU_MESH_PARTITION_SCOTCH_HH__
#define __AKANTU_MESH_PARTITION_SCOTCH_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh_partition.hh"

#if ! defined(AKANTU_USE_PTSCOTCH)
# ifndef AKANTU_SCOTCH_NO_EXTERN
extern "C" {
# endif //AKANTU_SCOTCH_NO_EXTERN
# include <scotch.h>
# ifndef AKANTU_SCOTCH_NO_EXTERN
}
# endif //AKANTU_SCOTCH_NO_EXTERN
#else //AKANTU_USE_PTSCOTCH
# include <ptscotch.h>
#endif //AKANTU_USE_PTSCOTCH
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class MeshPartitionScotch : public MeshPartition {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MeshPartitionScotch(const Mesh & mesh, UInt spatial_dimension,
		      const MemoryID & memory_id = 0);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void partitionate(UInt nb_part);

  virtual void reorder();

  /// function to print the contain of the class
  //virtual void printself(std::ostream & stream, int indent = 0) const;

private:
  SCOTCH_Mesh * createMesh();

  void destroyMesh(SCOTCH_Mesh * meshptr);

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

// #include "mesh_partition_scotch_inline_impl.cc"

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const MeshPartitionScotch & _this)
// {
//   _this.printself(stream);
//   return stream;
// }


__END_AKANTU__

#endif /* __AKANTU_MESH_PARTITION_SCOTCH_HH__ */
