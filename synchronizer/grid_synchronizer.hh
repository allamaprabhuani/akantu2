/**
 * @file   grid_synchronizer.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Sep 16 15:05:52 2011
 *
 * @brief  synchronizer based in RegularGrid
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

#ifndef __AKANTU_GRID_SYNCHRONIZER_HH__
#define __AKANTU_GRID_SYNCHRONIZER_HH__

__BEGIN_AKANTU__

class GridSynchronizer : protected DistributedSynchronizer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  GridSynchronizer(ID & id = "grid_synchronizer", MemoryID memory_id = 0);

  virtual ~GridSynchronizer();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  static GridSynchronizer *
  createGridSynchronizer(Mesh & mesh,
			 const RegularGrid * grid,
			 UInt root = 0,
			 SynchronizerID id = "grid_synchronizer",
			 MemoryID memory_id = 0);

  /// function to print the contain of the class
  // virtual void printself(std::ostream & stream, int indent = 0) const;

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

#include "grid_synchronizer_inline_impl.cc"

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const GridSynchronizer & _this)
// {
//   _this.printself(stream);
//   return stream;
// }


__END_AKANTU__

#endif /* __AKANTU_GRID_SYNCHRONIZER_HH__ */
