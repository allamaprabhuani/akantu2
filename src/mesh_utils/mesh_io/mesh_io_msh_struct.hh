/**
 * @file   mesh_io_msh_struct.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Jul 15 19:41:58 2011
 *
 * @brief  Read/Write for MSH files
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

#ifndef __AKANTU_MESH_IO_MSH_STRUCT_HH__
#define __AKANTU_MESH_IO_MSH_STRUCT_HH__

/* -------------------------------------------------------------------------- */
#include "mesh_io_msh.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__


class MeshIOMSHStruct : public MeshIOMSH {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MeshIOMSHStruct();

};


__END_AKANTU__

#endif /* __AKANTU_MESH_IO_MSH_STRUCT_HH__ */
