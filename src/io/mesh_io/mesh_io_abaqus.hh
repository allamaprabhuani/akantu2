/**
 * @file   mesh_io_abaqus.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Sep 26 2010
 * @date last modification: Tue Jun 30 2015
 *
 * @brief  read a mesh from an abaqus input file
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "mesh_io.hh"
#include "mesh_accessor.hh"

#ifndef __AKANTU_MESH_IO_ABAQUS_HH__
#define __AKANTU_MESH_IO_ABAQUS_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
class MeshIOAbaqus : public MeshIO, private MeshAccessor {
public:
  MeshIOAbaqus();
  virtual ~MeshIOAbaqus();

  /// read a mesh from the file
  virtual void read(const std::string & filename, Mesh & mesh);

  /// write a mesh to a file
  //  virtual void write(const std::string & filename, const Mesh & mesh);

private:
  /// correspondence between msh element types and akantu element types
  std::map<std::string, ElementType> _abaqus_to_akantu_element_types;
};

__END_AKANTU__

#endif /* __AKANTU_MESH_IO_ABAQUS_HH__ */
