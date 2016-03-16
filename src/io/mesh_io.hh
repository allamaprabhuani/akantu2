/**
 * @file   mesh_io.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Mon Jun 01 2015
 *
 * @brief  interface of a mesh io class, reader and writer
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
#ifndef __AKANTU_MESH_IO_HH__
#define __AKANTU_MESH_IO_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class MeshIO {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MeshIO();

  virtual ~MeshIO();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void read(const std::string & filename, Mesh & mesh, const MeshIOType & type);
  void write(const std::string & filename, Mesh & mesh,
             const MeshIOType & type);

  /// read a mesh from the file
  virtual void read(__attribute__((unused)) const std::string & filename,
                    __attribute__((unused)) Mesh & mesh) {}

  /// write a mesh to a file
  virtual void write(__attribute__((unused)) const std::string & filename,
                     __attribute__((unused)) const Mesh & mesh) {}

  /// function to request the manual construction of the physical names maps
  virtual void constructPhysicalNames(const std::string & tag_name,
                                      Mesh & mesh);

  /// method to permit to be printed to a generic stream
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /// static contruction of a meshio object
  static MeshIO * getMeshIO(const std::string & filename,
                            const MeshIOType & type);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  std::map<UInt, std::string> & getPhysicalNameMap() { return phys_name_map; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  bool canReadSurface;

  bool canReadExtendedData;

  /// correspondance between a tag and physical names (if applicable)
  std::map<UInt, std::string> phys_name_map;
};

/* -------------------------------------------------------------------------- */

inline std::ostream & operator<<(std::ostream & stream, const MeshIO & _this) {
  _this.printself(stream);
  return stream;
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__

#include "mesh_io_msh.hh"
#include "mesh_io_diana.hh"
#include "mesh_io_abaqus.hh"

#if defined(AKANTU_STRUCTURAL_MECHANICS)
#include "mesh_io_msh_struct.hh"
#endif

#endif /* __AKANTU_MESH_IO_HH__ */
