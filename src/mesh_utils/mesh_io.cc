/**
 * @file   mesh_io.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Thu Jul 15 00:41:12 2010
 *
 * @brief  common part for all mesh io classes
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
#include "mesh_io.hh"

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MeshIO::MeshIO() {
  canReadSurface      = false;
  canReadExtendedData = false;
}

/* -------------------------------------------------------------------------- */
MeshIO::~MeshIO() {

}

/* -------------------------------------------------------------------------- */
MeshIO * MeshIO::getMeshIO(const std::string & filename, const MeshIOType & type) {
  MeshIOType t = type;
  if(type == _miot_auto) {
    std::string::size_type idx = filename.rfind('.');
    std::string ext;
    if(idx != std::string::npos) {
      ext = filename.substr(idx+1);
    }

    if(ext == "msh") t = _miot_gmsh;
    else if(ext == "diana") t = _miot_diana;
    else AKANTU_EXCEPTION("Cannot guess the type of file of "
			  << filename << " (ext "<< ext <<"). "
			  << "Please provide the MeshIOType to the read function");
  }

  switch(t) {
  case _miot_gmsh: return new MeshIOMSH();
  case _miot_diana: return new MeshIODiana();
  default:
    return NULL;
  }
}

/* -------------------------------------------------------------------------- */
void MeshIO::read(const std::string & filename, Mesh & mesh, const MeshIOType & type) {
  MeshIO * mesh_io = getMeshIO(filename, type);
  mesh_io->read(filename, mesh);
  delete mesh_io;
}

/* -------------------------------------------------------------------------- */
void MeshIO::write(const std::string & filename, Mesh & mesh, const MeshIOType & type) {
  MeshIO * mesh_io = getMeshIO(filename, type);
  mesh_io->write(filename, mesh);
  delete mesh_io;
}

__END_AKANTU__
