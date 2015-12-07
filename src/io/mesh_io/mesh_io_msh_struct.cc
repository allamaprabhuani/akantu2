/**
 * @file   mesh_io_msh_struct.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jul 15 2011
 * @date last modification: Fri Jul 04 2014
 *
 * @brief  Read/Write for MSH files generated by gmsh
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "mesh_io.hh"


/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

MeshIOMSHStruct::MeshIOMSHStruct() : MeshIOMSH() {
  canReadSurface      = true;
  canReadExtendedData = true;

  _msh_to_akantu_element_types.clear();
  _msh_to_akantu_element_types[_msh_not_defined   ] = _not_defined;
  _msh_to_akantu_element_types[_msh_segment_2     ] = _bernoulli_beam_2;
  _msh_to_akantu_element_types[_msh_triangle_3    ] = _kirchhoff_shell;

  _akantu_to_msh_element_types.clear();
  _akantu_to_msh_element_types[_not_defined     ] = _msh_not_defined;
  _akantu_to_msh_element_types[_bernoulli_beam_2] = _msh_segment_2;
  _akantu_to_msh_element_types[_bernoulli_beam_3] = _msh_segment_2;
  _akantu_to_msh_element_types[_kirchhoff_shell] = _msh_triangle_3;

  std::map<ElementType, MSHElementType>::iterator it;
  for(it = _akantu_to_msh_element_types.begin();
      it != _akantu_to_msh_element_types.end(); ++it) {
    UInt nb_nodes = _msh_nodes_per_elem[it->second];
    std::vector<UInt> tmp = new UInt[nb_nodes];
    for (UInt i = 0; i < nb_nodes; ++i) tmp[i] = i;
    _read_order[it->first] = tmp;
  }
}


/* -------------------------------------------------------------------------- */
void MeshIOMSHStruct::read(const std::string & filename, Mesh & mesh) {
  if(mesh.getSpatialDimension() == 2) {
    _msh_to_akantu_element_types[_msh_segment_2     ] = _bernoulli_beam_2;
  } else if (mesh.getSpatialDimension() == 3) {
    _msh_to_akantu_element_types[_msh_segment_2     ] = _bernoulli_beam_3;
    AKANTU_DEBUG_WARNING("The MeshIOMSHStruct is reading bernoulli beam 3D be sure to provide the missing normals");
  }

  MeshIOMSH::read(filename, mesh);
}


__END_AKANTU__

