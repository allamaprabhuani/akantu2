/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "mesh_io_msh_struct.hh"
/* -------------------------------------------------------------------------- */
#include <numeric>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
MeshIOMSHStruct::MeshIOMSHStruct() {
  canReadSurface = true;
  canReadExtendedData = true;

  _msh_to_akantu_element_types.clear();
  _msh_to_akantu_element_types[_msh_not_defined] = _not_defined;
  _msh_to_akantu_element_types[_msh_segment_2] = _bernoulli_beam_2;
  _msh_to_akantu_element_types[_msh_triangle_3] =
      _discrete_kirchhoff_triangle_18;

  _akantu_to_msh_element_types.clear();
  _akantu_to_msh_element_types[_not_defined] = _msh_not_defined;
  _akantu_to_msh_element_types[_bernoulli_beam_2] = _msh_segment_2;
  _akantu_to_msh_element_types[_bernoulli_beam_3] = _msh_segment_2;
  _akantu_to_msh_element_types[_discrete_kirchhoff_triangle_18] =
      _msh_triangle_3;

  for (auto & kv_pair : _akantu_to_msh_element_types) {
    Int nb_nodes = _msh_nodes_per_elem[kv_pair.second];
    std::vector<Idx> tmp(nb_nodes);
    std::iota(tmp.begin(), tmp.end(), 0);
    _read_order[kv_pair.first] = tmp;
  }
}

/* -------------------------------------------------------------------------- */
void MeshIOMSHStruct::read(const std::string & filename, Mesh & mesh) {
  if (mesh.getSpatialDimension() == 2) {
    _msh_to_akantu_element_types[_msh_segment_2] = _bernoulli_beam_2;
  } else if (mesh.getSpatialDimension() == 3) {
    _msh_to_akantu_element_types[_msh_segment_2] = _bernoulli_beam_3;
    AKANTU_DEBUG_WARNING("The MeshIOMSHStruct is reading bernoulli beam 3D be "
                         "sure to provide the missing normals with the element "
                         "data \"extra_normal\"");
  }

  MeshIOMSH::read(filename, mesh);
}

} // namespace akantu
