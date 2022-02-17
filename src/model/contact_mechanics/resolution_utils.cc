/**
 * @file   resolution_utils.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon May 20 2019
 * @date last modification: Sun Jun 06 2021
 *
 * @brief  Implementation of various utilities neede for resolution class
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "resolution_utils.hh"
#include "element_class_helper.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void ResolutionUtils::computeShapeFunctionMatric(
    const ContactElement & element, const Vector<Real> & projection,
    Matrix<Real> & shape_matric) {

  shape_matric.zero();

  const ElementType & type = element.master.type;

  auto surface_dimension = Mesh::getSpatialDimension(type);
  auto spatial_dimension = surface_dimension + 1;
  UInt nb_nodes_per_contact = element.getNbNodes();

  AKANTU_DEBUG_ASSERT(spatial_dimension == shape_matric.rows() &&
                          spatial_dimension * nb_nodes_per_contact ==
                              shape_matric.cols(),
                      "Shape Matric dimensions are not correct");

  auto && shapes = ElementClassHelper<_ek_regular>::getN(projection, type);

  for (auto i : arange(nb_nodes_per_contact)) {
    for (auto j : arange(spatial_dimension)) {
      if (i == 0) {
        shape_matric(j, i * spatial_dimension + j) = 1;
        continue;
      }
      shape_matric(j, i * spatial_dimension + j) = -shapes[i - 1];
    }
  }
}

} // namespace akantu
