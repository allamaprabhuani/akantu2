/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
//#include "mesh.hh"
#include "shape_igfem.hh"
/* -------------------------------------------------------------------------- */

#if defined(AKANTU_IGFEM)

namespace akantu {

/* -------------------------------------------------------------------------- */
ShapeLagrange<_ek_igfem>::ShapeLagrange(const Mesh & mesh, const ID & id)
    : ShapeFunctions(mesh, id), shapes("shapes_generic", id),
      shapes_derivatives("shapes_derivatives_generic", id),
      igfem_integration_points("igfem_integration_points", id),
      shapes_at_enrichments("shapes_at_enrichments", id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/*-------------------------------------------------------------------------- */
void ShapeLagrange<_ek_igfem>::extractValuesAtStandardNodes(
    const Array<Real> & nodal_values, Array<Real> & extracted_values,
    GhostType ghost_type) const {

  AKANTU_DEBUG_ASSERT(nodal_values.getNbComponent() ==
                          extracted_values.getNbComponent(),
                      "The arrays are not of the same size!!!!!");
  extracted_values.zero();
  Int spatial_dimension = mesh.getSpatialDimension();
  Mesh::type_iterator it =
      mesh.firstType(spatial_dimension, ghost_type, _ek_igfem);
  Mesh::type_iterator end =
      mesh.lastType(spatial_dimension, ghost_type, _ek_igfem);
  for (; it != end; ++it) {
    ElementType type = *it;
    UInt nb_elements = mesh.getNbElement(type, ghost_type);
    UInt nb_parent_nodes = 0;
    UInt nb_nodes_per_element = 0;
#define GET_NODES_INFO(type)                                                   \
  const ElementType parent_type =                                              \
      ElementClassProperty<type>::parent_element_type;                         \
  nb_parent_nodes =                                                            \
      ElementClass<parent_type>::getNbNodesPerInterpolationElement();          \
  nb_nodes_per_element =                                                       \
      ElementClass<type>::getNbNodesPerInterpolationElement();

    AKANTU_BOOST_IGFEM_ELEMENT_SWITCH(GET_NODES_INFO);
#undef GET_NODES_INFO

    UInt * conn_val = mesh.getConnectivity(type, ghost_type).data();
    for (Int e = 0; e < nb_elements; ++e) {
      /// copy the value at standard nodes
      UInt offset = e * nb_nodes_per_element;
      for (Int n = 0; n < nb_parent_nodes; ++n) {
        UInt node = conn_val[offset + n];
        for (Int i = 0; i < nodal_values.getNbComponent(); ++i)
          extracted_values(node, i) = nodal_values(node, i);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void ShapeLagrange<_ek_igfem>::printself(std::ostream & stream,
                                         int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "Shapes Lagrange [" << std::endl;
  ShapeFunctions::printself(stream, indent + 1);
  shapes.printself(stream, indent + 1);
  shapes_derivatives.printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}
} // namespace akantu

#endif
