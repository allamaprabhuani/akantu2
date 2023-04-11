/**
 * Copyright (©) 2017-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "shape_lagrange_base.hh"
#include "mesh_iterators.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

ShapeLagrangeBase::ShapeLagrangeBase(const Mesh & mesh, Int spatial_dimension,
                                     ElementKind kind, const ID & id)
    : ShapeFunctions(mesh, spatial_dimension, id), _kind(kind) {}

/* -------------------------------------------------------------------------- */
ShapeLagrangeBase::~ShapeLagrangeBase() = default;

/* -------------------------------------------------------------------------- */
void ShapeLagrangeBase::computeShapesOnIntegrationPoints(
    const Array<Real> & nodes, const Ref<const MatrixXr> integration_points,
    Array<Real> & shapes, ElementType type, GhostType ghost_type,
    const Array<Idx> & filter_elements) const {
  auto && call = [&](auto && enum_type) {
    constexpr ElementType type = std ::decay_t<decltype(enum_type)>::value;
    this->computeShapesOnIntegrationPoints<type>(
        nodes, integration_points, shapes, ghost_type, filter_elements);
  };
  tuple_dispatch<
      tuple::cat_t<ElementTypes_t<_ek_regular>, ElementTypes_t<_ek_cohesive>>>(
      call, type);
}

/* -------------------------------------------------------------------------- */
void ShapeLagrangeBase::onElementsAdded(const Array<Element> & new_elements) {
  AKANTU_DEBUG_IN();
  const auto & nodes = mesh.getNodes();

  for (auto elements_range : MeshElementsByTypes(new_elements)) {
    auto type = elements_range.getType();
    auto ghost_type = elements_range.getGhostType();

    if (mesh.getSpatialDimension(type) != _spatial_dimension) {
      continue;
    }

    if (Mesh::getKind(type) != _kind) {
      continue;
    }

    const auto & elements = elements_range.getElements();

    auto itp_type = FEEngine::getInterpolationType(type);

    if (not shapes.exists(itp_type, ghost_type)) {
      auto size_of_shapes = this->getShapeSize(type);
      this->shapes.alloc(0, size_of_shapes, itp_type, ghost_type);
    }

    const auto & natural_coords = integration_points(type, ghost_type);
    computeShapesOnIntegrationPoints(nodes, natural_coords,
                                     shapes(itp_type, ghost_type), type,
                                     ghost_type, elements);

    if (_spatial_dimension != mesh.getSpatialDimension()) {
      continue;
    }

    if (not this->shapes_derivatives.exists(itp_type, ghost_type)) {
      auto size_of_shapesd = this->getShapeDerivativesSize(type);
      this->shapes_derivatives.alloc(0, size_of_shapesd, itp_type, ghost_type);
    }

    computeShapeDerivativesOnIntegrationPoints(
        nodes, natural_coords, shapes_derivatives(itp_type, ghost_type), type,
        ghost_type, elements);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ShapeLagrangeBase::onElementsRemoved(
    const Array<Element> &, const ElementTypeMapArray<Idx> & new_numbering) {
  this->shapes.onElementsRemoved(new_numbering);
  this->shapes_derivatives.onElementsRemoved(new_numbering);
}

/* -------------------------------------------------------------------------- */
void ShapeLagrangeBase::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "Shapes Lagrange [" << std::endl;
  ShapeFunctions::printself(stream, indent + 1);
  shapes.printself(stream, indent + 1);
  shapes_derivatives.printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}

} // namespace akantu
