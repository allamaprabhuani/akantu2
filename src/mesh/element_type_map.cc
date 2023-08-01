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
#include "fe_engine.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

FEEngineElementTypeMapArrayInitializer::FEEngineElementTypeMapArrayInitializer(
    const FEEngine & fe_engine, Int nb_component, Int spatial_dimension,
    GhostType ghost_type, ElementKind element_kind,
    const ElementTypeMapArray<Idx> * filter)
    : FEEngineElementTypeMapArrayInitializer(
          fe_engine,
          [nb_component](auto && /*type*/, auto && /*ghost_type*/) {
            return nb_component;
          },
          spatial_dimension, ghost_type, element_kind, filter) {}

FEEngineElementTypeMapArrayInitializer::FEEngineElementTypeMapArrayInitializer(
    const FEEngine & fe_engine,
    const ElementTypeMapArrayInitializer::CompFunc & nb_component,
    Int /*spatial_dimension*/, GhostType ghost_type, ElementKind element_kind,
    const ElementTypeMapArray<Idx> * filter)
    : MeshElementTypeMapArrayInitializer(
          fe_engine.getMesh(), nb_component, fe_engine.getElementDimension(),
          ghost_type, element_kind, true, false, filter),
      fe_engine(fe_engine) {}

Int FEEngineElementTypeMapArrayInitializer::size(ElementType type) const {
  auto size = MeshElementTypeMapArrayInitializer::size(type);
  if (size != 0) {
    return size * fe_engine.getNbIntegrationPoints(type, this->ghost_type);
  }

  return 0;
}

auto FEEngineElementTypeMapArrayInitializer::elementTypes() const
    -> std::vector<ElementType> {
  std::vector<ElementType> types;
  if (this->filter != nullptr) {
    for (auto type : this->filter->elementTypes(spatial_dimension, ghost_type,
                                                element_kind)) {
      types.emplace_back(type);
    }
    return types;
  }
  for (auto type : this->fe_engine.elementTypes(spatial_dimension, ghost_type,
                                                element_kind)) {
    types.emplace_back(type);
  }
  return types;
}

} // namespace akantu
