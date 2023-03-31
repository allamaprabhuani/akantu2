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

#ifndef AKANTU_DUMPER_IGFEM_MATERIAL_INTERNAL_FIELD_HH_
#define AKANTU_DUMPER_IGFEM_MATERIAL_INTERNAL_FIELD_HH_
/* -------------------------------------------------------------------------- */
#include "dumper_igfem_quadrature_points_field.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {
namespace dumpers {
  /* --------------------------------------------------------------------------
   */

  template <typename T, bool filtered = false>
  class IGFEMInternalMaterialField
      : public IGFEMGenericElementalField<SingleType<T, Vector, filtered>,
                                          igfem_quadrature_point_iterator> {

    /* ------------------------------------------------------------------------
     */
    /* Typedefs */
    /* ------------------------------------------------------------------------
     */

  public:
    typedef SingleType<T, Vector, filtered> types;
    typedef IGFEMGenericElementalField<types, igfem_quadrature_point_iterator>
        parent;
    typedef typename types::field_type field_type;

    /* ------------------------------------------------------------------------
     */
    /* Constructors/Destructors */
    /* ------------------------------------------------------------------------
     */

    IGFEMInternalMaterialField(const field_type & field,
                               Int spatial_dimension = _all_dimensions,
                               GhostType ghost_type = _not_ghost,
                               ElementKind kind = _ek_igfem)
        : parent(field, spatial_dimension, ghost_type, kind) {}
  };

} // namespace dumpers
} // namespace akantu

#endif /* AKANTU_DUMPER_IGFEM_MATERIAL_INTERNAL_FIELD_HH_ */
