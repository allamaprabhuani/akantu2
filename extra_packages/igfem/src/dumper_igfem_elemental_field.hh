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

#ifndef AKANTU_DUMPER_IGFEM_ELEMENTAL_FIELD_HH_
#define AKANTU_DUMPER_IGFEM_ELEMENTAL_FIELD_HH_
/* -------------------------------------------------------------------------- */
#include "dumper_field.hh"
#include "dumper_igfem_generic_elemental_field.hh"
#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {
namespace dumpers {
  /* --------------------------------------------------------------------------
   */

  template <typename T, template <class> class ret = Vector,
            bool filtered = false>
  class IGFEMElementalField
      : public IGFEMGenericElementalField<SingleType<T, ret, filtered>,
                                          igfem_elemental_field_iterator> {

  public:
    /* ------------------------------------------------------------------------
     */
    /* Typedefs */
    /* ------------------------------------------------------------------------
     */

    typedef SingleType<T, ret, filtered> types;
    typedef typename types::field_type field_type;
    typedef elemental_field_iterator<types> iterator;

    /* ------------------------------------------------------------------------
     */
    /* Constructors/Destructors */
    /* ------------------------------------------------------------------------
     */

    IGFEMElementalField(const field_type & field,
                        Int spatial_dimension = _all_dimensions,
                        GhostType ghost_type = _not_ghost,
                        ElementKind element_kind = _ek_igfem)
        : IGFEMGenericElementalField<types, igfem_elemental_field_iterator>(
              field, spatial_dimension, ghost_type, element_kind) {}
  };

} // namespace dumpers
} // namespace akantu

#endif /* AKANTU_DUMPER_IGFEM_ELEMENTAL_FIELD_HH_ */
