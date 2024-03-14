/**
 * Copyright (©) 2014-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
namespace akantu {
namespace dumpers {
#ifdef AKANTU_IGFEM
#include "dumper_igfem_element_partition.hh"
#endif
  /* ------------------------------------------------------------------------ */
  template <class types>
  class element_partition_field_iterator
      : public element_iterator<types, element_partition_field_iterator> {
  public:
    using parent =
        element_iterator<types, dumpers::element_partition_field_iterator>;
    using return_type =
        typename SingleType<int, Vector<int>, true>::return_type;
    using array_iterator = typename types::array_iterator;
    using field_type = typename types::field_type;

    element_partition_field_iterator(
        const field_type & field,
        const typename field_type::type_iterator & t_it,
        const typename field_type::type_iterator & t_it_end,
        const array_iterator & array_it, const array_iterator & array_it_end,
        const GhostType ghost_type = _not_ghost)
        : parent(field, t_it, t_it_end, array_it, array_it_end, ghost_type) {
      prank = Communicator::getWorldCommunicator().whoAmI();
    }

    return_type operator*() {
      return_type ret(1);
      ret.fill(prank);
      return ret;
    }

  protected:
    Int prank;
  };

  /* ------------------------------------------------------------------------ */
  template <bool filtered = false>
  class ElementPartitionField
      : public GenericElementalField<SingleType<Int, Vector<Int>, filtered>,
                                     element_partition_field_iterator> {
  public:
    using types = SingleType<Int, Vector<Int>, filtered>;
    using iterator = element_partition_field_iterator<types>;
    using parent =
        GenericElementalField<types, element_partition_field_iterator>;
    using field_type = typename types::field_type;

    ElementPartitionField(const field_type & field,
                          Int spatial_dimension = _all_dimensions,
                          GhostType ghost_type = _not_ghost,
                          ElementKind element_kind = _ek_not_defined)
        : parent(field, spatial_dimension, ghost_type, element_kind) {
      this->homogeneous = true;
    }

    Int getDim() override { return 1; }
  };

  /* ------------------------------------------------------------------------ */

} // namespace dumpers
} // namespace akantu
