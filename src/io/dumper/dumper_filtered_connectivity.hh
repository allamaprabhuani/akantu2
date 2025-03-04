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

#include "dumper_generic_elemental_field.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {
namespace dumpers {
  /* ------------------------------------------------------------------------ */

  template <class types>
  class filtered_connectivity_field_iterator
      : public element_iterator<types, filtered_connectivity_field_iterator> {
    /* ---------------------------------------------------------------------- */
    /* Typedefs                                                               */
    /* ---------------------------------------------------------------------- */
  public:
    using parent =
        element_iterator<types, dumpers::filtered_connectivity_field_iterator>;
    using return_type = typename types::return_type;
    using field_type = typename types::field_type;
    using array_iterator = typename types::array_iterator;

    /* ---------------------------------------------------------------------- */
    /* Constructors/Destructors                                               */
    /* ---------------------------------------------------------------------- */
  public:
    filtered_connectivity_field_iterator(
        const field_type & field,
        const typename field_type::type_iterator & t_it,
        const typename field_type::type_iterator & t_it_end,
        const array_iterator & array_it, const array_iterator & array_it_end,
        const GhostType ghost_type = _not_ghost)
        : parent(field, t_it, t_it_end, array_it, array_it_end, ghost_type) {}

    /* ---------------------------------------------------------------------- */
    /* Methods                                                                */
    /* ---------------------------------------------------------------------- */
  public:
    return_type operator*() {
      const auto & old_connect = *this->array_it;
      Vector<Int> new_connect(old_connect.size());
      auto nodes_begin = nodal_filter->begin();
      auto nodes_end = nodal_filter->end();
      for (Int i(0); i < old_connect.size(); ++i) {
        auto new_id = std::find(nodes_begin, nodes_end, old_connect(i));
        if (new_id == nodes_end) {
          AKANTU_EXCEPTION("Node not found in the filter!");
        }
        new_connect(i) = new_id - nodes_begin;
      }
      return new_connect;
    }

    void setNodalFilter(const Array<Idx> & new_nodal_filter) {
      nodal_filter = &new_nodal_filter;
    }

    /* ---------------------------------------------------------------------- */
    /* Class Members                                                          */
    /* ---------------------------------------------------------------------- */
  private:
    const Array<Idx> * nodal_filter;
  };

  /* ------------------------------------------------------------------------ */

  class FilteredConnectivityField
      : public GenericElementalField<SingleType<Idx, Vector<Idx>, true>,
                                     filtered_connectivity_field_iterator> {
    /* ---------------------------------------------------------------------- */
    /* Typedefs                                                               */
    /* ---------------------------------------------------------------------- */
  public:
    using types = SingleType<Idx, Vector<Idx>, true>;
    using iterator = filtered_connectivity_field_iterator<types>;
    using field_type = types::field_type;
    using parent =
        GenericElementalField<types, filtered_connectivity_field_iterator>;

    /* ---------------------------------------------------------------------- */
    /* Constructors/Destructors                                               */
    /* ---------------------------------------------------------------------- */
  public:
    FilteredConnectivityField(const field_type & field,
                              const Array<Idx> & nodal_filter,
                              Int spatial_dimension = _all_dimensions,
                              GhostType ghost_type = _not_ghost,
                              ElementKind element_kind = _ek_not_defined)
        : parent(field, spatial_dimension, ghost_type, element_kind),
          nodal_filter(nodal_filter) {}

    ~FilteredConnectivityField() override {
      // since the field is created in registerFilteredMesh it is destroyed here
      delete const_cast<field_type *>(&this->field);
    }

    /* ---------------------------------------------------------------------- */
    /* Methods                                                                */
    /* ---------------------------------------------------------------------- */
  public:
    iterator begin() override {
      auto it = parent::begin();
      it.setNodalFilter(nodal_filter);
      return it;
    }

    iterator end() override {
      auto it = parent::end();
      it.setNodalFilter(nodal_filter);
      return it;
    }

    /* ---------------------------------------------------------------------- */
    /* Class Members                                                          */
    /* ---------------------------------------------------------------------- */
  private:
    const Array<Idx> & nodal_filter;
  };

  /* ------------------------------------------------------------------------ */

} // namespace dumpers
} // namespace akantu

/* -------------------------------------------------------------------------- */
