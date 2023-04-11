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
#ifndef AKANTU_DUMPER_IGFEM_GENERIC_ELEMENTAL_FIELD_HH_
#define AKANTU_DUMPER_IGFEM_GENERIC_ELEMENTAL_FIELD_HH_
/* -------------------------------------------------------------------------- */
#include "dumper_generic_elemental_field.hh"
#include "dumper_igfem_element_iterator.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {
namespace dumpers {
  /* --------------------------------------------------------------------------
   */
  template <class _types, template <class> class iterator_type>
  class IGFEMGenericElementalField
      : public GenericElementalField<_types, iterator_type> {

    /* ------------------------------------------------------------------------
     */
    /* Typedefs */
    /* ------------------------------------------------------------------------
     */

  public:
    typedef _types types;
    typedef typename types::data_type data_type;
    typedef typename types::it_type it_type;
    typedef typename types::field_type field_type;
    typedef typename types::array_type array_type;
    typedef typename types::array_iterator array_iterator;
    typedef typename field_type::type_iterator field_type_iterator;
    typedef iterator_type<types> iterator;

    /* ------------------------------------------------------------------------
     */
    /* Constructors/Destructors */
    /* ------------------------------------------------------------------------
     */

  public:
    IGFEMGenericElementalField(const field_type & field,
                               Int spatial_dimension = _all_dimensions,
                               GhostType ghost_type = _not_ghost,
                               ElementKind kind = _ek_igfem)
        :

          GenericElementalField<types, iterator_type>(field, spatial_dimension,
                                                      ghost_type, kind) {
      this->checkHomogeneity();
    }

    /* ------------------------------------------------------------------------
     */
    /* Methods */
    /* ------------------------------------------------------------------------
     */

  public:
    /// return the size of the contained data: i.e. the number of elements ?
    virtual UInt size() {
      this->checkHomogeneity();
      return ((this->nb_total_element) * 2);
    }

    virtual iterator begin() {
      field_type_iterator tit;
      field_type_iterator end;
      UInt sub_element = 0;

      /// type iterators on the elemental field
      tit = this->field.firstType(this->spatial_dimension, this->ghost_type,
                                  this->element_kind);
      end = this->field.lastType(this->spatial_dimension, this->ghost_type,
                                 this->element_kind);

      /// skip all types without data
      ElementType type = *tit;
      for (; tit != end && this->field(*tit, this->ghost_type).getSize() == 0;
           ++tit) {
      }
      type = *tit;

      if (tit == end)
        return this->end();

      /// getting information for the field of the given type
      const array_type & vect = this->field(type, this->ghost_type);
      UInt nb_data_per_elem = this->getNbDataPerElem(type);
      UInt nb_component = vect.getNbComponent();
      UInt size = (vect.getSize() * nb_component) / nb_data_per_elem;

      /// define element-wise iterator
      array_iterator it = vect.begin_reinterpret(nb_data_per_elem, size);
      array_iterator it_end = vect.end_reinterpret(nb_data_per_elem, size);
      /// define data iterator
      iterator rit = iterator(this->field, tit, end, it, it_end,
                              this->ghost_type, sub_element);
      rit.setNbDataPerElem(this->nb_data_per_elem);
      return rit;
    }

    virtual iterator end() {
      field_type_iterator tit;
      field_type_iterator end;
      UInt sub_element = 0;

      tit = this->field.firstType(this->spatial_dimension, this->ghost_type,
                                  this->element_kind);
      end = this->field.lastType(this->spatial_dimension, this->ghost_type,
                                 this->element_kind);

      ElementType type = *tit;
      for (; tit != end; ++tit)
        type = *tit;

      const array_type & vect = this->field(type, this->ghost_type);
      UInt nb_data = this->getNbDataPerElem(type);
      UInt nb_component = vect.getNbComponent();
      UInt size = (vect.getSize() * nb_component) / nb_data;
      array_iterator it = vect.end_reinterpret(nb_data, size);

      iterator rit = iterator(this->field, end, end, it, it, this->ghost_type,
                              sub_element);
      rit.setNbDataPerElem(this->nb_data_per_elem);
      return rit;
    }

    /* ------------------------------------------------------------------------
     */
    /* Class Members */
    /* ------------------------------------------------------------------------
     */

  protected:
  };

} // namespace dumpers
} // namespace akantu

#endif /* AKANTU_DUMPER_IGFEM_GENERIC_ELEMENTAL_FIELD_HH_ */
