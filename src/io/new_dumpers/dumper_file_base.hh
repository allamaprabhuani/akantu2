/**
 * @file   dumper_file_base.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Thu Oct 12 2017
 *
 * @brief A Documented file.
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#include "dumper_types.hh"
/* -------------------------------------------------------------------------- */
#include "dumper_field.hh"
#include "dumper_variable.hh"
#include "support.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_DUMPER_FILE_BASE_HH
#define AKANTU_DUMPER_FILE_BASE_HH

namespace akantu {
namespace dumper {

  /* ------------------------------------------------------------------------ */
  class FileBase {
  public:
    /* ---------------------------------------------------------------------- */
    FileBase(SupportBase & support) : support(support) {}
    FileBase(const FileBase & other) = default;
    FileBase & operator=(const FileBase & other) = delete;
    virtual ~FileBase() = default;

  protected:
    /* ---------------------------------------------------------------------- */
    virtual void dump(FieldArrayBase & field) = 0;
    virtual void dump(FieldElementMapArrayBase & field) = 0;

    virtual void dump(FieldNodalArrayBase & field) {
      dump(aka::as_type<FieldArrayBase>(field));
    }

    virtual void dump(FieldElementalArrayBase & field) {
      dump(aka::as_type<FieldArrayBase>(field));
    }

    /* ---------------------------------------------------------------------- */
    virtual void dump(Support<Mesh> & support) = 0;
    virtual void dump(Support<ElementGroup> & /*support*/) {
      AKANTU_TO_IMPLEMENT();
    }

    /* ---------------------------------------------------------------------- */
    virtual void dump(FieldBase & field) {
      tuple_dispatch<AllFieldTypes>(
          [&](auto field_type) { this->dump(field_cast(field, field_type)); },
          field.getFieldType());
    }

    /* ---------------------------------------------------------------------- */
    virtual void dump(SupportBase & support) {
      tuple_dispatch<AllSupportTypes>(
          [&](auto type) { this->dump(support_cast(support, type)); },
          support.getType());
    }

    /* ---------------------------------------------------------------------- */
    virtual void read(FieldNodalArrayBase & /*field*/) {
      AKANTU_TO_IMPLEMENT();
    }
    /* ---------------------------------------------------------------------- */
    virtual void read(FieldArrayBase & /*field*/) { AKANTU_TO_IMPLEMENT(); }
    /* ---------------------------------------------------------------------- */
    virtual void read(FieldElementMapArrayBase & /*field*/) {
      AKANTU_TO_IMPLEMENT();
    }
    /* ---------------------------------------------------------------------- */
    virtual void read(Support<Mesh> & /*support*/) { AKANTU_TO_IMPLEMENT(); }
    /* ---------------------------------------------------------------------- */
    virtual void read(Support<ElementGroup> & /*support*/) {
      AKANTU_TO_IMPLEMENT();
    }

    /* ---------------------------------------------------------------------- */
    virtual void read(FieldBase & field) {
      tuple_dispatch<AllFieldTypes>(
          [&](auto field_type) { this->read(field_cast(field, field_type)); },
          field.getFieldType());
    }

    /* ---------------------------------------------------------------------- */
    virtual void read(SupportBase & support) {
      tuple_dispatch<AllSupportTypes>(
          [&](auto type) { this->read(support_cast(support, type)); },
          support.getType());
    }

  public:
    virtual void dump() { dump(support); }
    virtual void read() { read(support); }

  protected:
    SupportBase & support;
  };

} // namespace dumper
} // namespace akantu

#endif /* AKANTU_DUMPER_FILE_BASE_HH */
