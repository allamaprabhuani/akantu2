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
#include "dumper_field.hh"
#include "dumper_variable.hh"
#include "support.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_DUMPER_FILE_BASE_HH__
#define __AKANTU_DUMPER_FILE_BASE_HH__

namespace akantu {
namespace dumper {

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
    /* ---------------------------------------------------------------------- */
    virtual void dump(Support<Mesh> & support) = 0;
    virtual void dump(Support<ElementGroup> & /*support*/) {
      AKANTU_TO_IMPLEMENT();
    }

    /* ---------------------------------------------------------------------- */
    virtual void dump(FieldBase & field) {
      using dumper::FieldType;
      switch (field.getFieldType()) {
      case FieldType::_node_array: /* FALLTHRU */
      case FieldType::_node_array_function:
        dump(aka::as_type<FieldArrayBase>(field));
        break;
      case FieldType::_element_map_array: /* FALLTHRU */
      case FieldType::_element_map_array_function:
        dump(aka::as_type<FieldElementMapArrayBase>(field));
        break;
      case FieldType::_internal_field: /* FALLTHRU */
      case FieldType::_internal_field_function:
        dump(aka::as_type<FieldElementMapArrayBase>(field));
        break;
      case FieldType::_not_defined: /* FALLTHRU */
      default:
        AKANTU_EXCEPTION("The field type is not properly defined");
        break;
      }
    }

  public:
    virtual void dump() { dump(support); }

    /* ---------------------------------------------------------------------- */
    virtual void dump(SupportBase & support) {
      using dumper::SupportType;
      switch (support.getType()) {
      case SupportType::_mesh:
        dump(aka::as_type<Support<Mesh>>(support));
        break;
      case SupportType::_element_group:
        dump(aka::as_type<Support<ElementGroup>>(support));
        break;
      }
    }

  protected:
    SupportBase & support;
  };

} // namespace dumper
} // namespace akantu

#endif /* __AKANTU_DUMPER_FILE_BASE_HH__ */
