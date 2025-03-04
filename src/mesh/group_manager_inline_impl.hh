/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "dumper_field.hh"
#include "element_group.hh"
#include "element_type_map_filter.hh"
/* -------------------------------------------------------------------------- */
#include "dumper_nodal_field.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T, template <bool> class dump_type>
std::shared_ptr<dumpers::Field>
GroupManager::createElementalField(const ElementTypeMapArray<T> & field,
                                   const std::string & group_name,
                                   Int spatial_dimension, ElementKind kind,
                                   ElementTypeMap<Int> nb_data_per_elem) {

  const auto * field_ptr = &field;
  if (field_ptr == nullptr) {
    return nullptr;
  }
  if (group_name == "all") {
    return this->createElementalField<dump_type<false>>(
        field, group_name, spatial_dimension, kind, nb_data_per_elem);
  }
  return this->createElementalFilteredField<dump_type<true>>(
      field, group_name, spatial_dimension, kind, nb_data_per_elem);
}

/* -------------------------------------------------------------------------- */
template <typename T, class ret_type,
          template <class, class, bool> class dump_type>
std::shared_ptr<dumpers::Field>
GroupManager::createElementalField(const ElementTypeMapArray<T> & field,
                                   const std::string & group_name,
                                   Int spatial_dimension, ElementKind kind,
                                   ElementTypeMap<Int> nb_data_per_elem) {

  const auto * field_ptr = &field;
  if (field_ptr == nullptr) {
    return nullptr;
  }
  if (group_name == "all") {
    return this->createElementalField<dump_type<T, ret_type, false>>(
        field, group_name, spatial_dimension, kind, nb_data_per_elem);
  }
  return this->createElementalFilteredField<dump_type<T, ret_type, true>>(
      field, group_name, spatial_dimension, kind, nb_data_per_elem);
}

/* -------------------------------------------------------------------------- */
template <typename T, template <typename T2, bool filtered>
                      class dump_type> ///< type of InternalMaterialField
std::shared_ptr<dumpers::Field>
GroupManager::createElementalField(const ElementTypeMapArray<T> & field,
                                   const std::string & group_name,
                                   Int spatial_dimension, ElementKind kind,
                                   ElementTypeMap<Int> nb_data_per_elem) {
  const auto * field_ptr = &field;

  if (field_ptr == nullptr) {
    return nullptr;
  }
  if (group_name == "all") {
    return this->createElementalField<dump_type<T, false>>(
        field, group_name, spatial_dimension, kind, nb_data_per_elem);
  }
  return this->createElementalFilteredField<dump_type<T, true>>(
      field, group_name, spatial_dimension, kind, nb_data_per_elem);
}

/* -------------------------------------------------------------------------- */
template <typename dump_type, typename field_type>
std::shared_ptr<dumpers::Field> GroupManager::createElementalField(
    const field_type & field, const std::string & group_name,
    Int spatial_dimension, ElementKind kind,
    const ElementTypeMap<Int> & nb_data_per_elem) {
  const field_type * field_ptr = &field;
  if (field_ptr == nullptr) {
    return nullptr;
  }

  if (group_name != "all") {
    throw;
  }

  auto dumper =
      std::make_shared<dump_type>(field, spatial_dimension, _not_ghost, kind);
  dumper->setNbDataPerElem(nb_data_per_elem);
  return dumper;
}

/* -------------------------------------------------------------------------- */
template <typename dump_type, typename field_type>
std::shared_ptr<dumpers::Field> GroupManager::createElementalFilteredField(
    const field_type & field, const std::string & group_name,
    Int spatial_dimension, ElementKind kind,
    ElementTypeMap<Int> nb_data_per_elem) {

  const auto * field_ptr = &field;
  if (field_ptr == nullptr) {
    return nullptr;
  }
  if (group_name == "all") {
    throw;
  }

  using T = typename field_type::value_type;
  auto & group = this->getElementGroup(group_name);
  auto dim = group.getDimension();
  if (dim != spatial_dimension) {
    throw;
  }
  const auto & elemental_filter = group.getElements();

  auto * filtered = new ElementTypeMapArrayFilter<T>(field, elemental_filter,
                                                     nb_data_per_elem);

  auto dumper = std::make_shared<dump_type>(*filtered, dim, _not_ghost, kind);
  dumper->setNbDataPerElem(nb_data_per_elem);

  return dumper;
}

/* -------------------------------------------------------------------------- */
template <typename type, bool flag, template <class, bool> class ftype>
std::shared_ptr<dumpers::Field>
GroupManager::createNodalField(const ftype<type, flag> * field,
                               const std::string & group_name,
                               Int padding_size) {
  return createStridedNodalField(field, group_name, 0, 0, padding_size);
}

/* -------------------------------------------------------------------------- */
template <typename type, bool flag, template <class, bool> class ftype>
std::shared_ptr<dumpers::Field>
GroupManager::createStridedNodalField(const ftype<type, flag> * field,
                                      const std::string & group_name, Int size,
                                      Int stride, Int padding_size) {
  if (not field) {
    return nullptr;
  }

  if (group_name == "all") {
    using DumpType = typename dumpers::NodalField<type, false>;
    auto dumper = std::make_shared<DumpType>(*field, size, stride);
    dumper->setPadding(padding_size);
    return dumper;
  }

  ElementGroup & group = this->getElementGroup(group_name);
  const auto & nodal_filter = group.getNodeGroup().getNodes();
  using DumpType = typename dumpers::NodalField<type, true>;
  auto dumper = std::make_shared<DumpType>(*field, size, stride, &nodal_filter);
  dumper->setPadding(padding_size);
  return dumper;
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
