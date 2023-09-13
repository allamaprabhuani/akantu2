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
#include "internal_field.hh"
#include "random_internal_field.hh"

#ifndef AKANTU_COHESIVE_INTERNAL_FIELD_HH_
#define AKANTU_COHESIVE_INTERNAL_FIELD_HH_

namespace akantu {

/// internal field class for cohesive materials
template <typename T> class CohesiveInternalField : public InternalField<T> {
public:
  CohesiveInternalField(const ID & id,
                        ConstitutiveLawInternalHandler & constitutive_law,
                        Int dim, const ID & fem_id,
                        const ElementTypeMapArray<Idx> & element_filter);

  /// initialize the field to a given number of component
  void initialize(Int nb_component) override;

private:
  friend class ConstitutiveLawInternalHandler;

  CohesiveInternalField operator=(const CohesiveInternalField & /*other*/){};
};

/* -------------------------------------------------------------------------- */
/* Facet Internal Field                                                       */
/* -------------------------------------------------------------------------- */
template <typename T> class FacetInternalField : public InternalField<T> {
public:
  FacetInternalField(const ID & id,
                     ConstitutiveLawInternalHandler & constitutive_law, Int dim,
                     const ID & fem_id,
                     const ElementTypeMapArray<Idx> & element_filter);

  std::shared_ptr<FacetInternalField> getPtr() {
    return aka::as_type<FacetInternalField>(this->shared_from_this());
  }

  /// initialize the field to a given number of component
  void initialize(Int nb_component) override;

  friend class ConstitutiveLawInternalHandler;
};

template <typename T>
using CohesiveRandomInternalField =
    RandomInternalField<T, CohesiveInternalField, RandomGenerator>;

template <typename T>
using FacetRandomInternalField =
    RandomInternalField<T, FacetInternalField, RandomGenerator>;

} // namespace akantu

#endif /* AKANTU_COHESIVE_INTERNAL_FIELD_HH_ */
