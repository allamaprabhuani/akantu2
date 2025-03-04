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
#include "aka_common.hh"
#include "aka_random_generator.hh"
#include "internal_field.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_RANDOM_INTERNAL_FIELD_HH_
#define AKANTU_RANDOM_INTERNAL_FIELD_HH_

namespace akantu {

/**
 * class for the internal fields of materials with a random
 * distribution
 */
template <typename T, template <typename> class BaseField = InternalField,
          template <typename> class Generator = RandomGenerator>
class RandomInternalField : public BaseField<T> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using BaseField<T>::BaseField;

  friend class ConstitutiveLawInternalHandler;
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  std::shared_ptr<RandomInternalField> getPtr() {
    return aka::as_type<RandomInternalField>(this->shared_from_this());
  }

  /// initialize the field to a given number of component
  void initialize(Int nb_component) override;

  /// set the field to a given value
  void setDefaultValue(const T & value) override;

  /// set the specified random distribution to a given parameter
  void setRandomDistribution(const RandomParameter<T> & param);

  /// print the content
  void printself(std::ostream & stream, int indent = 0) const override;

protected:
  void setArrayValues(T * begin, T * end) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  inline operator Real() const;

  AKANTU_GET_MACRO(RandomParameter, random_parameter,
                   const RandomParameter<T> &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// random parameter containing the distribution and base value
  RandomParameter<T> random_parameter{T()};
};

/// standard output stream operator
template <typename T, template <typename> class BaseField = InternalField,
          template <typename> class Generator = RandomGenerator>
inline std::ostream &
operator<<(std::ostream & stream,
           const RandomInternalField<T, BaseField, Generator> & _this) {
  _this.printself(stream);
  return stream;
}

template <typename T>
using DefaultRandomInternalField =
    RandomInternalField<T, InternalField, RandomGenerator>;

} // namespace akantu

#endif /* AKANTU_RANDOM_INTERNAL_FIELD_HH_ */
