/**
 * @file   random_internal_field_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Nov  1 10:24:11 2013
 *
 * @brief  Random internal material parameter implementation
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
#include "aka_common.hh"
#include "aka_random_generator.hh"
#include "internal_field_tmpl.hh"

/* -------------------------------------------------------------------------- */


#ifndef __AKANTU_RANDOM_INTERNAL_FIELD_TMPL_HH__
#define __AKANTU_RANDOM_INTERNAL_FIELD_TMPL_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<typename T,
	 template<typename> class BaseField,
	 template<typename> class Generator>
RandomInternalField<T, BaseField, Generator>::RandomInternalField(const ID & id, Material & material) :
  BaseField<T>(id, material), random_parameter(T()) {
}

/* -------------------------------------------------------------------------- */
template<typename T,
	 template<typename> class BaseField,
	 template<typename> class Generator>
RandomInternalField<T, BaseField, Generator>::~RandomInternalField() {  }

/* -------------------------------------------------------------------------- */
template<typename T,
	 template<typename> class BaseField,
	 template<typename> class Generator>
void RandomInternalField<T, BaseField, Generator>::initialize(UInt nb_component) {
  this->internalInitialize(nb_component);
}

/* ------------------------------------------------------------------------ */
template<typename T,
	 template<typename> class BaseField,
	 template<typename> class Generator>
void RandomInternalField<T, BaseField, Generator>::setDefaultValue(const T &value) {
  random_parameter.setBaseValue(value);
  this->reset();
}

/* ------------------------------------------------------------------------ */
template<typename T,
	 template<typename> class BaseField,
	 template<typename> class Generator>
void RandomInternalField<T, BaseField, Generator>::setRandomDistribution(const RandomParameter<T> & param) {
  random_parameter = param;
  this->reset();
}

/* ------------------------------------------------------------------------ */
template<typename T,
	 template<typename> class BaseField,
	 template<typename> class Generator>
void RandomInternalField<T, BaseField, Generator>::printself(std::ostream & stream, unsigned int indent) const {
  stream << "RandomInternalField [ ";
  random_parameter.printself(stream);
  stream << " ]";
#if !defined(AKANTU_NDEBUG)
  if(AKANTU_DEBUG_TEST(dblDump)) {
    stream << std::endl;
    InternalField<T>::printself(stream, indent);
  }
#endif
}

/* -------------------------------------------------------------------------- */
template<typename T,
	 template<typename> class BaseField,
	 template<typename> class Generator>
void RandomInternalField<T, BaseField, Generator>::setArrayValues(T * begin, T * end) {
  random_parameter.template setValues<Generator>(begin, end);
}

/* -------------------------------------------------------------------------- */
template<typename T,
	 template<typename> class BaseField,
	 template<typename> class Generator>
inline RandomInternalField<T, BaseField, Generator>::operator Real() const {
  return random_parameter.getBaseValue();
}

/* -------------------------------------------------------------------------- */
template<>
inline void ParsableParamTyped< RandomInternalField<Real> >::parseParam(const ParserParameter & in_param) {
  ParsableParam::parseParam(in_param);
  RandomParameter<Real> r = in_param;
  param.setRandomDistribution(r);
}

/* -------------------------------------------------------------------------- */


__END_AKANTU__

#endif /* __AKANTU_RANDOM_INTERNAL_FIELD_TMPL_HH__ */
