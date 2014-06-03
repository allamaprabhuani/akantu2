/**
 * @file   random_internal_field.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Nov  1 10:24:11 2013
 *
 * @brief  Random internal material parameter
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
#include "internal_field.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_RANDOM_INTERNAL_FIELD_HH__
#define __AKANTU_RANDOM_INTERNAL_FIELD_HH__

__BEGIN_AKANTU__

template<typename T,
	 template<typename> class BaseField = InternalField,
	 template<typename> class Generator = RandGenerator>
class RandomInternalField : public BaseField<T> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  RandomInternalField(const ID & id, Material & material);

  virtual ~RandomInternalField();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:
  virtual void initialize(UInt nb_component);

  void setDefaultValue(const T & value);

  void setRandomDistribution(const RandomParameter<T> & param);

  virtual void printself(std::ostream & stream, unsigned int indent = 0) const;

protected:
  virtual void setArrayValues(T * begin, T * end);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  inline operator Real() const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// random parameter containing the distribution and base value
  RandomParameter<T> random_parameter;
};


/// standard output stream operator
template<typename T>
inline std::ostream & operator <<(std::ostream & stream, const RandomInternalField<T> & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AKANTU_RANDOM_INTERNAL_FIELD_HH__ */
