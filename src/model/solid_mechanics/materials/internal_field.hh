/**
 * @file   internal_field.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Nov  5 14:25:18 2013
 *
 * @brief  Material internal properties
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
#include "by_element_type.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_INTERNAL_FIELD_HH__
#define __AKANTU_INTERNAL_FIELD_HH__

__BEGIN_AKANTU__

class Material;
class FEM;

template<typename T>
class InternalField : public ByElementTypeArray<T> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  InternalField(const ID & id, Material & material);
  virtual ~InternalField();

protected:
  InternalField(const ID & id, Material & material, FEM & fem,
		const ByElementTypeArray<UInt> & element_filter);
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the field to a given number of component
  virtual void initialize(UInt nb_component);

  /// resize the arrays and set the new element to 0
  virtual void resize();

  /// set the field to a given value v
  virtual void setDefaultValue(const T & v);

  /// reset all the fields to the default value
  virtual void reset();

  /// remove the quadrature points corresponding to suppressed elements
  virtual void removeQuadraturePoints(const ByElementTypeUInt & new_numbering);

  /// print the content
  void printself(std::ostream & stream, UInt indent = 0) const;

protected:
  /// initialize the arrays in the ByElementTypeArray<T>
  void internalInitialize(UInt nb_component);

  /// set the values for new internals
  virtual void setArrayValues(T * begin, T * end);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// the material for which this is an internal parameter
  Material & material;

  /// the fem containing the mesh and the element informations
  FEM & fem;

  /// Element filter if needed
  const ByElementTypeArray<UInt> & element_filter;

  /// default value
  T default_value;

  /// spatial dimension of the element to consider
  UInt spatial_dimension;

  /// ElementKind of the element to consider
  ElementKind element_kind;

  /// Number of component of the internal field
  UInt nb_component;

  /// Is the field initialized
  bool is_init;
};

/// standard output stream operator
template<typename T>
inline std::ostream & operator <<(std::ostream & stream, const InternalField<T> & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AKANTU_INTERNAL_FIELD_HH__ */
