/**
 * @file   internal_field.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Mar 26 2021
 *
 * @brief  Constitutive law internal properties
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "element_type_map.hh"
/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_INTERNAL_FIELD_HH_
#define AKANTU_INTERNAL_FIELD_HH_

namespace akantu {
class ConstitutiveLawInternalHandler;
class FEEngine;
} // namespace akantu

namespace akantu {

class InternalFieldBase
    : public std::enable_shared_from_this<InternalFieldBase> {
public:
  InternalFieldBase(const ID & id) : id(id) {}

  /// activate the history of this field
  virtual void initializeHistory() = 0;

  /// resize the arrays and set the new element to 0
  virtual void resize() = 0;

  /// save the current values in the history
  virtual void saveCurrentValues() = 0;

  /// restore the previous values from the history
  virtual void restorePreviousValues() = 0;

  /// remove the quadrature points corresponding to suppressed elements
  virtual void
  removeIntegrationPoints(const ElementTypeMapArray<UInt> & new_numbering) = 0;

  virtual bool hasHistory() const = 0;

  AKANTU_GET_MACRO(ID, id, const ID &);

protected:
  ID id;
};

/**
 * class for the internal fields of constitutive law
 * to store values for each quadrature
 */
template <class ConstitutiveLaw_, typename T>
class InternalFieldTmpl : public InternalFieldBase,
                          public ElementTypeMapArray<T> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using ConstitutiveLaw = ConstitutiveLaw_;

  InternalFieldTmpl(const ID & id, ConstitutiveLaw & constitutive_law);
  ~InternalFieldTmpl() override;

  /// This constructor is only here to let cohesive elements compile
  InternalFieldTmpl(const ID & id, ConstitutiveLaw & constitutive_law,
                    const ID & fem_id,
                    const ElementTypeMapArray<UInt> & element_filter);

  /// More general constructor
  InternalFieldTmpl(const ID & id, ConstitutiveLaw & constitutive_law, UInt dim,
                    const ID & fem_id,
                    const ElementTypeMapArray<UInt> & element_filter);

  InternalFieldTmpl(const ID & id,
                    const InternalFieldTmpl<ConstitutiveLaw_, T> & other);

  InternalFieldTmpl operator=(const InternalFieldTmpl &) = delete;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// function to reset the FEEngine for the internal fieldx
  //  virtual void setFEEngine(FEEngine & fe_engine);

  /// function to reset the element kind for the internal
  virtual void setElementKind(ElementKind element_kind);

  /// initialize the field to a given number of component
  virtual void initialize(UInt nb_component);

  /// activate the history of this field
  void initializeHistory() override;

  /// resize the arrays and set the new element to 0
  void resize() override;

  /// set the field to a given value v
  virtual void setDefaultValue(const T & v);

  /// reset all the fields to the default value
  virtual void reset();

  /// save the current values in the history
  void saveCurrentValues() override;

  /// restore the previous values from the history
  void restorePreviousValues() override;

  /// remove the quadrature points corresponding to suppressed elements
  void removeIntegrationPoints(
      const ElementTypeMapArray<UInt> & new_numbering) override;

  /// print the content
  void printself(std::ostream & stream, int /*indent*/ = 0) const override;

  /// get the default value
  inline operator T() const;

  virtual FEEngine & getFEEngine() { return fem; }

  virtual const FEEngine & getFEEngine() const { return fem; }

protected:
  /// initialize the arrays in the ElementTypeMapArray<T>
  void internalInitialize(UInt nb_component);

  /// set the values for new internals
  virtual void setArrayValues(T * begin, T * end);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get filter types for range loop
  decltype(auto) elementTypes(GhostType ghost_type = _not_ghost) const {
    return ElementTypeMapArray<T>::elementTypes(
        _spatial_dimension = this->spatial_dimension,
        _element_kind = this->element_kind, _ghost_type = ghost_type);
  }

  /// get filter types for range loop
  decltype(auto) filterTypes(GhostType ghost_type = _not_ghost) const {
    return this->element_filter.elementTypes(
        _spatial_dimension = this->spatial_dimension,
        _element_kind = this->element_kind, _ghost_type = ghost_type);
  }

  /// get the array for a given type of the element_filter
  const Array<UInt> & getFilter(ElementType type,
                                GhostType ghost_type = _not_ghost) const {
    return this->element_filter(type, ghost_type);
  }

  /// get the Array corresponding to the type en ghost_type specified
  virtual Array<T> & operator()(ElementType type,
                                GhostType ghost_type = _not_ghost) {
    return ElementTypeMapArray<T>::operator()(type, ghost_type);
  }

  virtual const Array<T> & operator()(ElementType type,
                                      GhostType ghost_type = _not_ghost) const {
    return ElementTypeMapArray<T>::operator()(type, ghost_type);
  }

  virtual Array<T> & previous(ElementType type,
                              GhostType ghost_type = _not_ghost) {
    AKANTU_DEBUG_ASSERT(previous_values != nullptr,
                        "The history of the internal "
                            << this->getID() << " has not been activated");
    return this->previous_values->operator()(type, ghost_type);
  }

  virtual const Array<T> & previous(ElementType type,
                                    GhostType ghost_type = _not_ghost) const {
    AKANTU_DEBUG_ASSERT(previous_values != nullptr,
                        "The history of the internal "
                            << this->getID() << " has not been activated");
    return this->previous_values->operator()(type, ghost_type);
  }

  virtual InternalFieldTmpl<ConstitutiveLaw_, T> & previous() {
    AKANTU_DEBUG_ASSERT(previous_values != nullptr,
                        "The history of the internal "
                            << this->getID() << " has not been activated");
    return *(this->previous_values);
  }

  virtual const InternalFieldTmpl<ConstitutiveLaw_, T> & previous() const {
    AKANTU_DEBUG_ASSERT(previous_values != nullptr,
                        "The history of the internal "
                            << this->getID() << " has not been activated");
    return *(this->previous_values);
  }

  /// check if the history is used or not
  bool hasHistory() const override { return (previous_values != nullptr); }

  /// get the kind treated by the internal
  ElementKind getElementKind() const { return element_kind; }

  /// return the number of components
  UInt getNbComponent() const { return nb_component; }

  /// return the spatial dimension corresponding to the internal element type
  /// loop filter
  UInt getSpatialDimension() const { return this->spatial_dimension; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// the constitutive_law for which this is an internal parameter
  ConstitutiveLaw & constitutive_law;

  /// the fem containing the mesh and the element informations
  FEEngine & fem;

  /// Element filter if needed
  const ElementTypeMapArray<UInt> & element_filter;

  /// default value
  T default_value{};

  /// spatial dimension of the element to consider
  UInt spatial_dimension{0};

  /// ElementKind of the element to consider
  ElementKind element_kind{_ek_regular};

  /// Number of component of the internal field
  UInt nb_component{0};

  /// Is the field initialized
  bool is_init{false};

  /// previous values
  std::unique_ptr<InternalFieldTmpl<ConstitutiveLaw_, T>> previous_values;
};

/// standard output stream operator
template <class ConstitutiveLaw_, typename T>
inline std::ostream &
operator<<(std::ostream & stream,
           const InternalFieldTmpl<ConstitutiveLaw_, T> & _this) {
  _this.printself(stream);
  return stream;
}

template <typename T>
using InternalField = InternalFieldTmpl<ConstitutiveLaw_, T>;

} // namespace akantu

#endif /* AKANTU_INTERNAL_FIELD_HH_ */
