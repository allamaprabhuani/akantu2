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

class InternalFieldBase {
public:
  InternalFieldBase(const ID & id) : id_(id) {}

  virtual ~InternalFieldBase() = default;

  /* ------------------------------------------------------------------------ */
  InternalFieldBase(const InternalFieldBase & /*other*/) = default;
  InternalFieldBase(InternalFieldBase && /*other*/) = default;
  InternalFieldBase & operator=(const InternalFieldBase & /*other*/) = default;
  InternalFieldBase & operator=(InternalFieldBase && /*other*/) = default;
  /* ------------------------------------------------------------------------ */

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
  removeIntegrationPoints(const ElementTypeMapArray<Idx> & new_numbering) = 0;

  [[nodiscard]] virtual bool hasHistory() const = 0;

  [[nodiscard]] auto getRegisterID() const { return id_; }

protected:
  ID id_;
};

/**
 * class for the internal fields of constitutive law
 * to store values for each quadrature
 */
template <typename T>
class InternalField : public InternalFieldBase, public ElementTypeMapArray<T> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  InternalField(const ID & id,
                ConstitutiveLawInternalHandler & constitutive_law);
  /// This constructor is only here to let cohesive elements compile
  InternalField(const ID & id,
                ConstitutiveLawInternalHandler & constitutive_law,
                const ID & fem_id,
                const ElementTypeMapArray<Idx> & element_filter);

  /// More general constructor
  InternalField(const ID & id,
                ConstitutiveLawInternalHandler & constitutive_law, Int dim,
                const ID & fem_id,
                const ElementTypeMapArray<Idx> & element_filter);

  InternalField(const ID & id, const InternalField<T> & other);

  friend class ConstitutiveLawInternalHandler;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// initialize the field to a given number of component
  virtual void initialize(Int nb_component);

public:
  /// function to reset the element kind for the internal
  virtual void setElementKind(ElementKind element_kind);

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
      const ElementTypeMapArray<Idx> & new_numbering) override;

  /// print the content
  void printself(std::ostream & stream, int /*indent*/ = 0) const override;

  /// get the default value
  inline operator T() const;

  virtual auto getFEEngine() -> FEEngine & { return fem; }

  [[nodiscard]] virtual auto getFEEngine() const -> const FEEngine & {
    return fem;
  }

protected:
  /// initialize the arrays in the ElementTypeMapArray<T>
  void internalInitialize(Int nb_component);

  /// set the values for new internals
  virtual void setArrayValues(T * begin, T * end);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
protected:
  /// get filter types for range loop
  auto elementTypesImpl(Int /*dim*/ = _all_dimensions,
                        GhostType ghost_type = _not_ghost,
                        ElementKind /*kind*/ = _ek_not_defined) const ->
      typename ElementTypeMapArray<T>::ElementTypesIteratorHelper override {
    return ElementTypeMapArray<T>::elementTypesImpl(
        this->spatial_dimension, ghost_type, this->element_kind);
  }

public:
  /// get filter types for range loop
  decltype(auto) filterTypes(GhostType ghost_type = _not_ghost) const {
    return this->element_filter.elementTypes(
        _spatial_dimension = this->spatial_dimension,
        _element_kind = this->element_kind, _ghost_type = ghost_type);
  }

  /// get the array for a given type of the element_filter
  decltype(auto) getFilter(ElementType type,
                           GhostType ghost_type = _not_ghost) const {
    return (this->element_filter(type, ghost_type));
  }

  virtual auto previous(ElementType type, GhostType ghost_type = _not_ghost)
      -> Array<T> & {
    AKANTU_DEBUG_ASSERT(previous_values != nullptr,
                        "The history of the internal "
                            << this->getID() << " has not been activated");
    return this->previous_values->operator()(type, ghost_type);
  }

  virtual auto previous(ElementType type,
                        GhostType ghost_type = _not_ghost) const
      -> const Array<T> & {
    AKANTU_DEBUG_ASSERT(previous_values != nullptr,
                        "The history of the internal "
                            << this->getID() << " has not been activated");
    return this->previous_values->operator()(type, ghost_type);
  }

  virtual InternalField<T> & previous() {
    AKANTU_DEBUG_ASSERT(previous_values != nullptr,
                        "The history of the internal "
                            << this->getID() << " has not been activated");
    return *(this->previous_values);
  }

  virtual const InternalField<T> & previous() const {
    AKANTU_DEBUG_ASSERT(previous_values != nullptr,
                        "The history of the internal "
                            << this->getID() << " has not been activated");
    return *(this->previous_values);
  }

  /// check if the history is used or not
  [[nodiscard]] auto hasHistory() const -> bool override {
    return (previous_values != nullptr);
  }

  /// get the kind treated by  the internal
  AKANTU_GET_MACRO_AUTO(ElementKind, element_kind);

  /// return the number of components
  AKANTU_GET_MACRO_AUTO(NbComponent, nb_component);

  /// return the spatial dimension corresponding to the internal element type
  /// loop filter
  AKANTU_GET_MACRO_AUTO(SpatialDimension, spatial_dimension);

  Int & getRelease(ElementType type, GhostType ghost_type) {
    return releases(type, ghost_type);
  }

  Int getRelease(ElementType type, GhostType ghost_type) const {
    return releases(type, ghost_type);
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// the constitutive_law for which this is an internal parameter
  ConstitutiveLawInternalHandler & constitutive_law;

  /// the fem containing the mesh and the element informations
  FEEngine & fem;

  /// Element filter if needed
  const ElementTypeMapArray<Int> & element_filter;

  /// default value
  T default_value{};

  /// spatial dimension of the element to consider
  Int spatial_dimension{0};

  /// ElementKind of the element to consider
  ElementKind element_kind{_ek_regular};

  /// Number of component of the internal field
  Int nb_component{0};

  /// Is the field initialized
  bool is_init{false};

  /// previous values
  std::shared_ptr<InternalField<T>> previous_values;

  ElementTypeMap<Int> releases;
};

/// standard output stream operator
template <typename T>
inline std::ostream & operator<<(std::ostream & stream,
                                 const InternalField<T> & _this) {
  _this.printself(stream);
  return stream;
}

// template <typename T> using InternalField = InternalFieldTmpl<Material, T>;

} // namespace akantu

#endif /* AKANTU_INTERNAL_FIELD_HH_ */
