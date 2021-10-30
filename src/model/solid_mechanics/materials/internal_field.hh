/**
 * @file   internal_field.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Mar 26 2021
 *
 * @brief  Material internal properties
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

#ifndef AKANTU_INTERNAL_FIELD_HH_
#define AKANTU_INTERNAL_FIELD_HH_

namespace akantu {

class Material;
class FEEngine;

/**
 * class for the internal fields of materials
 * to store values for each quadrature
 */
template <class Material, typename T>
class InternalFieldTmpl : public ElementTypeMapArray<T> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  InternalFieldTmpl(const ID & id, Material & material);
  ~InternalFieldTmpl() override;

  /// This constructor is only here to let cohesive elements compile
  InternalFieldTmpl(const ID & id, Material & material, FEEngine & fem,
                    const ElementTypeMapArray<Idx> & element_filter);

  /// More general constructor
  InternalFieldTmpl(const ID & id, Material & material, Int dim, FEEngine & fem,
                    const ElementTypeMapArray<Idx> & element_filter);

  InternalFieldTmpl(const ID & id,
                    const InternalFieldTmpl<Material, T> & other);

  auto operator=(const InternalFieldTmpl &) -> InternalFieldTmpl = delete;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// function to reset the FEEngine for the internal field
  virtual void setFEEngine(FEEngine & fe_engine);

  /// function to reset the element kind for the internal
  virtual void setElementKind(ElementKind element_kind);

  /// initialize the field to a given number of component
  virtual void initialize(UInt nb_component);

  /// activate the history of this field
  virtual void initializeHistory();

  /// resize the arrays and set the new element to 0
  virtual void resize();

  /// set the field to a given value v
  virtual void setDefaultValue(const T & v);

  /// reset all the fields to the default value
  virtual void reset();

  /// save the current values in the history
  virtual void saveCurrentValues();

  /// restore the previous values from the history
  virtual void restorePreviousValues();

  /// remove the quadrature points corresponding to suppressed elements
  virtual void
  removeIntegrationPoints(const ElementTypeMapArray<Int> & new_numbering);

  /// print the content
  void printself(std::ostream & stream, int /*indent*/ = 0) const override;

  /// get the default value
  inline operator T() const;

  virtual auto getFEEngine() -> FEEngine & { return *fem; }

  virtual auto getFEEngine() const -> const FEEngine & { return *fem; }

protected:
  /// initialize the arrays in the ElementTypeMapArray<T>
  void internalInitialize(Int nb_component);

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
  decltype(auto) getFilter(ElementType type,
                           GhostType ghost_type = _not_ghost) const {
    return (this->element_filter(type, ghost_type));
  }

  /// get the Array corresponding to the type en ghost_type specified
  virtual auto operator()(ElementType type, GhostType ghost_type = _not_ghost)
      -> Array<T> & {
    return ElementTypeMapArray<T>::operator()(type, ghost_type);
  }

  virtual auto operator()(ElementType type,
                          GhostType ghost_type = _not_ghost) const
      -> const Array<T> & {
    return ElementTypeMapArray<T>::operator()(type, ghost_type);
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

  virtual auto previous() -> InternalFieldTmpl & {
    AKANTU_DEBUG_ASSERT(previous_values != nullptr,
                        "The history of the internal "
                            << this->getID() << " has not been activated");
    return *(this->previous_values);
  }

  virtual auto previous() const -> const InternalFieldTmpl & {
    AKANTU_DEBUG_ASSERT(previous_values != nullptr,
                        "The history of the internal "
                            << this->getID() << " has not been activated");
    return *(this->previous_values);
  }

  /// check if the history is used or not
  auto hasHistory() const -> bool { return (previous_values != nullptr); }

  /// get the kind treated by  the internal
  AKANTU_GET_MACRO_AUTO(ElementKind, element_kind);

  /// return the number of components
  AKANTU_GET_MACRO_AUTO(NbComponent, nb_component);

  /// return the spatial dimension corresponding to the internal element type
  /// loop filter
  AKANTU_GET_MACRO_AUTO(SpatialDimension, spatial_dimension);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// the material for which this is an internal parameter
  Material & material;

  /// the fem containing the mesh and the element informations
  FEEngine * fem{nullptr};

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
  std::unique_ptr<InternalFieldTmpl<Material, T>> previous_values;
};

/// standard output stream operator
template <class Material, typename T>
inline auto operator<<(std::ostream & stream, const InternalFieldTmpl<Material, T> & _this)
    -> std::ostream & {
  _this.printself(stream);
  return stream;
}

template <typename T> using InternalField = InternalFieldTmpl<Material, T>;

} // namespace akantu

#endif /* AKANTU_INTERNAL_FIELD_HH_ */
