/**
 * Copyright (©) 2020-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_factory.hh"
#include "data_accessor.hh"
#include "integration_point.hh"
#include "parsable.hh"
#include "parser.hh"
/* -------------------------------------------------------------------------- */
#include "internal_field.hh"
#include "random_internal_field.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PHASEFIELD_HH__
#define __AKANTU_PHASEFIELD_HH__

/* -------------------------------------------------------------------------- */
namespace akantu {
class Model;
class PhaseFieldModel;
class PhaseField;
} // namespace akantu

namespace akantu {

template <typename T>
using InternalPhaseField = InternalFieldTmpl<PhaseField, T>;

template <>
inline void
ParameterTyped<RandomInternalField<Real, InternalPhaseField>>::setAuto(
    const ParserParameter & in_param) {
  Parameter::setAuto(in_param);
  RandomParameter<Real> random_param = in_param;
  param.setRandomDistribution(random_param);
}

using PhaseFieldFactory =
    Factory<PhaseField, ID, const ID &, PhaseFieldModel &, const ID &>;

class PhaseField : public DataAccessor<Element>, public Parsable {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  PhaseField(const PhaseField & phase) = delete;
  PhaseField & operator=(const PhaseField & phase) = delete;

  /// Initialize phasefield with defaults
  PhaseField(PhaseFieldModel & model, const ID & id = "");

  /// Initialize phasefield with custom mesh & fe_engine
  PhaseField(PhaseFieldModel & model, Int dim, const Mesh & mesh,
             FEEngine & fe_engine, const ID & id = "");

  /// Destructor
  ~PhaseField() override;

protected:
  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  template <typename T>
  void registerInternal(InternalPhaseField<T> & /*vect*/) {
    AKANTU_TO_IMPLEMENT();
  }

  template <typename T>
  void unregisterInternal(InternalPhaseField<T> & /*vect*/) {
    AKANTU_TO_IMPLEMENT();
  }

  /// initialize the phasefield computed parameter
  virtual void initPhaseField();

  ///
  virtual void beforeSolveStep();

  ///
  virtual void afterSolveStep();

  /// assemble the residual for this phasefield
  virtual void assembleInternalForces(GhostType ghost_type);

  /// assemble the stiffness matrix for this phasefield
  virtual void assembleStiffnessMatrix(GhostType ghost_type);

  /// compute the driving force for this phasefield
  virtual void computeAllDrivingForces(GhostType ghost_type = _not_ghost);

  /// save the phi in the phi internal field if needed
  virtual void savePreviousState();

  /// add an element to the local mesh filter
  inline Int addElement(const Element & element);

  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override;

protected:
  /// compute the dissipated energy by element
  void computeDissipatedEnergyByElements();

  /// add an element to the local mesh filter
  inline Int addElement(const ElementType & type, Idx element,
                        const GhostType & ghost_type);

  /// resize the internals arrrays
  virtual void resizeInternals();

  /// function called to updatet the internal parameters when the
  /// modifiable parameters are modified
  virtual void updateInternalParameters();

  // constitutive law for driving force
  virtual void computeDrivingForce(ElementType /* el_type */,
                                   GhostType /* ghost_type */ = _not_ghost) {
    AKANTU_TO_IMPLEMENT();
  }

  /// compute the dissiapted energy
  virtual void computeDissipatedEnergy(ElementType el_type);

  /// compute the dissipated energy for an element
  virtual void
  computeDissipatedEnergyByElement(const Element & /*element*/,
                                   Vector<Real> & /*edis_on_quad_points*/) {
    AKANTU_TO_IMPLEMENT();
  }

  /// compute the dissipated energy for an element
  virtual void
  computeDissipatedEnergyByElement(ElementType /*type*/, Idx /*index*/,
                                   Vector<Real> & /*edis_on_quad_points*/) {
    AKANTU_TO_IMPLEMENT();
  }

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:
  inline Int getNbData(const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

  template <typename T>
  inline void packElementDataHelper(const ElementTypeMapArray<T> & data_to_pack,
                                    CommunicationBuffer & buffer,
                                    const Array<Element> & elements,
                                    const ID & fem_id = ID()) const;

  template <typename T>
  inline void unpackElementDataHelper(ElementTypeMapArray<T> & data_to_unpack,
                                      CommunicationBuffer & buffer,
                                      const Array<Element> & elements,
                                      const ID & fem_id = ID());

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
protected:
  /// return the damage energyfor the provided element
  virtual Real getEnergy(ElementType type, Idx index);

public:
  /// return the damage energyfor the subset of elements contained
  /// by the phasefield
  virtual Real getEnergy();

  /// Compute dissipated energy for an individual element
  Real getEnergy(const Element & element);

  AKANTU_GET_MACRO(Name, name, const std::string &);

  AKANTU_GET_MACRO(Model, model, const PhaseFieldModel &)

  AKANTU_GET_MACRO(ID, id, const ID &);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Strain, strain, Real);

  AKANTU_GET_MACRO(Strain, strain, const ElementTypeMapArray<Real> &);

  AKANTU_GET_MACRO_NOT_CONST(Strain, strain, ElementTypeMapArray<Real> &);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Damage, damage_on_qpoints, Real);

  AKANTU_GET_MACRO_NOT_CONST(Damage, damage_on_qpoints,
                             ElementTypeMapArray<Real> &);
  AKANTU_GET_MACRO(Damage, damage_on_qpoints,
                   const ElementTypeMapArray<Real> &);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ElementFilter, element_filter, Idx);

  AKANTU_GET_MACRO(ElementFilter, element_filter,
                   const ElementTypeMapArray<Idx> &);

  template <typename T>
  const Array<T> & getArray(const ID & id, ElementType type,
                            GhostType ghost_type = _not_ghost) const;
  template <typename T>
  Array<T> & getArray(const ID & id, ElementType type,
                      GhostType ghost_type = _not_ghost);

  template <typename T>
  const InternalPhaseField<T> & getInternal(const ID & id) const;
  template <typename T> InternalPhaseField<T> & getInternal(const ID & id);

  template <typename T>
  inline bool isInternal(const ID & id, const ElementKind & element_kind) const;

  template <typename T> inline void setParam(const ID & param, T value);
  inline const Parameter & getParam(const ID & param) const;

  template <typename T>
  void flattenInternal(const std::string & field_id,
                       ElementTypeMapArray<T> & internal_flat,
                       GhostType ghost_type = _not_ghost,
                       ElementKind element_kind = _ek_not_defined) const;

  template <typename T>
  void inflateInternal(const std::string & field_id,
                       const ElementTypeMapArray<T> & field,
                       GhostType ghost_type = _not_ghost,
                       ElementKind element_kind = _ek_not_defined);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// boolean to know if the material has been initialized
  bool is_init;

  std::map<ID, InternalPhaseField<Real> *> internal_vectors_real;
  std::map<ID, InternalPhaseField<Int> *> internal_vectors_int;
  std::map<ID, InternalPhaseField<bool> *> internal_vectors_bool;

protected:
  ID id;

  /// Link to the fem object in the model
  FEEngine & fem;

  /// phasefield name
  std::string name;

  /// The model to whch the phasefield belong
  PhaseFieldModel & model;

  /// length scale parameter
  Real l0;

  /// critical energy release rate
  // Real g_c;
  RandomInternalField<Real, InternalPhaseField> g_c;

  /// Young's modulus
  Real E;

  /// Poisson ratio
  Real nu;

  /// Isotropic formulation
  bool isotropic{true};

  /// Lame's first parameter
  Real lambda;

  /// Lame's second paramter
  Real mu;

  /// spatial dimension
  Int spatial_dimension;

  /// list of element handled by the phasefield
  ElementTypeMapArray<Idx> element_filter;

  /// damage arrays ordered by element types
  InternalPhaseField<Real> damage_on_qpoints;

  /// grad_d arrays ordered by element types
  InternalPhaseField<Real> gradd;

  /// phi arrays ordered by element types
  InternalPhaseField<Real> phi;

  /// strain arrays ordered by element types
  InternalPhaseField<Real> strain;

  /// driving force ordered by element types
  InternalPhaseField<Real> driving_force;

  /// driving energy ordered by element types
  InternalPhaseField<Real> driving_energy;

  /// damage energy ordered by element types
  InternalPhaseField<Real> damage_energy;

  /// damage energy density ordered by element types
  InternalPhaseField<Real> damage_energy_density;

  /// dissipated energy by element
  InternalPhaseField<Real> dissipated_energy;
};

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const PhaseField & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "phasefield_inline_impl.hh"

#include "internal_field_tmpl.hh"
#include "random_internal_field_tmpl.hh"

/* -------------------------------------------------------------------------- */
#define PHASEFIELD_DEFAULT_ALLOCATOR(id, phase_name)                           \
  [](const ID &, PhaseFieldModel & model,                                      \
     const ID & id) -> std::unique_ptr<PhaseField> {                           \
    return std::make_unique<phase_name>(model, id);                            \
  }

#define INSTANTIATE_PHASEFIELD(id, phase_name)                                 \
  static bool phasefield_is_allocated_##id [[gnu::unused]] =                   \
      PhaseFieldFactory::getInstance().registerAllocator(                      \
          #id, PHASEFIELD_DEFAULT_ALLOCATOR(id, phase_name))

#endif
