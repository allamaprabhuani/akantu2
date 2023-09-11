/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "boundary_condition.hh"
#include "data_accessor.hh"
#include "fe_engine.hh"
#include "model.hh"
/* -------------------------------------------------------------------------- */
#include <array>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PHASE_FIELD_MODEL_HH__
#define __AKANTU_PHASE_FIELD_MODEL_HH__

namespace akantu {
class PhaseField;
class PhaseFieldSelector;
template <ElementKind kind, class IntegrationOrderFuntor> class IntegratorGauss;
template <ElementKind kind> class ShapeLagrange;
} // namespace akantu

/* -------------------------------------------------------------------------- */
namespace akantu {

/* -------------------------------------------------------------------------- */
class PhaseFieldModel : public Model,
                        public DataAccessor<Element>,
                        public DataAccessor<Idx>,
                        public BoundaryCondition<PhaseFieldModel> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using FEEngineType = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;

  PhaseFieldModel(Mesh & mesh, Int dim = _all_dimensions,
                  const ID & id = "phase_field_model",
                  std::shared_ptr<DOFManager> dof_manager = nullptr,
                  ModelType model_type = ModelType::_phase_field_model);

  ~PhaseFieldModel() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// generic function to initialize everything ready for explicit dynamics
  void initFullImpl(const ModelOptions & options) override;

  /// initialize all internal array for phasefields
  void initPhaseFields();

  /// allocate all vectors
  void initSolver(TimeStepSolverType /*unused*/,
                  NonLinearSolverType /*unused*/) override;

  /// initialize the model
  void initModel() override;

  /// predictor
  void predictor() override;

  /// corrector
  void corrector() override;

  /// compute the heat flux
  void assembleResidual() override;

  /// get the type of matrix needed
  MatrixType getMatrixType(const ID & /*unused*/) const override;

  /// callback to assemble a Matrix
  void assembleMatrix(const ID & /*unused*/) override;

  /// callback to assemble a lumped Matrix
  void assembleLumpedMatrix(const ID & /*unused*/) override;

  /// function to print the containt of the class
  void printself(std::ostream & stream, int indent = 0) const override;

protected:
  /* ------------------------------------------------------------------------ */
  TimeStepSolverType getDefaultSolverType() const override;

  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const override;

  /* ------------------------------------------------------------------------ */
  /* Materials (phase_field_model.cc)                                         */
  /* ------------------------------------------------------------------------ */
public:
  /// register an empty phasefield of a given type
  PhaseField & registerNewPhaseField(const ID & phase_name,
                                     const ID & phase_type,
                                     const ID & opt_param);

  /// reassigns phasefields depending on the phasefield selector
  void reassignPhaseField();

protected:
  /// register a phasefield in the dynamic database
  PhaseField & registerNewPhaseField(const ParserSection & phase_section);

  /// read the phasefield files to instantiate all the phasefields
  void instantiatePhaseFields();

  /// set the element_id_by_phasefield and add the elements to the good
  /// phasefields
  void
  assignPhaseFieldToElements(const ElementTypeMapArray<Idx> * filter = nullptr);

  /* ------------------------------------------------------------------------ */
  /* Methods for static                                                       */
  /* ------------------------------------------------------------------------ */
public:
  /// assembles the phasefield stiffness matrix
  virtual void assembleStiffnessMatrix();

  /// compute the internal forces
  virtual void assembleInternalForces();

  // compute the internal forces
  void assembleInternalForces(GhostType ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Methods for dynamic                                                      */
  /* ------------------------------------------------------------------------ */
public:
  /// set the stable timestep
  void setTimeStep(Real time_step, const ID & solver_id = "") override;

protected:
  /// callback for the solver, this is called at beginning of solve
  void beforeSolveStep() override;
  /// callback for the solver, this is called at end of solve
  void afterSolveStep(bool converged = true) override;

  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:
  Int getNbData(const Array<Element> & elements,
                const SynchronizationTag & tag) const override;

  void packData(CommunicationBuffer & buffer, const Array<Element> & elements,
                const SynchronizationTag & tag) const override;

  void unpackData(CommunicationBuffer & buffer, const Array<Element> & elements,
                  const SynchronizationTag & tag) override;

  Int getNbData(const Array<Idx> & indexes,
                const SynchronizationTag & tag) const override;

  void packData(CommunicationBuffer & buffer, const Array<Idx> & indexes,
                const SynchronizationTag & tag) const override;

  void unpackData(CommunicationBuffer & buffer, const Array<Idx> & indexes,
                  const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
protected:
  /// Compute dissipated energy for an element type and phasefield index
  Real getEnergy(ElementType type, Idx index);

public:
  /// return the damage array
  AKANTU_GET_MACRO_DEREF_PTR(Damage, damage);

  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(Damage, damage);

  /// get the PhaseFieldModel::internal_force vector (internal forces)
  AKANTU_GET_MACRO_DEREF_PTR(InternalForce, internal_force);

  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(InternalForce, internal_force);

  /// get the PhaseFieldModel::external_force vector (external forces)
  AKANTU_GET_MACRO_DEREF_PTR(ExternalForce, external_force);

  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(ExternalForce, external_force);

  /// get the PhaseFieldModel::force vector (external forces)
  Array<Real> & getForce() {
    AKANTU_DEBUG_WARNING("getForce was maintained for backward compatibility, "
                         "use getExternalForce instead");
    return *external_force;
  }

  /// get the SolidMechanicsModel::blocked_dofs array
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(BlockedDOFs, blocked_dofs);

  /// get the PhaseFieldModel::blocked_dofs vector
  AKANTU_GET_MACRO_DEREF_PTR(BlockedDOFs, blocked_dofs);

  /// get an iterable on the phasefields
  inline decltype(auto) getPhaseFields();

  /// get an iterable on the phasefields
  inline decltype(auto) getPhaseFields() const;

  /// get a particular phasefield (by phasefield index)
  inline PhaseField & getPhaseField(Idx mat_index);

  /// get a particular phasefield (by phasefield index)
  inline const PhaseField & getPhaseField(Idx mat_index) const;

  /// get a particular phasefield (by phasefield name)
  inline PhaseField & getPhaseField(const std::string & name);

  /// get a particular phasefield (by phasefield name)
  inline const PhaseField & getPhaseField(const std::string & name) const;

  /// get a particular phasefield id from is name
  inline Idx getPhaseFieldIndex(const std::string & name) const;

  /// give the number of phasefields
  inline Idx getNbPhaseFields() const { return phasefields.size(); }

  /// give the phasefield internal index from its id
  Idx getInternalIndexFromID(const ID & id) const;

  /**
   * @brief Returns the total dissipated energy
   *
   */
  Real getEnergy();

  /// Compute dissipated energy for an individual element
  Real getEnergy(const Element & element) {
    return getEnergy(element.type, element.element);
  }

  /// Compute dissipated energy for an element group
  Real getEnergy(const ID & group_id);

  AKANTU_GET_MACRO(PhaseFieldByElement, phasefield_index,
                   const ElementTypeMapArray<Idx> &);
  AKANTU_GET_MACRO(PhaseFieldLocalNumbering, phasefield_local_numbering,
                   const ElementTypeMapArray<Idx> &);

  /// vectors containing local material element index for each global element
  /// index
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(PhaseFieldByElement, phasefield_index,
                                         Idx);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(PhaseFieldByElement, phasefield_index, Idx);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(PhaseFieldLocalNumbering,
                                         phasefield_local_numbering, Idx);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(PhaseFieldLocalNumbering,
                                   phasefield_local_numbering, Idx);

  AKANTU_GET_MACRO_NOT_CONST(PhaseFieldSelector, *phasefield_selector,
                             PhaseFieldSelector &);

  AKANTU_SET_MACRO(PhaseFieldSelector, phasefield_selector,
                   std::shared_ptr<PhaseFieldSelector>);

  FEEngine & getFEEngineBoundary(const ID & name = "") override;

  /* ------------------------------------------------------------------------ */
  /* Dumpable Interface                                                       */
  /* ------------------------------------------------------------------------ */
public:
  std::shared_ptr<dumpers::Field>
  createNodalFieldReal(const std::string & field_name,
                       const std::string & group_name,
                       bool padding_flag) override;

  std::shared_ptr<dumpers::Field>
  createNodalFieldBool(const std::string & field_name,
                       const std::string & group_name,
                       bool padding_flag) override;

  std::shared_ptr<dumpers::Field>
  createElementalField(const std::string & field_name,
                       const std::string & group_name, bool padding_flag,
                       Int spatial_dimension, ElementKind kind) override;

  //! flatten a given phasefield internal field
  ElementTypeMapArray<Real> &
  flattenInternal(const std::string & field_name, ElementKind kind,
                  GhostType ghost_type = _not_ghost);

  //! inverse operation of the flatten
  void inflateInternal(const std::string & field_name,
                       const ElementTypeMapArray<Real> & field,
                       ElementKind kind, GhostType ghost_type = _not_ghost);

  /// Compute strain tensors from grad_u tensors
  void computeStrain(GhostType ghost_type);

  /// Save previous damage
  void savePreviousDamage();

  /// Save previous damage
  void savePreviousState();

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// number of iterations
  Int n_iter;

  /// damage array
  std::unique_ptr<Array<Real>> damage;

  /// damage array at the previous time step
  std::unique_ptr<Array<Real>> previous_damage;

  /// boundary vector
  std::unique_ptr<Array<bool>> blocked_dofs;

  /// external force vector
  std::unique_ptr<Array<Real>> external_force;

  /// residuals array
  std::unique_ptr<Array<Real>> internal_force;

  /// Arrays containing the phasefield index for each element
  ElementTypeMapArray<Idx> phasefield_index;

  /// Arrays containing the position in the element filter of the phasefield
  /// (phasefield's local numbering)
  ElementTypeMapArray<Idx> phasefield_local_numbering;

  /// class defining of to choose a phasefield
  std::shared_ptr<PhaseFieldSelector> phasefield_selector;

  /// mapping between phasefield name and phasefield internal id
  std::map<std::string, Idx> phasefields_names_to_id;

  /// list of used phasefields
  std::vector<std::unique_ptr<PhaseField>> phasefields;

  using flatten_internal_map =
      std::map<std::pair<std::string, ElementKind>,
               std::unique_ptr<ElementTypeMapArray<Real>>>;

  /// tells if the phasefields are instantiated
  flatten_internal_map registered_internals;

  /// tells if the phasefield are instantiated
  bool are_phasefields_instantiated{false};
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "parser.hh"
#include "phasefield.hh"

#include "phase_field_model_inline_impl.hh"
/* -------------------------------------------------------------------------- */

#endif
