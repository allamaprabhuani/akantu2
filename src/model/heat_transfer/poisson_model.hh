/**
 * @file   poisson_model.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Srinivasa Babu Ramisetti <srinivasa.ramisetti@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Rui Wang <rui.wang@epfl.ch>
 *
 * @date creation: Sun May 01 2011
 * @date last modification: Mon Mar 15 2021
 *
 * @brief  Model of Generic Poisson Equation
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
#include "boundary_condition.hh"
#include "data_accessor.hh"
#include "fe_engine.hh"
#include "model.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_POISSON_MODEL_HH_
#define AKANTU_POISSON_MODEL_HH_

namespace akantu {
class ConstitutiveLaw;
class ConstitutiveLawSelector;  
template <ElementKind kind, class IntegrationOrderFunctor> class IntegratorGauss;
template <ElementKind kind> class ShapeLagrange;
} // namespace akantu

namespace akantu {

class PoissonModel : public Model,
		     public DataAccessor<Element>,
		     public DataAccessor<UInt>,
		     public BoundaryCondition<PoissonModel> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using FEEngineType = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;

  PoissonModel(Mesh & mesh, UInt dim = _all_dimensions,
	       const ID & id = "poisson_model",
	       std::shared_ptr<DOFManager> dof_manager = nullptr,
	       ModelType model_type = ModelType::_poisson_model);

  ~PoissonModel() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// generic function to initialize everything ready for explicit dynamics
  void initFullImpl(const ModelOptions & options) override;

  /// read one material file to instantiate all the materials
  void initConstitutiveLaws();

  /// allocate all vectors
  void initSolver(TimeStepSolverType time_step_solver_type,
                  NonLinearSolverType non_linear_solver_type) override;

  /// initialize the model
  void initModel() override;

  void predictor() override;

  /// callback for the solver, this is called at end of solve
  void corrector() override;

  /// compute the heat flux
  void assembleResidual() override;

  /// get the type of matrix needed
  MatrixType getMatrixType(const ID & matrix_id) const override;

  /// callback to assemble a Matrix
  void assembleMatrix(const ID & matrix_id) override;

  /// callback to assemble a lumped Matrix
  void assembleLumpedMatrix(const ID & matrix_id) override;

  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const override;
  
  /// function to print the containt of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  /// 
  TimeStepSolverType getDefaultSolverType() const override;

  /* ------------------------------------------------------------------------ */
  /* Constitutive laws                                                        */
  /* ------------------------------------------------------------------------ */
public:
  /// register an empty constitutive law of a given type
  ConstitutiveLaw & registerNewConstitutiveLaw(const ID & phase_name,
					       const ID & phase_type,
					       const ID & opt_param);

  /// reassigns constitutive Laws depending on the constitutive law selector
  void reassignConstitutiveLaw();

protected:
  /// register a constitutive law in the dynamic database
  ConstitutiveLaw & registerNewConstitutiveLaw(const ParserSection & phase_section);

  /// read the constitutive laws to instantiate all the constitutive_laws
  void instantiateConstitutiveLaws();

  /// set the element_id_by_constitutive_law and add the elements to the good
  /// constitutive law
  void assignConstitutiveLawToElements(
      const ElementTypeMapArray<UInt> * filter = nullptr);
  
  /* ------------------------------------------------------------------------ */
  /* Methods for explicit                                                     */
  /* ------------------------------------------------------------------------ */
public:
  /// compute and get the stable time step
  Real getStableTimeStep();

  /// set the stable timestep
  void setTimeStep(Real time_step, const ID & solver_id = "") override;

protected:
  /// callback for the solver, this is called at beginning of solve
  void beforeSolveStep() override;
  /// callback for the solver, this is called at end of solve
  void afterSolveStep(bool converged = true) override;
  
  /// compute the stable time step
  Real getStableTimeStep(GhostType ghost_type);

  
public:
  /// compute the internal dof rate \todo Need code review: currently not
  /// public method
  void assembleInternalDofRate();

  /// assemble the stiffness matrix
  void assembleStiffnessMatrix(bool need_to_reassemble = false);


public:
  /// calculate the lumped capacity vector for heat transfer problem
  void assembleCapacityLumped();

  /// assemble the capacity matrix
  void assembleCapacity();


public:
  /// assemble the lumped capacity matrix for local and ghost elements
  void assembleCapacityLumped(GhostType ghost_type);

  /// assemble the capacity matrix for either _ghost or _not_ghost elements
  void assembleCapacity(GhostType ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:
  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const override;
  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;
  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

  inline UInt getNbData(const Array<UInt> & indexes,
                        const SynchronizationTag & tag) const override;
  inline void packData(CommunicationBuffer & buffer,
                       const Array<UInt> & indexes,
                       const SynchronizationTag & tag) const override;
  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<UInt> & indexes,
                         const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
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
                       UInt spatial_dimension, ElementKind kind) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the current value of the time step
  AKANTU_GET_MACRO(TimeStep, time_step, Real);
  /// get the assembled heat flux
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(InternalDofRate, internal_dof_rate);
  /// get the assembled heat flux
  AKANTU_GET_MACRO_DEREF_PTR(InternalDofRate, internal_dof_rate); 
  /// get the external dof rate vector
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(ExternalDofRate, external_dof_rate);
  /// get the external dof rate vector
  AKANTU_GET_MACRO_DEREF_PTR(ExternalDofRate, external_dof_rate);
  /// get the degree of freedom
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(Dof, dof);
  /// get the degree of freedom
  AKANTU_GET_MACRO_DEREF_PTR(Dof, dof); 
  /// get the dof rate
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(DofRate, dof_rate);
  /// get the dof rate
  AKANTU_GET_MACRO_DEREF_PTR(DofRate, dof_rate);
  /// get the PoissonModel::blocked_dofs array
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(BlockedDOFs, blocked_dofs);
  /// get the blocked_dofs vector
  AKANTU_GET_MACRO_DEREF_PTR(BlockedDOFs, blocked_dofs);

  /// get an iterable on the constitutive_laws
  inline decltype(auto) getConstitutiveLaws();
  
  /// get an iterable on the constitutive_laws
  inline decltype(auto) getConstitutiveLaws() const;

  /// get a particular constitutive law (by constitutive law index)
  inline ConstitutiveLaw & getConstitutiveLaw(UInt mat_index);

  /// get a particular constitutive law (by constitutive law index)
  inline const ConstitutiveLaw & getConstitutiveLaw(UInt mat_index) const;

  /// get a particular constitutive law (by constitutive law name)
  inline ConstitutiveLaw & getConstitutiveLaw(const std::string & name);

  /// get a particular constitutive law (by constitutive law name)
  inline const ConstitutiveLaw & getConstitutiveLaw(const std::string & name) const;

  /// get a particular constitutive law id from is name
  inline UInt getConstitutiveLawIndex(const std::string & name) const;

  /// give the number of constitutive_laws
  inline UInt getNbConstitutiveLaws() const { return constitutive_laws.size(); }

  /// give the constitutive law internal index from its id
  Int getInternalIndexFromID(const ID & id) const;


  AKANTU_GET_MACRO(ConstitutiveLawByElement, constitutive_law_index,
                   const ElementTypeMapArray<UInt> &);
  AKANTU_GET_MACRO(ConstitutiveLawLocalNumbering, constitutive_law_local_numbering,
                   const ElementTypeMapArray<UInt> &);

  /// vectors containing local material element index for each global element
  /// index
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ConstitutiveLawByElement, constitutive_law_index,
                                         UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(ConstitutiveLawByElement, constitutive_law_index, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ConstitutiveLawLocalNumbering,
                                         constitutive_law_local_numbering, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(ConstitutiveLawLocalNumbering,
                                   constitutive_law_local_numbering, UInt);

  AKANTU_GET_MACRO_NOT_CONST(ConstitutiveLawSelector, *constitutive_law_selector,
                             ConstitutiveLawSelector &);

  AKANTU_SET_MACRO(ConstitutiveLawSelector, constitutive_law_selector,
                   std::shared_ptr<ConstitutiveLawSelector>);

  FEEngine & getFEEngineBoundary(const ID & name = "") override;


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// time step
  Real time_step;

  /// temperatures array
  std::unique_ptr<Array<Real>> dof;

  /// temperatures derivatives array
  std::unique_ptr<Array<Real>> dof_rate;

  /// increment array (@f$\delta \dot T@f$ or @f$\delta T@f$)
  std::unique_ptr<Array<Real>> increment;

  /// external flux vector
  std::unique_ptr<Array<Real>> external_dof_rate;

  /// residuals array
  std::unique_ptr<Array<Real>> internal_dof_rate;

  /// boundary vector
  std::unique_ptr<Array<bool>> blocked_dofs;

    /// Arrays containing the constitutive law index for each element
  ElementTypeMapArray<UInt> constitutive_law_index;

  /// Arrays containing the position in the element filter of the constitutive law
  /// (constitutive law's local numbering)
  ElementTypeMapArray<UInt> constitutive_law_local_numbering;

  /// class defining of to choose a constitutive law
  std::shared_ptr<ConstitutiveLawSelector> constitutive_law_selector;

  /// mapping between constitutive Law name and law internal id
  std::map<std::string, UInt> constitutive_laws_names_to_id;

  /// list of used constitutive laws
  std::vector<std::unique_ptr<ConstitutiveLaw>> constitutive_laws;

  /// tells if the constitutive law are instantiated
  bool are_constitutive_laws_instantiated{false};


  
  bool need_to_reassemble_capacity{true};

  bool need_to_reassemble_capacity_lumped{true};

  UInt dof_release{0};

  UInt conductivity_matrix_release{UInt(-1)};

};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "parser.hh"
#include "constitutive_law.hh"

#include "poisson_model_inline_impl.hh"

#endif /* AKANTU_POISSON_MODEL_HH_ */
