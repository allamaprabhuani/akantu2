/**
 * @file   coupler_solid_cohesive_contact.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Thu Jan 17 2019
 * @date last modification: Thu Jan 17 2019
 *
 * @brief  class for coupling of solid mechanics and conatct mechanics
 * model in explicit
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "boundary_condition.hh"
#include "contact_mechanics_model.hh"
#include "data_accessor.hh"
#include "model.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "sparse_matrix.hh"
#include "time_step_solver.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COUPLER_SOLID_COHESIVE_CONTACT_HH__
#define __AKANTU_COUPLER_SOLID_COHESIVE_CONTACT_HH__

/* ------------------------------------------------------------------------ */
/* Coupling : Solid Mechanics / Contact Mechanics                           */
/* ------------------------------------------------------------------------ */
namespace akantu {
template <ElementKind kind, class IntegrationOrderFunctor>
class IntegratorGauss;
template <ElementKind kind> class ShapeLagrange;
class DOFManager;
} // namespace akantu

/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
class CouplerSolidCohesiveContact
    : public Model,
      public DataAccessor<Element>,
      public DataAccessor<UInt>,
      public BoundaryCondition<CouplerSolidCohesiveContact> {

  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructor                                                   */
  /* ------------------------------------------------------------------------ */

  using MyFEEngineCohesiveType =
      FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_cohesive>;
  using MyFEEngineFacetType =
      FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_regular,
                       FacetsCohesiveIntegrationOrderFunctor>;

public:
  CouplerSolidCohesiveContact(
      Mesh & mesh, UInt spatial_dimension = _all_dimensions,
      const ID & id = "coupler_solid_cohesive_contact",
      std::shared_ptr<DOFManager> dof_manager = nullptr,
      const ModelType model_type = ModelType::_coupler_solid_cohesive_contact);

  ~CouplerSolidCohesiveContact() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// initialize completely the model
  void initFullImpl(const ModelOptions & options) override;

  /// initialize the modelType
  void initModel() override;

  /// get some default values for derived classes
  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  /* ------------------------------------------------------------------------ */
  /* Solver Interface                                                         */
  /* ------------------------------------------------------------------------ */
public:
  /// assembles the contact stiffness matrix
  virtual void assembleStiffnessMatrix();

  /// assembles the contant internal forces
  virtual void assembleInternalForces();

protected:
  /// callback for the solver, this adds f_{ext} - f_{int} to the residual
  void assembleResidual() override;

  /// callback for the solver, this adds f_{ext} or  f_{int} to the residual
  void assembleResidual(const ID & residual_part) override;
  bool canSplitResidual() override { return true; }

  /// get the type of matrix needed
  MatrixType getMatrixType(const ID & matrix_id) override;

  /// callback for the solver, this assembles different matrices
  void assembleMatrix(const ID & matrix_id) override;

  /// callback for the solver, this assembles the stiffness matrix
  void assembleLumpedMatrix(const ID & matrix_id) override;

  /// callback for the solver, this is called at beginning of solve
  void predictor() override;

  /// callback for the solver, this is called at end of solve
  void corrector() override;

  /// callback for the solver, this is called at beginning of solve
  void beforeSolveStep() override;

  /// callback for the solver, this is called at end of solve
  void afterSolveStep(bool converged = true) override;

  /// callback for the model to instantiate the matricess when needed
  void initSolver(TimeStepSolverType, NonLinearSolverType) override;

  /* ------------------------------------------------------------------------ */
  /* Mass matrix for solid mechanics model                                    */
  /* ------------------------------------------------------------------------ */
public:
  /// assemble the lumped mass matrix
  void assembleMassLumped();

  /// assemble the mass matrix for consistent mass resolutions
  void assembleMass();

protected:
  /// assemble the lumped mass matrix for local and ghost elements
  void assembleMassLumped(GhostType ghost_type);

  /// assemble the mass matrix for either _ghost or _not_ghost elements
  void assembleMass(GhostType ghost_type);

protected:
  /* --------------------------------------------------------------------------
   */
  TimeStepSolverType getDefaultSolverType() const override;
  /* --------------------------------------------------------------------------
   */
  ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const;

public:
  bool isDefaultSolverExplicit() { return method == _explicit_dynamic_contact; }

  /* ------------------------------------------------------------------------ */
public:
  // DataAccessor<Element>
  UInt getNbData(const Array<Element> &,
                 const SynchronizationTag &) const override {
    return 0;
  }
  void packData(CommunicationBuffer &, const Array<Element> &,
                const SynchronizationTag &) const override {}
  void unpackData(CommunicationBuffer &, const Array<Element> &,
                  const SynchronizationTag &) override {}

  // DataAccessor<UInt> nodes
  UInt getNbData(const Array<UInt> &,
                 const SynchronizationTag &) const override {
    return 0;
  }
  void packData(CommunicationBuffer &, const Array<UInt> &,
                const SynchronizationTag &) const override {}
  void unpackData(CommunicationBuffer &, const Array<UInt> &,
                  const SynchronizationTag &) override {}

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  FEEngine & getFEEngineBoundary(const ID & name = "") override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// return the dimension of the system space
  AKANTU_GET_MACRO(SpatialDimension, Model::spatial_dimension, UInt);

  /// get the ContactMechanicsModel::displacement vector
  AKANTU_GET_MACRO(Displacement, *displacement, Array<Real> &);

  /// get  the ContactMechanicsModel::increment  vector \warn  only  consistent
  /// if ContactMechanicsModel::setIncrementFlagOn has been called before
  AKANTU_GET_MACRO(Increment, *displacement_increment, Array<Real> &);

  /// get the ContactMechanicsModel::external_force vector (external forces)
  AKANTU_GET_MACRO(ExternalForce, *external_force, Array<Real> &);

  /// get the ContactMechanicsModel::force vector (external forces)
  Array<Real> & getForce() {
    AKANTU_DEBUG_WARNING("getForce was maintained for backward compatibility, "
                         "use getExternalForce instead");
    return *external_force;
  }

  /// get the solid mechanics model
  AKANTU_GET_MACRO(SolidMechanicsModelCohesive, *solid,
                   SolidMechanicsModelCohesive &);

  /// get the contact mechanics model
  AKANTU_GET_MACRO(ContactMechanicsModel, *contact, ContactMechanicsModel &);

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

  virtual void dump(const std::string & dumper_name);

  virtual void dump(const std::string & dumper_name, UInt step);

  virtual void dump(const std::string & dumper_name, Real time, UInt step);

  void dump() override;

  virtual void dump(UInt step);

  virtual void dump(Real time, UInt step);

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
private:
  /// solid mechanics model
  SolidMechanicsModelCohesive * solid{nullptr};

  /// contact mechanics model
  ContactMechanicsModel * contact{nullptr};

  ///
  Array<Real> * displacement{nullptr};

  ///
  Array<Real> * displacement_increment{nullptr};

  /// external forces array
  Array<Real> * external_force{nullptr};

  bool search_for_contact;

  UInt step;
};

} // namespace akantu

#endif /* __COUPLER_SOLID_CONTACT_HH__  */
