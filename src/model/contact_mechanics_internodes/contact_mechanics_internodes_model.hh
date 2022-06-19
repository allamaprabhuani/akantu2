/**
 * @file   contact_mechanics_internodes_model.cc
 *
 * @author Moritz Waldleben <moritz.waldleben@epfl.ch>
 *
 * @date creation: Thu Jul 09 2022
 * @date last modification: Thu Jul 17 2022
 *
 * @brief Model for Contact Mechanics Internodes
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
#include "contact_detector_internodes.hh"
#include "solid_mechanics_model.hh"
#include "model.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONTACT_MECHANICS_INTERNODES_MODEL_HH__
#define __AKANTU_CONTACT_MECHANICS_INTERNODES_MODEL_HH__

/* -------------------------------------------------------------------------- */
namespace akantu {

/* -------------------------------------------------------------------------- */
class ContactMechanicsInternodesModel : public Model {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ContactMechanicsInternodesModel(
      Mesh & mesh, UInt dim = _all_dimensions,
      const ID & id = "contact_mechanics_internodes_model",
      std::shared_ptr<DOFManager> dof_manager = nullptr,
      ModelType model_type = ModelType::_solid_mechanics_model);

  ~ContactMechanicsInternodesModel() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// initialize completely the model
  void initFullImpl(const ModelOptions & options) override;

  /// function to print the containt of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  /* ------------------------------------------------------------------------ */
  /* Solver interface                                                         */
  /* ------------------------------------------------------------------------ */
public:
  /// costum solve step for Internodes
  void solveStep(SolverCallback & callback, const ID & solver_id = "") override;
  void solveStep(const ID & solver_id = "") override;

  // assemble extended matrix K
  void assembleInternodesMatrix();

  /// assemble an interface matrix, eather master or slave
  Matrix<Real> assembleInterfaceMass(const NodeGroup & contact_node_group);

protected:
  // call back for the solver, computes the force residual
  void assembleResidual() override;

  // get the type of matrix needed
  MatrixType getMatrixType(const ID & matrix_id) const override;

  // callback for the solver, this assembles different matrices
  void assembleMatrix(const ID & matrix_id) override;

  // callback for the solver, this assembles the stiffness matrix
  void assembleLumpedMatrix(const ID & matrix_id) override;

  // get some default values for derived classes
  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  /// callback for the solver, this is called at beginning of solve
  void predictor() override;
  /// callback for the solver, this is called at end of solve
  void corrector() override;

  /// callback for the solver, this is called at beginning of solve
  void beforeSolveStep() override;

  /// callback for the solver, this is called at end of solve
  void afterSolveStep(bool converged = true) override;

  /// allocate all vectors
  void initSolver(TimeStepSolverType time_step_solver_type,
      NonLinearSolverType non_linear_solver_type) override;

protected:
  TimeStepSolverType getDefaultSolverType() const override;

  ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get contact detector
  AKANTU_GET_MACRO(ContactDetectorInternodes, *detector,
                   ContactDetectorInternodes &);

  /// get solid mechanics model
  AKANTU_GET_MACRO(SolidMechanicsModel, *solid, SolidMechanicsModel &);

  /// get lambdas from internodes formulation
  AKANTU_GET_MACRO(Lambdas, *lambdas, Array<Real>);

  /// get constraints 
  AKANTU_GET_MACRO(Contstraints, *constraints, Array<Real>);
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// contact detector
  std::unique_ptr<ContactDetectorInternodes> detector;

  /// solid mechanics model
  std::unique_ptr<SolidMechanicsModel> solid;

  /// lambdas array
  std::unique_ptr<Array<Real>> lambdas;

  /// contstraints for lambdas array
  std::unique_ptr<Array<Real>> constraints;

  /// blocked dofs for lambdas
  std::unique_ptr<Array<bool>> blocked_dofs;
};

} // namespace akantu

#endif /* __AKANTU_CONTACT_MECHANICS_INTERNODES_MODEL_HH__ */
