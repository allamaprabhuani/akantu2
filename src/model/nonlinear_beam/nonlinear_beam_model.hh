#include "data_accessor.hh"
#include "boundary_condition.hh"
#include "fe_engine.hh"
#include "model.hh"
#include <array>

#ifndef AKANTU_NONLINEAR_BEAM_MODEL_HH_
#define AKANTU_NONLINEAR_BEAM_MODEL_HH_

namespace akantu {
class Material;
class MaterialSelector;
class DumperIOHelper;
template <ElementKind kind, class IntegrationOrderFunctor>
class IntegratorGauss;
template <ElementKind kind> class ShapeLagrange;
} // namespace akantu

namespace akantu {

class NonlinearBeamModel : public Model,
                           public DataAccessor<Element>,
                           public DataAccessor<UInt> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using MyFEEngineType = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;

  NonlinearBeamModel(Mesh & mesh, UInt dim = _all_dimensions,
                     const ID & id = "nonlinear_beam_model",
                     std::shared_ptr<DOFManager> dof_manager = nullptr);

  ~NonlinearBeamModel() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// generic function to initialize everything ready for explicit dynamics
  void initFullImpl(const ModelOptions & options) override;

  /// read one material file to instantiate all the materials
  void readMaterials();

  /// allocate all vectors
  void initSolver(TimeStepSolverType, NonLinearSolverType) override;

  /// initialize the model
  void initModel() override;

  //void predictor() override;
  
  void assembleResidual() override;

  /// get the type of matrix needed
  MatrixType getMatrixType(const ID & matrix_id) const override;

  /// callback to assemble a lumped Matrix
  void assembleLumpedMatrix(const ID & matrix_id) override;
  
  /// callback to assemble a Matrix
  void assembleMatrix(const ID & matrix_id) override;

  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const override;

  TimeStepSolverType
  getDefaultSolverType() const override;

  Matrix<Real> N_matrix();

  Matrix<Real> N_rotator_matrix();

  Matrix<Real> N_grad_matrix();

  void assembleInternalForces();

  void assembleMass();

  void assembleStiffnessMatrix();

  void assembleMassLumped();
  /* ------------------------------------------------------------------------ */
  /* Methods for explicit                                                     */
  /* ------------------------------------------------------------------------ */
public:
  /// compute and get the stable time step
  //Real getStableTimeStep();

  /// set the stable timestep
  void setTimeStep(Real time_step, const ID & solver_id = "") override;

  ///RAJOUTER DES FUNCTION POUR LE CALCUL EXPLICIT

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
  /* ------------------------------------------------------------------------ */
public:
  /*
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
  */

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the current value of the time step
  AKANTU_GET_MACRO(TimeStep, time_step, Real);

  /// get the density value 
  AKANTU_GET_MACRO(Density, density, Real);

  /// get the value of the Young Modulus
  AKANTU_GET_MACRO(E, E, Real);

  /// get the value of the Poisson ratio
  AKANTU_GET_MACRO(Nu, nu, Real);

  /// get the value of the section area
  AKANTU_GET_MACRO(A, A, Real);

  /// get the Inertia
  AKANTU_GET_MACRO(J_11, J_11, Real);
  AKANTU_GET_MACRO(J_22, J_22, Real);
  AKANTU_GET_MACRO(J_33, J_33, Real);
  AKANTU_GET_MACRO(J_12, J_12, Real);
  AKANTU_GET_MACRO(J_13, J_13, Real);
  AKANTU_GET_MACRO(J_23, J_23, Real);
  
  /// get the NonlinearBeamModel::linear_angular_displacement array
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(Linear_Angular_Displacement, linear_angular_displacement);
  /// get the NonlinearBeamModel::linear_angular_displacement array
  AKANTU_GET_MACRO_DEREF_PTR(Linear_Angular_Displacement, linear_angular_displacement);

  /// get the NonlinearBeamModel::linear_angular_velocity array
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(Linear_Angular_Velocity, linear_angular_velocity);
  /// get the NonlinearBeamModel::linear_angular_velocity array
  AKANTU_GET_MACRO_DEREF_PTR(Linear_Angular_Velocity, linear_angular_velocity);

  /// get the NonlinearBeamModel::linear_angular_acceleration array
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(Linear_Angular_Acceleration, linear_angular_acceleration);
  /// get the NonlinearBeamModel::linear_angular_acceleration array
  AKANTU_GET_MACRO_DEREF_PTR(Linear_Angular_Acceleration, linear_angular_acceleration);

  /// get the NonlinearBeamModel::external_force_torque array
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(ExternalForceTorque, external_force_torque);
  /// get the NonlinearBeamModel::external_force_torque array
  AKANTU_GET_MACRO_DEREF_PTR(ExternalForceTorque, external_force_torque);

  /// get the NonlinearBeamModel::internal_force_torque array
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(InternalForceTorque, internal_force_torque);
  /// get the NonlinearBeamModel::internal_force_torque array
  AKANTU_GET_MACRO_DEREF_PTR(InternalForceTorque, internal_force_torque);

  /// get the NonlinearBeamModel::blocked_dofs array
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(BlockedDOFs, blocked_dofs);
  /// get the NonlinearBeamModel::blocked_dofs array
  AKANTU_GET_MACRO_DEREF_PTR(BlockedDOFs, blocked_dofs);

protected:
  /* ------------------------------------------------------------------------ */
  FEEngine & getFEEngineBoundary(const ID & name = "") override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  ElementType type;
  /// time step
  Real time_step;

  /// the density
  Real density;

  /// Young modulus
  Real E;

  /// poisson Ratio
  Real nu;

  /// Section Area
  Real A;

  /// Inertia
  Real J_11;
  Real J_22;
  Real J_33;
  Real J_12;
  Real J_13;
  Real J_23;

  /// linear and angular displacement array
  std::unique_ptr<Array<Real>> linear_angular_displacement;
  
  /// linear and angular velocity array
  std::unique_ptr<Array<Real>> linear_angular_velocity;
  
  /// linear and angular acceleration array
  std::unique_ptr<Array<Real>> linear_angular_acceleration;

  /// external force and torque array
  std::unique_ptr<Array<Real>> external_force_torque;
  
  /// internal force and torque array
  std::unique_ptr<Array<Real>> internal_force_torque;

  /// blocked dofs array
  std::unique_ptr<Array<bool>> blocked_dofs;



};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#endif /* AKANTU_NONLINEAR_BEAM_MODEL_HH_ */
  
  
  
