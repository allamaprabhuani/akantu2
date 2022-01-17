#include "data_accessor.hh"
#include "fe_engine.hh"
#include "model.hh"
#include <array>

#ifndef AKANTU_NONLINEAR_BEAM_MODEL_HH_
#define AKANTU_NONLINEAR_BEAM_MODEL_HH_

namespace akantu {
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
  using FEEngine = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;

  NonlinearBeamModel(Mesh & mesh, UInt dim = _all_dimensions,
                     const ID & id = "nonlinear_beam_model",
                     std::shred_ptr<DOFManager> dof_manager = nullptr);

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

  void predictor() override;
  
  void assembleResidual() override;

  /// get the type of matrix needed
  MatrixType getMatrixType(const ID & matrix_id) override;

  /// callback to assemble a lumped Matrix
  void assembleLumpedMatrix(const ID & matrix_id) override;
  
  /// callback to assemble a Matrix
  void assembleMatrix(const ID & matrix_id) override;

  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const override;
  /* ------------------------------------------------------------------------ */
  /* Methods for explicit                                                     */
  /* ------------------------------------------------------------------------ */
public:
  /// compute and get the stable time step
  Real getStableTimeStep();

  /// set the stable timestep
  void setTimeStep(Real time_step, const ID & solver_id = "") override;

  ///RAJOUTER DES FUNCTION POUR LE CALCUL EXPLICIT

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
  AKANTU_GET_MACRO(Density, density, Real);
  AKANTU_GET_MACRO(Capacity, capacity, Real);
  /// get the current value of the time step
  AKANTU_GET_MACRO(TimeStep, time_step, Real);
  /// RAJOUTER DES ACCESSORS

protected:
  /* ------------------------------------------------------------------------ */
  FEEngine & getFEEngineBoundary(const ID & name = "") override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// time step
  Real time_step;

  /// the density
  Real density;

  /// displacement array
  std::unique_ptr<Array<Real>> displacement;
  /// velocity array
  std::unique_ptr<Array<Real>> velocity;
  /// acceleration array
  std::unique_ptr<Array<Real>> acceleration;
  

  /// RAJOUTER


};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#endif /* AKANTU_NONLINEAR_BEAM_MODEL_HH_ */
  
  
  
