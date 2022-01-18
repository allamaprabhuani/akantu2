#include "nonlinear_beam_model.hh"
#include "dumpable_inline_impl.hh"
#include "element_synchronizer.hh"
#include "fe_engine_template.hh"
#include "generalized_trapezoidal.hh"
#include "group_manager_inline_impl.hh"
#include "integrator_gauss.hh"
#include "mesh.hh"
#include "parser.hh"
#include "shape_lagrange.hh"

#ifdef AKANTU_USE_IOHELPER
#include "dumper_element_partition.hh"
#include "dumper_elemental_field.hh"
#include "dumper_internal_material_field.hh"
#include "dumper_iohelper_paraview.hh"
#endif

/* -------------------------------------------------------------------------- */
namespace akantu {

namespace nonlinear_beam {
  namespace details {
    class ComputeRhoFunctor {
    public:
      ComputeRhoFunctor(const NonlinearBeamModel & model) : model(model){};

      void operator()(Matrix<Real> & rho, const Element & /*unused*/) {
        rho.set(model.getCapacity() * model.getDensity());
      }

    private:
      const NonlinearBeamModel & model;
    };
  } // namespace details
} // namespace nonlinear_beam

/* -------------------------------------------------------------------------- */
NonlinearBeamModel::NonlinearBeamModel(Mesh & mesh, UInt dim, const ID & id,
                                       std::shared_ptr<DOFManager> dof_manager)
    : Model(mesh, ModelType::_nonlinear_neam_model, dof_manager, dim, id){
  AKANTU_DEBUG_IN();

  this->registerDataAccessor(*this);

  /*
  if (this->mesh.isDistributed()) {
    auto & synchronizer = this->mesh.getElementSynchronizer();
    /// Syncronisation
    this->registerSynchronizer(synchronizer,
                               SynchronizationTag::_htm_temperature);
    this->registerSynchronizer(synchronizer,
                               SynchronizationTag::_htm_gradient_temperature);
  }
  */

  registerFEEngineObject<FEEngineType>(id + ":fem", mesh, spatial_dimension);

#ifdef AKANTU_USE_IOHELPER
  this->mesh.registerDumper<DumperParaview>("nonlinear_beam", id, true);
  this->mesh.addDumpMesh(mesh, spatial_dimension, _not_ghost, _ek_regular);
#endif
  this->registerParam("E", E, _pat_parsmod);
  this->registerParam("nu", nu, _pat_parsmod);
  this->registerParam("A", A, _pat_parsmod);
  this->registerParam("J_11", J_11, _pat_parsmod);
  this->registerParam("J_22", J_22, _pat_parsmod);
  this->registerParam("J_33", J_33, _pat_parsmod);
  this->registerParam("J_12", J_12, _pat_parsmod);
  this->registerParam("J_13", J_13, _pat_parsmod);
  this->registerParam("J_23", J_23, _pat_parsmod);
  this->registerParam("density", density, _pat_parsmod);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
NonlinearBeamModel::~NonlinearBeamModel() = default;

/* -------------------------------------------------------------------------- */
void NonlinearBeamModel::setTimeStep(Real time_step, const ID & solver_id) {
  Model::setTimeStep(time_step, solver_id);

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.getDumper().setTimeStep(time_step);
#endif
}
/* -------------------------------------------------------------------------- */
/* Initialization                                                             */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void NonlinearBeamModel::initFullImpl(const ModelOptions & options) {
  Model::initFullImpl(options);

  readMaterials();
}

/* -------------------------------------------------------------------------- */
void NonlinearBeamModel::initModel() {
  auto & fem = this->getFEEngine();
  fem.initShapeFunctions(_not_ghost);
  fem.initShapeFunctions(_ghost);
  
}

/* -------------------------------------------------------------------------- */
TimeStepSolverType NonlinearBeamModel::getDefaultSolverType() const {
  return TimeStepSolverType::_dynamic_lumped;
}

/* -------------------------------------------------------------------------- */
ModelSolverOptions NonlinearBeamModel::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  case TimeStepSolverType::_dynamic_lumped: {
    options.non_linear_solver_type = NonLinearSolverType::_lumped;
    options.integration_scheme_type["linear_angular_displacement"] =
        IntegrationSchemeType::_central_difference;
    options.solution_type["linear_angular_displacement"] = IntegrationScheme::_acceleration;
    break;
  }
  case TimeStepSolverType::_static: {
    options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
    options.integration_scheme_type["linear_angular_displacement"] =
        IntegrationSchemeType::_pseudo_time;
    options.solution_type["linear_angular_displacement"] = IntegrationScheme::_not_defined;
    break;
  }
  default:
    AKANTU_EXCEPTION(type << " is not a valid time step solver type");
  }

  return options;
}

/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
NonlinearBeamModel::getDefaultSolverID(const AnalysisMethod & method) {
  switch (method) {
  case _explicit_lumped_mass: {
    return std::make_tuple("explicit_lumped",
                           TimeStepSolverType::_dynamic_lumped);
  }
  case _static: {
    return std::make_tuple("static", TimeStepSolverType::_static);
  }
  default:
    return std::make_tuple("unknown", TimeStepSolverType::_not_defined);
  }
}

/* -------------------------------------------------------------------------- */
void NonlinearBeamModel::initSolver(TimeStepSolverType time_step_solver_type,
                                     NonLinearSolverType /*unused*/) {
  auto & dof_manager = this->getDOFManager();
  /* ------------------------------------------------------------------------ */
  // for alloc type of solvers
  this->allocNodalField(this->displacement, spatial_dimension, "displacement");
  this->allocNodalField(this->angle, spatial_dimension, "angle");
  this->allocNodalField(this->linear_angular_displacement, 2*spatial_dimension, "linear_angular_displacement");
  //
  this->allocNodalField(this->internal_force, spatial_dimension,
                        "internal_force");
  this->allocNodalField(this->internal_torque, spatial_dimension,
                        "internal_torque");
  this->allocNodalField(this->internal_force_torque, 2*spatial_dimension,
                        "internal_force_torque");
  //
  this->allocNodalField(this->external_force, spatial_dimension,
                        "external_force");
  this->allocNodalField(this->external_torque, spatial_dimension,
                        "external_torque");
  this->allocNodalField(this->external_force_torque, 2*spatial_dimension,
                        "external_force_torque");
  //
  this->allocNodalField(this->blocked_dofs, 2*spatial_dimension, "blocked_dofs");

  /* ------------------------------------------------------------------------ */
  if (!dof_manager.hasDOFs("linear_angular_displacement")) {
    dof_manager.registerDOFs("linear_angular_displacement", *this->linear_angular_displacement, _dst_nodal);
    dof_manager.registerBlockedDOFs("linear_angular_displacement", *this->blocked_dofs);
    
  }

  /* ------------------------------------------------------------------------ */
  // for dynamic
  if (time_step_solver_type == TimeStepSolverType::_dynamic_lumped) {
    this->allocNodalField(this->velocity, spatial_dimension, "velocity");
    this->allocNodalField(this->angular_velocity, spatial_dimension, "angular_velocity");
    this->allocNodalField(this->linear_angular_velocity, 2*spatial_dimension, "linear_angular_velocity");
    //
    this->allocNodalField(this->acceleration, spatial_dimension,
                          "acceleration");
    this->allocNodalField(this->angular_acceleration, spatial_dimension, "angular_acceleration");
    this->allocNodalField(this->linear_angular_acceleration, 2*spatial_dimension, "linear_angular_acceleration");

    if (!dof_manager.hasDOFsDerivatives("linear_angular_displacement", 1)) {
      dof_manager.registerDOFsDerivative("linear_angular_displacement", 1, *this->linear_angular_velocity);
      dof_manager.registerDOFsDerivative("linear_angular_displacement", 2, *this->linear_angular_acceleration);
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonlinearBeamModel::assembleResidual() {
  AKANTU_DEBUG_IN();

  /* ------------------------------------------------------------------------ */
  // computes the internal forces
  this->assembleInternalForces();

  /* ------------------------------------------------------------------------ */
  this->getDOFManager().assembleToResidual("displacement",
                                           *this->external_force_torque, 1);
  this->getDOFManager().assembleToResidual("displacement",
                                           *this->internal_force_torque, 1);

  AKANTU_DEBUG_OUT();
}


MatrixType NonlinearBeamModel::getMatrixType(const ID & matrix_id) const {
  // \TODO check the materials to know what is the correct answer
  if (matrix_id == "C") {
    return _mt_not_defined;
  }

  if (matrix_id == "K") {
    auto matrix_type = _unsymmetric;

    for (auto & material : materials) {
      matrix_type = std::max(matrix_type, material->getMatrixType(matrix_id));
    }
  }
  return _symmetric;
}

/* -------------------------------------------------------------------------- */
void NonlinearBeamModel::assembleMatrix(const ID & matrix_id) {
  if (matrix_id == "K") {
    this->assembleStiffnessMatrix();
  } else if (matrix_id == "M") {
    this->assembleMass();
  }
}

/* -------------------------------------------------------------------------- */
void NonlinearBeamModel::assembleLumpedMatrix(const ID & matrix_id) {
  if (matrix_id == "M") {
    this->assembleMassLumped();
  }
}

void NonlinearBeamModel::performStep(Real time_step) {
  AKANTU_DEBUG_IN();
  
  const Array<bool> & blocked_dofs =
      this->dof_manager.getBlockedDOFs(this->dof_id);
  
  this->prediction_Velocity(time_step, blocked_dofs);

  AKANTU_DEBUG_OUT();
}

void NonlinearBeamModel::prediction_velocity(Real time_step,const Array<bool> & blocked_dofs) {
  this->linear_angular_velocity += time_step / 2 * this->linear_angular_acceleration;
}
