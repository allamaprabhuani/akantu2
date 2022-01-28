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

/* -------------------------------------------------------------------------- */
NonlinearBeamModel::NonlinearBeamModel(Mesh & mesh, UInt dim, const ID & id,
                                       std::shared_ptr<DOFManager> dof_manager)
    : Model(mesh, model_type, dof_manager, dim, id),Ns(0,2*spatial_dimension*2*spatial_dimension*Mesh::getNbNodesPerElement(_segment_3), "shape function matrix"), dNs(0,2*spatial_dimension*2*spatial_dimension*Mesh::getNbNodesPerElement(_segment_3), "derivative shape function matrix"), M(0,2*spatial_dimension*2*spatial_dimension*Mesh::getNbNodesPerElement(_segment_3)*Mesh::getNbNodesPerElement(_segment_3), "Mass matrix") {
  AKANTU_DEBUG_IN();

  this->registerDataAccessor(*this);
  auto nb_element = mesh.getNbElement(_segment_3);
  if (nb_element ==0){
    AKANTU_EXCEPTION("Model only works with _segment_3 elements");
  }
  this->type = _segment_3;
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

  registerFEEngineObject<FEEngineType>(id + ":fem", mesh, 1);

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
  this->registerParam("J_0", J_0, _pat_parsmod);
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
  Ns.resize(mesh.getNbElement(type) * fem.getNbIntegrationPoints(type));
  dNs.resize(mesh.getNbElement(type) * fem.getNbIntegrationPoints(type));
  M.resize(mesh.getNbElement(type) * fem.getNbIntegrationPoints(type));
  this->N_matrix(Ns);
  this->N_grad_matrix(dNs);
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
  
  // for alloc type of solvers
  this->allocNodalField(this->linear_angular_displacement, 2*spatial_dimension, "linear_angular_displacement");
  //
  this->allocNodalField(this->internal_force_torque, 2*spatial_dimension,
                        "internal_force_torque");
  this->allocNodalField(this->inertial_force_torque, 2*spatial_dimension,
                        "inertial_force_torque");
  //
  this->allocNodalField(this->external_force_torque, 2*spatial_dimension,
                        "external_force_torque");
  //
  this->allocNodalField(this->blocked_dofs, 2*spatial_dimension, "blocked_dofs");
  //
  this->allocNodalField(this->initial_angle, spatial_dimension, "initial_angle");
  

  /* ------------------------------------------------------------------------ */
  
  if (!dof_manager.hasDOFs("linear_angular_displacement")) {
    dof_manager.registerDOFs("linear_angular_displacement", *this->linear_angular_displacement, _dst_nodal);
    dof_manager.registerBlockedDOFs("linear_angular_displacement", *this->blocked_dofs);
    
  }
  

  /* ------------------------------------------------------------------------ */
  // for dynamic
  
  if (time_step_solver_type == TimeStepSolverType::_dynamic_lumped) {
    this->allocNodalField(this->linear_angular_velocity, 2*spatial_dimension, "linear_angular_velocity");
    //
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
  
  // computes the internal forces
  this->assembleInternalForces();
  this->computeInertialForces();

  this->getDOFManager().assembleToResidual("displacement",
                                           *this->external_force_torque, 1);
  this->getDOFManager().assembleToResidual("displacement",
                                           *this->internal_force_torque, -1);
  this->getDOFManager().assembleToResidual("displacement",
                                           *this->inertial_force_torque, -1);
  
  AKANTU_DEBUG_OUT();
}


MatrixType NonlinearBeamModel::getMatrixType(const ID & matrix_id) const {
  // \TODO check the materials to know what is the correct answer
  /*
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
  */
}


/* -------------------------------------------------------------------------- */
void NonlinearBeamModel::assembleMatrix(const ID & matrix_id) {
  
  if (matrix_id == "K") {
    //this->assembleStiffnessMatrix();
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

/* -------------------------------------------------------------------------- */
void NonlinearBeamModel::N_matrix(Array<Real> &Ns) {
  
  AKANTU_DEBUG_IN();
  Ns.set(0.);
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  const auto & shape_vector = getFEEngine().getShapes(type);
  
  for(auto && data : zip(make_view(shape_vector, nb_nodes_per_element), make_view(Ns, 2* spatial_dimension, 2*spatial_dimension*nb_nodes_per_element))) {
    auto & shapes = std::get<0>(data);
    auto & N = std::get<1>(data);
    for (UInt nd=0; nd < nb_nodes_per_element; ++nd) {
      N.block(Matrix<Real>::eye(6,1.) * shapes(nd),0, nd * 2*spatial_dimension);
    }
  }
   
  AKANTU_DEBUG_OUT();
  
}

void NonlinearBeamModel::N_grad_matrix(Array<Real> &dNs) {
  
  AKANTU_DEBUG_IN();
  dNs.set(0.);
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  const auto & derivative_shape_vector = getFEEngine().getShapesDerivatives(type);
  
  for(auto && data : zip(make_view(derivative_shape_vector, nb_nodes_per_element), make_view(dNs, 2* spatial_dimension, 2*spatial_dimension*nb_nodes_per_element))) {
    auto & deri_shapes = std::get<0>(data);
    auto & dN = std::get<1>(data);
    for (UInt nd=0; nd < nb_nodes_per_element; ++nd) {
      dN.block(Matrix<Real>::eye(6,1.) * deri_shapes(nd),0, nd * 2*spatial_dimension);
    }
  }
   
  AKANTU_DEBUG_OUT();
  
}

  
void NonlinearBeamModel::N_rotator_matrix(Array<Real> & N_rot_mat) {
  AKANTU_DEBUG_IN();
  
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  const auto & derivative_shape_vector = getFEEngine().getShapesDerivatives(type);
  const auto & shape_vector = getFEEngine().getShapes(type);
  UInt dim = spatial_dimension;
  //Array<Real> N_rot_mat(derivative_shape_vector.size(),2*spatial_dimension*2*spatial_dimension*nb_nodes_per_element);
  N_rot_mat.set(0.);
  Vector<Real> phip(3);

  auto conn_it = make_view(mesh.getConnectivity(type), nb_nodes_per_element).begin();
  auto disp_it = make_view(*(this->linear_angular_displacement), dim, 2).begin();
  auto init_pos_it = make_view(mesh.getNodes(), dim).begin();
  auto nb_quad_points = getFEEngine().getNbIntegrationPoints(type);
  
  for(auto && data : enumerate(make_view(derivative_shape_vector, nb_nodes_per_element), make_view(shape_vector, nb_nodes_per_element), make_view(N_rot_mat, 2* spatial_dimension, 2*spatial_dimension*nb_nodes_per_element))) {
    auto & deri_shapes = std::get<1>(data);
    auto & shapes = std::get<2>(data);
    auto & N_rot = std::get<3>(data);
    auto elem = std::get<0>(data)/nb_quad_points;

    Vector<UInt> conn = conn_it[elem];
    
    phip.set(0.);
    for (UInt nd = 0; nd<nb_nodes_per_element; ++nd) {
      Vector<Real> disp = Matrix<Real> (disp_it[conn(nd)])(0);
      Vector<Real> init = init_pos_it[conn(nd)];
      phip += (disp + init) * deri_shapes(nd);
    }
    Matrix<Real> phip_mat = skew(phip);
    for (UInt nd=0; nd < nb_nodes_per_element; ++nd) {
      N_rot.block(phip_mat * shapes(nd), 0, nd * 2*spatial_dimension + 3);
    }
  }

  AKANTU_DEBUG_OUT();

}


void NonlinearBeamModel::B_matrix(Array<Real> B_mat) {
  
  AKANTU_DEBUG_IN();
  UInt nb_node_per_element = Mesh::getNbNodesPerElement(type);
  const auto & derivative_shape_vector = getFEEngine().getShapesDerivatives(type);
  
  Array<Real> L(derivative_shape_vector.size(), spatial_dimension*nb_node_per_element);
  Array<Real> N_rot_mat(derivative_shape_vector.size(),2*spatial_dimension*2*spatial_dimension*nb_node_per_element);
  this->get_rotation_matrix(L);
  this->N_rotator_matrix(N_rot_mat);
  Matrix<Real> Rot(6,6);

  for (auto && data : zip(make_view(L, spatial_dimension, nb_node_per_element), make_view(N_rot_mat, 2* spatial_dimension, 2*spatial_dimension*nb_node_per_element), make_view(dNs, 2* spatial_dimension, 2*spatial_dimension*nb_node_per_element), make_view(B_mat, 2* spatial_dimension, 2*spatial_dimension*nb_node_per_element))) {
    auto & LL = std::get<0>(data);
    auto & N_rot = std::get<1>(data);
    auto & dN = std::get<2>(data);
    auto & B = std::get<3>(data);

    Rot.set(0.);
    Rot.block(LL, 0, 0);
    Rot.block(LL, 3, 3);

    B = Rot.transpose() * (dN + N_rot);
  }
  AKANTU_DEBUG_OUT();
    
}

void NonlinearBeamModel::get_rotation_matrix(Array<Real> L, bool origin) {
  
  AKANTU_DEBUG_IN();
  const auto & shape_vector = getFEEngine().getShapes(type);
  UInt nb_node_per_element = Mesh::getNbNodesPerElement(type);
  
  Array<Real> inter_angle(shape_vector.size(), nb_node_per_element);
  Array<Real> inter_initAngle(shape_vector.size(), nb_node_per_element);
  Array<Real> angle(this->linear_angular_displacement->size(), spatial_dimension);
  
  for (auto && data : zip(make_view(angle, spatial_dimension), make_view(*(this->linear_angular_displacement), spatial_dimension, 2))) {
    std::get<0>(data) = std::get<1>(data)(1);
  }
    
  this->interpolate(angle, inter_angle);
  this->interpolate(*(this->initial_angle), inter_initAngle);

  for (auto && data : zip(make_view(L, spatial_dimension,nb_node_per_element), make_view(inter_angle, spatial_dimension), make_view(inter_initAngle, spatial_dimension))) {
    auto & LL = std::get<0>(data);
;
    LL = expMap(std::get<2>(data));
    if (not origin) {
      LL = expMap(std::get<1>(data)) * LL;
    }
  }
  AKANTU_DEBUG_OUT();
  
}

void NonlinearBeamModel::interpolate(Array<Real> &field, Array<Real> &interField) {
  AKANTU_DEBUG_IN();
  const auto & shape_vector = getFEEngine().getShapes(type);
  UInt nb_node_per_element = Mesh::getNbNodesPerElement(type);

  auto conn_it = make_view(mesh.getConnectivity(type), nb_node_per_element).begin();
  auto field_it = make_view(field, spatial_dimension, 1).begin();
  auto nb_quad_points = getFEEngine().getNbIntegrationPoints(type);
  
  for (auto && data : enumerate(make_view(shape_vector, nb_node_per_element), make_view(interField, nb_node_per_element))) {
    auto & N = std::get<1>(data);
    auto & intF = std::get<2>(data);
    auto elem = std::get<0>(data)/nb_quad_points;

    Vector<UInt> conn = conn_it[elem];
    for (UInt nd = 0; nd<nb_node_per_element; ++nd) {
      Vector<Real> F = Matrix<Real> (field_it[conn(nd)])(0);
      intF += F * N(nd);
    }   
  }
  AKANTU_DEBUG_OUT();  
}

void NonlinearBeamModel::grad_interpolate(Array<Real> &field, Array<Real> &gradinterField) {
  AKANTU_DEBUG_IN();
  const auto & derivative_shape_vector = getFEEngine().getShapesDerivatives(type);
  UInt nb_node_per_element = Mesh::getNbNodesPerElement(type);

  auto conn_it = make_view(mesh.getConnectivity(type), nb_node_per_element).begin();
  auto field_it = make_view(field, spatial_dimension, 1).begin();
  auto nb_quad_points = getFEEngine().getNbIntegrationPoints(type);
  
  for (auto && data : enumerate(make_view(derivative_shape_vector, nb_node_per_element), make_view(gradinterField, nb_node_per_element))) {
    auto & dN = std::get<1>(data);
    auto & gintF = std::get<2>(data);
    auto elem = std::get<0>(data)/nb_quad_points;

    Vector<UInt> conn = conn_it[elem];
    for (UInt nd = 0; nd<nb_node_per_element; ++nd) {
      Vector<Real> F = Matrix<Real> (field_it[conn(nd)])(0);
      gintF += F * dN(nd);
    }   
  }
  AKANTU_DEBUG_OUT();  
}


void NonlinearBeamModel::computeStrains(Array<Real> strains, bool origin) {

  AKANTU_DEBUG_IN();
  UInt nb_node_per_element = Mesh::getNbNodesPerElement(type);
  const auto & shape_vector = getFEEngine().getShapes(type);

  Array<Real> pos(this->linear_angular_displacement->size(), spatial_dimension);

  for (auto && data : zip(make_view(pos, spatial_dimension), make_view(mesh.getNodes(), spatial_dimension), make_view(*(this->linear_angular_displacement), spatial_dimension, 2))) {
    auto P = std::get<0>(data);
    auto init_P = std::get<1>(data);
    auto lin_disp = std::get<2>(data)(0);

    P = init_P;
    if (not origin) {
      P = init_P + Vector<Real> (lin_disp);
    }
  }

  Array<Real> L(shape_vector.size(), spatial_dimension*nb_node_per_element);
  Array<Real> L0(shape_vector.size(), spatial_dimension*nb_node_per_element);
  this->get_rotation_matrix(L, origin);
  this->get_rotation_matrix(L0, true);

  Array<Real> theta(this->linear_angular_displacement->size(), spatial_dimension);
  for (auto && data : zip(make_view(theta, spatial_dimension), make_view(*(this->linear_angular_displacement), spatial_dimension, 2))) {
    std::get<0>(data) = std::get<1>(data)(1);
  }

  Array<Real> Gamma(shape_vector.size(), spatial_dimension);
  Array<Real> B(shape_vector.size(), spatial_dimension);
  Array<Real> Bp(shape_vector.size(), spatial_dimension);
  this->grad_interpolate(pos, Gamma);
  this->interpolate(theta, B);
  this->grad_interpolate(theta, Bp);

  Matrix<Real> Omega_hat(3,3);
  Vector<Real> e3 = {0., 0., 1.};

  for (auto && data : zip(make_view(L,spatial_dimension,nb_node_per_element), make_view(L0, spatial_dimension,nb_node_per_element), make_view(Gamma, nb_node_per_element), make_view(B, nb_node_per_element), make_view(Bp, nb_node_per_element), make_view(strains, spatial_dimension, 2))) {
    auto & LL = std::get<0>(data);
    auto & LL0 = std::get<1>(data);
    auto & Gam = std::get<2>(data);
    auto & BB = std::get<3>(data);
    auto & BBp = std::get<4>(data);
    auto & strain = std::get<5>(data);
    
    Omega_hat = LL.transpose() * expDerivative(BB, BBp) * LL0;
    strain(1) = skew2vec(Omega_hat);
    strain(0) = LL.transpose() * Gam - e3;
  }
  AKANTU_DEBUG_OUT();
}

  

void NonlinearBeamModel::constitutive_C(Array<Real> C) {
  AKANTU_DEBUG_IN();
  Real E = this->E;
  Real G = this->E/(2*(1+this->nu));
  Real A = this->A;
  Real J_22 = this->J_22;
  Real J_12 = this->J_12;
  Real J_11 = this->J_11;
  Real J_0 = this->J_0;
  for (auto && data : zip(make_view(C, 2*spatial_dimension, 2*spatial_dimension))) {
    auto & CC = std::get<0>(data);
    CC(0, 0) = G * A;
    CC(1, 1) = G * A;
    CC(2, 2) = E * A;
    CC(3, 3) = J_22 * E;
    CC(3, 4) = -J_12 * E;
    CC(4, 3) = -J_12 * E;
    CC(4, 4) = J_11 * E;
    CC(5, 5) = J_0 * G;
  }
  AKANTU_DEBUG_OUT();

}

void NonlinearBeamModel::computeJbar(Array<Real> Jbar) {
  AKANTU_DEBUG_IN();
  UInt nb_node_per_element = Mesh::getNbNodesPerElement(type);
  const auto & shape_vector = getFEEngine().getShapes(type);
  
  Array<Real> strains(shape_vector.size(), 2*spatial_dimension);
  Array<Real> strains0(shape_vector.size(), 2*spatial_dimension);
  this->computeStrains(strains, false);
  this->computeStrains(strains0, true);
  
  Array<Real> L(shape_vector.size(), spatial_dimension*nb_node_per_element);
  this->get_rotation_matrix(L);

  
  Array<Real> ID(shape_vector.size(), spatial_dimension*nb_node_per_element);
  for (auto && data : zip(make_view(ID, spatial_dimension, nb_node_per_element))){
    std::get<0>(data) = Matrix<Real>::eye(spatial_dimension,1);
  }
  
  Array<Real> DX(shape_vector.size(), spatial_dimension*spatial_dimension);
  Array<Real> DX0(shape_vector.size(), spatial_dimension*spatial_dimension);
  this->computeDX(DX, strains, L, 0., 0.);
  this->computeDX(DX0, strains0, ID, 0., 0.);

  Array<Real> F(shape_vector.size(), spatial_dimension*spatial_dimension);
  for (auto && data : zip(make_view(F, spatial_dimension, spatial_dimension), make_view(DX, spatial_dimension, spatial_dimension), make_view(DX0, spatial_dimension, spatial_dimension), make_view(Jbar))) {
    auto & FF = std::get<0>(data);
    auto & D = std::get<1>(data);
    auto & D0 = std::get<2>(data);
    FF = D * D0.inverse();
    std::get<3>(data) = FF.det();
  }
  AKANTU_DEBUG_OUT();
  
}

void NonlinearBeamModel::computeDX(Array<Real> DX, Array<Real> strains, Array<Real> Lambda, Real x1, Real x2) {
  AKANTU_DEBUG_IN();
  Vector<Real> e1 = {1., 0., 0.};
  Vector<Real> e2 = {0., 1., 0.};
  Vector<Real> e3 = {0., 0., 1.};
  Vector<Real> Gamma(3);
  Vector<Real> Omega(3);
  Vector<Real> t1(3);
  Vector<Real> t2(3);
  Vector<Real> vs(3);
  Vector<Real> vec(3);
  Vector<Real> cp(3);
  Matrix<Real> op(3,3);
  Matrix<Real> res(3,3);

  for (auto && data : zip(make_view(strains, spatial_dimension, 2), make_view(Lambda, spatial_dimension, spatial_dimension), make_view(DX, spatial_dimension, spatial_dimension))) {
    auto & strain = std::get<0>(data);

    Gamma = strain(0);
    Omega = strain(1);

    t1 = std::get<1>(data) * e1;
    t2 = std::get<1>(data) * e2;
    vs = t1 * x1 + t2 * x2;
    vec = Gamma + cp.crossProduct(Omega, std::get<1>(data) * vs);
    op.outerProduct(vec, e3);
    res = Matrix<Real>::eye(3, 1.) + op;
    std::get<2>(data) = std::get<1>(data) * res;
  }   

  AKANTU_DEBUG_OUT();
}

void NonlinearBeamModel::computeStresses(Array<Real> stresses) {

  AKANTU_DEBUG_IN();
  UInt nb_node_per_element = Mesh::getNbNodesPerElement(type);
  const auto & shape_vector = getFEEngine().getShapes(type);
  
  Array<Real> strains(shape_vector.size(), 2*spatial_dimension);
  this->computeStrains(strains, false);
  Array<Real> C(shape_vector.size(), 2*spatial_dimension*2*spatial_dimension);
  this->constitutive_C(C);

  for (auto && data :zip(make_view(stresses, 2*spatial_dimension), make_view(strains, 2*spatial_dimension), make_view(C, 2*spatial_dimension, 2*nb_node_per_element))){
    auto & stress = std::get<0>(data);
    auto & strain = std::get<1>(data);
    auto & CC = std::get<2>(data);

    stress = CC * strain;
  }

  AKANTU_DEBUG_OUT();
}


void NonlinearBeamModel::computeInertiaTensor(Array<Real> J) {
  
  AKANTU_DEBUG_IN();

  Real J_22 = this->J_22;
  Real J_12 = this->J_12;
  Real J_21 = J_12;
  Real J_11 = this->J_11;
  Real rho = this->density;
  Matrix<Real> J1(spatial_dimension, spatial_dimension);
  J1 = {{J_11, J_12, 0}, {J_21, J_22, 0}, {0, 0, 0}};
  for (auto && data : zip(make_view(J, spatial_dimension, spatial_dimension))) {
    std::get<0>(data) = rho * Matrix<Real>::eye(spatial_dimension, 1.) * (J_11+J_22) - rho * J1;
  }

  AKANTU_DEBUG_OUT();
  
}



void NonlinearBeamModel::assembleMass() {
  AKANTU_DEBUG_IN();
  UInt nb_node_per_element = Mesh::getNbNodesPerElement(type);
  const auto & shape_vector = getFEEngine().getShapes(type);
  
  Array<Real> J(shape_vector.size(), spatial_dimension*spatial_dimension);
  this->computeInertiaTensor(J);
  Array<Real> L(shape_vector.size(), spatial_dimension*nb_node_per_element);
  this->get_rotation_matrix(L);

  Matrix<Real> R(2*spatial_dimension, 2*spatial_dimension);
  R.set(0.);
  Matrix<Real> Imat(2*spatial_dimension, 2*spatial_dimension);
  Imat.set(0.);
  Imat.block(Matrix<Real>::eye(spatial_dimension, 1.) * this->density * this->A, 0, 0);
  
  for (auto && data : zip(make_view(M, 2*spatial_dimension*nb_node_per_element, 2*spatial_dimension*nb_node_per_element), make_view(J, spatial_dimension, spatial_dimension), make_view(L, spatial_dimension, spatial_dimension), make_view(Ns, 2* spatial_dimension, 2*spatial_dimension*nb_node_per_element))) {
    auto & Mass = std::get<0>(data);
    auto & LL = std::get<2>(data);
    auto & N = std::get<3>(data);
      
    R.block(LL, 0, 0);
    R.block(LL, 3, 3);
    Imat.block(std::get<1>(data), 3, 3);
      
    Mass = N.transpose() * R * Imat * R.transpose() * N;
    //Mass /= L/2;
  }
  AKANTU_DEBUG_OUT();

}

void NonlinearBeamModel::convected_angular_velocity(Array<Real> conv_rv) {
  AKANTU_DEBUG_IN();
  UInt nb_node_per_element = Mesh::getNbNodesPerElement(type);
  const auto & shape_vector = getFEEngine().getShapes(type);

  Array<Real> L(shape_vector.size(), spatial_dimension*nb_node_per_element);
  Array<Real> L0(shape_vector.size(), spatial_dimension*nb_node_per_element);
  this->get_rotation_matrix(L, false);
  this->get_rotation_matrix(L0, true);

  Array<Real> inter_rv(shape_vector.size(), nb_node_per_element);
  Array<Real> rv(this->linear_angular_displacement->size(), spatial_dimension);
  Array<Real> inter_theta(shape_vector.size(), nb_node_per_element);
  Array<Real> theta(this->linear_angular_displacement->size(), spatial_dimension);
  
  for (auto && data : zip(make_view(rv, spatial_dimension), make_view(*(this->linear_angular_velocity), spatial_dimension, 2), make_view(theta, spatial_dimension), make_view(*(this->linear_angular_displacement), spatial_dimension, 2))) {
    std::get<0>(data) = std::get<1>(data)(1);
    std::get<2>(data) = std::get<3>(data)(1);
  }    
  this->interpolate(rv, inter_rv);
  this->interpolate(theta, inter_theta);

  Matrix<Real> W_hat(spatial_dimension, spatial_dimension);
  for (auto && data : zip(make_view(L, spatial_dimension, nb_node_per_element), make_view(L0, spatial_dimension, nb_node_per_element), make_view(inter_theta, spatial_dimension), make_view(inter_rv, spatial_dimension), make_view(conv_rv, spatial_dimension))) {
    auto & LL = std::get<0>(data);
    auto & LL0 = std::get<1>(data);

    W_hat = LL.transpose() * expDerivative(std::get<2>(data), std::get<3>(data)) * LL0;
    std::get<4>(data) = skew2vec(W_hat);
  }
  AKANTU_DEBUG_OUT();
}
  


void NonlinearBeamModel::computeInertialBodyForce(Array<Real> IBF) {
  
  AKANTU_DEBUG_IN();
  UInt nb_node_per_element = Mesh::getNbNodesPerElement(type);
  const auto & shape_vector = getFEEngine().getShapes(type);

  Array<Real> J(shape_vector.size(), spatial_dimension*spatial_dimension);
  this->computeInertiaTensor(J);

  Array<Real> w(shape_vector.size(), spatial_dimension);
  this->convected_angular_velocity(w);

  Array<Real> L(shape_vector.size(), spatial_dimension*nb_node_per_element);
  this->get_rotation_matrix(L);

  Vector<Real> f(spatial_dimension);

  for (auto && data : zip(make_view(IBF, spatial_dimension, 2), make_view(J, spatial_dimension, spatial_dimension), make_view(w, spatial_dimension), make_view(L, spatial_dimension, nb_node_per_element))) {
    auto & LL = std::get<3>(data);

    f = skew(std::get<2>(data)) * LL * std::get<1>(data) * LL.transpose() * std::get<2>(data);
    std::get<0>(data)(0) = Vector<Real> {0., 0., 0.};
    std::get<0>(data)(1) = f;
  }

  AKANTU_DEBUG_OUT();
}
    
  
void NonlinearBeamModel::computeInertialForces() {
  AKANTU_DEBUG_IN();
  UInt nb_node_per_element = Mesh::getNbNodesPerElement(type);
  const auto & shape_vector = getFEEngine().getShapes(type);

  Array<Real> IBF(shape_vector.size(), 2*spatial_dimension);
  this->computeInertialBodyForce(IBF);

  for (auto && data : zip(make_view(*(this->inertial_force_torque), 2*spatial_dimension*nb_node_per_element), make_view(Ns, 2* spatial_dimension, 2*spatial_dimension*nb_node_per_element), make_view(IBF, 2*spatial_dimension))) {
    auto & N = std::get<1>(data);

    std::get<0>(data) = N.transpose() * std::get<2>(data);
    //std::get<0>(data) /= L/2.;
  }
  AKANTU_DEBUG_OUT();
}


void NonlinearBeamModel::assembleInternalForces() {
  AKANTU_DEBUG_IN();
  UInt nb_node_per_element = Mesh::getNbNodesPerElement(type);
  const auto & shape_vector = getFEEngine().getShapes(type);

  Array<Real> B(shape_vector.size(), 2*spatial_dimension*2*spatial_dimension*nb_node_per_element);
  this->B_matrix(B);
  Array<Real> stresses(shape_vector.size(), 2*spatial_dimension);
  this->computeStresses(stresses);

  for (auto && data : zip(make_view(*(this->internal_force_torque), 2*spatial_dimension*nb_node_per_element), make_view(B, 2*spatial_dimension, 2*spatial_dimension*nb_node_per_element), make_view(stresses, 2*spatial_dimension))) {
    std::get<0>(data) = std::get<1>(data).transpose() * std::get<2>(data);
  }
  AKANTU_DEBUG_OUT();
}




void NonlinearBeamModel::computeAcceleration() {
  AKANTU_DEBUG_IN();
  UInt nb_node_per_element = Mesh::getNbNodesPerElement(type);
  const auto & shape_vector = getFEEngine().getShapes(type);
  
  this->assembleMass();

  for (auto && data : zip(make_view(M, 2*spatial_dimension*nb_node_per_element, 2*spatial_dimension*nb_node_per_element), make_view(*(this->blocked_dofs), 2*spatial_dimension*nb_node_per_element))){

  }

  AKANTU_DEBUG_OUT();
}

  



































} // namespace akantu

