/**
 * @file   heat_transfer_interface_model.cc
 *
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 *
 * @date creation: Thu Apr 13 2023
 * @date last modification: Wed May 22 2024
 *
 * @brief  Implementation of HeatTransferInterfaceModel class
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
#include "heat_transfer_interface_model.hh"
#include "dumpable_inline_impl.hh"
#include "element_class.hh"
#include "element_synchronizer.hh"
#include "fe_engine_template.hh"
#include "generalized_trapezoidal.hh"
#include "group_manager_inline_impl.hh"
#include "integrator_gauss.hh"
#include "mesh.hh"
#include "parser.hh"
#include "shape_cohesive.hh"
#include "shape_lagrange.hh"
/* -------------------------------------------------------------------------- */
#include "dumper_element_partition.hh"
#include "dumper_elemental_field.hh"
#include "dumper_internal_material_field.hh"
#include "dumper_iohelper_paraview.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

namespace heat_transfer_interface {
  namespace details {
    class ComputeRhoFunctor {
    public:
      ComputeRhoFunctor(const HeatTransferInterfaceModel & model)
          : model(model){};

      void operator()(Matrix<Real> & rho, const Element & el) {
        // auto opening_it = model.getOpening(el.type, el.ghost_type).begin();

        // rho.set(model.getCapacityInCrack() * model.getDensityInCrack() *
        //         opening_it[el.element]);
        rho.set(model.getCapacityInCrack() * model.getDensityInCrack() *
                model.getDefaultOpening());
      }

    private:
      const HeatTransferInterfaceModel & model;
    };
  } // namespace details
} // namespace heat_transfer_interface

/* -------------------------------------------------------------------------- */
HeatTransferInterfaceModel::HeatTransferInterfaceModel(
    Mesh & mesh, UInt dim, const ID & id,
    std::shared_ptr<DOFManager> dof_manager)
    : HeatTransferModel(mesh, dim, id,
                        ModelType::_heat_transfer_interface_model, dof_manager),
      temperature_on_qpoints_coh("temperature_on_qpoints_coh", id),
      opening_on_qpoints("opening_on_qpoints", id),
      opening_rate("opening_rate", id), k_long_w("k_long_w", id) {
  AKANTU_DEBUG_IN();

  registerFEEngineObject<MyFEEngineCohesiveType>("InterfacesFEEngine", mesh,
                                                 Model::spatial_dimension);

  this->mesh.registerDumper<DumperParaview>("heat_interfaces", id);
  this->mesh.addDumpMeshToDumper("heat_interfaces", mesh, spatial_dimension,
                                 _not_ghost, _ek_cohesive);

  if (this->mesh.isDistributed()) {
    /// create the distributed synchronizer for cohesive elements
    this->cohesive_synchronizer = std::make_unique<ElementSynchronizer>(
        mesh, "cohesive_distributed_synchronizer");

    auto & synchronizer = mesh.getElementSynchronizer();
    this->cohesive_synchronizer->split(synchronizer, [](auto && el) {
      return Mesh::getKind(el.type) == _ek_cohesive;
    });

    this->registerSynchronizer(*cohesive_synchronizer,
                               SynchronizationTag::_htm_gradient_temperature);
    this->registerSynchronizer(*cohesive_synchronizer,
                               SynchronizationTag::_htm_temperature_on_qpoints);

    auto & mesh_facets = this->mesh.getMeshFacets();
    auto & facet_synchronizer = mesh_facets.getElementSynchronizer();
    const auto & cfacet_synchronizer = facet_synchronizer;
    // update the cohesive element synchronizer
    cohesive_synchronizer->updateSchemes([&](auto && scheme, auto && proc,
                                             auto && direction) {
      auto & facet_scheme =
          cfacet_synchronizer.getCommunications().getScheme(proc, direction);

      for (auto && facet : facet_scheme) {
        const auto & cohesive_element = const_cast<const Mesh &>(mesh_facets)
                                            .getElementToSubelement(facet)[1];

        if (cohesive_element == ElementNull or
            cohesive_element.kind() != _ek_cohesive) {
          continue;
        }
        scheme.push_back(cohesive_element);
      }
    });
  }

  this->registerParam("longitudinal_conductivity", longitudinal_conductivity,
                      _pat_parsmod);
  this->registerParam("transversal_conductivity", transversal_conductivity,
                      _pat_parsmod);
  this->registerParam("default_opening", default_opening, 1.e-5, _pat_parsmod);
  this->registerParam("capacity_in_crack", capacity_in_crack, _pat_parsmod);
  this->registerParam("density_in_crack", density_in_crack, _pat_parsmod);
  this->registerParam("use_opening_rate", use_opening_rate, _pat_parsmod);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferInterfaceModel::initModel() {

  Parent::initModel();

  auto & fem = this->getFEEngine("InterfacesFEEngine");
  fem.initShapeFunctions(_not_ghost);
  fem.initShapeFunctions(_ghost);

  this->temperature_gradient.initialize(fem, _nb_component = spatial_dimension);
  this->k_long_w.initialize(fem, _nb_component = (spatial_dimension - 1) *
                                                 (spatial_dimension - 1));
  this->temperature_on_qpoints_coh.initialize(fem, _nb_component = 1);
  this->opening_on_qpoints.initialize(fem, _nb_component = 1);
  if (use_opening_rate) {
    opening_rate.initialize(fem, _nb_component = 1);
  }
}

/* -------------------------------------------------------------------------- */
HeatTransferInterfaceModel::~HeatTransferInterfaceModel() = default;

/* -------------------------------------------------------------------------- */
void HeatTransferInterfaceModel::readMaterials() {

  auto sect = this->getParserSection();

  if (not std::get<1>(sect)) {
    const auto & section = std::get<0>(sect);
    this->parseSection(section);
  }

  this->conductivity_on_qpoints.set(conductivity);
  this->initial_conductivity_array.set(conductivity);
  this->density_array.set(density);
  this->capacity_array.set(capacity);
  this->conductivity_variation_array.set(conductivity_variation);
  opening_on_qpoints.set(default_opening);
}
/* --------------------------------------------------------------------------
 */
void HeatTransferInterfaceModel::assembleCapacity() {
  AKANTU_DEBUG_IN();
  auto ghost_type = _not_ghost;

  Parent::assembleCapacity();

  auto & fem_interface = this->getFEEngine("InterfacesFEEngine");

  heat_transfer_interface::details::ComputeRhoFunctor rho_functor(*this);

  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_cohesive)) {
    fem_interface.assembleFieldMatrix(rho_functor, "M", "temperature",
                                      this->getDOFManager(), type, ghost_type);
  }

  need_to_reassemble_capacity = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferInterfaceModel::assembleConductivityMatrix() {
  AKANTU_DEBUG_IN();

  Parent::assembleConductivityMatrix();

  this->computeKLongOnQuadPoints(_not_ghost);

  if ((crack_conductivity_release[_not_ghost] ==
       crack_conductivity_matrix_release) and
      (crack_conductivity_matrix_release ==
       this->conductivity_matrix_release)) {
    return;
  }

  auto & fem_interface = this->getFEEngine("InterfacesFEEngine");

  for (auto && type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_cohesive)) {
    auto nb_element = mesh.getNbElement(type);
    if (nb_element == 0)
      continue;
    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points = fem_interface.getNbIntegrationPoints(type);
    auto && shape_derivatives = fem_interface.getShapesDerivatives(type);
    auto && shapes = fem_interface.getShapes(type);

    auto A = ExtendingOperators::getAveragingOperator(type);
    auto C = ExtendingOperators::getDifferencingOperator(type);

    auto && long_cond = this->k_long_w(type);
    Array<Real> at_bt_k_long_w_b_a(nb_element * nb_quadrature_points,
                                   nb_nodes_per_element * nb_nodes_per_element,
                                   "A^t*B^t*k_long*B*A");
    Array<Real> ct_nt_k_perp_n_c(nb_element * nb_quadrature_points,
                                 nb_nodes_per_element * nb_nodes_per_element,
                                 "C^t*N^t*k_perp*N*C");
    Array<Real> tangent_conductivity(
        nb_element * nb_quadrature_points,
        nb_nodes_per_element * nb_nodes_per_element, "tangent_conductivity");

    Matrix<Real> B_A(spatial_dimension - 1, nb_nodes_per_element);
    Matrix<Real> k_long_w_B_A(spatial_dimension - 1, nb_nodes_per_element);
    Matrix<Real> N_C(1, nb_nodes_per_element);

    for (auto && data :
         zip(make_view(long_cond, spatial_dimension - 1, spatial_dimension - 1),
             make_view(shape_derivatives, spatial_dimension - 1,
                       nb_nodes_per_element / 2),
             make_view(at_bt_k_long_w_b_a, nb_nodes_per_element,
                       nb_nodes_per_element),
             make_view(shapes, 1, nb_nodes_per_element / 2),
             make_view(ct_nt_k_perp_n_c, nb_nodes_per_element,
                       nb_nodes_per_element),
             make_view(tangent_conductivity, nb_nodes_per_element,
                       nb_nodes_per_element))) {
      // assemble conductivity contribution along crack
      const auto & _k_long_w = std::get<0>(data);
      const auto & B = std::get<1>(data);
      auto & At_Bt_k_long_w_B_A = std::get<2>(data);

      B_A.mul<false, false>(B, A);
      k_long_w_B_A.mul<false, false>(_k_long_w, B_A);
      At_Bt_k_long_w_B_A.mul<true, false>(k_long_w_B_A, B_A);

      // assemble conductivity contribution perpendicular to the crack
      const auto & N = std::get<3>(data);
      auto & Ct_Nt_k_perp_N_C = std::get<4>(data);
      auto & tangent = std::get<5>(data);

      N_C.mul<false, false>(N, C);
      Ct_Nt_k_perp_N_C.mul<true, false>(N_C, N_C,
                                        this->transversal_conductivity);

      // summing two contributions of tangent operator
      tangent = At_Bt_k_long_w_B_A + Ct_Nt_k_perp_N_C;
    }

    auto K_e = std::make_unique<Array<Real>>(
        nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_e");

    fem_interface.integrate(tangent_conductivity, *K_e,
                            nb_nodes_per_element * nb_nodes_per_element, type);

    this->getDOFManager().assembleElementalMatricesToMatrix(
        "K", "temperature", *K_e, type, _unsymmetric);
  }

  crack_conductivity_matrix_release = this->conductivity_matrix_release;
  crack_conductivity_release[_not_ghost] = this->conductivity_matrix_release;

  AKANTU_DEBUG_OUT();
}
/* --------------------------------------------------------------------------
 */
void HeatTransferInterfaceModel::computeGradAndDeltaT(GhostType ghost_type) {
  auto & fem_interface =
      this->getFEEngineClass<MyFEEngineCohesiveType>("InterfacesFEEngine");
  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_cohesive)) {
    auto nb_quadrature_points = fem_interface.getNbIntegrationPoints(type);
    auto nb_element = mesh.getNbElement(type, ghost_type);
    auto & gradient = temperature_gradient(type, ghost_type);
    auto surface_gradient = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points, spatial_dimension - 1,
        "surface_gradient");
    auto jump = std::make_unique<Array<Real>>(nb_element * nb_quadrature_points,
                                              1, "temperature_jump");

    this->getFEEngine("InterfacesFEEngine")
        .gradientOnIntegrationPoints(*temperature, *surface_gradient, 1, type,
                                     ghost_type);

#define COMPUTE_JUMP(type)                                                     \
  fem_interface.getShapeFunctions()                                            \
      .interpolateOnIntegrationPoints<type, CohesiveReduceFunctionDifference>( \
          *temperature, *jump, 1, ghost_type);

    AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(COMPUTE_JUMP);
#undef COMPUTE_JUMP
    for (auto && data :
         zip(make_view(gradient, spatial_dimension),
             make_view(*surface_gradient, spatial_dimension - 1), *jump)) {
      auto & grad = std::get<0>(data);
      auto & surf_grad = std::get<1>(data);
      auto & delta = std::get<2>(data);
      for (auto i : arange(spatial_dimension - 1)) {
        grad(i) = surf_grad(i);
      }
      grad(spatial_dimension - 1) = delta;
    }
  }
}
/* ------------------------------------------------------------------- */
void HeatTransferInterfaceModel::computeTempOnQpoints(GhostType ghost_type) {

  Parent::computeTempOnQpoints(ghost_type);

  auto & fem_interface =
      this->getFEEngineClass<MyFEEngineCohesiveType>("InterfacesFEEngine");
  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_cohesive)) {
    auto & t_on_qpoints = temperature_on_qpoints_coh(type, ghost_type);

#define COMPUTE_T(type)                                                        \
  fem_interface.getShapeFunctions()                                            \
      .interpolateOnIntegrationPoints<type, CohesiveReduceFunctionMean>(       \
          *temperature, t_on_qpoints, 1, ghost_type);

    AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(COMPUTE_T);
#undef COMPUTE_T
  }
  this->synchronizeField(SynchronizationTag::_htm_temperature_on_qpoints);
}
/* -------------------------------------------------------------------------- */
void HeatTransferInterfaceModel::updateNormalOpeningAtQuadraturePoints(
    Array<Real> positions, GhostType ghost_type) {
  auto & fem_interface =
      this->getFEEngineClass<MyFEEngineCohesiveType>("InterfacesFEEngine");
  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_cohesive)) {
    auto nb_quadrature_points = fem_interface.getNbIntegrationPoints(type);
    auto nb_element = mesh.getNbElement(type, ghost_type);

    auto & normal_openings = opening_on_qpoints(type, ghost_type);
    auto opening_vectors =
        std::make_unique<Array<Real>>(nb_element * nb_quadrature_points,
                                      spatial_dimension, "opening_vectors");
    auto normal_vectors = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points, spatial_dimension, "normal_vectors");

#define COMPUTE_JUMP(type)                                                     \
  fem_interface.getShapeFunctions()                                            \
      .interpolateOnIntegrationPoints<type, CohesiveReduceFunctionDifference>( \
          positions, *opening_vectors, spatial_dimension, ghost_type);

    AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(COMPUTE_JUMP);
#undef COMPUTE_JUMP
#define COMPUTE_NORMALS(type)                                                  \
  fem_interface.getShapeFunctions()                                            \
      .computeNormalsOnIntegrationPoints<type, CohesiveReduceFunctionMean>(    \
          mesh.getNodes(), *normal_vectors, ghost_type);
    AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(COMPUTE_NORMALS);
#undef COMPUTE_NORMALS

    for (auto && data :
         zip(make_view(*opening_vectors, spatial_dimension),
             make_view(*normal_vectors, spatial_dimension), normal_openings)) {
      auto & opening = std::get<0>(data);
      auto & normal = std::get<1>(data);
      auto & normal_opening = std::get<2>(data);
      normal_opening = opening.dot(normal);
    }
  }
}
/* -------------------------------------------------------------------------- */
void HeatTransferInterfaceModel::computeNormalsAtQuadraturePoints(
    Array<Real> positions, Array<Real> normals, GhostType ghost_type) {
  auto & fem_interface =
      this->getFEEngineClass<MyFEEngineCohesiveType>("InterfacesFEEngine");
  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_cohesive)) {

#define COMPUTE_NORMALS(type)                                                  \
  fem_interface.getShapeFunctions()                                            \
      .computeNormalsOnIntegrationPoints<type, CohesiveReduceFunctionMean>(    \
          positions, normals, ghost_type);
    AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(COMPUTE_NORMALS);
#undef COMPUTE_NORMALS
  }
}

/* -------------------------------------------------------------------------- */
void HeatTransferInterfaceModel::computeKLongOnQuadPoints(
    GhostType ghost_type) {

  // if opening did not change, longitudinal conductivity will neither
  if (opening_release == crack_conductivity_release[ghost_type]) {
    return;
  }

  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_cohesive)) {
    auto & opening = opening_on_qpoints(type, ghost_type);
    auto & long_cond_w = k_long_w(type, ghost_type);
    UInt dim = spatial_dimension - 1;

    Matrix<Real> long_cond_vect(dim, dim);
    for (auto i : arange(dim)) {
      long_cond_vect(i, i) = this->longitudinal_conductivity;
    }

    for (auto && tuple : zip(make_view(long_cond_w, dim, dim), opening)) {
      auto & k = std::get<0>(tuple);
      auto w = std::get<1>(tuple);
      // limiting the lower limit of w to avoid numeric instability
      if (w < this->default_opening) {
        w = this->default_opening;
      }
      k = long_cond_vect * w;
    }
  }

  crack_conductivity_release[ghost_type] = opening_release;

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------- */
void HeatTransferInterfaceModel::assembleInternalHeatRate() {
  AKANTU_DEBUG_IN();

  Parent::assembleInternalHeatRate();

  computeTempOnQpoints(_not_ghost);
  computeGradAndDeltaT(_not_ghost);
  // communicate the stresses
  AKANTU_DEBUG_INFO("Send data for residual assembly");
  this->asynchronousSynchronize(SynchronizationTag::_htm_gradient_temperature);

  computeLongHeatRate(_not_ghost);
  computeTransHeatRate(_not_ghost);
  computeInertialHeatRate(_not_ghost);

  // finalize communications
  AKANTU_DEBUG_INFO("Wait distant stresses");
  this->waitEndSynchronize(SynchronizationTag::_htm_gradient_temperature);

  computeLongHeatRate(_ghost);
  computeTransHeatRate(_ghost);
  computeInertialHeatRate(_ghost);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------*/
void HeatTransferInterfaceModel::computeLongHeatRate(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  computeKLongOnQuadPoints(ghost_type);
  auto & internal_rate = const_cast<Array<Real> &>(this->getInternalHeatRate());

  auto & fem = this->getFEEngine("InterfacesFEEngine");
  for (auto type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_cohesive)) {

    const auto & shapes_derivatives =
        fem.getShapesDerivatives(type, ghost_type);

    auto & gradient = temperature_gradient(type, ghost_type);

    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);
    auto nb_element = mesh.getNbElement(type, ghost_type);

    auto A = ExtendingOperators::getAveragingOperator(type);
    Matrix<Real> B_A(spatial_dimension - 1, nb_nodes_per_element);
    Array<Real> k_long_w_gradT(nb_element * nb_quadrature_points,
                               spatial_dimension - 1);
    Array<Real> bt_k_w_gradT(nb_element * nb_quadrature_points,
                             nb_nodes_per_element);
    for (auto && values :
         zip(make_view(k_long_w(type, ghost_type), spatial_dimension - 1,
                       spatial_dimension - 1),
             make_view(gradient, spatial_dimension),
             make_view(k_long_w_gradT, spatial_dimension - 1),
             make_view(shapes_derivatives, spatial_dimension - 1,
                       nb_nodes_per_element / 2),
             make_view(bt_k_w_gradT, nb_nodes_per_element))) {
      const auto & k_w = std::get<0>(values);
      const auto & full_gradT = std::get<1>(values);
      Vector<Real> surf_gradT(spatial_dimension - 1);
      for (auto i : arange(spatial_dimension - 1)) {
        surf_gradT(i) = full_gradT(i);
      }

      auto & k_w_gradT = std::get<2>(values);
      k_w_gradT.mul<false>(k_w, surf_gradT);
      /// compute @f$B_a t_i@f$
      const auto & B = std::get<3>(values);
      auto & At_Bt_vector = std::get<4>(values);
      B_A.mul<false, false>(B, A);
      At_Bt_vector.template mul<true>(B_A, k_w_gradT);
    }

    Array<Real> long_heat_rate(nb_element, nb_nodes_per_element);

    fem.integrate(bt_k_w_gradT, long_heat_rate, nb_nodes_per_element, type,
                  ghost_type);

    /// assemble
    this->getDOFManager().assembleElementalArrayLocalArray(
        long_heat_rate, internal_rate, type, ghost_type, 1);
  }

  AKANTU_DEBUG_OUT();
}
/* --------------------------------------------------------------------------
 */
void HeatTransferInterfaceModel::computeTransHeatRate(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & internal_rate = const_cast<Array<Real> &>(this->getInternalHeatRate());

  auto & fem = this->getFEEngine("InterfacesFEEngine");
  for (auto type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_cohesive)) {

    const auto & shapes = fem.getShapes(type, ghost_type);
    auto size_of_shapes = shapes.getNbComponent();

    auto & gradient = temperature_gradient(type, ghost_type);

    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);
    auto nb_element = mesh.getNbElement(type, ghost_type);
    auto natural_dimension = spatial_dimension - 1;

    // difference computing operator
    auto C = ExtendingOperators::getDifferencingOperator(type);

    auto nt_k_deltaT = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points, nb_nodes_per_element);
    auto k_perp_deltaT = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points, 1, "k_perp_deltaT");

    for (auto && values :
         zip(*k_perp_deltaT, /*k_perp_over_w(type, ghost_type),*/
             make_view(gradient, spatial_dimension),
             make_view(shapes, size_of_shapes),
             make_view(*nt_k_deltaT, nb_nodes_per_element))) {
      auto & _k_perp_deltaT = std::get<0>(values);
      const auto & full_gradT = std::get<1>(values);
      const auto & N = std::get<2>(values);
      auto & Nt_k_deltaT = std::get<3>(values);

      _k_perp_deltaT =
          this->transversal_conductivity * full_gradT(spatial_dimension - 1);

      Nt_k_deltaT.mul<true>(C, N, _k_perp_deltaT);
    }

    auto perp_heat_rate = std::make_unique<Array<Real>>(
        nb_element, nb_nodes_per_element, "perp_heat_rate");

    fem.integrate(*nt_k_deltaT, *perp_heat_rate, nb_nodes_per_element, type,
                  ghost_type);

    /// assemble
    this->getDOFManager().assembleElementalArrayLocalArray(
        *perp_heat_rate, internal_rate, type, ghost_type, 1);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------*/
void HeatTransferInterfaceModel::computeInertialHeatRate(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  if (not this->use_opening_rate)
    return;

  auto & internal_rate = const_cast<Array<Real> &>(this->getInternalHeatRate());

  auto & fem_interface =
      this->getFEEngineClass<MyFEEngineCohesiveType>("InterfacesFEEngine");

  for (auto type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_cohesive)) {

    const auto & shapes = fem_interface.getShapes(type, ghost_type);
    UInt size_of_shapes = shapes.getNbComponent();

    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points =
        fem_interface.getNbIntegrationPoints(type, ghost_type);
    auto nb_element = mesh.getNbElement(type, ghost_type);
    auto natural_dimension = spatial_dimension - 1;

    // averaging operator
    auto A = ExtendingOperators::getAveragingOperator(type);

    auto && opening_rates = this->opening_rate(type, ghost_type);
    auto * nt_opening_rate = new Array<Real>(nb_element * nb_quadrature_points,
                                             nb_nodes_per_element);

    for (auto && values : zip(make_view(shapes, size_of_shapes),
                              make_view(*nt_opening_rate, nb_nodes_per_element),
                              opening_rates)) {
      const auto & N = std::get<0>(values);
      auto & Nt_opening_rate = std::get<1>(values);
      auto & open_rate = std::get<2>(values);
      Nt_opening_rate.mul<true>(A, N, open_rate);
    }

    auto * inertial_heat_rate =
        new Array<Real>(nb_element, nb_nodes_per_element, "inertial_heat_rate");

    fem_interface.integrate(*nt_opening_rate, *inertial_heat_rate,
                            nb_nodes_per_element, type, ghost_type);
    delete nt_opening_rate;

    /// assemble
    this->getDOFManager().assembleElementalArrayLocalArray(
        *inertial_heat_rate, internal_rate, type, ghost_type, 1);

    delete inertial_heat_rate;
  }

  AKANTU_DEBUG_OUT();
}
/* --------------------------------------------------------------------------
 */
Real HeatTransferInterfaceModel::getStableTimeStep() {
  AKANTU_DEBUG_IN();

  auto HTM_time_step = Parent::getStableTimeStep();

  Real el_size;
  Real min_el_size = std::numeric_limits<Real>::max();
  Real max_opening_on_qpoint = std::numeric_limits<Real>::min();

  for (const auto & type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_cohesive)) {

    auto nb_nodes_per_element = mesh.getNbNodesPerElement(type);
    auto & opening = opening_on_qpoints(type, _not_ghost);
    auto nb_quad_per_element = opening.getNbComponent();
    auto mesh_dim = this->mesh.getSpatialDimension();
    const auto facet_type = Mesh::getFacetType(type);

    Array<Real> coord(0, nb_nodes_per_element * mesh_dim);
    FEEngine::extractNodalToElementField(mesh, mesh.getNodes(), coord, type,
                                         _not_ghost);

    for (auto && data : zip(make_view(coord, mesh_dim, nb_nodes_per_element),
                            make_view(opening, nb_quad_per_element))) {
      Matrix<Real> & el_coord = std::get<0>(data);
      Vector<Real> & el_opening = std::get<1>(data);
      el_size = getFEEngine().getElementInradius(el_coord, facet_type);
      min_el_size = std::min(min_el_size, el_size);
      max_opening_on_qpoint =
          std::max(max_opening_on_qpoint, el_opening.norm<L_inf>());
    }

    AKANTU_DEBUG_INFO("The minimum element size : "
                      << min_el_size
                      << " and the max opening is : " << max_opening_on_qpoint);
  }

  Real min_dt = min_el_size * min_el_size / 4 * this->density_in_crack *
                this->capacity_in_crack / this->longitudinal_conductivity;

  mesh.getCommunicator().allReduce(min_dt, SynchronizerOperation::_min);

  min_dt = std::min(min_dt, HTM_time_step);

  AKANTU_DEBUG_OUT();

  return min_dt;
}
/* --------------------------------------------------------------------------
 */

void HeatTransferInterfaceModel::setTimeStep(Real time_step,
                                             const ID & solver_id) {
  Parent::setTimeStep(time_step, solver_id);

  this->mesh.getDumper("heat_interfaces").setTimeStep(time_step);
}

/* -------------------------------------------------------------------------*/
std::shared_ptr<dumpers::Field>
HeatTransferInterfaceModel::createElementalField(const std::string & field_name,
                                                 const std::string & group_name,
                                                 bool padding_flag,
                                                 UInt spatial_dimension,
                                                 ElementKind element_kind) {
  std::shared_ptr<dumpers::Field> field;

  auto getNbDataPerElem = [&](auto & field) {
    const auto & fe_engine = getFEEngine("InterfacesFEEngine");
    ElementTypeMap<UInt> res;
    for (auto ghost_type : ghost_types) {
      for (auto & type :
           field.elementTypes(spatial_dimension, ghost_type, element_kind)) {
        auto nb_quadrature_points =
            fe_engine.getNbIntegrationPoints(type, ghost_type);
        res(type, ghost_type) =
            field(type, ghost_type).getNbComponent() * nb_quadrature_points;
      }
    }

    return res;
  };

  if (element_kind == _ek_regular) {
    field = Parent::createElementalField(field_name, group_name, padding_flag,
                                         spatial_dimension, element_kind);
  } else if (element_kind == _ek_cohesive) {
    if (field_name == "partitions") {
      field = mesh.createElementalField<UInt, dumpers::ElementPartitionField>(
          mesh.getConnectivities(), group_name, this->spatial_dimension,
          element_kind);
    } else if (field_name == "opening_on_qpoints") {
      auto nb_data_per_elem = getNbDataPerElem(opening_on_qpoints);
      field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
          opening_on_qpoints, group_name, this->spatial_dimension, element_kind,
          nb_data_per_elem);
    } else if (field_name == "temperature_on_qpoints_coh") {
      auto nb_data_per_elem = getNbDataPerElem(temperature_on_qpoints_coh);
      field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
          temperature_on_qpoints_coh, group_name, this->spatial_dimension,
          element_kind, nb_data_per_elem);
    } else if (field_name == "temperature_gradient") {
      auto nb_data_per_elem = getNbDataPerElem(temperature_gradient);

      field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
          temperature_gradient, group_name, this->spatial_dimension,
          element_kind, nb_data_per_elem);
    }
  }
  return field;
}
} // namespace akantu
