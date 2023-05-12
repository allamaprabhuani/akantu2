/**
 * @file   heat_transfer_interface_model.hh
 *
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 *
 * @date creation: Thu Apr 13 2023
 * @date last modification: Thu Apr 13 2023
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
        auto opening_it = model.getOpening(el.type, el.ghost_type).begin();

        rho.set(model.getCapacityInCrack() * model.getDensityInCrack() *
                opening_it[el.element]);
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
      k_long_w_gradT_on_qpoints("k_long_w_gradT_on_qpoints", id),
      k_perp_deltaT_over_w_on_qpoints("k_perp_deltaT_over_w_on_qpoints", id),
      temperature_gradient_on_surface("temperature_gradient_on_surface", id),
      temperature_jump("temperature_jump", id),
      opening_on_qpoints("opening_on_qpoints", id),
      opening_rate("opening_rate", id), k_long_w("k_long_w", id),
      k_perp_over_w("k_perp_over_w", id) {
  AKANTU_DEBUG_IN();

  registerFEEngineObject<MyFEEngineCohesiveType>("InterfacesFEEngine", mesh,
                                                 Model::spatial_dimension);

  this->mesh.registerDumper<DumperParaview>("heat_interfaces", id);
  this->mesh.addDumpMeshToDumper("heat_interfaces", mesh, spatial_dimension,
                                 _not_ghost, _ek_cohesive);
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

  this->temperature_gradient_on_surface.initialize(
      fem, _nb_component = spatial_dimension - 1);
  this->temperature_jump.initialize(fem, _nb_component = 1);
  this->k_long_w.initialize(fem, _nb_component = (spatial_dimension - 1) *
                                                 (spatial_dimension - 1));
  this->k_perp_over_w.initialize(fem, _nb_component = 1);
  this->k_long_w_gradT_on_qpoints.initialize(fem, _nb_component =
                                                      spatial_dimension - 1);
  this->k_perp_deltaT_over_w_on_qpoints.initialize(fem, _nb_component = 1);
  this->opening_on_qpoints.initialize(fem, _nb_component = 1);
  if (use_opening_rate) {
    opening_rate.initialize(fem, _nb_component = 1);
  }
}

/* -------------------------------------------------------------------------- */
HeatTransferInterfaceModel::~HeatTransferInterfaceModel() = default;
/* -------------------------------------------------------------------------- */
void HeatTransferInterfaceModel::initFullImpl(const ModelOptions & options) {
  Model::initFullImpl(options);

  readMaterials();
}
/* -------------------------------------------------------------------------- */
void HeatTransferInterfaceModel::predictor() {
  Parent::predictor();
  ++opening_release;
}

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
// void HeatTransferInterfaceModel::assembleCapacityLumped(GhostType ghost_type)
// {
//   AKANTU_DEBUG_IN();

//   Parent::assembleCapacityLumped(ghost_type);

//   auto & fem_interface =
//       this->getFEEngineClass<MyFEEngineCohesiveType>("InterfacesFEEngine");
//   heat_transfer_interface::details::ComputeRhoFunctor compute_rho(*this);

//   for (auto && type :
//        mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
//     fem_interface.assembleFieldLumped(compute_rho, "M", "temperature",
//                                       this->getDOFManager(), type,
//                                       ghost_type);
//   }

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
void HeatTransferInterfaceModel::assembleConductivityMatrix() {
  AKANTU_DEBUG_IN();

  Parent::assembleConductivityMatrix();

  this->computeKLongOnQuadPoints(_not_ghost);
  this->computeKTransOnQuadPoints(_not_ghost);

  if ((long_conductivity_release[_not_ghost] ==
       this->crack_conductivity_matrix_release) and
      (perp_conductivity_release[_not_ghost] ==
       this->crack_conductivity_matrix_release)) {
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

    auto A = getAveragingOperator(type);
    auto C = getDifferencingOperator(type);

    // tangent_conductivity_matrix.zero();
    auto && long_cond = this->k_long_w(type);
    auto && perp_cond = this->k_perp_over_w(type);
    Array<Real> at_bt_k_long_b_a(nb_element * nb_quadrature_points,
                                 nb_nodes_per_element * nb_nodes_per_element,
                                 "A^t*B^t*k_long*B*A");
    Array<Real> ct_nt_k_perp_n_c(nb_element * nb_quadrature_points,
                                 nb_nodes_per_element * nb_nodes_per_element,
                                 "C^t*N^t*k_perp*N*C");
    Array<Real> tangent_conductivity(
        nb_element * nb_quadrature_points,
        nb_nodes_per_element * nb_nodes_per_element, "tangent_conductivity");

    Matrix<Real> B_A(spatial_dimension - 1, nb_nodes_per_element);
    Matrix<Real> k_long_B_A(spatial_dimension - 1, nb_nodes_per_element);
    // Matrix<Real> N(spatial_dimension, nb_nodes_per_element);
    Matrix<Real> N_C(1, nb_nodes_per_element);

    for (auto && data :
         zip(make_view(long_cond, spatial_dimension - 1, spatial_dimension - 1),
             make_view(shape_derivatives, spatial_dimension - 1,
                       nb_nodes_per_element / 2),
             make_view(at_bt_k_long_b_a, nb_nodes_per_element,
                       nb_nodes_per_element),
             perp_cond, make_view(shapes, 1, nb_nodes_per_element / 2),
             make_view(ct_nt_k_perp_n_c, nb_nodes_per_element,
                       nb_nodes_per_element),
             make_view(tangent_conductivity, nb_nodes_per_element,
                       nb_nodes_per_element))) {
      // assemble conductivity contribution along crack
      const auto & k_long = std::get<0>(data);
      const auto & B = std::get<1>(data);
      auto & At_Bt_k_long_B_A = std::get<2>(data);

      B_A.mul<false, false>(B, A);
      k_long_B_A.mul<false, false>(k_long, B_A);
      At_Bt_k_long_B_A.mul<true, false>(k_long_B_A, B_A);

      // assemble conductivity contribution perpendicular to the crack
      const auto & k_perp = std::get<3>(data);
      const auto & N = std::get<4>(data);
      auto & Ct_Nt_k_perp_N_C = std::get<5>(data);
      auto & tangent = std::get<6>(data);

      N_C.mul<false, false>(N, C);
      Ct_Nt_k_perp_N_C.mul<true, false>(N_C, N_C, k_perp);

      // summing two contributions of tangent operator
      tangent = At_Bt_k_long_B_A + Ct_Nt_k_perp_N_C;
    }

    auto K_e = std::make_unique<Array<Real>>(
        nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_e");

    fem_interface.integrate(tangent_conductivity, *K_e,
                            nb_nodes_per_element * nb_nodes_per_element, type);

    this->getDOFManager().assembleElementalMatricesToMatrix(
        "K", "temperature", *K_e, type, _not_ghost, _unsymmetric);
  }

  this->crack_conductivity_matrix_release =
      long_conductivity_release[_not_ghost];

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
// void HeatTransferInterfaceModel::computeNormal(const Array<Real> & position,
//                                                Array<Real> & normal,
//                                                ElementType type,
//                                                GhostType ghost_type) {
//   AKANTU_DEBUG_IN();

//   auto & fem_interface =
//       this->getFEEngineClass<MyFEEngineCohesiveType>("InterfacesFEEngine");

//   normal.zero();

// #define COMPUTE_NORMAL(type)                                                   \
//   fem_interface.getShapeFunctions()                                            \
//       .computeNormalsOnIntegrationPoints<type, CohesiveReduceFunctionMean>(    \
//           position, normal, ghost_type);

//   AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(COMPUTE_NORMAL);
// #undef COMPUTE_NORMAL

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
// void HeatTransferInterfaceModel::computeBasis(const Array<Real> & position,
//                                               Array<Real> & basis,
//                                               ElementType type,
//                                               GhostType ghost_type) {
//   AKANTU_DEBUG_IN();

//   auto & fem_interface =
//       this->getFEEngineClass<MyFEEngineCohesiveType>("InterfacesFEEngine");

//   basis.zero();

// #define COMPUTE_BASIS(type)                                                    \
//   fem_interface.getShapeFunctions()                                            \
//       .computeAllignedBasisOnIntegrationPoints<type,                           \
//                                                CohesiveReduceFunctionMean>(    \
//           position, basis, ghost_type);

//   AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(COMPUTE_BASIS);
// #undef COMPUTE_BASIS

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
void HeatTransferInterfaceModel::computeKLongOnQuadPoints(
    GhostType ghost_type) {

  // if opening did not change, longitudinal conductivity will neither
  if (opening_release == long_conductivity_release[ghost_type]) {
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
      auto & w = std::get<1>(tuple);

      k = long_cond_vect * w;
    }
  }

  long_conductivity_release[ghost_type] = opening_release;

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
void HeatTransferInterfaceModel::computeKTransOnQuadPoints(
    GhostType ghost_type) {

  // if opening did not change, longitudinal conductivity will neither
  if (opening_release == perp_conductivity_release[ghost_type]) {
    return;
  }

  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_cohesive)) {
    auto & opening = opening_on_qpoints(type, ghost_type);
    auto & trans_cond_over_w = k_perp_over_w(type, ghost_type);

    for (auto && tuple : zip(trans_cond_over_w, opening)) {
      auto & k_perp = std::get<0>(tuple);
      auto & w = std::get<1>(tuple);

      k_perp = transversal_conductivity / w;
    }
  }
  perp_conductivity_release[ghost_type] = opening_release;

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
void HeatTransferInterfaceModel::computeKLongGradT(GhostType ghost_type) {

  computeKLongOnQuadPoints(ghost_type);
  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_cohesive)) {
    UInt dim = spatial_dimension - 1;
    auto & gradient = temperature_gradient_on_surface(type, ghost_type);
    auto & opening = opening_on_qpoints(type, ghost_type);

    this->getFEEngine("InterfacesFEEngine")
        .gradientOnIntegrationPoints(*temperature, gradient, 1, type,
                                     ghost_type);
    for (auto && values :
         zip(make_view(k_long_w(type, ghost_type), dim, dim),
             make_view(gradient, dim),
             make_view(k_long_w_gradT_on_qpoints(type, ghost_type), dim))) {
      const auto & k_w = std::get<0>(values);
      const auto & gradT = std::get<1>(values);
      auto & k_w_gradT = std::get<2>(values);

      k_w_gradT.mul<false>(k_w, gradT);
    }
  }
  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
void HeatTransferInterfaceModel::computeKTransDeltaT(GhostType ghost_type) {

  computeKTransOnQuadPoints(ghost_type);

  auto & fem_interface =
      this->getFEEngineClass<MyFEEngineCohesiveType>("InterfacesFEEngine");

  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_cohesive)) {
    auto & jump = temperature_jump(type, ghost_type);

#define COMPUTE_JUMP(type)                                                     \
  fem_interface.getShapeFunctions()                                            \
      .interpolateOnIntegrationPoints<type, CohesiveReduceFunctionOpening>(    \
          *temperature, jump, 1, ghost_type);

    AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(COMPUTE_JUMP);
#undef COMPUTE_JUMP

    for (auto && values : zip(k_perp_deltaT_over_w_on_qpoints(type, ghost_type),
                              k_perp_over_w(type, ghost_type), jump)) {
      auto & k_perp_deltaT_over_w = std::get<0>(values);
      auto & k_perp_over_w = std::get<1>(values);
      const auto & deltaT = std::get<2>(values);
      k_perp_deltaT_over_w = k_perp_over_w * deltaT;
    }
  }
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
void HeatTransferInterfaceModel::assembleInternalHeatRate() {
  AKANTU_DEBUG_IN();

  Parent::assembleInternalHeatRate();

  for (auto ghost_type : ghost_types) {
    computeKLongGradT(ghost_type);
    computeKTransDeltaT(ghost_type);

    computeLongHeatRate(ghost_type);
    computeTransHeatRate(ghost_type);

    computeInertialHeatRate(ghost_type);
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
void HeatTransferInterfaceModel::computeLongHeatRate(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & internal_rate = const_cast<Array<Real> &>(this->getInternalHeatRate());

  auto & fem = this->getFEEngine("InterfacesFEEngine");
  for (auto type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_cohesive)) {

    const auto & shapes_derivatives =
        fem.getShapesDerivatives(type, ghost_type);

    auto & k_w_gradT = k_long_w_gradT_on_qpoints(type, ghost_type);

    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);
    auto nb_element = mesh.getNbElement(type, ghost_type);
    auto natural_dimension = spatial_dimension - 1;

    auto A = getAveragingOperator(type);
    Matrix<Real> B_A(natural_dimension, nb_nodes_per_element);

    /// compute @f$B_a t_i@f$
    auto * bt_k_w_gradT = new Array<Real>(nb_element * nb_quadrature_points,
                                          nb_nodes_per_element);

    for (auto && values : zip(make_view(shapes_derivatives, natural_dimension,
                                        nb_nodes_per_element / 2),
                              make_view(k_w_gradT, natural_dimension),
                              make_view(*bt_k_w_gradT, nb_nodes_per_element))) {
      const auto & B = std::get<0>(values);
      const auto & vector = std::get<1>(values);
      auto & At_Bt_vector = std::get<2>(values);
      B_A.mul<false, false>(B, A);
      At_Bt_vector.template mul<true>(B_A, vector);
    }

    auto * long_heat_rate =
        new Array<Real>(nb_element, nb_nodes_per_element, "long_heat_rate");

    fem.integrate(*bt_k_w_gradT, *long_heat_rate, nb_nodes_per_element, type,
                  ghost_type);

    delete bt_k_w_gradT;

    /// assemble
    this->getDOFManager().assembleElementalArrayLocalArray(
        *long_heat_rate, internal_rate, type, ghost_type, 1);

    delete long_heat_rate;
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

    auto & k_perp_deltaT = k_perp_deltaT_over_w_on_qpoints(type, ghost_type);

    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);
    auto nb_element = mesh.getNbElement(type, ghost_type);
    auto natural_dimension = spatial_dimension - 1;

    // difference computing operator
    auto C = getDifferencingOperator(type);

    auto * nt_k_deltaT = new Array<Real>(nb_element * nb_quadrature_points,
                                         nb_nodes_per_element);

    for (auto && values :
         zip(make_view(shapes, size_of_shapes),
             make_view(*nt_k_deltaT, nb_nodes_per_element), k_perp_deltaT)) {
      const auto & N = std::get<0>(values);
      auto & Nt_k_deltaT = std::get<1>(values);
      auto & k_deltaT = std::get<2>(values);
      Nt_k_deltaT.mul<true>(C, N, k_deltaT);
    }

    auto * perp_heat_rate =
        new Array<Real>(nb_element, nb_nodes_per_element, "perp_heat_rate");

    fem.integrate(*nt_k_deltaT, *perp_heat_rate, nb_nodes_per_element, type,
                  ghost_type);

    delete nt_k_deltaT;

    /// assemble
    this->getDOFManager().assembleElementalArrayLocalArray(
        *perp_heat_rate, internal_rate, type, ghost_type, 1);

    delete perp_heat_rate;
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
    auto A = getAveragingOperator(type);

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
  }

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------*/
Matrix<Real> HeatTransferInterfaceModel::getAveragingOperator(
    const ElementType & type, const UInt & nb_degree_of_freedom) {
  AKANTU_DEBUG_IN();
  auto kind = Mesh::getKind(type);
  AKANTU_DEBUG_ASSERT(kind == _ek_cohesive,
                      "Extending operators work only for cohesive elements");
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

  // averaging operator
  Matrix<Real> A(nb_degree_of_freedom * nb_nodes_per_element / 2,
                 nb_degree_of_freedom * nb_nodes_per_element);

  for (UInt i = 0; i < nb_degree_of_freedom * nb_nodes_per_element / 2; ++i) {
    A(i, i) = 0.5;
    A(i, i + nb_degree_of_freedom * nb_nodes_per_element / 2) = 0.5;
  }
  AKANTU_DEBUG_OUT();
  return A;
}

/* --------------------------------------------------------------------------*/
Matrix<Real> HeatTransferInterfaceModel::getDifferencingOperator(
    const ElementType & type, const UInt & nb_degree_of_freedom) {
  AKANTU_DEBUG_IN();
  auto kind = Mesh::getKind(type);
  AKANTU_DEBUG_ASSERT(kind == _ek_cohesive,
                      "Extending operators work only for cohesive elements");
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

  // averaging operator
  Matrix<Real> A(nb_degree_of_freedom * nb_nodes_per_element / 2,
                 nb_degree_of_freedom * nb_nodes_per_element);

  for (UInt i = 0; i < nb_degree_of_freedom * nb_nodes_per_element / 2; ++i) {
    A(i, i) = 1;
    A(i, i + nb_degree_of_freedom * nb_nodes_per_element / 2) = -1;
  }
  AKANTU_DEBUG_OUT();
  return A;
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

  Real min_dt = 2. * min_el_size * min_el_size / 4 * density_in_crack *
                capacity_in_crack / longitudinal_conductivity;

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

// /* --------------------------------------------------------------------------
// */ void HeatTransferInterfaceModel::assignPropertyToPhysicalGroup(
//     const std::string & property_name, const std::string & group_name,
//     Real value) {
//   AKANTU_DEBUG_ASSERT(property_name != "conductivity",
//                       "Scalar was provided instead of a conductivity
//                       matrix");
//   auto && el_group = mesh.getElementGroup(group_name);
//   auto & fem = this->getFEEngine();

//   for (auto ghost_type : ghost_types) {
//     for (auto && type :
//          mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
//       auto nb_quadrature_points = fem.getNbIntegrationPoints(type);
//       auto && elements = el_group.getElements(type, ghost_type);

//       Array<Real>::scalar_iterator field_it;
//       if (property_name == "density") {
//         field_it = density_array(type, ghost_type).begin();
//       } else if (property_name == "capacity") {
//         field_it = capacity_array(type, ghost_type).begin();
//       } else {
//         AKANTU_EXCEPTION(property_name +
//                          " is not a valid material property name.");
//       }
//       for (auto && el : elements) {
//         for (auto && qpoint : arange(nb_quadrature_points)) {
//           field_it[el * nb_quadrature_points + qpoint] = value;
//         }
//       }
//       need_to_reassemble_capacity = true;
//       need_to_reassemble_capacity_lumped = true;
//     }
//   }
// }
// /* --------------------------------------------------------------------------
// */ void HeatTransferInterfaceModel::assignPropertyToPhysicalGroup(
//     const std::string & property_name, const std::string & group_name,
//     Matrix<Real> cond_matrix) {
//   AKANTU_DEBUG_ASSERT(property_name == "conductivity",
//                       "When updating material parameters, only conductivity "
//                       "accepts matrix as an input");
//   auto && el_group = mesh.getElementGroup(group_name);
//   auto & fem = this->getFEEngine();
//   auto dim = this->getMesh().getSpatialDimension();

//   for (auto ghost_type : ghost_types) {
//     for (auto && type :
//          mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
//       auto nb_quadrature_points = fem.getNbIntegrationPoints(type);
//       auto && elements = el_group.getElements(type, ghost_type);

//       auto init_cond_it =
//           make_view(initial_conductivity_array(type, ghost_type), dim, dim)
//               .begin();
//       auto cond_on_quad_it =
//           make_view(conductivity_on_qpoints(type, ghost_type), dim, dim)
//               .begin();

//       for (auto && el : elements) {
//         for (auto && qpoint : arange(nb_quadrature_points)) {
//           init_cond_it[el * nb_quadrature_points + qpoint] = cond_matrix;
//           cond_on_quad_it[el * nb_quadrature_points + qpoint] = cond_matrix;
//         }
//       }
//     }
//     conductivity_release[ghost_type] += 1;
//   }
// }

/* -------------------------------------------------------------------------*/
std::shared_ptr<dumpers::Field>
HeatTransferInterfaceModel::createElementalField(const std::string & field_name,
                                                 const std::string & group_name,
                                                 bool padding_flag,
                                                 UInt spatial_dimension,
                                                 ElementKind element_kind) {
  std::shared_ptr<dumpers::Field> field;

  if (element_kind == _ek_regular) {
    field = Parent::createElementalField(field_name, group_name, padding_flag,
                                         spatial_dimension, element_kind);
  } else if (element_kind == _ek_cohesive) {
    if (field_name == "partitions") {
      field = mesh.createElementalField<UInt, dumpers::ElementPartitionField>(
          mesh.getConnectivities(), group_name, this->spatial_dimension,
          element_kind);
    } else if (field_name == "opening_on_qpoints") {
      ElementTypeMap<UInt> nb_data_per_elem =
          this->mesh.getNbDataPerElem(opening_on_qpoints);

      field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
          opening_on_qpoints, group_name, this->spatial_dimension, element_kind,
          nb_data_per_elem);
    } else if (field_name == "temperature_gradient_on_surface") {
      ElementTypeMap<UInt> nb_data_per_elem =
          this->mesh.getNbDataPerElem(temperature_gradient_on_surface);

      field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
          temperature_gradient_on_surface, group_name, this->spatial_dimension,
          element_kind, nb_data_per_elem);
    } else if (field_name == "temperature_jump") {
      ElementTypeMap<UInt> nb_data_per_elem =
          this->mesh.getNbDataPerElem(temperature_jump);
      if (!field) {
        std::cout << "toto" << std::endl;
      }
      field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
          temperature_jump, group_name, this->spatial_dimension, element_kind,
          nb_data_per_elem);
      if (!field) {
        std::cout << "toto" << std::endl;
      }
    }
  }
  return field;
}
/* -------------------------------------------------------------------------- */
} // namespace akantu
