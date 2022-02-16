/**
 * @file   material.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue Jul 27 2010
 * @date last modification: Fri Apr 09 2021
 *
 * @brief  Implementation of the common part of the material class
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
#include "material.hh"
#include "mesh_iterators.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
Material::Material(SolidMechanicsModel & model, const ID & id)
    : Parsable(ParserType::_material, id), id(id), fem(model.getFEEngine()),
      model(model), spatial_dimension(this->model.getSpatialDimension()),
      element_filter("element_filter", id), stress("stress", *this),
      eigengradu("eigen_grad_u", *this), gradu("grad_u", *this),
      green_strain("green_strain", *this),
      piola_kirchhoff_2("piola_kirchhoff_2", *this),
      potential_energy("potential_energy", *this),
      interpolation_inverse_coordinates("interpolation inverse coordinates",
                                        *this),
      interpolation_points_matrices("interpolation points matrices", *this),
      eigen_grad_u(model.getSpatialDimension(), model.getSpatialDimension()) {
  eigen_grad_u.fill(0.);

  this->registerParam("eigen_grad_u", eigen_grad_u, _pat_parsable,
                      "EigenGradU");

  /// for each connectivity types allocate the element filer array of the
  /// material
  element_filter.initialize(model.getMesh(),
                            _spatial_dimension = spatial_dimension,
                            _element_kind = _ek_regular);
  this->initialize();
}

/* -------------------------------------------------------------------------- */
Material::Material(SolidMechanicsModel & model, Int dim, const Mesh & mesh,
                   FEEngine & fe_engine, const ID & id)
    : Parsable(ParserType::_material, id), id(id), fem(fe_engine), model(model),
      spatial_dimension(dim), element_filter("element_filter", id),
      stress("stress", *this, dim, fe_engine, this->element_filter),
      eigengradu("eigen_grad_u", *this, dim, fe_engine, this->element_filter),
      gradu("gradu", *this, dim, fe_engine, this->element_filter),
      green_strain("green_strain", *this, dim, fe_engine, this->element_filter),
      piola_kirchhoff_2("piola_kirchhoff_2", *this, dim, fe_engine,
                        this->element_filter),
      potential_energy("potential_energy", *this, dim, fe_engine,
                       this->element_filter),
      interpolation_inverse_coordinates("interpolation inverse_coordinates",
                                        *this, dim, fe_engine,
                                        this->element_filter),
      interpolation_points_matrices("interpolation points matrices", *this, dim,
                                    fe_engine, this->element_filter),
      eigen_grad_u(dim, dim) {
  eigen_grad_u.fill(0.);

  element_filter.initialize(mesh, _spatial_dimension = spatial_dimension,
                            _element_kind = _ek_regular);

  this->initialize();
}

/* -------------------------------------------------------------------------- */
Material::~Material() = default;

/* -------------------------------------------------------------------------- */
void Material::initialize() {
  registerParam("rho", rho, Real(0.), _pat_parsable | _pat_modifiable,
                "Density");
  registerParam("name", name, std::string(), _pat_parsable | _pat_readable);
  registerParam("finite_deformation", finite_deformation, false,
                _pat_parsable | _pat_readable, "Is finite deformation");
  registerParam("inelastic_deformation", inelastic_deformation, false,
                _pat_internal, "Is inelastic deformation");

  /// allocate gradu stress for local elements
  eigengradu.initialize(spatial_dimension * spatial_dimension);
  gradu.initialize(spatial_dimension * spatial_dimension);
  stress.initialize(spatial_dimension * spatial_dimension);

  potential_energy.initialize(1);

  this->model.registerEventHandler(*this);
}

/* -------------------------------------------------------------------------- */
void Material::initMaterial() {
  AKANTU_DEBUG_IN();

  if (finite_deformation) {
    this->piola_kirchhoff_2.initialize(spatial_dimension * spatial_dimension);
    this->piola_kirchhoff_2.initializeHistory();
    this->green_strain.initialize(spatial_dimension * spatial_dimension);
  }

  this->stress.initializeHistory();
  this->gradu.initializeHistory();

  this->resizeInternals();

  auto dim = spatial_dimension;
  for (const auto & type :
       element_filter.elementTypes(_element_kind = _ek_regular)) {
    for (auto & eigen_gradu : make_view(eigengradu(type), dim, dim)) {
      eigen_gradu = eigen_grad_u;
    }
  }

  is_init = true;

  updateInternalParameters();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::savePreviousState() {
  AKANTU_DEBUG_IN();

  for (auto pair : internal_vectors_real) {
    if (pair.second->hasHistory()) {
      pair.second->saveCurrentValues();
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::restorePreviousState() {
  AKANTU_DEBUG_IN();

  for (auto pair : internal_vectors_real) {
    if (pair.second->hasHistory()) {
      pair.second->restorePreviousValues();
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Compute the internal forces by assembling @f$\int_{e} \sigma_e \frac{\partial
 * \varphi}{\partial X} dX @f$
 *
 * @param[in] ghost_type compute the internal forces for _ghost or _not_ghost
 * element
 */
void Material::assembleInternalForces(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Int spatial_dimension = model.getSpatialDimension();

  if (not finite_deformation) {

    auto & internal_force = model.getInternalForce();

    // Mesh & mesh = fem.getMesh();
    for (auto && type :
         element_filter.elementTypes(spatial_dimension, ghost_type)) {
      auto && elem_filter = element_filter(type, ghost_type);
      auto nb_element = elem_filter.size();

      if (nb_element == 0) {
        continue;
      }

      const Array<Real> & shapes_derivatives =
          fem.getShapesDerivatives(type, ghost_type);

      UInt size_of_shapes_derivatives = shapes_derivatives.getNbComponent();
      UInt nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);
      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

      /// compute @f$\sigma \frac{\partial \varphi}{\partial X}@f$ by
      /// @f$\mathbf{B}^t \mathbf{\sigma}_q@f$
      auto * sigma_dphi_dx =
          new Array<Real>(nb_element * nb_quadrature_points,
                          size_of_shapes_derivatives, "sigma_x_dphi_/_dX");

      fem.computeBtD(stress(type, ghost_type), *sigma_dphi_dx, type, ghost_type,
                     elem_filter);

      /**
       * compute @f$\int \sigma  * \frac{\partial \varphi}{\partial X}dX@f$ by
       * @f$ \sum_q \mathbf{B}^t
       * \mathbf{\sigma}_q \overline w_q J_q@f$
       */
      auto * int_sigma_dphi_dx =
          new Array<Real>(nb_element, nb_nodes_per_element * spatial_dimension,
                          "int_sigma_x_dphi_/_dX");

      fem.integrate(*sigma_dphi_dx, *int_sigma_dphi_dx,
                    size_of_shapes_derivatives, type, ghost_type, elem_filter);
      delete sigma_dphi_dx;

      /// assemble
      model.getDOFManager().assembleElementalArrayLocalArray(
          *int_sigma_dphi_dx, internal_force, type, ghost_type, -1,
          elem_filter);
      delete int_sigma_dphi_dx;
    }
  } else {
    switch (spatial_dimension) {
    case 1:
      this->assembleInternalForces<1>(ghost_type);
      break;
    case 2:
      this->assembleInternalForces<2>(ghost_type);
      break;
    case 3:
      this->assembleInternalForces<3>(ghost_type);
      break;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Compute  the  stress from the gradu
 *
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */
void Material::computeAllStresses(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto spatial_dimension = model.getSpatialDimension();

  for (const auto & type :
       element_filter.elementTypes(spatial_dimension, ghost_type)) {
    auto & elem_filter = element_filter(type, ghost_type);

    if (elem_filter.empty()) {
      continue;
    }
    auto & gradu_vect = gradu(type, ghost_type);

    /// compute @f$\nabla u@f$
    fem.gradientOnIntegrationPoints(model.getDisplacement(), gradu_vect,
                                    spatial_dimension, type, ghost_type,
                                    elem_filter);

    gradu_vect -= eigengradu(type, ghost_type);

    /// compute @f$\mathbf{\sigma}_q@f$ from @f$\nabla u@f$
    computeStress(type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::computeAllCauchyStresses(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(finite_deformation, "The Cauchy stress can only be "
                                          "computed if you are working in "
                                          "finite deformation.");

  for (auto type : element_filter.elementTypes(spatial_dimension, ghost_type)) {
    switch (spatial_dimension) {
    case 1:
      this->StoCauchy<1>(type, ghost_type);
      break;
    case 2:
      this->StoCauchy<2>(type, ghost_type);
      break;
    case 3:
      this->StoCauchy<3>(type, ghost_type);
      break;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void Material::StoCauchy(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  for (auto && data :
       zip(make_view<dim, dim>(this->gradu(el_type, ghost_type)),
           make_view<dim, dim>(this->piola_kirchhoff_2(el_type, ghost_type)),
           make_view<dim, dim>(this->stress(el_type, ghost_type)))) {
    auto && grad_u = std::get<0>(data);
    auto && piola = std::get<1>(data);
    auto && sigma = std::get<2>(data);

    this->StoCauchy<dim>(gradUToF<dim>(grad_u), piola, sigma);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::setToSteadyState(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const auto & displacement = model.getDisplacement();

  // resizeInternalArray(gradu);

  auto spatial_dimension = model.getSpatialDimension();

  for (auto type : element_filter.elementTypes(spatial_dimension, ghost_type)) {
    auto & elem_filter = element_filter(type, ghost_type);
    auto & gradu_vect = gradu(type, ghost_type);

    /// compute @f$\nabla u@f$
    fem.gradientOnIntegrationPoints(displacement, gradu_vect, spatial_dimension,
                                    type, ghost_type, elem_filter);

    setToSteadyState(type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Compute  the stiffness  matrix by  assembling @f$\int_{\omega}  B^t  \times D
 * \times B d\omega @f$
 *
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */
void Material::assembleStiffnessMatrix(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto spatial_dimension = model.getSpatialDimension();

  for (auto type : element_filter.elementTypes(spatial_dimension, ghost_type)) {
    if (finite_deformation) {
      switch (spatial_dimension) {
      case 1: {
        assembleStiffnessMatrixNL<1>(type, ghost_type);
        assembleStiffnessMatrixL2<1>(type, ghost_type);
        break;
      }
      case 2: {
        assembleStiffnessMatrixNL<2>(type, ghost_type);
        assembleStiffnessMatrixL2<2>(type, ghost_type);
        break;
      }
      case 3: {
        assembleStiffnessMatrixNL<3>(type, ghost_type);
        assembleStiffnessMatrixL2<3>(type, ghost_type);
        break;
      }
      }
    } else {
      switch (spatial_dimension) {
      case 1: {
        assembleStiffnessMatrix<1>(type, ghost_type);
        break;
      }
      case 2: {
        assembleStiffnessMatrix<2>(type, ghost_type);
        break;
      }
      case 3: {
        assembleStiffnessMatrix<3>(type, ghost_type);
        break;
      }
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void Material::assembleStiffnessMatrix(ElementType type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const auto & elem_filter = element_filter(type, ghost_type);
  if (elem_filter.empty()) {
    AKANTU_DEBUG_OUT();
    return;
  }

  // const Array<Real> & shapes_derivatives =
  //     fem.getShapesDerivatives(type, ghost_type);

  auto & gradu_vect = gradu(type, ghost_type);

  auto nb_element = elem_filter.size();
  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

  gradu_vect.resize(nb_quadrature_points * nb_element);

  fem.gradientOnIntegrationPoints(model.getDisplacement(), gradu_vect, dim,
                                  type, ghost_type, elem_filter);

  auto tangent_size = getTangentStiffnessVoigtSize(dim);

  auto tangent_stiffness_matrix = std::make_unique<Array<Real>>(
      nb_element * nb_quadrature_points, tangent_size * tangent_size,
      "tangent_stiffness_matrix");

  tangent_stiffness_matrix->zero();

  computeTangentModuli(type, *tangent_stiffness_matrix, ghost_type);

  /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  auto bt_d_b_size = dim * nb_nodes_per_element;

  auto bt_d_b = std::make_unique<Array<Real>>(
      nb_element * nb_quadrature_points, bt_d_b_size * bt_d_b_size, "B^t*D*B");

  fem.computeBtDB(*tangent_stiffness_matrix, *bt_d_b, 4, type, ghost_type,
                  elem_filter);

  /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  auto K_e = std::make_unique<Array<Real>>(nb_element,
                                           bt_d_b_size * bt_d_b_size, "K_e");

  fem.integrate(*bt_d_b, *K_e, bt_d_b_size * bt_d_b_size, type, ghost_type,
                elem_filter);

  model.getDOFManager().assembleElementalMatricesToMatrix(
      "K", "displacement", *K_e, type, ghost_type, _symmetric, elem_filter);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void Material::assembleStiffnessMatrixNL(ElementType type,
                                         GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const auto & shapes_derivatives = fem.getShapesDerivatives(type, ghost_type);

  auto & elem_filter = element_filter(type, ghost_type);
  // Array<Real> & gradu_vect = delta_gradu(type, ghost_type);

  auto nb_element = elem_filter.size();
  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

  auto shapes_derivatives_filtered = std::make_unique<Array<Real>>(
      nb_element * nb_quadrature_points, dim * nb_nodes_per_element,
      "shapes derivatives filtered");

  FEEngine::filterElementalData(fem.getMesh(), shapes_derivatives,
                                *shapes_derivatives_filtered, type, ghost_type,
                                elem_filter);

  /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  auto bt_s_b_size = dim * nb_nodes_per_element;

  auto bt_s_b = std::make_unique<Array<Real>>(
      nb_element * nb_quadrature_points, bt_s_b_size * bt_s_b_size, "B^t*D*B");

  auto piola_matrix_size = getCauchyStressMatrixSize(dim);

  Matrix<Real> B(piola_matrix_size, bt_s_b_size);
  Matrix<Real> S(piola_matrix_size, piola_matrix_size);

  for (auto && data :
       zip(make_view(*bt_s_b, bt_s_b_size, bt_s_b_size),
           make_view(piola_kirchhoff_2(type, ghost_type), dim, dim),
           make_view(*shapes_derivatives_filtered, spatial_dimension,
                     nb_nodes_per_element))) {
    auto && Bt_S_B = std::get<0>(data);
    auto && Piola_kirchhoff_matrix = std::get<1>(data);
    auto && shapes_derivatives = std::get<2>(data);

    setCauchyStressMatrix<dim>(Piola_kirchhoff_matrix, S);
    VoigtHelper<dim>::transferBMatrixToBNL(shapes_derivatives, B,
                                           nb_nodes_per_element);
    Bt_S_B = B.transpose() * S * B;
  }

  /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  auto K_e = std::make_unique<Array<Real>>(nb_element,
                                           bt_s_b_size * bt_s_b_size, "K_e");

  fem.integrate(*bt_s_b, *K_e, bt_s_b_size * bt_s_b_size, type, ghost_type,
                elem_filter);

  model.getDOFManager().assembleElementalMatricesToMatrix(
      "K", "displacement", *K_e, type, ghost_type, _symmetric, elem_filter);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void Material::assembleStiffnessMatrixL2(ElementType type,
                                         GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const auto & shapes_derivatives = fem.getShapesDerivatives(type, ghost_type);

  auto & elem_filter = element_filter(type, ghost_type);
  auto & gradu_vect = gradu(type, ghost_type);

  auto nb_element = elem_filter.size();
  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

  gradu_vect.resize(nb_quadrature_points * nb_element);

  fem.gradientOnIntegrationPoints(model.getDisplacement(), gradu_vect, dim,
                                  type, ghost_type, elem_filter);

  constexpr auto tangent_size = getTangentStiffnessVoigtSize(dim);

  auto tangent_stiffness_matrix = std::make_unique<Array<Real>>(
      nb_element * nb_quadrature_points, tangent_size * tangent_size,
      "tangent_stiffness_matrix");

  tangent_stiffness_matrix->zero();

  computeTangentModuli(type, *tangent_stiffness_matrix, ghost_type);

  auto shapes_derivatives_filtered = std::make_unique<Array<Real>>(
      nb_element * nb_quadrature_points, dim * nb_nodes_per_element,
      "shapes derivatives filtered");
  FEEngine::filterElementalData(fem.getMesh(), shapes_derivatives,
                                *shapes_derivatives_filtered, type, ghost_type,
                                elem_filter);

  /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  auto bt_d_b_size = dim * nb_nodes_per_element;
  auto bt_d_b = std::make_unique<Array<Real>>(
      nb_element * nb_quadrature_points, bt_d_b_size * bt_d_b_size, "B^t*D*B");

  Matrix<Real, tangent_size> B(tangent_size, bt_d_b_size);
  Matrix<Real, tangent_size> B2(tangent_size, bt_d_b_size);

  for (auto && data :
       zip(make_view(*bt_d_b, bt_d_b_size, bt_d_b_size),
           make_view<dim, dim>(gradu_vect),
           make_view<tangent_size, tangent_size>(*tangent_stiffness_matrix),
           make_view<dim, Eigen::Dynamic>(*shapes_derivatives_filtered, dim,
                                          nb_nodes_per_element))) {
    auto && Bt_D_B = std::get<0>(data);
    auto && grad_u = std::get<1>(data);
    auto && D = std::get<2>(data);
    auto && shapes_derivative = std::get<3>(data);

    // transferBMatrixToBL1<dim > (*shapes_derivatives_filtered_it, B,
    // nb_nodes_per_element);
    VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(shapes_derivative, B,
                                                       nb_nodes_per_element);
    VoigtHelper<dim>::transferBMatrixToBL2(shapes_derivative, grad_u, B2,
                                           nb_nodes_per_element);
    B += B2;
    Bt_D_B = B.transpose() * D * B;
  }

  /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  auto K_e = std::make_unique<Array<Real>>(nb_element,
                                           bt_d_b_size * bt_d_b_size, "K_e");

  fem.integrate(*bt_d_b, *K_e, bt_d_b_size * bt_d_b_size, type, ghost_type,
                elem_filter);

  model.getDOFManager().assembleElementalMatricesToMatrix(
      "K", "displacement", *K_e, type, ghost_type, _symmetric, elem_filter);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim> void Material::assembleInternalForces(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & internal_force = model.getInternalForce();
  auto & mesh = fem.getMesh();
  for (auto type : element_filter.elementTypes(_ghost_type = ghost_type)) {
    const auto & shapes_derivatives =
        fem.getShapesDerivatives(type, ghost_type);

    auto & elem_filter = element_filter(type, ghost_type);
    if (elem_filter.size() == 0)
      continue;
    auto size_of_shapes_derivatives = shapes_derivatives.getNbComponent();
    auto nb_element = elem_filter.size();
    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

    auto shapesd_filtered = std::make_unique<Array<Real>>(
        nb_element, size_of_shapes_derivatives, "filtered shapesd");

    FEEngine::filterElementalData(mesh, shapes_derivatives, *shapesd_filtered,
                                  type, ghost_type, elem_filter);

    // Set stress vectors
    auto stress_size = getTangentStiffnessVoigtSize(dim);
    auto bt_s_size = dim * nb_nodes_per_element;

    auto bt_s = std::make_unique<Array<Real>>(nb_element * nb_quadrature_points,
                                              bt_s_size, "B^t*S");

    Matrix<Real> B_tensor(stress_size, bt_s_size);
    Matrix<Real> B2_tensor(stress_size, bt_s_size);

    for (auto && data :
         zip(make_view<dim, dim>(this->gradu(type, ghost_type)),
             make_view<dim, dim>(this->piola_kirchhoff_2(type, ghost_type)),
             make_view(*bt_s, bt_s_size),
             make_view<dim, Eigen::Dynamic>(*shapesd_filtered, dim,
                                            nb_nodes_per_element))) {
      auto && grad_u = std::get<0>(data);
      auto && S = std::get<1>(data);
      auto && r = std::get<2>(data);
      auto && shapes_derivative = std::get<3>(data);

      VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(
          shapes_derivative, B_tensor, nb_nodes_per_element);

      VoigtHelper<dim>::transferBMatrixToBL2(shapes_derivative, grad_u,
                                             B2_tensor, nb_nodes_per_element);

      B_tensor += B2_tensor;
      r = B_tensor.transpose() * Material::stressToVoigt<dim>(S);
    }

    /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
    auto r_e = std::make_unique<Array<Real>>(nb_element, bt_s_size, "r_e");
    fem.integrate(*bt_s, *r_e, bt_s_size, type, ghost_type, elem_filter);

    model.getDOFManager().assembleElementalArrayLocalArray(
        *r_e, internal_force, type, ghost_type, -1., elem_filter);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::computePotentialEnergyByElements() {
  AKANTU_DEBUG_IN();

  for (auto type : element_filter.elementTypes(spatial_dimension, _not_ghost)) {
    computePotentialEnergy(type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::computePotentialEnergy(ElementType) { AKANTU_TO_IMPLEMENT(); }

/* -------------------------------------------------------------------------- */
Real Material::getPotentialEnergy() {
  AKANTU_DEBUG_IN();
  Real epot = 0.;

  computePotentialEnergyByElements();

  /// integrate the potential energy for each type of elements
  for (auto type : element_filter.elementTypes(spatial_dimension, _not_ghost)) {
    epot += fem.integrate(potential_energy(type, _not_ghost), type, _not_ghost,
                          element_filter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return epot;
}

/* -------------------------------------------------------------------------- */
Real Material::getPotentialEnergy(ElementType type, Int index) {
  return getPotentialEnergy({type, index, _not_ghost});
}

/* -------------------------------------------------------------------------- */
Real Material::getPotentialEnergy(const Element & element) {
  AKANTU_DEBUG_IN();
  Real epot = 0.;

  Vector<Real> epot_on_quad_points(fem.getNbIntegrationPoints(element.type));
  computePotentialEnergyByElement(element, epot_on_quad_points);

  epot = fem.integrate(epot_on_quad_points,
                       {element.type,
                        element_filter(element.type)(element.element),
                        _not_ghost});

  AKANTU_DEBUG_OUT();
  return epot;
}

/* -------------------------------------------------------------------------- */
Real Material::getEnergy(const std::string & type) {
  AKANTU_DEBUG_IN();
  if (type == "potential") {
    return getPotentialEnergy();
  }
  AKANTU_DEBUG_OUT();
  return 0.;
}

/* -------------------------------------------------------------------------- */
Real Material::getEnergy(const std::string & energy_id,
                         const Element & element) {
  AKANTU_DEBUG_IN();
  if (energy_id == "potential") {
    return getPotentialEnergy(element);
  }
  AKANTU_DEBUG_OUT();
  return 0.;
}

/* -------------------------------------------------------------------------- */
void Material::initElementalFieldInterpolation(
    const ElementTypeMapArray<Real> & interpolation_points_coordinates) {
  AKANTU_DEBUG_IN();

  this->fem.initElementalFieldInterpolationFromIntegrationPoints(
      interpolation_points_coordinates, this->interpolation_points_matrices,
      this->interpolation_inverse_coordinates, &(this->element_filter));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::interpolateStress(ElementTypeMapArray<Real> & result,
                                 const GhostType ghost_type) {

  this->fem.interpolateElementalFieldFromIntegrationPoints(
      this->stress, this->interpolation_points_matrices,
      this->interpolation_inverse_coordinates, result, ghost_type,
      &(this->element_filter));
}

/* -------------------------------------------------------------------------- */
void Material::interpolateStressOnFacets(
    ElementTypeMapArray<Real> & result,
    ElementTypeMapArray<Real> & by_elem_result, const GhostType ghost_type) {

  interpolateStress(by_elem_result, ghost_type);

  auto stress_size = this->stress.getNbComponent();

  const auto & mesh = this->model.getMesh();
  const auto & mesh_facets = mesh.getMeshFacets();

  for (auto type : element_filter.elementTypes(spatial_dimension, ghost_type)) {
    auto & elem_fil = element_filter(type, ghost_type);
    auto & by_elem_res = by_elem_result(type, ghost_type);
    auto nb_element_full = mesh.getNbElement(type, ghost_type);
    auto nb_interpolation_points_per_elem =
        by_elem_res.size() / nb_element_full;

    const auto & facet_to_element =
        mesh_facets.getSubelementToElement(type, ghost_type);
    auto type_facet = Mesh::getFacetType(type);
    auto nb_facet_per_elem = facet_to_element.getNbComponent();
    auto nb_quad_per_facet =
        nb_interpolation_points_per_elem / nb_facet_per_elem;
    Element element_for_comparison{type, 0, ghost_type};
    const Array<std::vector<Element>> * element_to_facet = nullptr;
    auto current_ghost_type = _not_ghost;
    auto result_vec_it =
        make_view(result(type_facet, current_ghost_type), stress_size, 2)
            .begin();

    auto result_it =
        make_view(by_elem_res, stress_size, nb_interpolation_points_per_elem)
            .begin();

    for (auto global_el : elem_fil) {
      element_for_comparison.element = global_el;

      for (Int f = 0; f < nb_facet_per_elem; ++f) {
        auto facet_elem = facet_to_element(global_el, f);
        auto global_facet = facet_elem.element;

        if (facet_elem.ghost_type != current_ghost_type) {
          current_ghost_type = facet_elem.ghost_type;
          element_to_facet = &mesh_facets.getElementToSubelement(
              type_facet, current_ghost_type);
          result_vec_it =
              make_view(result(type_facet, current_ghost_type), stress_size, 2)
                  .begin();
        }

        auto is_second_element =
            ((*element_to_facet)(global_facet)[0] != element_for_comparison);

        for (Int q = 0; q < nb_quad_per_facet; ++q) {
          auto && result_local =
              result_vec_it[global_facet * nb_quad_per_facet + q](
                  is_second_element);
          result_local = result_it[global_el](f * nb_quad_per_facet + q);
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void Material::addElements(const Array<Element> & elements_to_add) {
  AKANTU_DEBUG_IN();

  UInt mat_id = model.getMaterialIndex(name);
  for (const auto & element : elements_to_add) {
    auto index = this->addElement(element);
    model.material_index(element) = mat_id;
    model.material_local_numbering(element) = index;
  }

  this->resizeInternals();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::removeElements(const Array<Element> & elements_to_remove) {
  AKANTU_DEBUG_IN();

  auto el_begin = elements_to_remove.begin();
  auto el_end = elements_to_remove.end();

  if (elements_to_remove.empty()) {
    return;
  }

  auto & mesh = this->model.getMesh();

  ElementTypeMapArray<Idx> material_local_new_numbering(
      "remove mat filter elem", id);

  material_local_new_numbering.initialize(
      mesh, _element_filter = &element_filter, _element_kind = _ek_not_defined,
      _with_nb_element = true);

  ElementTypeMapArray<Idx> element_filter_tmp("element_filter_tmp", id);

  element_filter_tmp.initialize(mesh, _element_filter = &element_filter,
                                _element_kind = _ek_not_defined);

  ElementTypeMap<Idx> new_ids, element_ids;

  for_each_element(
      mesh,
      [&](auto && el) {
        if (not new_ids(el.type, el.ghost_type)) {
          element_ids(el.type, el.ghost_type) = 0;
        }

        auto & element_id = element_ids(el.type, el.ghost_type);
        auto l_el = Element{el.type, element_id, el.ghost_type};
        if (std::find(el_begin, el_end, el) != el_end) {
          material_local_new_numbering(l_el) = UInt(-1);
          return;
        }

        element_filter_tmp(el.type, el.ghost_type).push_back(el.element);
        if (not new_ids(el.type, el.ghost_type)) {
          new_ids(el.type, el.ghost_type) = 0;
        }

        auto & new_id = new_ids(el.type, el.ghost_type);

        material_local_new_numbering(l_el) = new_id;
        model.material_local_numbering(el) = new_id;

        ++new_id;
        ++element_id;
      },
      _element_filter = &element_filter, _element_kind = _ek_not_defined);

  for (auto ghost_type : ghost_types) {
    for (const auto & type : element_filter.elementTypes(
             _ghost_type = ghost_type, _element_kind = _ek_not_defined)) {
      element_filter(type, ghost_type)
          .copy(element_filter_tmp(type, ghost_type));
    }
  }

  for (auto it = internal_vectors_real.begin();
       it != internal_vectors_real.end(); ++it) {
    it->second->removeIntegrationPoints(material_local_new_numbering);
  }

  for (auto it = internal_vectors_int.begin(); it != internal_vectors_int.end();
       ++it) {
    it->second->removeIntegrationPoints(material_local_new_numbering);
  }

  for (auto it = internal_vectors_bool.begin();
       it != internal_vectors_bool.end(); ++it) {
    it->second->removeIntegrationPoints(material_local_new_numbering);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::resizeInternals() {
  AKANTU_DEBUG_IN();
  for (auto it = internal_vectors_real.begin();
       it != internal_vectors_real.end(); ++it) {
    it->second->resize();
  }

  for (auto it = internal_vectors_int.begin(); it != internal_vectors_int.end();
       ++it) {
    it->second->resize();
  }

  for (auto it = internal_vectors_bool.begin();
       it != internal_vectors_bool.end(); ++it) {
    it->second->resize();
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::onElementsAdded(const Array<Element> & /*unused*/,
                               const NewElementsEvent & /*unused*/) {
  this->resizeInternals();
}

/* -------------------------------------------------------------------------- */
void Material::onElementsRemoved(
    const Array<Element> & element_list,
    const ElementTypeMapArray<Idx> & new_numbering,
    [[gnu::unused]] const RemovedElementsEvent & event) {
  auto my_num = model.getInternalIndexFromID(getID());

  ElementTypeMapArray<Idx> material_local_new_numbering(
      "remove mat filter elem", getID());

  auto el_begin = element_list.begin();
  auto el_end = element_list.end();

  for (auto && gt : ghost_types) {
    for (auto && type :
         new_numbering.elementTypes(_all_dimensions, gt, _ek_not_defined)) {

      if (not element_filter.exists(type, gt) ||
          element_filter(type, gt).empty()) {
        continue;
      }

      auto & elem_filter = element_filter(type, gt);
      auto & mat_indexes = this->model.material_index(type, gt);
      auto & mat_loc_num = this->model.material_local_numbering(type, gt);
      auto nb_element = this->model.getMesh().getNbElement(type, gt);

      // all materials will resize of the same size...
      mat_indexes.resize(nb_element);
      mat_loc_num.resize(nb_element);

      if (!material_local_new_numbering.exists(type, gt)) {
        material_local_new_numbering.alloc(elem_filter.size(), 1, type, gt);
      }

      auto & mat_renumbering = material_local_new_numbering(type, gt);
      const auto & renumbering = new_numbering(type, gt);
      Array<Idx> elem_filter_tmp;
      UInt ni = 0;
      Element el{type, 0, gt};

      for (auto && data : enumerate(elem_filter)) {
        el.element = std::get<1>(data);
        if (std::find(el_begin, el_end, el) == el_end) {
          auto new_el = renumbering(el.element);
          AKANTU_DEBUG_ASSERT(new_el != -1,
                              "A not removed element as been badly renumbered");
          elem_filter_tmp.push_back(new_el);
          mat_renumbering(std::get<0>(data)) = ni;

          mat_indexes(new_el) = my_num;
          mat_loc_num(new_el) = ni;
          ++ni;
        } else {
          mat_renumbering(std::get<0>(data)) = -1;
        }
      }

      elem_filter.resize(elem_filter_tmp.size());
      elem_filter.copy(elem_filter_tmp);
    }
  }

  for (auto it = internal_vectors_real.begin();
       it != internal_vectors_real.end(); ++it) {
    it->second->removeIntegrationPoints(material_local_new_numbering);
  }

  for (auto it = internal_vectors_int.begin(); it != internal_vectors_int.end();
       ++it) {
    it->second->removeIntegrationPoints(material_local_new_numbering);
  }

  for (auto it = internal_vectors_bool.begin();
       it != internal_vectors_bool.end(); ++it) {
    it->second->removeIntegrationPoints(material_local_new_numbering);
  }
}

/* -------------------------------------------------------------------------- */
void Material::beforeSolveStep() { this->savePreviousState(); }

/* -------------------------------------------------------------------------- */
void Material::afterSolveStep(bool converged) {
  if (not converged) {
    this->restorePreviousState();
    return;
  }

  for (const auto & type : element_filter.elementTypes(
           _all_dimensions, _not_ghost, _ek_not_defined)) {
    this->updateEnergies(type);
  }
}
/* -------------------------------------------------------------------------- */
void Material::onDamageIteration() { this->savePreviousState(); }

/* -------------------------------------------------------------------------- */
void Material::onDamageUpdate() {
  for (const auto & type : element_filter.elementTypes(
           _all_dimensions, _not_ghost, _ek_not_defined)) {
    this->updateEnergiesAfterDamage(type);
  }
}

/* -------------------------------------------------------------------------- */
void Material::onDump() {
  if (this->isFiniteDeformation()) {
    this->computeAllCauchyStresses(_not_ghost);
  }
}

/* -------------------------------------------------------------------------- */
void Material::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);
  std::string type = getID().substr(getID().find_last_of(':') + 1);

  stream << space << "Material " << type << " [" << std::endl;
  Parsable::printself(stream, indent);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
/// extrapolate internal values
void Material::extrapolateInternal(const ID & id, const Element & element,
                                   const Matrix<Real> & /*point*/,
                                   Matrix<Real> & extrapolated) {
  if (this->isInternal<Real>(id, element.kind())) {
    auto name = this->getID() + ":" + id;
    auto nb_quads =
        this->internal_vectors_real[name]->getFEEngine().getNbIntegrationPoints(
            element.type, element.ghost_type);
    const auto & internal =
        this->getArray<Real>(id, element.type, element.ghost_type);
    auto nb_component = internal.getNbComponent();
    auto internal_it = make_view(internal, nb_component, nb_quads).begin();
    auto local_element = this->convertToLocalElement(element);

    /// instead of really extrapolating, here the value of the first GP
    /// is copied into the result vector. This works only for linear
    /// elements
    /// @todo extrapolate!!!!
    AKANTU_DEBUG_WARNING("This is a fix, values are not truly extrapolated");

    auto && values = internal_it[local_element.element];
    Int index = 0;
    Vector<Real> tmp(nb_component);
    for (Int j = 0; j < values.cols(); ++j) {
      tmp = values(j);
      if (tmp.norm() > 0) {
        index = j;
        break;
      }
    }

    for (Int i = 0; i < extrapolated.size(); ++i) {
      extrapolated(i) = values(index);
    }
  } else {
    extrapolated.fill(0.);
  }
}

/* -------------------------------------------------------------------------- */
void Material::applyEigenGradU(const Matrix<Real> & prescribed_eigen_grad_u,
                               const GhostType ghost_type) {

  for (auto && type : element_filter.elementTypes(_all_dimensions, _not_ghost,
                                                  _ek_not_defined)) {
    if (element_filter(type, ghost_type).empty()) {
      continue;
    }

    for (auto & eigengradu : make_view(this->eigengradu(type, ghost_type),
                                       spatial_dimension, spatial_dimension)) {
      eigengradu = prescribed_eigen_grad_u;
    }
  }
}

/* -------------------------------------------------------------------------- */
MaterialFactory & Material::getFactory() {
  return MaterialFactory::getInstance();
}

} // namespace akantu
