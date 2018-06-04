/**
 * @file   material_FE2.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @brief Material for multi-scale simulations. It stores an
 * underlying RVE on each integration point of the material.
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_FE2.hh"
#include "communicator.hh"
#include "solid_mechanics_model_RVE.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialFE2<spatial_dimension>::MaterialFE2(SolidMechanicsModel & model,
                                            const ID & id)
    : Parent(model, id), C("material_stiffness", *this),
      gelstrain("gelstrain", *this), non_reacted_gel("non_reacted_gel", *this) {
  AKANTU_DEBUG_IN();

  this->C.initialize(voigt_h::size * voigt_h::size);
  this->gelstrain.initialize(spatial_dimension * spatial_dimension);
  this->non_reacted_gel.initialize(1);
  this->initialize();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialFE2<spatial_dimension>::~MaterialFE2() = default;

/* -------------------------------------------------------------------------- */
template <UInt dim> void MaterialFE2<dim>::initialize() {
  this->registerParam("element_type", el_type, _triangle_3,
                      _pat_parsable | _pat_modifiable,
                      "element type in RVE mesh");
  this->registerParam("mesh_file", mesh_file, _pat_parsable | _pat_modifiable,
                      "the mesh file for the RVE");
  this->registerParam("nb_gel_pockets", nb_gel_pockets,
                      _pat_parsable | _pat_modifiable,
                      "the number of gel pockets in each RVE");
  this->registerParam("k", k, _pat_parsable | _pat_modifiable,
                      "pre-exponential factor of Arrhenius law");
  this->registerParam("activ_energy", activ_energy,
                      _pat_parsable | _pat_modifiable,
                      "activation energy of ASR in Arrhenius law");
  this->registerParam("R", R, _pat_parsable | _pat_modifiable,
                      "universal gas constant R in Arrhenius law");
  this->registerParam("saturation_constant", sat_const, Real(0.0),
                      _pat_parsable | _pat_modifiable, "saturation constant");
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  Parent::initMaterial();

  /// create a Mesh and SolidMechanicsModel on each integration point of the
  /// material
  auto C_it = this->C(this->el_type).begin(voigt_h::size, voigt_h::size);
  const auto & comm = this->model.getMesh().getCommunicator();
  UInt prank = comm.whoAmI();
  UInt C_size = this->C(this->el_type).size();

  for (auto && data :
       enumerate(make_view(C(this->el_type), voigt_h::size, voigt_h::size))) {
    auto q = std::get<0>(data);
    auto & C = std::get<1>(data);

    meshes.emplace_back(std::make_unique<Mesh>(
        spatial_dimension, "RVE_mesh_" + std::to_string(prank * C_size + q),
        q + 1));

    auto & mesh = *meshes.back();
    mesh.read(mesh_file);

    RVEs.emplace_back(std::make_unique<SolidMechanicsModelRVE>(
        mesh, true, this->nb_gel_pockets, _all_dimensions,
        "SMM_RVE_" + std::to_string(prank * C_size + q), q + 1));

    auto & RVE = *RVEs.back();
    RVE.initFull(_analysis_method = _static);

    /// compute intial stiffness of the RVE
    RVE.homogenizeStiffness(C);
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::computeStress(ElementType el_type,
                                                   GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  // Compute thermal stresses first

  Parent::computeStress(el_type, ghost_type);
  Array<Real>::const_scalar_iterator sigma_th_it =
      this->sigma_th(el_type, ghost_type).begin();

  // Wikipedia convention:
  // 2*eps_ij (i!=j) = voigt_eps_I
  // http://en.wikipedia.org/wiki/Voigt_notation

  Array<Real>::const_matrix_iterator C_it =
      this->C(el_type, ghost_type).begin(voigt_h::size, voigt_h::size);

  // create vectors to store stress and strain in Voigt notation
  // for efficient computation of stress
  Vector<Real> voigt_strain(voigt_h::size);
  Vector<Real> voigt_stress(voigt_h::size);

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  const Matrix<Real> & C_mat = *C_it;
  const Real & sigma_th = *sigma_th_it;

  /// copy strains in Voigt notation
  for (UInt I = 0; I < voigt_h::size; ++I) {
    /// copy stress in
    Real voigt_factor = voigt_h::factors[I];
    UInt i = voigt_h::vec[I][0];
    UInt j = voigt_h::vec[I][1];

    voigt_strain(I) = voigt_factor * (grad_u(i, j) + grad_u(j, i)) / 2.;
  }

  // compute stresses in Voigt notation
  voigt_stress.mul<false>(C_mat, voigt_strain);

  /// copy stresses back in full vectorised notation
  for (UInt I = 0; I < voigt_h::size; ++I) {
    UInt i = voigt_h::vec[I][0];
    UInt j = voigt_h::vec[I][1];
    sigma(i, j) = sigma(j, i) = voigt_stress(I) + (i == j) * sigma_th;
  }

  ++C_it;
  ++sigma_th_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::computeTangentModuli(
    const ElementType & el_type, Array<Real> & tangent_matrix,
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Array<Real>::const_matrix_iterator C_it =
      this->C(el_type, ghost_type).begin(voigt_h::size, voigt_h::size);

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
  tangent.copy(*C_it);
  ++C_it;
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::advanceASR(
    const Matrix<Real> & prestrain) {
  AKANTU_DEBUG_IN();

  for (auto && data :
       zip(RVEs, make_view(this->gradu(this->el_type), spatial_dimension,
                           spatial_dimension),
           make_view(this->eigengradu(this->el_type), spatial_dimension,
                     spatial_dimension),
           make_view(this->C(this->el_type), voigt_h::size, voigt_h::size),
           this->delta_T(this->el_type))) {
    auto & RVE = *(std::get<0>(data));

    /// apply boundary conditions based on the current macroscopic displ.
    /// gradient
    RVE.applyBoundaryConditions(std::get<1>(data));

    /// apply homogeneous temperature field to each RVE to obtain thermoelastic
    /// effect
    RVE.applyHomogeneousTemperature(std::get<4>(data));

    /// advance the ASR in every RVE
    RVE.advanceASR(prestrain);

    /// remove temperature field - not to mess up with the stiffness
    /// homogenization
    /// further
    RVE.removeTemperature();

    /// compute the average eigen_grad_u
    RVE.homogenizeEigenGradU(std::get<2>(data));

    /// compute the new effective stiffness of the RVE
    RVE.homogenizeStiffness(std::get<3>(data));
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::advanceASR(const Real & delta_time) {
  AKANTU_DEBUG_IN();

  for (auto && data :
       zip(RVEs, make_view(this->gradu(this->el_type), spatial_dimension,
                           spatial_dimension),
           make_view(this->eigengradu(this->el_type), spatial_dimension,
                     spatial_dimension),
           make_view(this->C(this->el_type), voigt_h::size, voigt_h::size),
           this->delta_T(this->el_type),
           make_view(this->gelstrain(this->el_type), spatial_dimension,
                     spatial_dimension))) {
    auto & RVE = *(std::get<0>(data));

    /// apply boundary conditions based on the current macroscopic displ.
    /// gradient
    RVE.applyBoundaryConditions(std::get<1>(data));

    /// apply homogeneous temperature field to each RVE to obtain thermoelastic
    /// effect
    RVE.applyHomogeneousTemperature(std::get<4>(data));

    /// compute new gel strain for every element
    if (this->sat_const)
      this->computeNewGelStrainTimeDependent(delta_time);
    else
      this->computeNewGelStrain(delta_time);

    /// advance the ASR in every RVE based on the new gel strain
    RVE.advanceASR(std::get<5>(data));

    /// remove temperature field - not to mess up with the stiffness
    /// homogenization
    /// further
    RVE.removeTemperature();

    /// compute the average eigen_grad_u
    RVE.homogenizeEigenGradU(std::get<2>(data));

    /// compute the new effective stiffness of the RVE
    RVE.homogenizeStiffness(std::get<3>(data));
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::computeNewGelStrain(
    const Real & delta_time) {
  AKANTU_DEBUG_IN();

  for (auto && data : zip(RVEs, make_view(this->gelstrain(this->el_type),
                                          spatial_dimension, spatial_dimension),
                          this->delta_T(this->el_type))) {

    auto & gelstrain = (std::get<1>(data));
    const auto & T = (std::get<2>(data));

    const auto & k = this->k;
    const auto & Ea = this->activ_energy;
    const auto & R = this->R;

    /// compute increase in gel strain value for interval of time delta_time
    /// as temperatures are stored in C, conversion to K is done
    Real delta_strain = k * std::exp(-Ea / (R * (T + 273.15))) * delta_time;

    for (UInt i = 0; i != spatial_dimension; ++i)
      gelstrain(i, i) += delta_strain;
  }
  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::computeNewGelStrainTimeDependent(
    const Real & delta_time) {
  AKANTU_DEBUG_IN();

  for (auto && data : zip(RVEs, make_view(this->gelstrain(this->el_type),
                                          spatial_dimension, spatial_dimension),
                          this->delta_T(this->el_type),
                          this->non_reacted_gel(this->el_type))) {

    auto & gelstrain = (std::get<1>(data));
    const auto & T = (std::get<2>(data));
    auto & non_reacted_gel = (std::get<3>(data));

    const auto & k = this->k;
    const auto & Ea = this->activ_energy;
    const auto & R = this->R;
    const auto & sat_const = this->sat_const;

    /// compute increase in gel strain value for interval of time delta_time
    /// as temperatures are stored in C, conversion to K is done
    Real delta_strain =
        non_reacted_gel * k * std::exp(-Ea / (R * (T + 273.15))) * delta_time;

    for (UInt i = 0; i != spatial_dimension; ++i)
      gelstrain(i, i) += delta_strain;

    non_reacted_gel -= std::exp(-Ea / (R * (T + 273.15))) *
                       (delta_time / 60 / 60 / 24) / sat_const;

    if (non_reacted_gel < 0.)
      non_reacted_gel = 0.;
  }

  AKANTU_DEBUG_OUT();
}

INSTANTIATE_MATERIAL(material_FE2, MaterialFE2);

} // namespace akantu
