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

/* ------------------------------------------------------------------- */
#include "material_FE2.hh"
#include "aka_iterators.hh"
#include "communicator.hh"
#include "solid_mechanics_model_RVE.hh"
#include <fstream>
/* ------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialFE2<spatial_dimension>::MaterialFE2(SolidMechanicsModel & model,
                                            const ID & id)
    : Parent(model, id), C("material_stiffness", *this),
      gelstrain("gelstrain",
                *this), /*non_reacted_gel("non_reacted_gel", *this),*/
      crack_volume_ratio("crack_volume_ratio", *this),
      crack_volume_ratio_paste("crack_volume_ratio_paste", *this),
      crack_volume_ratio_agg("crack_volume_ratio_agg", *this) {
  AKANTU_DEBUG_IN();

  this->C.initialize(voigt_h::size * voigt_h::size);
  this->gelstrain.initialize(spatial_dimension * spatial_dimension);
  // this->non_reacted_gel.initialize(1);
  // this->non_reacted_gel.setDefaultValue(1.0);
  this->crack_volume_ratio.initialize(1);
  this->crack_volume_ratio_paste.initialize(1);
  this->crack_volume_ratio_agg.initialize(1);
  this->initialize();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialFE2<spatial_dimension>::~MaterialFE2() = default;

/* -------------------------------------------------------------------------- */
template <UInt dim> void MaterialFE2<dim>::initialize() {
  this->registerParam("element_type", el_type, _triangle_3, _pat_parsmod,
                      "element type in RVE mesh");
  this->registerParam("mesh_file", mesh_file, _pat_parsmod,
                      "the mesh file for the RVE");
  this->registerParam("nb_gel_pockets", nb_gel_pockets, _pat_parsmod,
                      "the number of gel pockets in each RVE");
  // this->registerParam("k", k, _pat_parsable | _pat_modifiable,
  //                     "pre-exponential factor of Arrhenius law");
  // this->registerParam("activ_energy", activ_energy,
  //                     _pat_parsable | _pat_modifiable,
  //                     "activation energy of ASR in Arrhenius law");
  // this->registerParam("R", R, _pat_parsable | _pat_modifiable,
  //                     "universal gas constant R in Arrhenius law");
  // this->registerParam("saturation_constant", sat_const, Real(0.0),
  //                     _pat_parsable | _pat_modifiable, "saturation
  //                     constant");
  this->registerParam("eps_inf", eps_inf, _pat_parsmod,
                      "asymptotic value of ASR expansion");
  this->registerParam("U_C", U_C, _pat_parsmod, "thermal activation energy C");
  this->registerParam("U_L", U_L, _pat_parsmod, "thermal activation energy L");
  this->registerParam("T_ref", T_ref, _pat_parsmod, "reference temperature");
  this->registerParam("time_lat_ref", time_lat_ref, _pat_parsmod,
                      "latency time at the reference temperature");
  this->registerParam("time_ch_ref", time_ch_ref, _pat_parsmod,
                      "characteristic time at the reference temperature");
  this->registerParam("reset_damage", reset_damage, false, _pat_parsmod,
                      "reset damage");
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  Parent::initMaterial();

  /// create a Mesh and SolidMechanicsModel on each integration point of the
  /// material
  auto & mesh = this->model.getMesh();
  auto & g_ids = mesh.template getData<UInt>("global_ids", this->el_type);
  // AKANTU_DEBUG_ASSERT(g_ids.size() > 0, "Global numbering array is empty");
  auto const & element_filter = this->getElementFilter()(this->el_type);

  for (auto && data :
       zip(arange(element_filter.size()), element_filter,
           make_view(C(this->el_type), voigt_h::size, voigt_h::size))) {
    UInt mat_el_id = std::get<0>(data);
    UInt proc_el_id = std::get<1>(data);
    UInt gl_el_id = g_ids(proc_el_id);
    auto & C = std::get<2>(data);

    meshes.emplace_back(std::make_unique<Mesh>(
        spatial_dimension, "RVE_mesh_" + std::to_string(gl_el_id),
        mat_el_id + 1));

    auto & mesh = *meshes.back();
    mesh.read(mesh_file);

    RVEs.emplace_back(std::make_unique<SolidMechanicsModelRVE>(
        mesh, true, this->nb_gel_pockets, _all_dimensions,
        "SMM_RVE_" + std::to_string(gl_el_id), mat_el_id + 1));

    auto & RVE = *RVEs.back();
    RVE.initFull(_analysis_method = _static);

    /// compute intial stiffness of the RVE
    RVE.homogenizeStiffness(C, false);
  }
  AKANTU_DEBUG_OUT();
}
/* ------------------------------------------------------------------*/
template <UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::computeStress(ElementType el_type,
                                                   GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  // Compute thermal stresses first
  Parent::computeStress(el_type, ghost_type);

  for (auto && data :
       zip(RVEs,
           make_view(this->gradu(el_type), spatial_dimension,
                     spatial_dimension),
           make_view(this->stress(el_type), spatial_dimension,
                     spatial_dimension),
           this->sigma_th(el_type),
           make_view(this->C(this->el_type), voigt_h::size, voigt_h::size),
           make_view(this->gelstrain(this->el_type), spatial_dimension,
                     spatial_dimension),
           this->delta_T(this->el_type))) {
    auto & RVE = *(std::get<0>(data));

    /// reset nodal and internal fields
    RVE.resetNodalFields();
    RVE.resetInternalFields();

    /// reset the damage to the previously converged state
    if (reset_damage)
      RVE.restoreDamageField();

    /// apply boundary conditions based on the current macroscopic displ.
    /// gradient
    RVE.applyBoundaryConditionsRve(std::get<1>(data));

    /// advance the ASR in every RVE based on the new gel strain
    RVE.advanceASR(std::get<5>(data));

    /// compute the average average rVE stress
    RVE.homogenizeStressField(std::get<2>(data));

    // /// compute the new effective stiffness of the RVE
    //   /// decide whether stiffness homogenization is done via tension
    //   RVE.setStiffHomogenDir(std::get<2>(data));
    //   /// compute the new effective stiffness of the RVE
    //   RVE.homogenizeStiffness(std::get<4>(data), RVE.isTensileHomogen());
  }
  AKANTU_DEBUG_OUT();
}

/* ------------------------------------------------------------------*/
template <UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::increaseGelStrain(Real & dt_day) {
  for (auto && data : zip(RVEs, this->delta_T(this->el_type),
                          make_view(this->gelstrain(this->el_type),
                                    spatial_dimension, spatial_dimension))) {
    auto & RVE = *(std::get<0>(data));
    auto & strain_matrix = std::get<2>(data);
    Real ASRStrain = strain_matrix(0, 0);
    RVE.computeASRStrainLarive(
        dt_day, std::get<1>(data), ASRStrain, this->eps_inf, this->time_ch_ref,
        this->time_lat_ref, this->U_C, this->U_L, this->T_ref);

    for (UInt i = 0; i != spatial_dimension; ++i)
      strain_matrix(i, i) = ASRStrain;
  }
}
/* ------------------------------------------------------------------*/
template <UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::increaseGelStrainArrhenius(
    Real & dt_day, const Real & k, const Real & Ea) {
  for (auto && data : zip(RVEs, this->delta_T(this->el_type),
                          make_view(this->gelstrain(this->el_type),
                                    spatial_dimension, spatial_dimension))) {
    auto & RVE = *(std::get<0>(data));
    auto & strain_matrix = std::get<2>(data);
    Real ASRStrain = strain_matrix(0, 0);
    RVE.computeASRStrainArrhenius(dt_day, std::get<1>(data), ASRStrain, k, Ea);

    for (UInt i = 0; i != spatial_dimension; ++i)
      strain_matrix(i, i) = ASRStrain;
  }
}
/* ------------------------------------------------------------------*/
template <UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::updateStiffness() {
  AKANTU_DEBUG_IN();
  for (const auto & type :
       this->element_filter.elementTypes(spatial_dimension, _not_ghost)) {
    Array<UInt> & elem_filter = this->element_filter(type, _not_ghost);

    if (elem_filter.size() == 0)
      return;
  }
  for (auto && data :
       zip(RVEs,
           make_view(this->C(this->el_type), voigt_h::size, voigt_h::size),
           make_view(this->stress(el_type), spatial_dimension,
                     spatial_dimension),
           this->delta_T(this->el_type),
           this->crack_volume_ratio(this->el_type),
           this->crack_volume_ratio_paste(this->el_type),
           this->crack_volume_ratio_agg(this->el_type))) {
    auto & RVE = *(std::get<0>(data));
    if (reset_damage)
      RVE.storeDamageField();

    /// decide whether stiffness homogenization is done via tension
    RVE.setStiffHomogenDir(std::get<2>(data));
    /// compute the new effective stiffness of the RVE
    RVE.homogenizeStiffness(std::get<1>(data), RVE.isTensileHomogen());

    /// compute crack volume ratio in each RVE
    RVE.computeCrackVolume(std::get<4>(data));
    /// compute crack volume ratio per material
    RVE.computeCrackVolumePerMaterial(std::get<5>(data), "paste");
    RVE.computeCrackVolumePerMaterial(std::get<6>(data), "aggregate");
  }
  AKANTU_DEBUG_OUT();
}

/* ---------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::computeTangentModuli(
    ElementType el_type, Array<Real> & tangent_matrix, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Array<Real>::const_matrix_iterator C_it =
      this->C(el_type, ghost_type).begin(voigt_h::size, voigt_h::size);

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
  tangent.copy(*C_it);
  ++C_it;
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* ----------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialFE2<spatial_dimension>::computeAverageGelStrain() {
  AKANTU_DEBUG_IN();

  Real av_gelstrain = 0;
  UInt nb_RVEs = 0;

  for (auto && data :
       enumerate(make_view(this->gelstrain(this->el_type), spatial_dimension,
                           spatial_dimension))) {
    av_gelstrain += std::get<1>(data)(0, 0);
    nb_RVEs = std::get<0>(data) + 1;
  }
  auto && comm = akantu::Communicator::getWorldCommunicator();
  comm.allReduce(av_gelstrain, SynchronizerOperation::_sum);
  comm.allReduce(nb_RVEs, SynchronizerOperation::_sum);

  return av_gelstrain /= nb_RVEs;

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::setDirectoryToRveDumper(
    const std::string & directory) {
  AKANTU_DEBUG_IN();

  for (auto && data : RVEs) {
    /// set default directory to all RVEs
    data->setDirectory(directory);
  }

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension>
UInt MaterialFE2<spatial_dimension>::getNbRVEs() {
  AKANTU_DEBUG_IN();
  UInt nb_RVEs = 0;
  for (auto && data : enumerate(RVEs)) {
    nb_RVEs = std::get<0>(data) + 1;
  }
  return nb_RVEs;
  AKANTU_DEBUG_OUT();
}

/* ------------------------------------------------------------------ */
template <UInt spatial_dimension> void MaterialFE2<spatial_dimension>::dump() {
  AKANTU_DEBUG_IN();

  for (auto && RVE : RVEs) {
    /// update stress field before dumping
    RVE->assembleInternalForces();
    /// dump all the RVEs
    RVE->dumpRve();
  }

  AKANTU_DEBUG_OUT();
}
/* ------------------------------------------------------------------ */
template <UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::dump(UInt dump_nb) {
  AKANTU_DEBUG_IN();

  for (auto && RVE : RVEs) {
    /// update stress field before dumping
    RVE->assembleInternalForces();
    /// dump all the RVEs
    RVE->dumpRve(dump_nb);
  }

  AKANTU_DEBUG_OUT();
}

/* ------------------------------------------------------------------*/
template <UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::setTimeStep(Real time_step) {
  AKANTU_DEBUG_IN();

  for (auto && RVE : RVEs) {
    /// set time step to all the RVEs
    RVE->setTimeStep(time_step);
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------*/
template <UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::saveRVEsState(std::string & restart_dir) {
  AKANTU_DEBUG_IN();

  auto & mesh = this->model.getMesh();
  auto & g_ids = mesh.template getData<UInt>("global_ids", this->el_type);
  GhostType gt = _not_ghost;

  for (auto && data : enumerate(RVEs)) {
    auto proc_RVE_id = std::get<0>(data);
    UInt gl_RVE_id = g_ids(proc_RVE_id);
    auto & RVE = *(std::get<1>(data));

    std::stringstream file_stream;
    file_stream << restart_dir << "/RVE_" << std::to_string(gl_RVE_id)
                << "_state.txt";
    std::string RVE_state_file = file_stream.str();
    std::ofstream file_output;
    file_output.open(RVE_state_file.c_str());

    if (!file_output.is_open())
      AKANTU_EXCEPTION("Could not create the file " + RVE_state_file +
                       ", does its folder exist?");

    UInt nb_materials = RVE.getNbMaterials();

    for (UInt m = 0; m < nb_materials; ++m) {
      const Material & mat = RVE.getMaterial(m);

      for (auto pair : mat.getInternalVectorsReal()) {
        // set up precision
        constexpr auto max_digits10 = std::numeric_limits<Real>::max_digits10;
        file_output << std::setprecision(max_digits10);
        /// get the internal field
        const InternalField<Real> & internal_real = *pair.second;
        /// loop over all types in that material
        for (auto element_type : internal_real.elementTypes(gt)) {
          const Array<UInt> & elem_filter =
              mat.getElementFilter(element_type, gt);
          if (!elem_filter.size())
            continue;
          auto && internal_field_array = internal_real(element_type, gt);
          for (UInt i = 0; i < internal_field_array.size(); ++i) {
            for (UInt j = 0; j < internal_field_array.getNbComponent(); ++j) {
              file_output << internal_field_array(i, j) << std::endl;

              if (file_output.fail())
                AKANTU_EXCEPTION("Error in writing data to file " +
                                 RVE_state_file);
            }
          }
        }
      }
      for (auto pair : mat.getInternalVectorsUInt()) {
        // set up precision
        constexpr auto max_digits10 = std::numeric_limits<UInt>::max_digits10;
        file_output << std::setprecision(max_digits10);
        /// get the internal field
        const InternalField<UInt> & internal_uint = *pair.second;
        /// loop over all types in that material
        for (auto element_type : internal_uint.elementTypes(gt)) {
          const Array<UInt> & elem_filter =
              mat.getElementFilter(element_type, gt);
          if (!elem_filter.size())
            continue;
          const Array<UInt> & internal_field_array =
              internal_uint(element_type, gt);
          for (UInt i = 0; i < internal_field_array.size(); ++i) {
            for (UInt j = 0; j < internal_field_array.getNbComponent(); ++j) {
              file_output << internal_field_array(i, j) << std::endl;

              if (file_output.fail())
                AKANTU_EXCEPTION("Error in writing data to file " +
                                 RVE_state_file);
            }
          }
        }
      }
      for (auto pair : mat.getInternalVectorsBool()) {
        // set up precision
        constexpr auto max_digits10 = std::numeric_limits<bool>::max_digits10;
        file_output << std::setprecision(max_digits10);
        /// get the internal field
        const InternalField<bool> & internal_bool = *pair.second;
        /// loop over all types in that material
        for (auto element_type : internal_bool.elementTypes(gt)) {
          const Array<UInt> & elem_filter =
              mat.getElementFilter(element_type, gt);
          if (!elem_filter.size())
            continue;
          const Array<bool> & internal_field_array =
              internal_bool(element_type, gt);
          for (UInt i = 0; i < internal_field_array.size(); ++i) {
            for (UInt j = 0; j < internal_field_array.getNbComponent(); ++j) {
              file_output << internal_field_array(i, j) << std::endl;

              if (file_output.fail())
                AKANTU_EXCEPTION("Error in writing data to file " +
                                 RVE_state_file);
            }
          }
        }
      }
    }
    file_output.close();
    AKANTU_DEBUG_OUT();
  }
}

/* ------------------------------------------------------------------ */
template <UInt spatial_dimension>
void MaterialFE2<spatial_dimension>::loadRVEsState(std::string & restart_dir) {
  /// variables for parallel execution
  auto & mesh = this->model.getMesh();
  auto & g_ids = mesh.template getData<UInt>("global_ids", this->el_type);
  GhostType gt = _not_ghost;

  for (auto && data : enumerate(RVEs)) {
    auto proc_RVE_id = std::get<0>(data);
    UInt gl_RVE_id = g_ids(proc_RVE_id);
    auto & RVE = *(std::get<1>(data));

    /// open the input file
    std::stringstream file_stream;
    file_stream << restart_dir << "/RVE_" << std::to_string(gl_RVE_id)
                << "_state.txt";
    std::string RVE_state_file = file_stream.str();
    std::ifstream file_input;
    file_input.open(RVE_state_file.c_str());

    if (!file_input.is_open())
      AKANTU_EXCEPTION("Could not open file " + RVE_state_file);

    auto nb_materials = RVE.getNbMaterials();

    for (UInt m = 0; m < nb_materials; ++m) {
      const Material & mat = RVE.getMaterial(m);

      /// get the internal field that need to be set
      for (auto pair : mat.getInternalVectorsReal()) {
        /// get the internal field
        InternalField<Real> & internal_real = *pair.second;
        /// loop over all types in that material
        for (auto element_type : internal_real.elementTypes(gt)) {
          const Array<UInt> & elem_filter =
              mat.getElementFilter(element_type, gt);
          if (!elem_filter.size())
            continue;
          Array<Real> & internal_field_array = internal_real(element_type, gt);
          std::string line;
          for (UInt i = 0; i < internal_field_array.size(); ++i) {
            for (UInt j = 0; j < internal_field_array.getNbComponent(); ++j) {
              std::getline(file_input, line);

              if (file_input.fail())
                AKANTU_EXCEPTION("Could not read data from file " +
                                 RVE_state_file);

              std::stringstream sstr(line);
              sstr >> internal_field_array(i, j);
            }
          }
        }
      }

      /// get the internal field that need to be set
      for (auto pair : mat.getInternalVectorsUInt()) {
        /// get the internal field
        InternalField<UInt> & internal_uint = *pair.second;
        /// loop over all types in that material
        for (auto element_type : internal_uint.elementTypes(gt)) {
          const Array<UInt> & elem_filter =
              mat.getElementFilter(element_type, gt);
          if (!elem_filter.size())
            continue;
          Array<UInt> & internal_field_array = internal_uint(element_type, gt);
          std::string line;
          for (UInt i = 0; i < internal_field_array.size(); ++i) {
            for (UInt j = 0; j < internal_field_array.getNbComponent(); ++j) {
              std::getline(file_input, line);

              if (file_input.fail())
                AKANTU_EXCEPTION("Could not read data from file " +
                                 RVE_state_file);

              std::stringstream sstr(line);
              sstr >> internal_field_array(i, j);
            }
          }
        }
      }

      /// get the internal field that need to be set
      for (auto pair : mat.getInternalVectorsBool()) {
        /// get the internal field
        InternalField<bool> & internal_bool = *pair.second;
        /// loop over all types in that material
        for (auto element_type : internal_bool.elementTypes(gt)) {
          const Array<UInt> & elem_filter =
              mat.getElementFilter(element_type, gt);
          if (!elem_filter.size())
            continue;
          Array<bool> & internal_field_array = internal_bool(element_type, gt);
          std::string line;
          for (UInt i = 0; i < internal_field_array.size(); ++i) {
            for (UInt j = 0; j < internal_field_array.getNbComponent(); ++j) {
              std::getline(file_input, line);

              if (file_input.fail())
                AKANTU_EXCEPTION("Could not read data from file " +
                                 RVE_state_file);

              std::stringstream sstr(line);
              sstr >> internal_field_array(i, j);
            }
          }
        }
      }
    }
    /// close the file
    file_input.close();
  }
}

INSTANTIATE_MATERIAL(material_FE2, MaterialFE2);

} // namespace akantu
