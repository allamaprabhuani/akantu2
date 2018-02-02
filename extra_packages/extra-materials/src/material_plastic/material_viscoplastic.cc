/**
 * @file   material_viscoplastic.cc
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 *
 *
 * @brief  Specialization of the material class for isotropic viscoplastic (small deformation)
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_viscoplastic.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt dim>
MaterialViscoPlastic<dim>::MaterialViscoPlastic(SolidMechanicsModel & model, const ID & id) :
  Material(model, id),
  MaterialPlastic<dim>(model, id) {
  AKANTU_DEBUG_IN();

  this->registerParam("rate", rate, 0., _pat_parsable | _pat_modifiable, "Rate sensitivity component");
  this->registerParam("edot0", edot0, 0., _pat_parsable | _pat_modifiable, "Reference strain rate");
  this->registerParam("ts", ts, 0., _pat_parsable | _pat_modifiable, "Time Step");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialViscoPlastic<spatial_dimension>::computeStress(ElementType el_type, GhostType ghost_type) {
    AKANTU_DEBUG_IN();

    Real * iso_hardening = this->iso_hardening(el_type, ghost_type).storage();

    Array<Real>::matrix_iterator previous_grad_u_it =
      this->gradu.previous(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);

    Array<Real>::matrix_iterator previous_sigma_it =
      this->stress.previous(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);

    Array<Real>::matrix_iterator inelastic_strain_it =
      this->inelastic_strain(el_type, ghost_type).begin(spatial_dimension,spatial_dimension);

    Array<Real>::matrix_iterator previous_inelastic_strain_it =
      this->inelastic_strain.previous(el_type, ghost_type).begin(spatial_dimension,spatial_dimension);

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

    computeStressOnQuad(grad_u, *previous_grad_u_it,
                        sigma, *previous_sigma_it,
                        *inelastic_strain_it, *previous_inelastic_strain_it,
                        *iso_hardening);

    ++inelastic_strain_it;
    ++iso_hardening;
    ++previous_grad_u_it;

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialViscoPlastic<spatial_dimension>::computeTangentModuli(__attribute__((unused)) const ElementType & el_type,
                                                                   Array<Real> & tangent_matrix,
                                                                   __attribute__((unused)) GhostType ghost_type) {
    AKANTU_DEBUG_IN();


    Array<Real>::matrix_iterator previous_sigma_it =
       this->stress.previous(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);

    Array<Real>::matrix_iterator previous_strain_it =
       this->gradu.previous(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);

    Real * iso_hardening= this->iso_hardening(el_type, ghost_type).storage();


    MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);

    Matrix<Real> & previous_grad_u = *previous_strain_it;
    Matrix<Real> & previous_sigma_tensor = *previous_sigma_it;

    computeTangentModuliOnQuad(tangent, grad_u, previous_grad_u, sigma_tensor, previous_sigma_tensor, *iso_hardening);
    ++previous_sigma_it;
    ++previous_strain_it;
    ++iso_hardening;
    MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

    AKANTU_DEBUG_OUT();
}

INSTANTIATE_MATERIAL(MaterialViscoPlastic);

__END_AKANTU__
