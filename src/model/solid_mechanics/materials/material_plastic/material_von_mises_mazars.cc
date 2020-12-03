/* -------------------------------------------------------------------------- */
#include "material_von_mises_mazars.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

 /* -------------------------------------------------------------------------- */
template <UInt dim>
MaterialVonMisesMazars<dim>::MaterialVonMisesMazars(
    SolidMechanicsModel & model, const ID & id)
  : MaterialLinearIsotropicHardening<dim>(model, id),
  MaterialMazars<dim>(model, id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialVonMisesMazars<spatial_dimension>::
   MaterialVonMisesMazars(SolidMechanicsModel & model, UInt dim,
			  const Mesh & mesh, FEEngine & fe_engine,
			  const ID & id) 
    : MaterialLinearIsotropicHardening<spatial_dimension>(model, id),
     MaterialMazars<spatial_dimension>(model, id) {}


/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialLinearIsotropicHardening<spatial_dimension>::computeStress(
    ElementType el_type, GhostType ghost_type) {

  AKANTU_DEBUG_IN();

  MaterialThermal<spatial_dimension>::computeStress(el_type, ghost_type);
  // infinitesimal and finite deformation
  auto sigma_th_it = this->sigma_th(el_type, ghost_type).begin();

  auto previous_sigma_th_it =
      this->sigma_th.previous(el_type, ghost_type).begin();

  auto previous_gradu_it = this->gradu.previous(el_type, ghost_type)
                               .begin(spatial_dimension, spatial_dimension);

  auto previous_stress_it = this->stress.previous(el_type, ghost_type)
                                .begin(spatial_dimension, spatial_dimension);

  auto inelastic_strain_it = this->inelastic_strain(el_type, ghost_type)
                                 .begin(spatial_dimension, spatial_dimension);

  auto previous_inelastic_strain_it =
      this->inelastic_strain.previous(el_type, ghost_type)
          .begin(spatial_dimension, spatial_dimension);

  auto iso_hardening_it = this->iso_hardening(el_type, ghost_type).begin();

  auto previous_iso_hardening_it =
      this->iso_hardening.previous(el_type, ghost_type).begin();

  Real * dam = this->damage(el_type, ghost_type).storage();
  
  //
  // Finite Deformations
  //
  if (this->finite_deformation) {
    auto previous_piola_kirchhoff_2_it =
        this->piola_kirchhoff_2.previous(el_type, ghost_type)
            .begin(spatial_dimension, spatial_dimension);

    auto green_strain_it = this->green_strain(el_type, ghost_type)
                               .begin(spatial_dimension, spatial_dimension);

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

    auto & inelastic_strain_tensor = *inelastic_strain_it;
    auto & previous_inelastic_strain_tensor = *previous_inelastic_strain_it;
    auto & previous_grad_u = *previous_gradu_it;
    auto & previous_sigma = *previous_piola_kirchhoff_2_it;

    auto & green_strain = *green_strain_it;
    this->template gradUToE<spatial_dimension>(grad_u, green_strain);
    Matrix<Real> previous_green_strain(spatial_dimension, spatial_dimension);
    this->template gradUToE<spatial_dimension>(previous_grad_u,
                                                         previous_green_strain);
    Matrix<Real> F_tensor(spatial_dimension, spatial_dimension);
    this->template gradUToF<spatial_dimension>(grad_u, F_tensor);

    MaterialLinearIsotropicHardening<spatial_dimension>::computeStressOnQuad(
			 green_strain, previous_green_strain, sigma,
			 previous_sigma, inelastic_strain_tensor,
                        previous_inelastic_strain_tensor, *iso_hardening_it,
                        *previous_iso_hardening_it, *sigma_th_it,
                        *previous_sigma_th_it, F_tensor);

    ++sigma_th_it;
    ++inelastic_strain_it;
    ++iso_hardening_it;
    ++previous_sigma_th_it;
    //++previous_stress_it;
    ++previous_gradu_it;
    ++green_strain_it;
    ++previous_inelastic_strain_it;
    ++previous_iso_hardening_it;
    ++previous_piola_kirchhoff_2_it;

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  }
  // Infinitesimal deformations
  else {
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

    Real Ehat = 0;
    
    auto & inelastic_strain_tensor = *inelastic_strain_it;
    auto & previous_inelastic_strain_tensor = *previous_inelastic_strain_it;
    auto & previous_grad_u = *previous_gradu_it;
    auto & previous_sigma = *previous_stress_it;

    MaterialLinearIsotropicHardening<spatial_dimension>::computeStressOnQuad(
        grad_u, previous_grad_u, sigma, previous_sigma, inelastic_strain_tensor,
        previous_inelastic_strain_tensor, *iso_hardening_it,
	*previous_iso_hardening_it, *sigma_th_it, *previous_sigma_th_it);

    Matrix<Real> grad_u_elastic(grad_u);
    grad_u_elastic -= inelastic_strain_tensor;
    
    MaterialMazars<spatial_dimension>::computeStressOnQuad(grad_u_elastic, sigma, *dam, Ehat);
    ++dam;
    
    ++sigma_th_it;
    ++inelastic_strain_it;
    ++iso_hardening_it;
    ++previous_sigma_th_it;
    ++previous_stress_it;
    ++previous_gradu_it;
    ++previous_inelastic_strain_it;
    ++previous_iso_hardening_it;

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  }

  AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(mazars, MaterialMazars);

}
