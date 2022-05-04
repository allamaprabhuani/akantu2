
/* -------------------------------------------------------------------------- */
#include "constitutive_law_diffusion.hh"
#include "poisson_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
ConstitutiveLawDiffusion::ConstitutiveLawDiffusion(PoissonModel & model,
					 const ID & id)
  : ConstitutiveLaw(model, id), was_stiffness_assembled(false) {
  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
ConstitutiveLawDiffusion::ConstitutiveLawDiffusion(PoissonModel & model,
					 UInt /*a_dim*/, const Mesh & mesh,
					 FEEngine & fe_engine, const ID & id)
    : ConstitutiveLaw(model, dim, mesh, fe_engine, id), was_stiffness_assembled(false) {
  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ConstitutiveLawDiffusion<dim>::initialize() {
  this->registerParam("diffusivity", diffusivity, _pat_readable, "Diffusion Coefficient");
}

/* -------------------------------------------------------------------------- */
void ConstitutiveLawDiffusion::initMaterial() {
  AKANTU_DEBUG_IN();
  ConstitutiveLaw::initMaterial();

  this->updateInternalParameters();
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void ConstitutiveLawDiffusion::computeFlux(ElementType el_type,
					   GhostType ghost_type) {


  computeDiffusivityOnQuad(el_type, ghost_type);

  // concentration gradient at quadrature points
  auto & gradient = gradient_dof(el_type, ghost_type);
  this->getFEEngine().gradientOnIntegrationPoints(model.getDof(), gradient, 1,
						  el_type, ghost_type);
  
  for (auto && values :
         zip(make_view(diffusivity_on_qpoints(el_type, ghost_type),
                       spatial_dimension, spatial_dimension),
             make_view(gradient, spatial_dimension),
             make_view(flux_dof(el_type, ghost_type),
                       spatial_dimension))) {
    const auto & D = std::get<0>(values);
    const auto & BC = std::get<1>(values);
    auto & d_BC = std::get<2>(values);
    
    d_BT.mul<false>(D, BC);
  }
  
}

/* -------------------------------------------------------------------------- */
void ConstitutiveLawDiffusion::computeTangentModuli(ElementType el_type,
						    Array<Real> & tangent_matrix,
						    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  for (auto && tuple :
	 make_view(tangent_matrix, spatial_dimension, spatial_dimension)) {
    auto & D = std::get<0>(tuple);
    D = this->diffusivity;
  }
  
  this->was_stiffness_assembled = true;

  AKANTU_DEBUG_OUT();
}

  
/* -------------------------------------------------------------------------- */
void ConstitutiveLawDiffusion::computeDiffusivityOnQuad(ElementType el_type,
							GhostType ghost_type) {
  
  AKANTU_DEBUG_IN();
  auto & diff_coef = diffusivity_on_qpoints(el_type, ghost_type);
  for (auto && tuple :
	 make_view(diff_coef, spatial_dimension, spatial_dimension)) {
    auto & D = std::get<0>(tuple);
    D = this->diffusivity;
  }
  AKANTU_DEBUG_OUT();
}
  
