
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
ConstitutiveLawDiffusion::ConstitutiveLawDiffusion(PoissonModel & model,
						   UInt dim, const Mesh & mesh,
						   FEEngine & fe_engine, const ID & id)
  : ConstitutiveLaw(model, dim, mesh, fe_engine, id), was_stiffness_assembled(false) {
  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
  void ConstitutiveLawDiffusion::initialize() {
  this->registerParam("diffusivity", diffusivity, _pat_readable, "Diffusion Coefficient");
}

/* -------------------------------------------------------------------------- */
void ConstitutiveLawDiffusion::initConstitutiveLaw() {
  AKANTU_DEBUG_IN();
  ConstitutiveLaw::initConstitutiveLaw();

  this->updateInternalParameters();
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void ConstitutiveLawDiffusion::computeFlux(ElementType el_type,
					   GhostType ghost_type) {


  Matrix<Real> identity(spatial_dimension, spatial_dimension);
  identity.eye();
  auto D = identity * this->diffusivity;
  
  // concentration gradient at quadrature points
  auto & gradient = gradient_dof(el_type, ghost_type);
  fem.gradientOnIntegrationPoints(model.getDof(), gradient, 1,
				  el_type, ghost_type);
  
  for (auto && values :
         zip(make_view(gradient, spatial_dimension),
             make_view(flux_dof(el_type, ghost_type), spatial_dimension))) {
    const auto & BC = std::get<0>(values);
    auto & d_BC = std::get<1>(values);
  
    d_BC.mul<false>(D, BC);
  }
  
}

/* -------------------------------------------------------------------------- */
void ConstitutiveLawDiffusion::computeTangentModuli(ElementType el_type,
						    Array<Real> & tangent_matrix,
						    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  
  Matrix<Real> identity(spatial_dimension, spatial_dimension);
  identity.eye();
  
  for (auto && tuple :
	 make_view(tangent_matrix, spatial_dimension, spatial_dimension)) {
    tuple = identity * this->diffusivity;
  }
  
  this->was_stiffness_assembled = true;

  AKANTU_DEBUG_OUT();
}

  
  
}
