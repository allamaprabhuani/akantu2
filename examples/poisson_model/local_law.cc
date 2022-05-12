
/* -------------------------------------------------------------------------- */
#include "local_law.hh"
#include "poisson_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
LocalLaw::LocalLaw(PoissonModel & model,
					 const ID & id)
  : ConstitutiveLaw(model, id), was_stiffness_assembled(false) {
  AKANTU_DEBUG_IN();
  this->registerParam("permitivity", epsilon, Real(0.),
		      _pat_parsable | _pat_modifiable, "Medium Permitivity");
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void LocalLaw::initConstitutiveLaw() {
  AKANTU_DEBUG_IN();
  ConstitutiveLaw::initConstitutiveLaw();

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void LocalLaw::computeFlux(ElementType el_type,
					   GhostType ghost_type) {


  Matrix<Real> identity(spatial_dimension, spatial_dimension);
  identity.eye();
  auto D = identity * this->epsilon;
  D *= -1;
  
  auto & potential_gradient = gradient_dof(el_type, ghost_type);
  
  for (auto && values :
         zip(make_view(potential_gradient, spatial_dimension),
             make_view(flux_dof(el_type, ghost_type), spatial_dimension))) {
    const auto & BC = std::get<0>(values);
    auto & d_BC = std::get<1>(values);

    d_BC.mul<false>(D, BC);
  }
}

/* -------------------------------------------------------------------------- */
void LocalLaw::computeTangentModuli(ElementType el_type,
						    Array<Real> & tangent_matrix,
						    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  
  Matrix<Real> identity(spatial_dimension, spatial_dimension);
  identity.eye();
  
  for (auto && tuple :
	 make_view(tangent_matrix, spatial_dimension, spatial_dimension)) {
    tuple = identity * this->epsilon;
  }
  
  this->was_stiffness_assembled = true;

  AKANTU_DEBUG_OUT();
}


static bool constiitutive_law_is_alocated_local_law [[gnu::unused]] =
    ConstitutiveLawFactory::getInstance().registerAllocator(
        "local_law",
        [](const ID &, PoissonModel & model,
           const ID & id) -> std::unique_ptr<ConstitutiveLaw> {
          return std::make_unique<LocalLaw>(model, id);
        });

} // namespace akantu
