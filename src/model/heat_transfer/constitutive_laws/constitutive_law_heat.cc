
/* -------------------------------------------------------------------------- */
#include "constitutive_law_heat.hh"
#include "poisson_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
ConstitutiveLawHeat::ConstitutiveLawHeat(PoissonModel & model,
					 const ID & id)
    : ConstitutiveLaw(model, id), was_stiffness_assembled(false) {
  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
ConstitutiveLawHeat::ConstitutiveLawHeat(PoissonModel & model,
					 UInt /*a_dim*/, const Mesh & mesh,
					 FEEngine & fe_engine, const ID & id)
    : ConstitutiveLaw(model, dim, mesh, fe_engine, id), was_stiffness_assembled(false) {
  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ConstitutiveLawHeat<dim>::initialize() {
  this->registerParam("density", density, _pat_readable,
                      "density");
  this->registerParam("capacity", capacity, _pat_readable, "Heat capacity");
  this->registerParam("conductivity", conductivity, _pat_readable, "Heat conductivity");
}

/* -------------------------------------------------------------------------- */
void ConstitutiveLawHeat::initMaterial() {
  AKANTU_DEBUG_IN();
  ConstitutiveLaw::initMaterial();

  this->updateInternalParameters();
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void ConstitutiveLawHeat::computeFlux(ElementType el_type,
				      GhostType ghost_type) {


  computeConductivityOnQuad(el_type, ghost_type);

  // temperature gradient at quadrature points
  auto & gradient = gradient_dof(el_type, ghost_type);
  this->getFEEngine().gradientOnIntegrationPoints(model.getDof(), gradient, 1,
						  el_type, ghost_type);
  
  for (auto && values :
         zip(make_view(conductivity_on_qpoints(el_type, ghost_type),
                       spatial_dimension, spatial_dimension),
             make_view(gradient, spatial_dimension),
             make_view(flux_dof(el_type, ghost_type),
                       spatial_dimension))) {
      const auto & C = std::get<0>(values);
      const auto & BT = std::get<1>(values);
      auto & k_BT = std::get<2>(values);

      k_BT.mul<false>(C, BT);
  }
    
}

/* -------------------------------------------------------------------------- */
  void ConstitutiveLawHeat::computeConductivityOnQuad(ElementType el_type,
						      GhostType ghost_type) {
  // if already computed once check if need to compute
  if (not initial_conductivity[ghost_type]) {
    // if temperature did not change, conductivity will not vary
    if (dof_release == conductivity_release[ghost_type]) {
      return;
    }
    // if conductivity_variation is 0 no need to recompute
    if (conductivity_variation == 0.) {
      return;
    }
  }

  auto & temperature_interpolated = temperature_on_qpoints(el_type, ghost_type);

  // compute the temperature on quadrature points
  this->getFEEngine().interpolateOnIntegrationPoints(
			model.getDof(), temperature_interpolated, 1, el_type, ghost_type);

  auto & cond = conductivity_on_qpoints(el_type, ghost_type);
  for (auto && tuple :
         zip(make_view(cond, spatial_dimension, spatial_dimension),
             temperature_interpolated)) {
    auto & C = std::get<0>(tuple);
    auto & T = std::get<1>(tuple);
    C = conductivity;

    Matrix<Real> variation(spatial_dimension, spatial_dimension,
			   conductivity_variation * (T - T_ref));
    // @TODO: Guillaume are you sure ? why due you compute variation then ?
    C += conductivity_variation;
  }
    

  conductivity_release[ghost_type] = dof_release;
  initial_conductivity[ghost_type] = false;

  AKANTU_DEBUG_OUT();
}
  
