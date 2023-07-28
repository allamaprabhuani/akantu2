#include "heat_diffusion.hh"
#include "heat_transfer_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
HeatDiffusion<dim>::HeatDiffusion(DiffusionModel & model, const ID & id,
                                  const ID & fe_engine_id)
    : DiffusionLaw(model, id, fe_engine_id) {
  this->registerParam("density", density, _pat_parsmod);
  this->registerParam("conductivity", conductivity, _pat_parsmod);
  this->registerParam("conductivity_variation", conductivity_variation, 0.,
                      _pat_parsmod);
  this->registerParam("temperature_reference", T_ref, 0., _pat_parsmod);
  this->registerParam("capacity", capacity, _pat_parsmod);
}

/* -------------------------------------------------------------------------- */
template <Int dim> void HeatDiffusion<dim>::updateInternalParameters() {
  Matrix<Real> tmp = conductivity.block<dim, dim>(0, 0);
  conductivity = tmp;
  this->diffusivity.set(conductivity);
  Parent::updateInternalParameters();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void HeatDiffusion<dim>::computeDiffusivityGradUOnQuadPoints(
    ElementType type, GhostType ghost_type) {
  this->computeDiffusivityOnQuadPoints(type, ghost_type);

  for (auto && args : getArguments(type, ghost_type)) {
    const auto & C = args["diffusivity"_n];
    const auto & BT = args["∇u"_n];
    auto & k_BT = args["D∇u"_n];

    k_BT = C * BT;
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void HeatDiffusion<dim>::computeDiffusivityOnQuadPoints(ElementType type,
                                                        GhostType ghost_type) {
  auto temperature_release = this->getHandler().getDiffusionRelease();
  auto & diffusivity_release = this->diffusivity.getRelease(type, ghost_type);
  if (diffusivity_release != -1 and
      diffusivity_release == temperature_release) {
    return;
  }

  for (auto && type : element_filter.elementTypes(dim, ghost_type)) {
    Array<Real> temperature_on_qpoints(0, 1);

    // compute the temperature on quadrature points
    this->getFEEngine().interpolateOnIntegrationPoints(
        getHandler().getDiffusion(), temperature_on_qpoints, 1, type,
        ghost_type, element_filter(type, ghost_type));

    for (auto && [C, T] :
         zip(make_view<dim, dim>(this->diffusivity(type, ghost_type)),
             temperature_on_qpoints)) {
      C = conductivity;

      C.array() += conductivity_variation * (T - T_ref);
    }
  }

  diffusivity_release = temperature_release;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real HeatDiffusion<dim>::getStableTimeStep(Real element_size) {
  Real conductivity_max{};
  Vector<Real> ce;
  conductivity.eig(ce);
  for (auto c : ce) {
    conductivity_max = std::max(c, conductivity_max);
  }

  Real min_dt = 2. * element_size * element_size / 4. * this->density *
                capacity / conductivity_max;

  return min_dt;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class iterator, class t_iterator>
void HeatDiffusion<dim>::getThermalEnergy(iterator Eth, t_iterator T_it,
                                          t_iterator T_end) const {
  for (; T_it != T_end; ++T_it, ++Eth) {
    *Eth = capacity * density * *T_it;
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
auto HeatDiffusion<dim>::getThermalEnergy(const Element & element) -> Real {
  AKANTU_DEBUG_IN();

  auto nb_quadrature_points =
      getFEEngine().getNbIntegrationPoints(element.type);
  Vector<Real> Eth_on_quarature_points(nb_quadrature_points);

  Array<Real> temperature_interpolated(0, 1);
  Array<Idx> filter(1, 1);
  filter(0) = element.element;

  this->getFEEngine().interpolateOnIntegrationPoints(
      getHandler().getDiffusion(), temperature_interpolated, 1, element.type,
      element.ghost_type, filter);

  auto T_it = temperature_interpolated.begin();
  auto T_end = T_it + nb_quadrature_points;

  getThermalEnergy(Eth_on_quarature_points.data(), T_it, T_end);

  return getFEEngine().integrate(Eth_on_quarature_points, element);
}

/* -------------------------------------------------------------------------- */
template <Int dim> auto HeatDiffusion<dim>::getThermalEnergy() -> Real {
  Real Eth = 0;

  auto & fem = getFEEngine();

  for (auto && type : element_filter.elementTypes(dim, _not_ghost)) {
    auto nb_element = element_filter(type).size();
    auto nb_quadrature_points = fem.getNbIntegrationPoints(type, _not_ghost);
    Array<Real> Eth_per_quad(nb_element * nb_quadrature_points, 1);

    Array<Real> temperature_interpolated(0, 1);
    // compute the temperature on quadrature points
    this->getFEEngine().interpolateOnIntegrationPoints(
        getHandler().getDiffusion(), temperature_interpolated, 1, type,
        _not_ghost, element_filter(type));

    auto T_it = temperature_interpolated.begin();
    auto T_end = temperature_interpolated.end();
    getThermalEnergy(Eth_per_quad.begin(), T_it, T_end);

    Eth += fem.integrate(Eth_per_quad, type);
  }

  return Eth;
}

/* -------------------------------------------------------------------------- */
template <Int dim> Real HeatDiffusion<dim>::getEnergy(const ID & energy_id) {
  if (energy_id == "thermal") {
    return getThermalEnergy();
  }

  return Parent::getEnergy(energy_id);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real HeatDiffusion<dim>::getEnergy(const ID & energy_id,
                                   const Element & element) {
  if (energy_id == "thermal") {
    return getThermalEnergy(element);
  }

  return Parent::getEnergy(energy_id, element);
}

/* -------------------------------------------------------------------------- */
template class HeatDiffusion<1>;
template class HeatDiffusion<2>;
template class HeatDiffusion<3>;

const bool diffusion_law_is_alocated_heat_diffusion [[maybe_unused]] =
    instantiateDiffusionLaw<HeatDiffusion, HeatTransferModel>("heat_diffusion");

} // namespace akantu
