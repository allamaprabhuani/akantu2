#include "diffusion_law.hh"

#ifndef AKANTU_HEAT_DIFFUSION_HH
#define AKANTU_HEAT_DIFFUSION_HH

namespace akantu {

template <Int dim> class HeatDiffusion : public DiffusionLaw {
  using Parent = DiffusionLaw;

public:
  HeatDiffusion(DiffusionModel & model, const ID & id,
                const ID & fe_engine_id = "");

  void updateInternalParameters() override;
  void computeDiffusivityGradUOnQuadPoints(ElementType type,
                                           GhostType ghost_type) override;
  void computeDiffusivityOnQuadPoints(ElementType type,
                                      GhostType ghost_type) override;

  decltype(auto) getArguments(ElementType type, GhostType ghost_type) {
    return DiffusionLaw::getArguments<dim>(type, ghost_type);
  }

  [[nodiscard]] Real getStableTimeStep(Real element_size) override;

  [[nodiscard]] Real getRho() const override { return density * capacity; }

  [[nodiscard]] Real getEnergy(const ID & energy_id) override;
  [[nodiscard]] Real getEnergy(const ID & energy_id,
                               const Element & element) override;

private:
  auto getThermalEnergy(const Element & element) -> Real;
  auto getThermalEnergy() -> Real;

  template <class iterator, class t_iterator>
  void getThermalEnergy(iterator Eth, t_iterator T_it, t_iterator T_end) const;

private:
  Real density{0.};

  /// capacity
  Real capacity{0.};

  // conductivity matrix
  Matrix<Real> conductivity;

  // linear variation of the conductivity (for temperature dependent
  // conductivity)
  Real conductivity_variation{0.};

  // reference temperature for the interpretation of temperature variation
  Real T_ref{0.};
};

} // namespace akantu

#endif // AKANTU_HEAT_DIFFUSION_HH
