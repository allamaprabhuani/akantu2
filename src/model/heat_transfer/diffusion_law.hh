#include "constitutive_law.hh"

#ifndef AKANTU_DIFFUSION_LAW_HH
#define AKANTU_DIFFUSION_LAW_HH

namespace akantu {
class DiffusionModel;
class DiffusionLaw;
} // namespace akantu

namespace akantu {

using DiffusionFactory =
    Factory<DiffusionLaw, ID, Int, const ID &, DiffusionModel &, const ID &>;

/* -------------------------------------------------------------------------- */
class DiffusionLaw : public ConstitutiveLaw<DiffusionModel> {
  using Parent = ConstitutiveLaw<DiffusionModel>;

public:
  DiffusionLaw(DiffusionModel & model, const ID & id = "diffusion_law",
               const ID & fe_engine_id = "");

  virtual void assembleInternalFlow(GhostType ghost_type);
  virtual void assembleDiffusivityMatrix();

  virtual void computeGradU(ElementType type, GhostType ghost_type);

  virtual void computeDiffusivityGradU(GhostType ghost_type = _not_ghost);
  virtual void computeDiffusivityGradUOnQuadPoints(ElementType /*type*/,
                                                   GhostType /*ghost_type*/) {
    AKANTU_TO_IMPLEMENT();
  }

  virtual void computeDiffusivityOnQuadPoints(ElementType /*type*/,
                                              GhostType /*ghost_type*/) {
    AKANTU_TO_IMPLEMENT();
  }

  template <Int dim>
  decltype(auto) getArguments(ElementType type, GhostType ghost_type) {
    return zip("∇u"_n = make_view<dim>(grad_u(type, ghost_type)),
               "D∇u"_n = make_view<dim>(d_grad_u(type, ghost_type)),
               "diffusivity"_n =
                   make_view<dim, dim>(diffusivity(type, ghost_type)));
  }

  [[nodiscard]] virtual Real getEnergy(const ID & energy_id);
  [[nodiscard]] virtual Real getEnergy(const ID & energy_id,
                                       const Element & element);

  [[nodiscard]] virtual Real getStableTimeStep(Real /*element_size*/) {
    return 0;
  }

  [[nodiscard]] virtual Real getRho() const { return 1; }

  /* ------------------------------------------------------------------------ */
  [[nodiscard]] Int getNbData(const Array<Element> & elements,
                              const SynchronizationTag & tag) const override;
  void packData(CommunicationBuffer & buffer, const Array<Element> & elements,
                const SynchronizationTag & tag) const override;
  void unpackData(CommunicationBuffer & buffer, const Array<Element> & elements,
                  const SynchronizationTag & tag) override;
  /* ------------------------------------------------------------------------ */

  /// static method to retrieve the diffusion factory
  static DiffusionFactory & getFactory();

protected:
  InternalField<Real> & grad_u;
  InternalField<Real> & d_grad_u;
  InternalField<Real> & diffusivity;
};

namespace {
  template <
      template <Int> class Law, class Model_,
      std::enable_if_t<std::is_base_of_v<DiffusionModel, Model_>> * = nullptr>
  bool instantiateDiffusionLaw(const ID & id) {
    return DiffusionFactory::getInstance().registerAllocator(
        id,
        [](Int dim, const ID & name, DiffusionModel & model, const ID & id) {
          if (not aka::is_of_type<Model_>(model)) {
            AKANTU_EXCEPTION("The diffusion law "
                             << name << " works only with model of type "
                             << debug::demangle(typeid(Model_).name()));
          }
          return tuple_dispatch<AllSpatialDimensions>(
              [&](auto && _) -> std::unique_ptr<DiffusionLaw> {
                constexpr auto && dim_ = aka::decay_v<decltype(_)>;
                return std::make_unique<Law<dim_>>(model, id);
              },
              dim);
        });
  }
} // namespace

} // namespace akantu
#endif // AKANTU_DIFFUSION_LAW_HH
