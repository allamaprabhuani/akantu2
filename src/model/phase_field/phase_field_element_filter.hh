#include "group_manager.hh"
#include "phase_field_model.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {

class PhaseFieldElementFilter : public GroupManager::ClusteringFilter {
public:
  PhaseFieldElementFilter(const PhaseFieldModel & model,
                          const Real max_damage = 1.)
      : model(model), is_unbroken(max_damage) {}

  bool operator()(const Element & el) const override {

    const Array<Idx> & mat_indexes =
        model.getPhaseFieldByElement(el.type, el.ghost_type);
    const Array<Idx> & mat_loc_num =
        model.getPhaseFieldLocalNumbering(el.type, el.ghost_type);

    const auto & mat = model.getPhaseField(mat_indexes(el.element));

    Idx el_index = mat_loc_num(el.element);
    Int nb_quad_per_element =
        model.getFEEngine("PhaseFieldFEEngine")
            .getNbIntegrationPoints(el.type, el.ghost_type);

    const Array<Real> & damage_array = mat.getDamage(el.type, el.ghost_type);

    AKANTU_DEBUG_ASSERT(nb_quad_per_element * el_index < damage_array.size(),
                        "This quadrature point is out of range");

    const Real * element_damage =
        damage_array.data() + nb_quad_per_element * el_index;

    UInt unbroken_quads = std::count_if(
        element_damage, element_damage + nb_quad_per_element, is_unbroken);

    return (unbroken_quads > 0);
  }

private:
  struct IsUnbrokenFunctor {
    IsUnbrokenFunctor(const Real & max_damage) : max_damage(max_damage) {}
    bool operator()(const Real & x) const { return x < max_damage; }
    const Real max_damage;
  };

  const PhaseFieldModel & model;
  const IsUnbrokenFunctor is_unbroken;
};

}
